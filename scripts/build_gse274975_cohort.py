"""Build the GSE274975 (Oncobox NSCLC ICI) cohort for the executed E1.

Data sources (all fetched 2026-09-19, URLs recorded in raw/):
  - Expression: GEO supplementary GSE274975_raw_counts.tsv.gz
    (58,233 Ensembl gene IDs x 61 samples, raw counts, no STAR
    special rows -- verified).
  - Clinical: Luppov et al., Front Immunol 2024 (PMC11669362)
    Supplementary Table 1 (paper_tables/Table1.xlsx): RECIST response
    (CR/PR/SD/PD), PFS time (months), plus histotype / therapy / PD-L1
    covariates.  The paper defines responders as CR+PR -- identical to
    the SPATBench RECIST convention.
  - Sample identity: the GEO sample titles (OB_pat_LuC_N), the
    Table-1 Sample_IDs (OB_pat_LuC_N) and the counts column names
    (LuC-N_..., LuC_N_..., Luc-N_..., LuC-Nc/CG-batch variants) share
    ONE patient numbering space.  Verified: perfect 61x61x61
    bijection on extracted patient numbers + histotype concordance
    57/58 hard-agreement (GEO tissue field vs Table-1 Histotype;
    the single discordance is patient 128 AD-vs-SQ, 7 further rows
    are Table-1 "Mixed phenotype" vs a dominant histotype -- not
    contradictions).  Mapping saved in raw/id_map.json.

Processing (pre-declared e1_gse274975_replication.py convention):
  1. Ensembl -> HGNC via Ensembl REST /lookup/id (batch 900, retry,
     JSON cache scripts/_gse274975_ensembl_map.json); unmapped genes
     dropped; duplicate symbols summed (raw counts are additive).
  2. log2(CPM + 1)  (FPKM not recoverable without gene lengths).
  3. obs: response (CR/PR = 1, SD/PD = 0), pfs_months, status01
     (Table-1 "Response status available" column, semantics undefined
     in the paper -- stored for provenance, NOT used for endpoints),
     histotype, therapy, pdl1_status, age, gender, tumor_site,
     patient_num, geo_title, t1_id, counts_column.
  4. Splits: stratified k = min(5, n_pos, n_neg) = 5, seed 42,
     saved in the standard {"folds": [...]} format with SHA-1.

Output:
  data/cohorts/GSE274975/processed/GSE274975_HGNC.h5ad
  data/cohorts/GSE274975/processed/GSE274975_splits.json
  data/cohorts/GSE274975/processed/build_qc.json
"""
import gzip
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data/cohorts/GSE274975/raw"
PROC = REPO / "data/cohorts/GSE274975/processed"
CACHE = REPO / "scripts/_gse274975_ensembl_map.json"
SEED = 42
HEADERS = {"Content-Type": "application/json", "Accept": "application/json",
           "User-Agent": "spatbench-ecohort-build"}


def ensembl_map(ids):
    """Ensembl REST /lookup/id in batches of 900 with retry + cache."""
    if CACHE.exists():
        emap = json.load(open(CACHE))
        print(f"  cache hit: {len(emap)} entries")
        return emap
    out = {}
    for i in range(0, len(ids), 900):
        batch = ids[i:i + 900]
        body = json.dumps({"ids": batch}).encode()
        for attempt in range(6):
            try:
                req = urllib.request.Request(
                    "https://rest.ensembl.org/lookup/id", data=body,
                    headers=HEADERS)
                with urllib.request.urlopen(req, timeout=90) as r:
                    res = json.load(r)
                for e, info in res.items():
                    if info and info.get("display_name"):
                        out[e] = info["display_name"]
                break
            except Exception as e:
                print(f"  batch {i} attempt {attempt} failed: {str(e)[:60]}")
                time.sleep(2 + attempt * 3)
        if (i // 900) % 10 == 9:
            print(f"  ...{i + len(batch)}/{len(ids)} "
                  f"(mapped {len(out)})", flush=True)
        time.sleep(0.35)
    json.dump(out, open(CACHE, "w"))
    print(f"  mapped {len(out)}/{len(ids)} -> {CACHE.name}")
    return out


def load_counts():
    with gzip.open(RAW / "GSE274975_raw_counts.tsv.gz", "rt") as f:
        header = f.readline().rstrip("\n").split("\t")[1:]
        cols = [c.replace("ReadsPerGene", "") for c in header]
        genes, mat = [], []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[0])
            mat.append(parts[1:])
    X = np.asarray(mat, dtype=np.float64)
    assert not np.isnan(X).any(), "NaN in raw counts"
    return genes, cols, X


def stratified_folds(y, k, seed):
    rng = np.random.RandomState(seed)
    fold_of = np.full(len(y), -1)
    for cls in np.unique(y):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        for pos, i in enumerate(idx):
            fold_of[i] = pos % k
    return [{"train": [int(i) for i in np.where(fold_of != f)[0]],
             "test": [int(i) for i in np.where(fold_of == f)[0]]}
            for f in range(k)]


def main():
    PROC.mkdir(parents=True, exist_ok=True)

    # ---- 1. raw counts ----
    genes, cols, X = load_counts()
    print(f"raw: {X.shape[0]} genes x {X.shape[1]} samples")

    # ---- 2. patient-number identity map ----
    idm = json.load(open(RAW / "id_map.json"))
    num2col = {int(k): v["col"] for k, v in idm.items()}
    num2t1 = {int(k): v["t1_id"] for k, v in idm.items()}
    assert set(num2col.values()) == set(cols), "counts columns != id_map"
    col_order = [num2col[n] for n in sorted(num2col)]  # by patient number
    col_idx = [cols.index(c) for c in col_order]
    X = X[:, col_idx]
    patient_nums = sorted(num2col)
    print(f"identity map: {len(patient_nums)} patients, columns ordered")

    # ---- 3. clinical from Table 1 ----
    t1 = pd.read_excel(RAW / "paper_tables/Table1.xlsx",
                       sheet_name="Supplementary_table_1")
    t1 = t1.set_index("Sample_ID")
    recist = [t1.loc[num2t1[n], "RECIST Response"] for n in patient_nums]
    response = [1 if r in ("CR", "PR") else 0 for r in recist]
    pfs = [float(t1.loc[num2t1[n], "PFS time, months"]) for n in patient_nums]
    status01 = [int(t1.loc[num2t1[n], "Response status available"])
                for n in patient_nums]
    n_res = int(sum(response))
    print(f"RECIST: responders(CR/PR)={n_res} non={len(response) - n_res}; "
          f"PFS<6mo={sum(p < 6 for p in pfs)}")

    # ---- 4. Ensembl -> HGNC ----
    emap = ensembl_map(genes)
    symbols = [emap.get(g) for g in genes]
    n_unmapped = sum(s is None for s in symbols)
    print(f"HGNC mapped: {len(genes) - n_unmapped}/{len(genes)} "
          f"({n_unmapped} unmapped dropped)")

    df = pd.DataFrame(X, index=[s if s else f"__unmapped_{g}"
                                for s, g in zip(symbols, genes)],
                      columns=[f"LuC_pat_{n}" for n in patient_nums])
    df = df[~df.index.str.startswith("__unmapped_")]
    dup = df.index.duplicated(keep=False)
    n_dup_rows = int(dup.sum())
    df = df.groupby(level=0).sum()  # sum duplicate symbols (additive counts)
    print(f"after aggregation: {df.shape[0]} unique symbols "
          f"({n_dup_rows} rows involved in duplicate groups)")
    assert (df.values < 0).sum() == 0

    # ---- 5. log2(CPM+1) ----
    lib = df.values.sum(axis=0)
    print(f"library sizes: min={lib.min():.0f} median={np.median(lib):.0f} "
          f"max={lib.max():.0f}")
    Xp = np.log2(df.values / np.maximum(lib, 1) * 1e6 + 1.0).T  # (61, G)
    Xp = np.asarray(Xp, dtype=np.float32)

    # ---- 6. AnnData ----
    import anndata as ad
    obs = pd.DataFrame({
        "response": np.array(response, dtype=int),
        "pfs_months": np.array(pfs, dtype=float),
        "status01": np.array(status01, dtype=int),
        "recist": recist,
        "histotype": [str(t1.loc[num2t1[n], "Histotype"])
                      for n in patient_nums],
        "therapy": [str(t1.loc[num2t1[n],
                              "Targeted therapeutics (next line)"])
                    for n in patient_nums],
        "pdl1_status": [str(t1.loc[num2t1[n], "PD-L1 status"])
                        for n in patient_nums],
        "age": [float(t1.loc[num2t1[n], "Age at biosampling"])
                for n in patient_nums],
        "gender": [str(t1.loc[num2t1[n], "Gender"]) for n in patient_nums],
        "tumor_site": [str(t1.loc[num2t1[n], "Tumor_site (biosampling)"])
                       for n in patient_nums],
        "primary_or_met": [str(t1.loc[num2t1[n],
                                      "Primary_tumor._metastasis"])
                           for n in patient_nums],
        "patient_num": patient_nums,
        "geo_title": [f"OB_pat_LuC_{n}" for n in patient_nums],
        "t1_id": [num2t1[n] for n in patient_nums],
        "counts_column": [num2col[n] for n in patient_nums],
    }, index=[f"LuC_pat_{n}" for n in patient_nums])
    adata = ad.AnnData(X=Xp, obs=obs)
    adata.var_names = df.index.tolist()
    adata.var["ensembl_ids"] = [
        ";".join(sorted({g for g, s in zip(genes, symbols) if s == sym}))
        for sym in df.index]
    adata.uns["provenance"] = (
        "GSE274975 raw counts (GEO) + Luppov 2024 Front Immunol "
        "Suppl Table 1 clinical; identity via patient-number bijection "
        "(id_map.json); log2(CPM+1); ENSG->HGNC via Ensembl REST")

    # ---- 7. splits ----
    y = adata.obs["response"].values.astype(int)
    k = min(5, int((y == 1).sum()), int((y == 0).sum()))
    folds = stratified_folds(y, k, SEED)
    payload = {"cohort": "GSE274975", "k": k, "seed": SEED,
               "endpoint": f"RECIST (responders {n_res}/{len(y)})",
               "folds": folds}
    blob = json.dumps(payload, sort_keys=True).encode()
    payload["splits_sha1_16"] = hashlib.sha1(blob).hexdigest()[:16]

    adata.write_h5ad(PROC / "GSE274975_HGNC.h5ad")
    json.dump(payload, open(PROC / "GSE274975_splits.json", "w"), indent=1)
    qc = {"n_raw_genes": int(len(genes)), "n_unmapped": int(n_unmapped),
          "n_unique_symbols": int(df.shape[0]),
          "n_dup_rows": n_dup_rows,
          "n_samples": int(Xp.shape[0]),
          "n_responders": n_res, "n_pfs_lt6": int(sum(p < 6 for p in pfs)),
          "lib_min": float(lib.min()), "lib_max": float(lib.max()),
          "splits_sha1_16": payload["splits_sha1_16"]}
    json.dump(qc, open(PROC / "build_qc.json", "w"), indent=1)
    print("\nQC:", json.dumps(qc, indent=1))
    print(f"\nWROTE {PROC / 'GSE274975_HGNC.h5ad'}")
    print(f"WROTE {PROC / 'GSE274975_splits.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
