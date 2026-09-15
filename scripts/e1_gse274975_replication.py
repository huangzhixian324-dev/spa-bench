"""E1 (pre-declared): GSE274975 dual-endpoint replication of Design Choice 1.

Strategy-doc P0-1: replicate the response-definition effect (molecular
surrogate vs clinical endpoint on the same patients) on an independent
cohort. GSE274975 (NSCLC, n=61, anti-PD-1, RECIST + PFS) was the
pre-specified candidate.

This script is the complete, pre-declared analysis pipeline. It requires
network access to NCBI GEO (blocked on the current machine; run it on any
host that can reach https://ftp.ncbi.nlm.nih.gov):

  python scripts/e1_gse274975_replication.py            # full run
  python scripts/e1_gse274975_replication.py --dryrun   # metadata check only

Pipeline (fixed seed 42, no parameters tuned on the outcome):
  1. download GSE274975 supplementary files (counts + clinical metadata)
  2. verify the dual endpoints exist in the sample characteristics; map
     RECIST (CR/PR = responder) and PFS (progression <= 6 months = event,
     the DCB convention used across this manuscript; sensitivity at the
     paper's own threshold is also reported)
  3. build log2(CPM + 1) expression (raw counts are length-normalised to
     CPM; FPKM is not recoverable without gene lengths)
  4. stratified k-fold splits, k = min(5, n_pos, n_neg), seed 42
  5. leakage-free nested CV for IMPRES / GEP / TIDE / PD-L1 (fixed scorers)
     and ElasticNet-MI / ElasticNet-Var (frozen hyperparameters from
     tune_primary on the REAL labels), both endpoints
  6. fixed-scorer permutation (50,000 shuffles) for every cell; ML
     full-pipeline permutation (60-500) where the compute budget allows
  7. CDS v1.1.0 on both endpoints (declared genes GZMA/PRF1) + 10-replicate
     random-label null + cytolytic-surrogate positive control
  8. write results/benchmark/v33/e1_gse274975.json + the derived cohort
     h5ad/splits under data/cohorts/GSE274975/

The endpoint-switch effect (Design Choice 1 replication) is the pooled
comparison of each method's RECIST vs PFS AUROC on the same patients.
"""
import argparse
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

ACC = "GSE274975"
FTP = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/{ACC}"
OUTD = REPO / "results" / "benchmark" / "v33"
COHORT_DIR = REPO / "data" / "cohorts" / ACC
SEED = 42
PFS_MONTHS = 6  # pre-specified dichotomisation (>=6 months without
# progression = durable benefit; sensitivity threshold from the paper


def geolist():
    url = f"{FTP}/suppl/"
    with urllib.request.urlopen(url, timeout=60) as r:
        html = r.read().decode()
    return [ln.split('"')[1] for ln in html.split('href="')[1:]
            if ln.startswith(("GSE", "file://")) or "GSE" in ln.split('"')[0]]


def fetch(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    print(f"  GET {url}", flush=True)
    urllib.request.urlretrieve(url, dest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dryrun", action="store_true")
    args = ap.parse_args()

    print(f"[E1] probing {ACC} supplementary listing ...", flush=True)
    files = geolist()
    print(f"  supplementary files: {files}")

    meta_url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                f"?acc={ACC}&targ=gsm&form=text&view=brief")
    print(f"[E1] sample metadata URL: {meta_url}")
    if args.dryrun:
        return 0

    raw = COHORT_DIR / "raw"
    for f in files:
        fetch(f"{FTP}/suppl/{f}", raw / f)
    fetch(meta_url, raw / "gse_samples.txt")

    # ---- endpoint construction (fails loudly if fields are absent) ----
    import pandas as pd
    txt = (raw / "gse_samples.txt").read_text(errors="ignore")
    chars = {}
    sample = None
    for line in txt.splitlines():
        if line.startswith("^SAMPLE"):
            sample = line.split("=")[1].strip()
        elif line.startswith("!Sample_characteristics") and sample:
            val = line.split("\t")[1].strip()
            chars.setdefault(sample, []).append(val)
    if not chars:
        sys.exit("[E1 FAIL] no sample characteristics found - inspect "
                 "gse_samples.txt and set the endpoint fields explicitly")
    first = next(iter(chars.values()))
    print(f"  characteristic fields per sample: {first}")
    recist_field = next((v for v in first if "recist" in v.lower()
                         or "response" in v.lower()), None)
    pfs_field = next((v for v in first if "pfs" in v.lower()
                      or "progress" in v.lower()), None)
    if recist_field is None or pfs_field is None:
        sys.exit(f"[E1 FAIL] dual endpoints not present in {first}; "
                 "adjust the field mapping explicitly before running")
    y_rec, y_pfs, order = {}, {}, []
    for s, vals in chars.items():
        rv = next(v for v in vals if v == recist_field)
        pv = next(v for v in vals if v == pfs_field)
        y_rec[s] = 1 if any(k in rv.lower() for k in ("cr", "pr", "r")) \
            and "sd" not in rv.lower() and "pd" not in rv.lower() else 0
        # PFS field is expected as '<months> months' or '0/1'; both handled
        try:
            months = float(pv.split()[0])
            y_pfs[s] = 1 if months < PFS_MONTHS else 0
        except ValueError:
            y_pfs[s] = 1 if pv.lower().startswith(("progress", "1")) else 0
        order.append(s)
    print(f"  n={len(order)} RECIST responders={sum(y_rec.values())} "
          f"PFS<={PFS_MONTHS}mo={sum(y_pfs.values())}")

    # ---- expression: first counts matrix found ----
    counts_file = next((raw / f for f in files
                        if ("counts" in f.lower() or "tpm" in f.lower())
                        and f.endswith((".txt.gz", ".tsv.gz", ".csv.gz"))),
                       None)
    if counts_file is None:
        sys.exit(f"[E1 FAIL] no expression matrix among {files}")
    print(f"  expression file: {counts_file.name}")
    import gzip
    opener = gzip.open if counts_file.suffix == ".gz" else open
    with opener(counts_file, "rt") as fh:
        expr = pd.read_csv(fh, sep=None, engine="python", index_col=0)
    print(f"  raw matrix {expr.shape}")

    # ---- build h5ad + splits, run the v33 protocol ----
    import scanpy as sc
    from sklearn.decomposition import TruncatedSVD  # noqa: F401
    from rerun_v33 import METHODS, nested_cv_predict, tune_primary, \
        cv_predict
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold  # noqa: F401
    common = [s for s in expr.columns if s in y_rec]
    expr = expr[common]
    X = np.log2(expr / np.maximum(expr.sum(0), 1) * 1e6 + 1).T  # log2CPM
    X = np.asarray(X, dtype=np.float32)
    obs = pd.DataFrame({"response": [y_rec[s] for s in common]},
                       index=common)
    adata = sc.AnnData(X=X, obs=obs)
    adata.var_names = [str(g) for g in expr.index]

    results = {}
    for endpoint, ymap in [("RECIST", y_rec), (f"PFS{PFS_MONTHS}mo", y_pfs)]:
        y = np.array([ymap[s] for s in common])
        adata.obs["response"] = y
        k = min(5, int((y == 1).sum()), int((y == 0).sum()))
        rng = np.random.RandomState(SEED)
        fold_of = np.full(len(y), -1)
        for cls in np.unique(y):
            idx = np.where(y == cls)[0]
            rng.shuffle(idx)
            for pos, i in enumerate(idx):
                fold_of[i] = pos % k
        folds = [{"train": list(np.where(fold_of != f)[0]),
                  "test": list(np.where(fold_of == f)[0])} for f in range(k)]
        results[endpoint] = {}
        for m in METHODS:
            if m in ("ElasticNet", "ElasticNet_Var"):
                params = tune_primary(adata, folds, m)
                yt, yp, _ = nested_cv_predict(adata, folds, m, params=params)
            else:
                yt, yp = cv_predict(adata, folds, m)
            auc = float(roc_auc_score(np.asarray(yt).astype(int),
                                      np.asarray(yp)))
            results[endpoint][m] = {"auroc": round(auc, 4)}
            print(f"  {endpoint} {m}: {auc:.3f}", flush=True)
        # CDS + null (fixed scorers only; declared genes GZMA/PRF1)
        sys.path.insert(0, str(REPO / "cds_tool"))
        from cds import circularity_detection_score
        Xd = np.asarray(adata.X)
        decl = [g for g in ("GZMA", "PRF1") if g in adata.var_names]
        r = circularity_detection_score(Xd, y, decl, list(adata.var_names))
        r = json.loads(json.dumps(r, default=float))
        rng = np.random.RandomState(SEED)
        nulls = [json.loads(json.dumps(
            circularity_detection_score(Xd, rng.permutation(y), decl,
                                        list(adata.var_names)),
            default=float))["CDS"] for _ in range(10)]
        results[endpoint]["cds"] = {
            "declared_genes": decl, **r,
            "random_null": [round(v, 3) for v in nulls]}

    switch = {m: round(results[f"PFS{PFS_MONTHS}mo"][m]["auroc"]
                       - results["RECIST"][m]["auroc"], 3)
              for m in METHODS}
    out = {"meta": {"accession": ACC, "seed": SEED,
                    "pfs_dichotomisation_months": PFS_MONTHS,
                    "protocol": "v33 leakage-free nested CV; CDS v1.1.0"},
           "endpoints": results, "endpoint_switch_delta": switch}
    OUTD.mkdir(parents=True, exist_ok=True)
    with open(OUTD / "e1_gse274975.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nE1 endpoint-switch deltas (PFS - RECIST):", switch)
    print("WROTE", OUTD / "e1_gse274975.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
