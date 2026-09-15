"""E3: literature endpoint CDS backtest + Hugo OS truth, via cBioPortal.

GEO/NCBI is unreachable from this machine, but cBioPortal (iAtlas
harmonizations) hosts the endpoint cohorts the manuscript's Table S17
calls 'not publicly accessible':

  blca_iatlas_imvigor210_2017   IMvigor210 (TIDE / GEP published cohorts)
  mel_iatlas_liu_2019           Liu 2019 DFCI anti-PD-1 melanoma (RECIST)
  skcm_vanderbilt_mskcc_2015    Van Allen 2015 samples (no mRNA profile in
                                the harmonization -> recorded as
                                expression-inaccessible)
  mel_iatlas_hugo_ucla_2016     Hugo 2016 (OS_STATUS truth backfill;
                                samples matched to the local h5ad rows by
                                expression-profile correlation because the
                                local matrix carries no barcodes)

For every cohort with a response endpoint this script computes CDS v1.1.0
(declared genes GZMA/PRF1 - the same convention as the internal Table S4b)
for (a) the clinical response endpoint, (b) the cytolytic surrogate median
split as a positive control, and (c) a 10-replicate random-label null.

Outputs:
  data/external/<studyId>/{clinical.tsv, expression.tsv.gz}
  results/benchmark/v33/e3_literature_cds_backtest.json
  results/benchmark/v33/hugo_os_truth.json
"""
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "cds_tool"))
from cds import circularity_detection_score  # noqa: E402

BASE = "https://www.cbioportal.org/api"
OUTD = REPO / "data" / "external"
RES = REPO / "results" / "benchmark" / "v33"

COHORTS = {
    "blca_iatlas_imvigor210_2017": {"endpoint": "RESPONSE",
                                    "responders": ["COMPLETE RESPONSE",
                                                   "PARTIAL RESPONSE"]},
    "mel_iatlas_liu_2019": {"endpoint": "RESPONSE",
                            "responders": ["COMPLETE RESPONSE",
                                           "PARTIAL RESPONSE"]},
    "skcm_vanderbilt_mskcc_2015": {"endpoint": None},
    "mel_iatlas_hugo_ucla_2016": {"endpoint": "RESPONSE",
                                  "responders": ["R", "CR", "PR"],
                                  "max_genes": 5000},
}


def sess():
    import requests
    from requests.adapters import HTTPAdapter, Retry
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=6, backoff_factor=3,
        status_forcelist=[429, 500, 502, 503, 504])))
    return s


S = sess()


def get_with_backoff(url, params=None, timeout=180):
    for attempt in range(6):
        try:
            r = S.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r
        except Exception as e:
            if attempt == 5:
                raise
            wait = 30 * (attempt + 1)
            print(f"    [retry in {wait}s] {e}", flush=True)
            time.sleep(wait)


def fetch_clinical(study):
    """Merge patient- and sample-level clinical attributes; index by
    SAMPLE ID (patient attributes replicated on each sample)."""
    merged = {}
    for dtype in ("PATIENT", "SAMPLE"):
        page = 0
        while True:
            r = get_with_backoff(
                f"{BASE}/studies/{study}/clinical-data",
                params={"clinicalDataType": dtype, "pageSize": 5000,
                        "pageNumber": page})
            recs = r.json()
            if not recs:
                break
            for rec in recs:
                key = rec.get("sampleId") or rec["patientId"]
                merged.setdefault(key, {})[rec["clinicalAttributeId"]] = \
                    rec["value"]
                merged[key]["PATIENT_ID"] = rec["patientId"]
            if len(recs) < 5000:
                break
            page += 1
            time.sleep(2)
    df = pd.DataFrame.from_dict(merged, orient="index")
    df.index = df.index.astype(str)
    df.index.name = "SAMPLE_ID"
    return df


def fetch_expression(study, max_genes=None):
    cached = OUTD / study / "expression.tsv.gz"
    if cached.exists() and cached.stat().st_size > 10000:
        print("  using cached expression.tsv.gz", flush=True)
        return pd.read_csv(cached, index_col=0), "cached"
    profs = get_with_backoff(
        f"{BASE}/studies/{study}/molecular-profiles",
        params={"pageSize": 100}).json()
    prof = next((p["molecularProfileId"] for p in profs
                 if p.get("molecularAlterationType") == "MRNA_EXPRESSION"
                 and p.get("datatype") == "CONTINUOUS"
                 and "rna" in p["molecularProfileId"]), None)
    if prof is None:
        prof = next((p["molecularProfileId"] for p in profs
                     if p.get("molecularAlterationType") == "MRNA_EXPRESSION"),
                    None)
    if prof is None:
        return None, None
    print(f"  profile: {prof}", flush=True)
    genes = get_with_backoff(f"{BASE}/genes",
                             params={"pageSize": 30000}).json()
    entrez = sorted({g["entrezGeneId"] for g in genes})
    symbols = {g["entrezGeneId"]: g["hugoGeneSymbol"] for g in genes}
    if max_genes and len(entrez) > max_genes:
        step = max(1, len(entrez) // max_genes)
        entrez = entrez[::step][:max_genes]
        print(f"  gene universe capped to {len(entrez)} (deterministic "
              "spread; CDS C1/C2 use only the declared genes, C3 uses "
              "top-500 variance vs 500 random)", flush=True)
    print(f"  gene universe: {len(entrez)}", flush=True)
    slist = None
    for sl in get_with_backoff(f"{BASE}/sample-lists",
                               params={"studyId": study,
                                       "pageSize": 1000}).json():
        if sl["sampleListId"] == f"{study}_all":
            slist = sl["sampleListId"]
            break
    sample_ids = None
    if not slist:
        sample_ids = [s["sampleId"] for s in get_with_backoff(
            f"{BASE}/studies/{study}/samples",
            params={"pageSize": 100000}).json()]

    import concurrent.futures as cf

    def fetch_chunk(i):
        chunk = entrez[i:i + 2000]
        body = ({"sampleListId": slist, "entrezGeneIds": chunk} if slist
                else {"sampleIds": sample_ids, "entrezGeneIds": chunk})
        for attempt in range(3):
            try:
                r = S.post(
                    f"{BASE}/molecular-profiles/{prof}/molecular-data/fetch",
                    json=body, timeout=600)
                r.raise_for_status()
                recs = r.json()
                if not recs:
                    return None
                df = pd.DataFrame(recs)
                piv = df.pivot_table(index="entrezGeneId",
                                     columns="sampleId", values="value",
                                     aggfunc="first")
                piv.index = [symbols.get(g, str(g)) for g in piv.index]
                return piv
            except Exception as e:
                if attempt == 2:
                    print(f"    [chunk {i} failed: {e}]", flush=True)
                    return None
                time.sleep(30 * (attempt + 1))

    frames = []
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        for j, piv in enumerate(ex.map(fetch_chunk,
                                       range(0, len(entrez), 2000))):
            if piv is not None:
                frames.append(piv)
            if (j + 1) % 5 == 0:
                print(f"    genes {min((j + 1) * 2000, len(entrez))}"
                      f"/{len(entrez)}", flush=True)
    if not frames:
        return None, None
    mat = pd.concat(frames).groupby(level=0).mean()
    return mat, prof


def cds_block(X, y, gene_list, declared=("GZMA", "PRF1")):
    Xd = np.asarray(X, dtype=float)
    gl = list(gene_list)

    def run(yy, genes):
        r = circularity_detection_score(Xd, yy, list(genes), gl)
        r = json.loads(json.dumps(r, default=float))
        return {"CDS": r["CDS"],
                "C1": r["components"]["C1_MI_endpoint_genes"],
                "C2": r["components"]["C2_rho_endpoint_response"],
                "C3": r["components"]["C3_MI_skew_ratio"],
                "risk": r["risk_level"].split("(")[0].strip()}

    missing = [g for g in declared if g not in gl]
    out = {"declared_genes_missing": missing}
    out["clinical_endpoint"] = run(y, declared)
    if not missing:
        cyt = Xd[:, [gl.index(g) for g in declared]].mean(axis=1)
        out["cytolytic_surrogate_positive_control"] = run(
            (cyt > np.median(cyt)).astype(int), declared)
    rng = np.random.RandomState(42)
    nulls = [run(rng.permutation(y), declared)["CDS"] for _ in range(10)]
    out["random_null"] = {"n": 10, "mean": round(float(np.mean(nulls)), 3),
                          "max": round(float(np.max(nulls)), 3),
                          "values": [round(v, 3) for v in nulls]}
    return out


def main():
    OUTD.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    backtest = {}

    for study, cfg in COHORTS.items():
        print(f"===== {study} =====", flush=True)
        time.sleep(10)
        d = OUTD / study
        d.mkdir(parents=True, exist_ok=True)
        clin = fetch_clinical(study)
        clin.to_csv(d / "clinical.tsv", sep="\t")
        mat, prof = fetch_expression(study, cfg.get("max_genes"))
        if mat is None:
            print("  [SKIP] no mRNA profile - recorded as "
                  "expression-inaccessible", flush=True)
            backtest[study] = {"expression_accessible": False,
                               "reason": "no mRNA molecular profile in "
                                         "the iAtlas harmonization"}
            continue
        backtest[study] = {"expression_accessible": True}
        mat.to_csv(d / "expression.tsv.gz", compression="gzip")
        if not cfg.get("endpoint"):
            continue
        resp = clin.get(cfg["endpoint"])
        if resp is None:
            print(f"  [SKIP] no {cfg['endpoint']} attribute", flush=True)
            continue
        m = mat.T.join(resp.rename("response"), how="inner")
        m = m[m["response"].notna()]
        y = m["response"].astype(str).str.strip().str.upper() \
            .isin(cfg["responders"]).astype(int).values
        if len(np.unique(y)) < 2:
            print(f"  [SKIP] endpoint degenerate: "
                  f"{m['response'].value_counts().to_dict()}", flush=True)
            continue
        print(f"  n={len(y)} responders={int(y.sum())}", flush=True)
        blk = cds_block(m.drop(columns=["response"]).values, y,
                        list(mat.index))
        blk["n"] = int(len(y))
        blk["n_responders"] = int(y.sum())
        blk["endpoint_attr"] = cfg["endpoint"]
        blk["endpoint_values"] = m["response"].value_counts().to_dict()
        blk["profile"] = prof
        blk["note"] = ("expression values are cBioPortal z-scores/published "
                       "values; CDS statistics (rank MI, Spearman) are "
                       "invariant to per-gene affine transforms")
        backtest[study] = blk
        cc = blk.get("cytolytic_surrogate_positive_control", {})
        print(f"  CDS clinical={blk['clinical_endpoint']['CDS']:.3f} "
              f"cytolytic_control={cc.get('CDS', float('nan')):.3f} "
              f"null_mean={blk['random_null']['mean']:.3f}", flush=True)

    def np_default(o):
        if hasattr(o, "item"):
            return o.item()
        return str(o)

    with open(RES / "e3_literature_cds_backtest.json", "w") as fh:
        json.dump(backtest, fh, indent=1, default=np_default)

    # ---------------- Hugo OS truth ----------------
    hugo = OUTD / "mel_iatlas_hugo_ucla_2016"
    if (hugo / "expression.tsv.gz").exists():
        mat = pd.read_csv(hugo / "expression.tsv.gz", index_col=0)
        clin = pd.read_csv(hugo / "clinical.tsv", sep="\t", index_col=0)
        clin.index = clin.index.astype(str)
        import scanpy as sc
        a = sc.read_h5ad(REPO / "data" / "cohorts" / "Hugo_2016" /
                         "processed" / "Hugo_2016_processed.h5ad")
        Xloc = (a.X.toarray() if hasattr(a.X, "toarray")
                else np.asarray(a.X))
        genes = list(a.var_names)
        common = [g for g in mat.index if g in genes]
        li = {g: i for i, g in enumerate(genes)}
        A = Xloc[:, [li[g] for g in common]]
        B = mat.loc[common].T.values
        Az = (A - A.mean(1, keepdims=True)) / A.std(1, keepdims=True)
        Bz = (B - B.mean(1, keepdims=True)) / B.std(1, keepdims=True)
        corr = Az @ Bz.T / len(common)
        n_loc, n_ext = corr.shape
        pairs = sorted(((corr[i, j], i, j) for i in range(n_loc)
                        for j in range(n_ext)), reverse=True)
        used_i, used_j, assign = set(), set(), {}
        for c, i, j in pairs:
            if i in used_i or j in used_j:
                continue
            used_i.add(i)
            used_j.add(j)
            assign[i] = (j, c)
            if len(assign) == min(n_loc, n_ext):
                break
        rows = []
        for i in range(n_loc):
            if i not in assign:
                rows.append({"local_row": i, "iatlas_sample": None,
                             "corr": None, "margin": None,
                             "os_status": None, "os_months": None,
                             "response_iatlas": None,
                             "local_os_days": float(a.obs["OS_days"]
                                                    .iloc[i]),
                             "local_response":
                                 int(a.obs["response"].iloc[i]),
                             "note": "no iAtlas mRNA sample available"})
                continue
            j, c = assign[i]
            pid = mat.columns[j]
            second = max(corr[i, jj] for jj in range(n_ext) if jj != j)
            row = clin.loc[pid] if pid in clin.index else {}
            osm = row.get("OS_MONTHS")
            rows.append({
                "local_row": i, "iatlas_sample": pid,
                "corr": round(float(c), 3),
                "margin": round(float(c - second), 3),
                "os_status": row.get("OS_STATUS"),
                "os_months": (None if osm in (None, "", "NA", "nan")
                              else float(osm)),
                "response_iatlas": row.get("RESPONDER",
                                           row.get("RESPONSE")),
                "local_os_days": float(a.obs["OS_days"].iloc[i]),
                "local_response": int(a.obs["response"].iloc[i])})
        out = {"n_matched": sum(1 for r in rows if r["iatlas_sample"]),
               "min_margin": min(r["margin"] for r in rows
                                 if r["margin"] is not None),
               "note": "greedy one-to-one matching by per-sample "
                       "expression-profile correlation (z-scored common "
                       "genes); the local h5ad carries no barcodes",
               "matches": rows}
        with open(RES / "hugo_os_truth.json", "w") as fh:
            json.dump(out, fh, indent=1, default=np_default)
        print(f"Hugo OS truth: {out['n_matched']} matched, min margin "
              f"{out['min_margin']:.3f}", flush=True)
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
