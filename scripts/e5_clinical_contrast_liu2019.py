"""E5 (new experiment, 2026-09-17 quality campaign): the first NON-circular
response-definition contrast — two clinical endpoints on the same patients.

Motivation. The manuscript states (DC3 relation paragraph) that "no
non-circular response-definition contrast (e.g., RECIST vs PFS in the same
patients) was available in any public cohort — a boundary of the claim".
The cBioPortal iAtlas harmonization of Liu 2019 carries both a RECIST
response annotation and PFS_MONTHS/PFS_STATUS, so that boundary can be
closed with data already in the archive.

Design (mirrors e1_replication_imvigor210.py exactly; seed 42):
  endpoint A (clinical): RECIST, CR/PR = responder
  endpoint B (clinical): DCB-style benefit from PFS —
      DCB  = PFS_MONTHS >= 6 (progression-free for at least 6 months,
             whether by event or censoring)
      NDB  = PFS_MONTHS < 6 AND PFS_STATUS == '1:PROGRESSED'
      excluded = PFS_MONTHS < 6 AND censored (event-free status unknown)
  Both endpoints on the identical sample set (inner join, both annotated).
  leakage-free nested CV, frozen hyperparameters from tune_primary,
  6 methods; exact 50,000-shuffle permutations for the four fixed scorers;
  CDS v1.1.0 on both endpoints + 10-replicate random-label null
  (declared genes GZMA/PRF1, the clinical-endpoint convention of Table 3).

Mechanistic prediction (the manuscript's own account): neither endpoint is
gene-defined, so no circularity channel exists — (i) no method should
approach the 0.95+ inflated regime, and (ii) the endpoint-switch delta
should shrink toward the between-method gap scale rather than the 0.4+
cytolytic collapse.

Output: results/benchmark/v33/e5_liu2019_clinical_contrast.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))
sys.path.insert(0, str(REPO / "cds_tool"))

from rerun_v33 import SEED, nested_cv_predict, tune_primary, cv_predict  # noqa: E402
from metrics import evaluate_all  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

EXT = REPO / "data" / "external" / "mel_iatlas_liu_2019"
OUT = REPO / "results" / "benchmark" / "v33" / \
    "e5_liu2019_clinical_contrast.json"
RESPONDERS = ["COMPLETE RESPONSE", "PARTIAL RESPONSE"]


def auroc_fixed(score, y):
    r = rankdata(score)
    n1 = int(y.sum())
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def stratified_folds(y, k, seed):
    rng = np.random.RandomState(seed)
    fold_of = np.full(len(y), -1)
    for cls in np.unique(y):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        for pos, i in enumerate(idx):
            fold_of[i] = pos % k
    return [{"train": list(np.where(fold_of != f)[0]),
             "test": list(np.where(fold_of == f)[0])} for f in range(k)]


def perm_fixed(score, y, n=50000, seed=SEED):
    rng = np.random.RandomState(seed)
    obs = auroc_fixed(score, y)
    ge = 0
    for _ in range(n):
        if auroc_fixed(score, rng.permutation(y)) >= obs:
            ge += 1
    return obs, (ge + 1) / (n + 1)


def build_dcb(clin):
    """DCB/NDB/excluded construction from PFS_MONTHS + PFS_STATUS."""
    pfs = pd.to_numeric(clin["PFS_MONTHS"], errors="coerce")
    status = clin["PFS_STATUS"].astype(str).str.strip()
    dcb = pd.Series(np.nan, index=clin.index)
    dcb[pfs >= 6] = 1                                   # >=6 mo free of progression
    dcb[(pfs < 6) & (status == "1:PROGRESSED")] = 0     # progressed before 6 mo
    # pfs < 6 & censored -> excluded (NaN retained)
    return dcb, pfs, status


def main():
    import scanpy as sc
    mat = pd.read_csv(EXT / "expression.tsv.gz", sep="\t", index_col=0)
    mat = mat.dropna(axis=0, how="any")  # drop genes not measured in every sample
    print(f"expression after NaN gene drop: {mat.shape}", flush=True)
    clin = pd.read_csv(EXT / "clinical.tsv", sep="\t", index_col=0)
    clin.index = clin.index.astype(str)

    # Clinical structure (iAtlas export): 244 sample rows for 122 patients —
    # tumour-biopsy rows (SAMPLE_ID = Liu_SampleN, carrying RESPONSE and
    # matching the expression columns) plus one row per patient keyed
    # 'PatientN' carrying the patient-level PFS annotation. Aggregate to
    # patient level (first non-null per column), then attach PFS to the
    # tumour samples via PATIENT_ID.
    resp = clin.get("RESPONSE")
    pat_level = clin.groupby("PATIENT_ID").first(numeric_only=False)
    samp2pat = clin["PATIENT_ID"]

    recist = resp.astype(str).str.strip().str.upper() \
        .isin(RESPONDERS).astype(float)
    recist[resp.isna()] = np.nan

    pfs = pd.to_numeric(pat_level.get("PFS_MONTHS"), errors="coerce")
    pfs_status = pat_level.get("PFS_STATUS").astype(str).str.strip()
    dcb_pat = pd.Series(np.nan, index=pat_level.index)
    dcb_pat[pfs >= 6] = 1                                   # >=6 mo free of progression
    dcb_pat[(pfs < 6) & (pfs_status == "1:PROGRESSED")] = 0  # progressed < 6 mo
    # pfs < 6 & censored -> excluded (NaN retained)

    # expression columns -> patients (1:1, verified)
    col_pat = samp2pat.reindex(mat.columns)          # column -> PATIENT_ID
    recist_col = recist.reindex(mat.columns)         # sample-indexed, direct
    dcb_col = pd.Series(
        dcb_pat.reindex(col_pat.values).values, index=mat.columns)
    keep = recist_col.notna() & dcb_col.notna() & col_pat.notna()
    cols = list(mat.columns[keep.values])
    print(f"Liu 2019: clinical rows={len(clin)}, RECIST-annotated="
          f"{int(recist.notna().sum())}, PFS-annotated patients="
          f"{int(pfs.notna().sum())}, usable tumour samples={len(cols)}",
          flush=True)

    m = mat[cols].T
    m.insert(0, "recist", recist_col.reindex(cols).values)
    m.insert(1, "dcb", dcb_col.reindex(cols).values)
    y_recist = m["recist"].astype(int).values
    y_dcb = m["dcb"].astype(int).values
    genes = list(mat.index)
    X = np.asarray(m.drop(columns=["recist", "dcb"]).values, dtype=float)

    # construction + concordance table
    tab = pd.crosstab(pd.Series(y_recist, name="RECIST"),
                      pd.Series(y_dcb, name="DCB"))
    n_tot = len(y_recist)
    expected = float((tab.sum(axis=1) * tab.sum(axis=0)).sum()) / n_tot
    denom = n_tot - expected
    kappa = (float(tab.values.diagonal().sum()) - expected) / denom \
        if denom > 0 else float("nan")
    print(f"endpoint construction: RECIST responders={int(y_recist.sum())}/"
          f"{len(y_recist)}, DCB={int(y_dcb.sum())}/{len(y_dcb)}", flush=True)
    print(f"cross-tab (rows=RECIST 0/1, cols=DCB 0/1):\n{tab}", flush=True)
    print(f"Cohen's kappa(RECIST, DCB) = {kappa:.3f}", flush=True)

    adata = sc.AnnData(X=X)
    adata.var_names = genes

    out = {"meta": {
        "cohort": "Liu 2019 (mel_iatlas_liu_2019, cBioPortal iAtlas)",
        "n": int(len(y_recist)),
        "endpoints": {
            "RECIST": {"definition": "CR/PR vs SD/PD",
                       "n_responders": int(y_recist.sum())},
            "DCB": {"definition": "PFS_MONTHS >= 6 vs progression < 6 mo "
                                  "(censored < 6 mo excluded)",
                    "n_benefit": int(y_dcb.sum())}},
        "concordance": {"kappa": round(kappa, 3),
                        "crosstab": {str(k): {str(c): int(v)
                                              for c, v in row.items()}
                                     for k, row in tab.iterrows()}},
        "seed": SEED,
        "expression": "cBioPortal z-scores (affine-invariant for the "
                      "benchmark metrics)",
        "note": "first non-circular (clinical-clinical) response-definition "
                "contrast; closes the boundary acknowledged in the DC3 "
                "relation paragraph",
    }}
    bench_methods = ["GEP", "TIDE", "IMPRES", "PD_L1_IHC", "ElasticNet",
                     "ElasticNet_Var"]
    label = {"IMPRES": "IMPRES", "GEP": "GEP", "TIDE": "TIDE",
             "PD_L1_IHC": "PD-L1 (CD274)", "ElasticNet": "ElasticNet (MI)",
             "ElasticNet_Var": "ElasticNet (Var)"}

    for ename, y in (("RECIST", y_recist), ("DCB", y_dcb)):
        adata.obs["response"] = y
        k = min(5, int((y == 1).sum()), int((y == 0).sum()))
        folds = stratified_folds(y, k, SEED)
        res = {"k": k}
        for meth in bench_methods:
            res[meth] = {}
            if meth in ("ElasticNet", "ElasticNet_Var"):
                params = tune_primary(adata, folds, meth)
                yt, yp, _ = nested_cv_predict(adata, folds, meth,
                                              params=params)
                res[meth]["params"] = params
            else:
                yt, yp = cv_predict(adata, folds, meth)
            r = evaluate_all(np.asarray(yt).astype(int),
                             np.asarray(yp, dtype=float))
            res[meth].update({
                "auroc": round(r["auroc"], 4),
                "auroc_ci": [round(r["auroc_ci_low"], 4),
                             round(r["auroc_ci_high"], 4)],
                "auprc": round(r["auprc"], 4)})
            print(f"  {ename} {label[meth]}: {r['auroc']:.3f} "
                  f"[{r['auroc_ci_low']:.3f}-{r['auroc_ci_high']:.3f}]",
                  flush=True)
        for meth in ("IMPRES", "GEP", "TIDE", "PD_L1_IHC"):
            yt_f, yp = cv_predict(adata, folds, meth)
            yt_f = np.asarray(yt_f).astype(int)
            obs, p = perm_fixed(np.asarray(yp, dtype=float), yt_f)
            res[meth]["perm_p"] = p
            res[meth]["n_shuffle"] = 50000
            print(f"  {ename} {label[meth]} perm p={p:.5f}", flush=True)
        out[ename] = res

    out["endpoint_switch_delta"] = {
        meth: round(out["RECIST"][meth]["auroc"] - out["DCB"][meth]["auroc"], 3)
        for meth in bench_methods}

    from cds import circularity_detection_score
    gi = {g: i for i, g in enumerate(genes)}
    cyt_genes = [g for g in ("GZMA", "PRF1") if g in gi]
    rng = np.random.RandomState(SEED)
    out.setdefault("cds", {})
    for ename, y in (("RECIST", y_recist), ("DCB", y_dcb)):
        r = circularity_detection_score(X, y, cyt_genes, genes)
        r = json.loads(json.dumps(r, default=float))
        nulls = [float(circularity_detection_score(
            X, rng.permutation(y), cyt_genes, genes)["CDS"])
            for _ in range(10)]
        out["cds"][ename] = {"CDS": r["CDS"], "risk": r["risk_level"],
                             "components": r["components"],
                             "null_mean": round(float(np.mean(nulls)), 3),
                             "null_max": round(float(np.max(nulls)), 3)}
        print(f"  CDS {ename}: {r['CDS']:.3f} "
              f"(null mean {np.mean(nulls):.3f})", flush=True)

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print("WROTE", OUT, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
