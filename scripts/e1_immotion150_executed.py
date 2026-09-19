"""E1 (executed): IMmotion150 (ccRCC, iAtlas) dual-endpoint replication.

Third executed E1 replication of Design Choice 1, extending the
endpoint-sensitivity result to renal cell carcinoma (4th cancer type in
the manuscript's evidence base after melanoma, NSCLC, UC).  Protocol
identical to e1_gse274975_executed.py / e1_replication_liu2019.py:
fixed seed 42, leakage-free nested CV with k = min(5, n_pos, n_neg),
frozen hyperparameters from tune_primary on real labels, 50,000-shuffle
exact permutation for the four fixed scorers, CDS v1.1.0 + 10-replicate
random-label null (declared genes GZMA/PRF1).

Endpoints (atezolizumab arm only, n = 174; 9 lack RECIST):
  RECIST   CR/PR vs SD/PD (RESPONDER field; trial's own coding) - n=165
  DCB6mo   PFS >= 6 months AND no progression event = 1 (the manuscript's
           Jung-2019 durable-benefit convention; this cohort HAS the
           event indicator) - primary PFS endpoint
  PFS6mo   PFS < 6 months, time-only (GSE274975-comparable sensitivity)
  cytolytic_surrogate  GZMA/PRF1 mean, median split (exploratory)

Caveat: iAtlas pools atezo monotherapy with atezo+bevacizumab arms.

Output: results/benchmark/v33/e1_immotion150.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))
sys.path.insert(0, str(REPO / "cds_tool"))

from rerun_v33 import (SEED, cv_predict,  # noqa: E402
                       nested_cv_predict, tune_primary)
from metrics import evaluate_all  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

COHORT = REPO / "data/cohorts/IMmotion150/processed/IMmotion150_HGNC.h5ad"
SPLITS = REPO / "data/cohorts/IMmotion150/processed/IMmotion150_splits.json"
OUT = REPO / "results/benchmark/v33/e1_immotion150.json"
DECLARED = ("GZMA", "PRF1")


def auroc_fixed(score, y):
    r = rankdata(score)
    n1 = int(y.sum())
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def perm_fixed(score, y, n=50000, seed=SEED):
    rng = np.random.RandomState(seed)
    obs = auroc_fixed(score, y)
    ge = 0
    for _ in range(n):
        if auroc_fixed(score, rng.permutation(y)) >= obs:
            ge += 1
    return obs, (ge + 1) / (n + 1)


def folds_from_subset(y, k, seed):
    rng = np.random.RandomState(seed)
    fold_of = np.full(len(y), -1)
    for cls in np.unique(y):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        for pos, i in enumerate(idx):
            fold_of[i] = pos % k
    return [{"train": list(np.where(fold_of != f)[0]),
             "test": list(np.where(fold_of == f)[0])} for f in range(k)]


def main():
    import scanpy as sc
    adata = sc.read_h5ad(COHORT)
    splits = json.load(open(SPLITS))
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else \
        np.asarray(adata.X)
    genes = list(adata.var_names)
    gi = {g: i for i, g in enumerate(genes)}
    cyt_genes = [g for g in DECLARED if g in gi]
    assert len(cyt_genes) == 2, f"declared genes missing: {DECLARED}"
    cyt_score = X[:, [gi[g] for g in cyt_genes]].mean(axis=1)
    y_cyt = (cyt_score > np.median(cyt_score)).astype(int)

    recist_ok = adata.obs["recist_available"].values.astype(int) == 1
    y_recist = adata.obs["responder_true"].values.astype(int)  # 1=TRUE,0=FALSE
    dcb6 = adata.obs["dcb6"].values.astype(int)
    pfs6 = adata.obs["pfs6_timeonly"].values.astype(int)
    n = len(y_recist)
    print(f"cohort: n={n} RECIST-ok={int(recist_ok.sum())} "
          f"(responders {int(y_recist[recist_ok].sum())}) "
          f"DCB6mo={int(dcb6.sum())} PFS<6mo={int(pfs6.sum())} "
          f"cytolytic-high={int(y_cyt.sum())}", flush=True)

    out = {"meta": {
        "cohort": "IMmotion150 (ccRCC; McDermott 2018 Nat Med; iAtlas "
                  "harmonized via cBioPortal API)",
        "n": int(n),
        "arm": "atezolizumab pooled (mono + atezo-bev arms not separable "
               "in iAtlas; sunitinib arm excluded)",
        "endpoints": {
            "RECIST": {"definition": "CR/PR vs SD/PD (RESPONDER field)",
                       "n_responders": int(y_recist[recist_ok].sum()),
                       "n_evaluable": int(recist_ok.sum())},
            "DCB6mo": {"definition": "PFS >= 6 months AND no progression "
                                     "event (Jung-2019 durable-benefit "
                                     "convention); primary PFS endpoint",
                       "n_benefit": int(dcb6.sum())},
            "PFS6mo": {"definition": "PFS < 6 months, time-only "
                                     "(GSE274975-comparable sensitivity)",
                       "n_events": int(pfs6.sum())},
            "cytolytic_surrogate": {"definition": "GZMA/PRF1 mean, median "
                                                  "split (as Riaz); "
                                                  "exploratory addition",
                                    "n_high": int(y_cyt.sum()),
                                    "declared_genes": cyt_genes}},
        "identity_mapping": "cBioPortal sampleIds join directly (1 sample "
                            "per patient)",
        "expression": "iAtlas TPM via API -> log2(TPM+1); duplicate symbols "
                      "mean-aggregated; GZMB/TIGIT/CD276 absent in profile "
                      "(GEP 7/8, exhaustion 5/6, IMPRES pairs degraded "
                      "automatically by the wrappers)",
        "seed": SEED,
        "protocol": "v33 leakage-free nested CV; CDS v1.1.0; identical to "
                    "executed e1_gse274975 / e1_liu2019",
    }}

    def run_endpoint(ename, y, folds_e):
        adata.obs["response"] = y
        res = {"k": int(min(5, int((y == 1).sum()), int((y == 0).sum())))}
        bench_methods = ["GEP", "TIDE", "IMPRES", "PD_L1_IHC",
                         "ElasticNet", "ElasticNet_Var"]
        for meth in bench_methods:
            res[meth] = {}
            if meth in ("ElasticNet", "ElasticNet_Var"):
                params = tune_primary(adata, folds_e, meth)
                yt, yp, _ = nested_cv_predict(adata, folds_e, meth,
                                              params=params)
                res[meth] = {"params": params}
            else:
                yt, yp = cv_predict(adata, folds_e, meth)
            r = evaluate_all(np.asarray(yt).astype(int),
                             np.asarray(yp, dtype=float))
            res[meth].update({
                "auroc": round(r["auroc"], 4),
                "auroc_ci": [round(r["auroc_ci_low"], 4),
                             round(r["auroc_ci_high"], 4)],
                "auprc": round(r["auprc"], 4)})
            print(f"  {ename} {meth}: {r['auroc']:.3f} "
                  f"[{r['auroc_ci_low']:.3f}-{r['auroc_ci_high']:.3f}]",
                  flush=True)
        for meth in ("IMPRES", "GEP", "TIDE", "PD_L1_IHC"):
            yt_f, yp = cv_predict(adata, folds_e, meth)
            yt_f = np.asarray(yt_f).astype(int)
            obs_, p = perm_fixed(np.asarray(yp, dtype=float), yt_f)
            res[meth]["perm_p"] = p
            res[meth]["n_shuffle"] = 50000
            print(f"  {ename} {meth} perm p={p:.5f}", flush=True)
        return res

    k_rec = splits["k"]
    out["RECIST"] = run_endpoint("RECIST", y_recist, splits["folds"])
    out["DCB6mo"] = run_endpoint("DCB6mo", dcb6,
                                 folds_from_subset(dcb6, k_rec, SEED))
    out["PFS6mo"] = run_endpoint("PFS6mo", pfs6,
                                 folds_from_subset(pfs6, k_rec, SEED))
    out["cytolytic_surrogate"] = run_endpoint("cytolytic_surrogate", y_cyt,
                                              folds_from_subset(y_cyt, k_rec,
                                                                SEED))

    bench = ["GEP", "TIDE", "IMPRES", "PD_L1_IHC", "ElasticNet",
             "ElasticNet_Var"]
    out["endpoint_switch_delta"] = {
        meth: {"DCB6mo_minus_RECIST": round(out["DCB6mo"][meth]["auroc"]
                                            - out["RECIST"][meth]["auroc"], 3),
               "PFS6mo_minus_RECIST": round(out["PFS6mo"][meth]["auroc"]
                                            - out["RECIST"][meth]["auroc"], 3),
               "surrogate_minus_RECIST": round(
                   out["cytolytic_surrogate"][meth]["auroc"]
                   - out["RECIST"][meth]["auroc"], 3)}
        for meth in bench}

    from cds import circularity_detection_score
    rng = np.random.RandomState(SEED)
    out.setdefault("cds", {})
    for ename, y in (("RECIST", y_recist), ("DCB6mo", dcb6),
                     ("PFS6mo", pfs6), ("cytolytic_surrogate", y_cyt)):
        r = circularity_detection_score(X, y, cyt_genes, genes)
        r = json.loads(json.dumps(r, default=float))
        nulls = [float(circularity_detection_score(
            X, rng.permutation(y), cyt_genes, genes)["CDS"])
            for _ in range(10)]
        out["cds"][ename] = {
            "CDS": r["CDS"], "risk": r["risk_level"],
            "components": r["components"],
            "null_mean": round(float(np.mean(nulls)), 3),
            "null_max": round(float(np.max(nulls)), 3)}
        print(f"  CDS {ename}: {r['CDS']:.3f} "
              f"(null mean {np.mean(nulls):.3f})", flush=True)

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print("WROTE", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
