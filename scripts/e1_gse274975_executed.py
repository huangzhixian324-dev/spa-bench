"""E1 (executed): GSE274975 dual-endpoint replication of Design Choice 1.

Executes the pre-declared pipeline (scripts/e1_gse274975_replication.py,
strategy-doc P0-1) on GSE274975 (Luppov et al., Front Immunol 2024,
NSCLC, n = 61, anti-PD(L)-1).  The pre-declared script expected RECIST
and PFS in the GEO sample characteristics; they are in fact published in
the paper's Supplementary Table 1, so the executed version reads
endpoints from Table1.xlsx via the verified patient-number identity map
(raw/id_map.json; perfect 61x61x61 bijection, histotype anchor
57/58).  Cohort construction: scripts/build_gse274975_cohort.py.

Design (fixed seed 42, no parameters tuned on the outcome) -- identical
to the executed Liu-2019 E1 (e1_replication_liu2019.py):
  endpoints: RECIST (CR/PR = responder; the paper's own definition,
             identical to the SPATBench convention) and PFS6mo
             (PFS < 6 months = early progressor; the pre-declared DCB
             dichotomisation, using PFS time only).
  - leakage-free nested CV, k = min(5, n_pos, n_neg) = 5, frozen
    hyperparameters from tune_primary on the real labels;
  - fixed scorers (IMPRES/GEP/TIDE/PD-L1): 50,000 prediction shuffles
    (exact null);
  - CDS v1.1.0 on both endpoints (declared genes GZMA/PRF1) +
    10-replicate random-label null + cytolytic positive control is
    inherent in the declared genes;
  - expression log2(CPM+1) (pre-declared; FPKM not recoverable).

Output: results/benchmark/v33/e1_gse274975.json
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

COHORT = REPO / "data/cohorts/GSE274975/processed/GSE274975_HGNC.h5ad"
SPLITS = REPO / "data/cohorts/GSE274975/processed/GSE274975_splits.json"
OUT = REPO / "results/benchmark/v33/e1_gse274975.json"
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


def main():
    import scanpy as sc
    adata = sc.read_h5ad(COHORT)
    folds = json.load(open(SPLITS))["folds"]
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else \
        np.asarray(adata.X)
    genes = list(adata.var_names)
    gi = {g: i for i, g in enumerate(genes)}
    y_recist = adata.obs["response"].values.astype(int)
    pfs = adata.obs["pfs_months"].values.astype(float)
    y_pfs = (pfs < 6).astype(int)  # pre-declared dichotomisation
    cyt_genes = [g for g in DECLARED if g in gi]
    assert len(cyt_genes) == 2, f"declared genes missing: {DECLARED}"
    cyt_score = X[:, [gi[g] for g in cyt_genes]].mean(axis=1)
    y_cyt = (cyt_score > np.median(cyt_score)).astype(int)

    n = len(y_recist)
    print(f"cohort: n={n} RECIST responders={int(y_recist.sum())} "
          f"PFS<6mo={int(y_pfs.sum())} cytolytic-high={int(y_cyt.sum())}",
          flush=True)

    out = {"meta": {
        "cohort": "GSE274975 (Oncobox NSCLC ICI; Luppov 2024 Front Immunol)",
        "n": int(n),
        "endpoints": {
            "RECIST": {"definition": "CR/PR vs SD/PD (paper's own "
                                     "definition; SPATBench convention)",
                       "n_responders": int(y_recist.sum())},
            "PFS6mo": {"definition": "PFS < 6 months = early progressor "
                                     "(pre-declared DCB dichotomisation; "
                                     "PFS time only)",
                       "n_events": int(y_pfs.sum())},
            "cytolytic_surrogate": {"definition": "GZMA/PRF1 mean, median "
                                                  "split (as Riaz); "
                                                  "exploratory addition",
                                    "n_high": int(y_cyt.sum()),
                                    "declared_genes": cyt_genes}},
        "identity_mapping": "patient-number bijection "
                            "(data/cohorts/GSE274975/raw/id_map.json); "
                            "histotype anchor 57/58",
        "expression": "log2(CPM+1) from GEO raw counts",
        "seed": SEED,
        "protocol": "v33 leakage-free nested CV; CDS v1.1.0; "
                    "identical to executed e1_liu2019_replication",
    }}
    y_map = (("RECIST", y_recist), ("PFS6mo", y_pfs),
             ("cytolytic_surrogate", y_cyt))

    for ename, y in y_map:
        k = min(5, int((y == 1).sum()), int((y == 0).sum()))
        folds_e = folds if ename == "RECIST" else None
        if folds_e is None:
            # same stratified construction as build script (seed 42)
            rng = np.random.RandomState(SEED)
            fold_of = np.full(len(y), -1)
            for cls in np.unique(y):
                idx = np.where(y == cls)[0]
                rng.shuffle(idx)
                for pos, i in enumerate(idx):
                    fold_of[i] = pos % k
            folds_e = [{"train": list(np.where(fold_of != f)[0]),
                        "test": list(np.where(fold_of == f)[0])}
                       for f in range(k)]
        adata.obs["response"] = y
        res = {"k": int(k)}
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
        # exact permutations for fixed scorers (shuffled labels must be in
        # the SAME fold-concatenated order as the OOF scores)
        for meth in ("IMPRES", "GEP", "TIDE", "PD_L1_IHC"):
            yt_f, yp = cv_predict(adata, folds_e, meth)
            yt_f = np.asarray(yt_f).astype(int)
            obs_, p = perm_fixed(np.asarray(yp, dtype=float), yt_f)
            res[meth]["perm_p"] = p
            res[meth]["n_shuffle"] = 50000
            print(f"  {ename} {meth} perm p={p:.5f}", flush=True)
        out[ename] = res

    # endpoint-switch deltas (Design Choice 1 replication)
    bench = ["GEP", "TIDE", "IMPRES", "PD_L1_IHC", "ElasticNet",
             "ElasticNet_Var"]
    out["endpoint_switch_delta"] = {
        meth: round(out["PFS6mo"][meth]["auroc"]
                    - out["RECIST"][meth]["auroc"], 3)
        for meth in bench}

    # CDS on all three endpoints
    from cds import circularity_detection_score
    rng = np.random.RandomState(SEED)
    out.setdefault("cds", {})
    for ename, y in y_map:
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
