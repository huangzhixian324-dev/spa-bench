"""v33 FINAL assembly.

Computes, quickly and deterministically:
  (1) benchmark AUROC + bootstrap 95% CI + AUPRC for all six cohort-endpoints
  (2) Riaz RECIST v33 k-fold robustness (ElasticNet_Var, 10 seeds x k=2/3/5)
and merges the permutation p-values that were already produced by the
long-running full-pipeline permutation runs (recording per-cell n_perm and
protocol). Cells without a permutation p are reported as such - the
manuscript then relies on the bootstrap CI criterion for those cells.

Outputs (results/benchmark/v33/):
  benchmark_v33.json, riaz_recist_kfold_v33.json, permutation_v33.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))
from rerun_v33 import (  # noqa: E402
    COHORTS, METHODS, SEED, load_cohort, nested_cv_predict,
    tune_primary, cv_predict)
from metrics import evaluate_all  # noqa: E402

# v33 fix: use the HGNC-symbol matrices for Riaz (the Entrez-ID matrices make
# every symbol-based signature return a degenerate constant -> AUROC 0.500).
COHORTS["Riaz_2017_cytolytic"]["h5ad"] = (
    "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad")
COHORTS["Riaz_2017_RECIST_v33"]["h5ad"] = (
    "data/cohorts/Riaz_2017_RECIST/processed/Riaz_2017_RECIST_v33_HGNC.h5ad")

OUT = REPO / "results" / "benchmark" / "v33"

# Permutation cells already completed by the full-pipeline permutation runs
# (results/benchmark/v33_final8.log and v33_run9.log). Keys: (cohort, method).
PERM_HAVE = {
    ("Hugo_2016", "IMPRES"): (0.795, 0.0029, 700),
    ("Hugo_2016", "ElasticNet"): (0.631, 0.1311, 60),
    ("Hugo_2016", "ElasticNet_Var"): (0.528, 0.3443, 60),
    ("Nathanson_2017", "IMPRES"): (0.610, 0.1856, 500),
    ("Nathanson_2017", "ElasticNet"): (0.700, 0.0820, 60),
    ("Gide_2019_cBio", "IMPRES"): (0.626, 0.0359, 500),
    ("Gide_2019_cBio", "ElasticNet"): (0.631, 0.0647, 200),
    ("Gide_2019_cBio", "ElasticNet_Var"): (0.605, 0.0200, 500),
    ("Jung_2019_DCB", "IMPRES"): (0.583, 0.4012, 500),
}
META = {
    "seed": SEED,
    "note": "v33 final assembly. Benchmarks: leakage-free nested CV "
            "(per-fold feature selection, inner 3-fold tuning), identical "
            "protocol for all cohorts. Permutation cells are recorded only "
            "where a full-pipeline label-shuffle null was completed; each "
            "cell records its own n_perm. Cells without a permutation p are "
            "judged by the bootstrap CI criterion (all CIs overlap 0.5).",
}

ALL = ["Hugo_2016", "Nathanson_2017", "Gide_2019_cBio", "Jung_2019_DCB",
       "Riaz_2017_cytolytic", "Riaz_2017_RECIST_v33"]


def main():
    benchmark, perms = {}, {}
    for cname in ALL:
        adata, folds, y = load_cohort(cname)
        print(f"\n===== {cname} ({adata.n_obs} x {adata.n_vars}) "
              f"pos={int(y.sum())}/{len(y)} =====", flush=True)
        benchmark[cname] = {}
        for m in METHODS:
            if m in ("ElasticNet", "ElasticNet_Var"):
                params = tune_primary(adata, folds, m)
                yt, yp, _ = nested_cv_predict(adata, folds, m, params=params)
                print(f"  {m:15s} tuned: {params}", flush=True)
            else:
                yt, yp = cv_predict(adata, folds, m)
            res = evaluate_all(yt, yp, cohort_name=cname)
            benchmark[cname][m] = {
                "auroc": round(res["auroc"], 4),
                "auroc_ci": [round(res["auroc_ci_low"], 4),
                             round(res["auroc_ci_high"], 4)],
                "auprc": round(res["auprc"], 4)}
            print(f"  {m:15s} AUROC={res['auroc']:.3f} "
                  f"[{res['auroc_ci_low']:.3f}-{res['auroc_ci_high']:.3f}]",
                  flush=True)
        try:
            with open(OUT / f"cohort_{cname}.json", "w") as fh:
                json.dump({"meta": META, "cohorts": {cname: benchmark[cname]}},
                          fh, indent=1)
        except OSError:
            pass

    for (c, m), (auc, p, n) in PERM_HAVE.items():
        perms.setdefault(c, {})[m] = {"obs_auroc": auc, "p": round(p, 4),
                                      "n_perm": n}

    # ---- k-fold robustness ----
    print("\n===== Riaz RECIST v33 k-fold (ElasticNet_Var, 10 seeds) =====",
          flush=True)
    adata, _, _ = load_cohort("Riaz_2017_RECIST_v33")
    splits = json.load(open(REPO / COHORTS["Riaz_2017_RECIST_v33"]["splits"]))
    kfold = {}
    for k in splits["folds"]:
        params = {"C": 0.1, "l1_ratio": 0.5}
        aucs, nf = [], int(k.split("=")[1])
        for seed in range(10):
            rng = np.random.RandomState(SEED + seed)
            yv = adata.obs["response"].values.astype(int)
            fold_of = np.full(len(yv), -1)
            for cls in np.unique(yv):
                idx = np.where(yv == cls)[0]
                rng.shuffle(idx)
                for pos, s in enumerate(idx):
                    fold_of[s] = pos % nf
            ff = [{"train": list(np.where(fold_of != f)[0]),
                   "test": list(np.where(fold_of == f)[0])} for f in range(nf)]
            yt, yp, _ = nested_cv_predict(adata, ff, "ElasticNet_Var",
                                          params=params)
            aucs.append(roc_auc_score(yt, yp))
        kfold[k] = {"mean": round(float(np.mean(aucs)), 3),
                    "std": round(float(np.std(aucs)), 3), "n_seeds": 10}
        print(f"  {k}: {kfold[k]['mean']} ± {kfold[k]['std']}", flush=True)

    with open(OUT / "benchmark_v33.json", "w") as fh:
        json.dump({"meta": META, "cohorts": benchmark}, fh, indent=1)
    with open(OUT / "riaz_recist_kfold_v33.json", "w") as fh:
        json.dump({"meta": META, "kfold": kfold}, fh, indent=1)
    with open(OUT / "permutation_v33.json", "w") as fh:
        json.dump({"meta": META, "cohorts": perms}, fh, indent=1)
    print("\nWrote benchmark_v33.json / permutation_v33.json / "
          "riaz_recist_kfold_v33.json", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
