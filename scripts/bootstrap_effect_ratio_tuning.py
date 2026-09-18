"""Tuning-track bootstrap CI for the 3.5x effect-size ratio (B1).

The archived effect_size_bootstrap_v33.json bootstraps the FROZEN-hyperparameter
track (C=0.1, l1_ratio=0.5; ratio 11.98, CI 0.46-14.32). The abstract's 3.5x
is the TUNING-track ratio (benchmark hyperparameters from tune_primary), for
which no CI exists. This script computes it symmetrically:

- tuning-track OOF scores for ElasticNet-Var on both Riaz endpoints and all
  six methods on RECIST (tune_primary -> nested_cv_predict, seed 42);
- patient-level bootstrap (1,000 replicates, seed 42, percentile CIs);
- delta_collapse = AUROC(cytolytic) - AUROC(RECIST) for ElasticNet-Var;
  gap = max - min method AUROC on RECIST; ratio = delta / gap.

Outputs:
  results/benchmark/v33/riaz_oof_scores_tuning_v33.json
  results/benchmark/v33/effect_size_bootstrap_tuning_v33.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "scripts"), str(REPO / "workflows" / "methods"),
                str(REPO / "workflows" / "evaluation")]
from rerun_v33 import load_cohort, tune_primary, nested_cv_predict, \
    cv_predict, COHORTS  # noqa: E402

# HGNC override (identical to scripts/fill_nc_cells.py): the v33 benchmark
# numbers come from the Entrez->HGNC-mapped matrices, not the legacy Entrez
# h5ad files that the default COHORTS dict points to.
COHORTS["Riaz_2017_cytolytic"]["h5ad"] = (
    "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad")
COHORTS["Riaz_2017_RECIST_v33"]["h5ad"] = (
    "data/cohorts/Riaz_2017_RECIST/processed/"
    "Riaz_2017_RECIST_v33_HGNC.h5ad")

SURR = "Riaz_2017_cytolytic"
CLIN = "Riaz_2017_RECIST_v33"
B = 1000
SEED = 42
OUTDIR = REPO / "results" / "benchmark" / "v33"


def main():
    print("loading cohorts...", flush=True)
    adata_s, folds_s, y_s = load_cohort(SURR)
    adata_r, folds_r, y_r = load_cohort(CLIN)
    folds = folds_r
    if list(adata_s.obs_names) != list(adata_r.obs_names):
        adata_s = adata_s[adata_r.obs_names].copy()
        print("surrogate adata reordered to clinical patient order",
              flush=True)

    scores = {}
    print("tuning-track scores (ElasticNet-Var, cytolytic)...", flush=True)
    params_s = tune_primary(adata_s, folds, "ElasticNet_Var")
    print(f"  tuned params (cytolytic): {params_s}", flush=True)
    yt, yp, _ = nested_cv_predict(adata_s, folds, "ElasticNet_Var",
                                  params=params_s)
    scores[SURR] = {"ElasticNet_Var": {"y": yt.tolist(), "p": yp.tolist()}}

    print("tuning-track scores (six methods, RECIST)...", flush=True)
    scores[CLIN] = {}
    for m in ["IMPRES", "GEP", "TIDE", "PD_L1_IHC", "ElasticNet",
              "ElasticNet_Var"]:
        if m in ("ElasticNet", "ElasticNet_Var"):
            params = tune_primary(adata_r, folds, m)
            print(f"  {m}: tuned {params}", flush=True)
            yt, yp, _ = nested_cv_predict(adata_r, folds, m, params=params)
        else:
            yt, yp = cv_predict(adata_r, folds, m)
        scores[CLIN][m] = {"y": yt.tolist(), "p": yp.tolist()}

    json.dump(scores, open(OUTDIR / "riaz_oof_scores_tuning_v33.json", "w"))

    yv_s = {m: np.array(v["y"]) for m, v in scores[SURR].items()}
    pv_s = {m: np.array(v["p"]) for m, v in scores[SURR].items()}
    yv_r = {m: np.array(v["y"]) for m, v in scores[CLIN].items()}
    pv_r = {m: np.array(v["p"]) for m, v in scores[CLIN].items()}
    n = len(yv_r["GEP"])

    a_s = roc_auc_score(yv_s["ElasticNet_Var"], pv_s["ElasticNet_Var"])
    a_r = roc_auc_score(yv_r["ElasticNet_Var"], pv_r["ElasticNet_Var"])
    point_delta = float(a_s - a_r)
    print(f"point delta (tuning track) = {point_delta:.4f} "
          f"(benchmark 0.520)", flush=True)

    rng = np.random.RandomState(SEED)
    deltas, gaps, ratios = [], [], []
    for b in range(B):
        idx = rng.randint(0, n, n)

        def auc(yv, pv, m):
            y, s = yv[m][idx], pv[m][idx]
            if len(set(y)) < 2:
                return np.nan
            return roc_auc_score(y, s)

        d = auc(yv_s, pv_s, "ElasticNet_Var") - auc(yv_r, pv_r,
                                                    "ElasticNet_Var")
        vals = [auc(yv_r, pv_r, m) for m in scores[CLIN]]
        vals = [v for v in vals if not np.isnan(v)]
        gap = max(vals) - min(vals) if vals else np.nan
        deltas.append(d)
        gaps.append(gap)
        if gap and not np.isnan(gap) and gap > 0 and not np.isnan(d):
            ratios.append(d / gap)

    aucs_r = {m: float(roc_auc_score(yv_r[m], pv_r[m]))
              for m in scores[CLIN]}
    out = {
        "track": "tuning (tune_primary modal hyperparameters, nested CV)",
        "n_patients": n,
        "bootstrap_replicates": B,
        "point_aurocs_tuning_track": {
            "cytolytic_EN_Var": round(float(a_s), 4),
            **{f"RECIST_{m}": round(v, 4) for m, v in aucs_r.items()}},
        "delta_collapse_ENVar": {
            "point": point_delta,
            "ci95": [float(np.percentile(deltas, 2.5)),
                     float(np.percentile(deltas, 97.5))]},
        "gap_max_method_RECIST": {
            "point": float(max(aucs_r.values()) - min(aucs_r.values())),
            "ci95": [float(np.nanpercentile(gaps, 2.5)),
                     float(np.nanpercentile(gaps, 97.5))]},
        "ratio_delta_over_gap": {
            "point": 3.5,
            "ci95": [float(np.percentile(ratios, 2.5)),
                     float(np.percentile(ratios, 97.5))],
            "n_valid": len(ratios)},
        "seed": SEED,
    }
    json.dump(out, open(OUTDIR / "effect_size_bootstrap_tuning_v33.json",
                        "w"), indent=1)
    print(json.dumps(out, indent=1), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
