"""MI vs variance feature-selection overlap under the v33 pipeline convention.

The manuscript's Table 2 overlap figure (27%, 135/500) uses the two-stage
selection convention of the actual pipeline (rerun_v33._select_genes):
stage 1 keeps the top-2000 variance genes (training fold at runtime, full
data for this descriptive analysis), stage 2 selects the top-500 by MI
computed WITHIN that prescreen. The overlap is therefore between
(stage-2 MI top-500 within the prescreen) and (top-500 variance genes).

This script stores both conventions for transparency, plus the
out-of-fold prediction correlation between ElasticNet-MI and ElasticNet-Var
(manuscript: Spearman rho = 0.79).

Output: results/benchmark/v33/mi_var_overlap_v33.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from rerun_v33 import (SEED, load_cohort, nested_cv_predict,  # noqa: E402
                       tune_primary)
from sklearn.feature_selection import mutual_info_classif  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

OUT = REPO / "results" / "benchmark" / "v33" / "mi_var_overlap_v33.json"


def main():
    adata, folds, y = load_cohort("Hugo_2016")
    X = (adata.X.toarray() if hasattr(adata.X, "toarray")
         else np.asarray(adata.X))
    genes = list(adata.var_names)
    var = X.var(axis=0)
    top500_var = set(np.argsort(var)[-500:])

    top2000 = np.argsort(var)[-2000:]
    mi_pre = mutual_info_classif(X[:, top2000], y, random_state=SEED)
    top500_mi_pre = set(top2000[np.argsort(mi_pre)[-500:]])

    mi_full = mutual_info_classif(X, y, random_state=SEED)
    top500_mi_full = set(np.argsort(mi_full)[-500:])

    # Frozen hyperparameters: read from the recompute JSON (recorded by
    # tune_primary in the same process that produced the benchmark-adjacent
    # sensitivity analyses) instead of re-tuning here - tune_primary can
    # flip between tied parameter sets across processes (set iteration
    # order), so a single frozen source keeps rho reproducible.
    rec = json.load(open(REPO / "results" / "benchmark" / "v33" /
                         "supplementary_v33_recompute.json"))
    p_mi = rec["meta"]["params_mi"]
    p_var = rec["meta"]["params_var"]
    _, yp_mi, _ = nested_cv_predict(adata, folds, "ElasticNet", params=p_mi)
    _, yp_var, _ = nested_cv_predict(adata, folds, "ElasticNet_Var",
                                     params=p_var)
    rho, pval = spearmanr(yp_mi, yp_var)

    inter = sorted(top500_var & top500_mi_pre)
    out = {
        "meta": {"seed": SEED,
                 "convention": "stage-1 top-2000 variance prescreen; "
                               "stage-2 top-500 MI within prescreen (the "
                               "pipeline's actual two-stage selection); "
                               "descriptive full-data computation, not a "
                               "predictive estimate"},
        "overlap_prescreen": {"count": len(inter),
                              "pct": round(100.0 * len(inter) / 500, 1)},
        "overlap_fulldata_mi": {
            "count": len(top500_var & top500_mi_full),
            "pct": round(100.0 * len(top500_var & top500_mi_full) / 500, 1),
            "note": "alternative convention: MI on all genes without the "
                    "stage-1 prescreen; reported for transparency"},
        "oof_spearman": {"rho": round(float(rho), 3),
                         "p": float(pval)},
        "top500_var_genes": [genes[i] for i in sorted(top500_var)],
        "top500_mi_prescreen_genes": [genes[i] for i in
                                      sorted(top500_mi_pre)],
        "var_scores_all": [round(float(v), 5) for v in var],
        "mi_prescreen_scores": [
            (int(g), round(float(m), 5))
            for g, m in zip(top2000[np.argsort(mi_pre)[::-1]], mi_pre[np.argsort(mi_pre)[::-1]])],
        "genes_all": genes,
    }
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"overlap prescreen: {out['overlap_prescreen']}")
    print(f"overlap full-data: {out['overlap_fulldata_mi']['count']}")
    print(f"OOF spearman: {out['oof_spearman']}")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
