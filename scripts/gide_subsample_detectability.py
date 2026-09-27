"""Gide within-cohort subsample detectability curve (reviewer R1-M2):
does detectability degrade with n inside the same cohort (isolating the
sample-size effect from cohort identity)?

For n in {25, 35, 45, 55, 65, 73}: 10 stratified subsample draws
(class proportion preserved, ~40/73), ElasticNet-MI with frozen modal
hyperparameters, 5-fold stratified CV AUROC per draw. Per n: mean/min AUROC,
detection floor = one-sample Hanley-McNeil minimum AUROC at 80% power
(one-sided alpha = 0.05, z = 1.645; identical formula to
scripts/power_onesample_v33.py), and the fraction of draws at or above the
floor.

Output: results/benchmark/v33/gide_subsample_detectability.json
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
spec = importlib.util.spec_from_file_location(
    "rerun_v33", REPO / "scripts" / "rerun_v33.py")
rerun = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rerun)
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402

OUT = rerun.OUT / "gide_subsample_detectability.json"
METHOD = "ElasticNet"
NS = [25, 35, 45, 55, 65, 73]
DRAWS = 10
SEED = 42


def se_auroc(A, n_pos, n_neg):
    Q1 = A / (2 - A)
    Q2 = 2 * A ** 2 / (1 + A)
    var = (A * (1 - A) + (n_pos - 1) * (Q1 - A ** 2) +
           (n_neg - 1) * (Q2 - A ** 2)) / (n_pos * n_neg)
    return float(np.sqrt(var))


def min_auroc(n_pos, n_neg, zcrit=1.645, target=0.80):
    lo, hi = 0.5, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        p = norm.cdf((mid - 0.5) / se_auroc(mid, n_pos, n_neg) - zcrit)
        if p < target:
            lo = mid
        else:
            hi = mid
    return float(hi)


def main():
    adata, folds, y = rerun.load_cohort("Gide_2019_cBio")
    adata.obs["response"] = y
    y = y.astype(int)
    params = rerun.tune_primary(adata, folds, METHOD)
    print("frozen params:", params, flush=True)
    n_total, n_pos_total = len(y), int(y.sum())
    out = {"meta": {"cohort": "Gide_2019_cBio", "method": METHOD,
                    "frozen_params": params, "draws_per_n": DRAWS,
                    "seed_base": SEED, "k": 5,
                    "floor": "one-sample Hanley-McNeil, 80% power, "
                             "one-sided alpha 0.05 (same formula as "
                             "power_onesample_v33.py)"},
           "curve": []}
    for n in NS:
        n_pos = int(round(n * n_pos_total / n_total))
        n_neg = n - n_pos
        floor = min_auroc(n_pos, n_neg)
        aucs = []
        for i in range(DRAWS):
            rng = np.random.RandomState(SEED + 1000 + i)
            pos_idx = np.where(y == 1)[0]
            neg_idx = np.where(y == 0)[0]
            sel_pos = rng.choice(pos_idx, size=min(n_pos, len(pos_idx)),
                                 replace=False)
            sel_neg = rng.choice(neg_idx, size=min(n_neg, len(neg_idx)),
                                 replace=False)
            idx = np.concatenate([sel_pos, sel_neg])
            rng.shuffle(idx)
            ysub = y[idx]
            skf = StratifiedKFold(n_splits=5, shuffle=True,
                                  random_state=SEED + 2000 + i)
            ff = [{"train": idx[tr].tolist(), "test": idx[te].tolist()}
                  for tr, te in skf.split(np.zeros(len(idx)), ysub)]
            yt, yp, _ = rerun.nested_cv_predict(adata, ff, METHOD,
                                                params=params)
            aucs.append(float(roc_auc_score(yt, yp)))
        rec = {"n": n, "n_pos": n_pos, "n_neg": n_neg,
               "floor_80pct": round(floor, 3),
               "mean_auroc": round(float(np.mean(aucs)), 3),
               "min_auroc": round(float(np.min(aucs)), 3),
               "max_auroc": round(float(np.max(aucs)), 3),
               "draws_at_or_above_floor": int(sum(
                   a >= floor for a in aucs)), "per_draw": [round(a, 4) for a in aucs]}
        out["curve"].append(rec)
        print(f"  n={n}: floor {rec['floor_80pct']}, mean AUROC "
              f"{rec['mean_auroc']}, draws>=floor "
              f"{rec['draws_at_or_above_floor']}/10", flush=True)
    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written:", OUT, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
