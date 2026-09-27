"""10,000-resample bootstrap upgrade for the Riaz Delta intervals
(reviewer minor R1-m7: 1,000 resamples give limited tail precision for the
95 percent CI of the primary inferential anchor).

Recomputes, from the archived per-patient OOF predictions
(riaz_oof_scores_v33.json / riaz_oof_scores_tuning_v33.json; same y vectors
across endpoints, verified), the within-patient collapse
Delta = AUROC(cytolytic EN-Var) - AUROC(RECIST EN-Var) on the 42
RECIST-evaluable patients, with 10,000 patient resamples (seed 42).

Output: results/benchmark/v33/effect_size_delta_10k_v33.json
"""
import json
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
V = REPO / "results" / "benchmark" / "v33"
N = 10000
SEED = 42

out = {}
for track, fn in (("frozen", "riaz_oof_scores_v33.json"),
                  ("tuning", "riaz_oof_scores_tuning_v33.json")):
    d = json.load(open(V / fn, encoding="utf-8"))
    cyt = d["Riaz_2017_cytolytic"]["ElasticNet_Var"]
    rec = d["Riaz_2017_RECIST_v33"]["ElasticNet_Var"]
    y_c, p_c = np.asarray(cyt["y"]), np.asarray(cyt["p"], dtype=float)
    y_r, p_r = np.asarray(rec["y"]), np.asarray(rec["p"], dtype=float)
    assert len(y_c) == len(y_r) == 42
    # rows are the same 42 RECIST-evaluable patients in the same order;
    # y vectors differ because the two endpoints have different labels
    auc_c = roc_auc_score(y_c, p_c)
    auc_r = roc_auc_score(y_r, p_r)
    delta_point = auc_c - auc_r
    rng = np.random.RandomState(SEED)
    n = len(y_c)
    deltas = np.empty(N)
    for i in range(N):
        idx = rng.randint(0, n, n)
        # resample with replacement; skip degenerate all-one-class draws
        if len(np.unique(y_c[idx])) < 2 or len(np.unique(y_r[idx])) < 2:
            deltas[i] = np.nan
            continue
        deltas[i] = (roc_auc_score(y_c[idx], p_c[idx]) -
                     roc_auc_score(y_r[idx], p_r[idx]))
    deltas = deltas[~np.isnan(deltas)]
    ci = (round(float(np.quantile(deltas, 0.025)), 4),
          round(float(np.quantile(deltas, 0.975)), 4))
    out[track] = {"n_resamples": int(len(deltas)), "seed": SEED,
                  "auc_cytolytic_full": round(float(auc_c), 4),
                  "auc_recist_full": round(float(auc_r), 4),
                  "delta_point": round(float(delta_point), 4),
                  "delta_ci95_10000": list(ci)}
    print(f"{track}: delta {delta_point:.4f}, 10k CI {ci} "
          f"(n_valid={len(deltas)})", flush=True)

(V / "effect_size_delta_10k_v33.json").write_text(
    json.dumps(out, indent=1), encoding="utf-8")
print("written:", V / "effect_size_delta_10k_v33.json", flush=True)
