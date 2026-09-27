"""XGBoost permutation-matrix completion (reviewer items R1-M5 / R3-M1):
full-pipeline label-shuffle permutation (1,000 shuffles, default XGBoost
hyperparameters, same patient-stratified folds as rerun_v33) for the five
cohort-endpoints not covered by the archived Gide re-run
(xgb_exclusion_evidence_v33.json / Methods exploratory paragraph).

Cohort-endpoints: Hugo_2016, Nathanson_2017, Jung_2019_DCB,
Riaz_2017_cytolytic, Riaz_2017_RECIST_v33. Gide_2019_cBio already archived.

BH within each cohort over the 7-cell family (6 original methods + XGBoost)
is reported for context; the analysis remains exploratory per the manuscript.

Checkpointed every 100 shuffles.
Output: results/benchmark/v33/xgb_matrix_completion.json
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
spec = importlib.util.spec_from_file_location(
    "rerun_v33", REPO / "scripts" / "rerun_v33.py")
rerun = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rerun)

from xgboost import XGBClassifier  # noqa: E402
from sklearn.model_selection import GridSearchCV  # noqa: E402

N_PERM = 1000
SEED_BASE = 50000
OUT = rerun.OUT / "xgb_matrix_completion.json"
COHORTS = ["Hugo_2016", "Nathanson_2017", "Jung_2019_DCB",
           "Riaz_2017_cytolytic", "Riaz_2017_RECIST_v33"]


def xgb_oof(adata, folds, y_override=None):
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    y = (adata.obs["response"].values.astype(int) if y_override is None
         else np.asarray(y_override).astype(int))
    yt, yp = [], []
    for f in folds:
        tr, te = np.array(f["train"]), np.array(f["test"])
        clf = XGBClassifier(max_depth=3, n_estimators=100,
                            learning_rate=0.1, subsample=1.0,
                            colsample_bytree=1.0,
                            eval_metric="logloss", verbosity=0,
                            n_jobs=4, random_state=42)
        clf.fit(X[tr], y[tr])
        yt.extend(y[te])
        yp.extend(clf.predict_proba(X[te])[:, 1])
    return np.array(yt).astype(int), np.array(yp)


def main():
    results = {}
    for cname in COHORTS:
        adata, folds, y = rerun.load_cohort(cname)
        adata.obs["response"] = y
        yt, yp = xgb_oof(adata, folds)
        obs = float(roc_auc_score(yt, yp))
        nulls = []
        ck = {"cohort": cname, "obs_auroc": obs, "n_perm_target": N_PERM,
              "done": 0, "nulls": []}
        ck_path = OUT
        for i in range(N_PERM):
            rng = np.random.RandomState(SEED_BASE + i)
            y_perm = rng.permutation(y)
            ytp, ypp = xgb_oof(adata, folds, y_override=y_perm)
            nulls.append(float(roc_auc_score(ytp, ypp)))
            ck["done"] = i + 1
            ck["nulls"] = nulls
            if (i + 1) % 100 == 0:
                na = np.array(nulls)
                p = (np.sum(na >= obs - 1e-12) + 1) / (len(na) + 1)
                ck_path.write_text(json.dumps(ck, indent=1), encoding="utf-8")
                print(f"  {cname} [{i+1}/{N_PERM}] null mean {na.mean():.4f} "
                      f"p={p:.4f}", flush=True)
        na = np.array(nulls)
        p = (np.sum(na >= obs - 1e-12) + 1) / (N_PERM + 1)
        results[cname] = {"obs_auroc": obs, "n_perm": N_PERM,
                          "null_mean": float(na.mean()),
                          "null_q95": float(np.quantile(na, 0.95)),
                          "p": p}
        ck_path.write_text(json.dumps({"completed": results,
                                       "current": ck}, indent=1),
                           encoding="utf-8")
        print(f"{cname}: obs {obs:.4f} p={p:.4f}", flush=True)
    final = {"meta": {"protocol": "full-pipeline label shuffles, default "
                                 "XGBoost hyperparameters, same folds as "
                                 "rerun_v33; exploratory (not part of the "
                                 "36-cell matrix)",
                      "seed_base": SEED_BASE},
             "results": results}
    OUT.write_text(json.dumps(final, indent=1), encoding="utf-8")
    print("written:", OUT, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
