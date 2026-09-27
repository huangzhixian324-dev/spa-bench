"""Frozen vs re-tuned permutation convergence for one negative cell
(reviewer items R1-M3 / R3-M4): Hugo 2016, ElasticNet-Var, RECIST.

For N_PERM label shuffles (separate stream, RandomState(90000+i)):
  frozen track: nested_cv_predict with modal hyperparameters frozen
  re-tune track: nested_cv_predict with per-fold inner 3-fold GridSearchCV
  (hyperparameters re-selected on the shuffled labels, standard practice)
Each track yields its own permutation null and p-value for the same
observed AUROC, testing whether the frozen protocol's negative conclusion
is an artifact of the frozen protocol.

Checkpointed every 25 shuffles.
Output: results/benchmark/v33/permutation_convergence_hugo_envar.json
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

COHORT = "Hugo_2016"
METHOD = "ElasticNet_Var"
N_PERM = 300
SEED_BASE = 90000
OUT = rerun.OUT / "permutation_convergence_hugo_envar.json"


def main():
    adata, folds, y = rerun.load_cohort(COHORT)
    adata.obs["response"] = y
    params = rerun.tune_primary(adata, folds, METHOD)
    yt, yp, _ = rerun.nested_cv_predict(adata, folds, METHOD, params=params)
    obs = float(roc_auc_score(yt, yp))
    print(f"Hugo EN-Var: obs AUROC {obs:.4f}, frozen params {params}",
          flush=True)

    ck = {"cohort": COHORT, "method": METHOD, "obs_auroc": obs,
          "frozen_params": params, "n_perm_target": N_PERM,
          "done": 0, "frozen_aucs": [], "retune_aucs": []}
    ck_path = OUT
    for i in range(N_PERM):
        rng = np.random.RandomState(SEED_BASE + i)
        y_perm = rng.permutation(y)
        yt1, yp1, _ = rerun.nested_cv_predict(adata, folds, METHOD,
                                              params=params,
                                              y_override=y_perm)
        a_f = float(roc_auc_score(yt1, yp1))
        yt2, yp2, _ = rerun.nested_cv_predict(adata, folds, METHOD,
                                              params=None,
                                              y_override=y_perm)
        a_r = float(roc_auc_score(yt2, yp2))
        ck["frozen_aucs"].append(a_f)
        ck["retune_aucs"].append(a_r)
        ck["done"] = i + 1
        if (i + 1) % 25 == 0:
            fa = np.array(ck["frozen_aucs"]); ra = np.array(ck["retune_aucs"])
            pf = (np.sum(fa >= obs - 1e-12) + 1) / (len(fa) + 1)
            pr = (np.sum(ra >= obs - 1e-12) + 1) / (len(ra) + 1)
            ck_path.write_text(json.dumps(ck, indent=1), encoding="utf-8")
            print(f"  [{i+1}/{N_PERM}] frozen p={pf:.4f} (mean null "
                  f"{fa.mean():.4f}) | retune p={pr:.4f} (mean null "
                  f"{ra.mean():.4f})", flush=True)
    fa = np.array(ck["frozen_aucs"]); ra = np.array(ck["retune_aucs"])
    ck["p_frozen"] = (np.sum(fa >= obs - 1e-12) + 1) / (len(fa) + 1)
    ck["p_retune"] = (np.sum(ra >= obs - 1e-12) + 1) / (len(ra) + 1)
    ck["null_mean_frozen"] = float(fa.mean())
    ck["null_mean_retune"] = float(ra.mean())
    ck["status"] = "complete"
    OUT.write_text(json.dumps(ck, indent=1), encoding="utf-8")
    print(f"FINAL: frozen p={ck['p_frozen']:.4f} | retune p={ck['p_retune']:.4f}",
          flush=True)
    print("written:", OUT, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
