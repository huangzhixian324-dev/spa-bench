"""Gide 2019 RECIST: 10-seed benchmark-level rerun (k-fold seed robustness).

Mirrors the archived Riaz RECIST k-fold robustness protocol exactly
(scripts/rerun_v33.py, "k-fold robustness" block; results in
riaz_recist_kfold_v33.json / Supplementary Table S16):

  - uses the pipeline's own functions (load_cohort / tune_primary /
    nested_cv_predict) with no reimplementation;
  - per seed, stratified k=5 folds are REGENERATED from scratch with
    RandomState(42 + seed) (class-wise shuffle, round-robin assignment);
  - hyperparameters are tuned once on the REAL labels via tune_primary
    (modal setting over the seed-42 primary folds) and frozen across the
    10 seed runs, as in the Riaz precedent (frozen modal setting);
  - per-seed AUROC for ElasticNet-MI ("ElasticNet") and ElasticNet-Var.

Permutation-level 10-seed reruns remain infeasible (60-90 s/shuffle x
5,000 x 10) and are out of scope, as disclosed in the revision.

Output -> results/benchmark/v33/gide_recist_10seed_v33.json
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "rerun_v33", REPO / "scripts" / "rerun_v33.py")
rerun = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rerun)

SEED = rerun.SEED
K = 5  # Gide primary design is k=5
OUT = rerun.OUT / "gide_recist_10seed_v33.json"


def stratified_folds(y, k, rng):
    """Verbatim regeneration logic from the Riaz k-fold precedent."""
    fold_of = np.full(len(y), -1)
    for cls in np.unique(y):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        for pos, s in enumerate(idx):
            fold_of[s] = pos % k
    return [{"train": list(np.where(fold_of != f)[0]),
             "test": list(np.where(fold_of == f)[0])}
            for f in range(k)]


def main():
    adata, primary_folds, y = rerun.load_cohort("Gide_2019_cBio")
    adata.obs["response"] = y  # Gide flip=False; keep obs consistent
    print(f"Gide_2019_cBio: {adata.n_obs} x {adata.n_vars} "
          f"pos={int(y.sum())}/{len(y)}", flush=True)

    # freeze hyperparameters: modal setting from tune_primary on the
    # seed-42 primary folds (same construction as the permutation runs)
    frozen = {}
    for m in ("ElasticNet", "ElasticNet_Var"):
        frozen[m] = rerun.tune_primary(adata, primary_folds, m)
        print(f"  frozen params {m}: {frozen[m]}", flush=True)

    out = {"meta": {"cohort": "Gide_2019_cBio", "endpoint": "RECIST",
                    "n": int(len(y)), "responders": int(y.sum()),
                    "k": K, "n_seeds": 10, "seed_base": SEED,
                    "frozen_params": frozen,
                    "protocol": "10 stratified seeds, folds regenerated "
                    "per seed (RandomState(42+seed)); pipeline-native "
                    "nested_cv_predict with frozen hyperparameters; "
                    "mirrors riaz_recist_kfold_v33.json"},
           "methods": {}}

    for m in ("ElasticNet", "ElasticNet_Var"):
        aucs = []
        per_seed = []
        for seed in range(10):
            rng = np.random.RandomState(SEED + seed)
            ff = stratified_folds(y.astype(int), K, rng)
            yt, yp, _ = rerun.nested_cv_predict(adata, ff, m,
                                                params=frozen[m])
            auc = float(roc_auc_score(yt, yp))
            aucs.append(auc)
            per_seed.append(round(auc, 4))
            print(f"  {m:15s} seed {SEED + seed}: AUROC = {auc:.4f}",
                  flush=True)
        out["methods"][m] = {
            "mean": round(float(np.mean(aucs)), 3),
            "std": round(float(np.std(aucs)), 3),
            "n_seeds": 10,
            "per_seed": per_seed,
        }
        print(f"  {m}: AUROC = {out['methods'][m]['mean']} ± "
              f"{out['methods'][m]['std']} (10 stratified seeds)",
              flush=True)

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"written: {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
