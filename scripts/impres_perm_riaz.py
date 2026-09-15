"""Permutation p for fixed (non-trainable) scorers via prediction shuffling.

For a scorer with no fitting step (IMPRES, GEP, TIDE, PD-L1) the full
pipeline permutation test (shuffle labels, re-score) is *identical* to the
prediction-permutation test: shuffle the labels against the fixed
out-of-fold scores and recompute the AUROC. We therefore estimate p with
50,000 label shuffles over the pooled out-of-fold predictions, which is
exact up to Monte Carlo error (~+/-0.002) and takes seconds.

Applies to every cohort/method already evaluated with cv_predict and lets us
report IMPRES/GEP/TIDE p-values (including the two Riaz cohorts on the
HGNC-mapped matrices) without re-running the slow per-fold MI path.
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
    COHORTS, SEED, load_cohort, cv_predict)

OUT = REPO / "results" / "benchmark" / "v33"
N_SHUFFLE = 50000
SIGNATURE_METHODS = ["IMPRES", "GEP", "TIDE"]


def prediction_perm_p(adata, folds, method, n=N_SHUFFLE, seed=SEED):
    yt, yp = cv_predict(adata, folds, method)
    yt = np.asarray(yt)
    yp = np.asarray(yp)
    obs = roc_auc_score(yt, yp)
    rng = np.random.RandomState(seed)
    ge = 0
    for _ in range(n):
        yp_ = rng.permutation(yt)
        try:
            auc = roc_auc_score(yp_, yp)
        except ValueError:
            auc = 0.5
        ge += int(auc >= obs - 1e-12)
    return float(obs), (ge + 1) / (n + 1)


def main():
    perms = {}
    targets = sys.argv[1:] or [
        "Riaz_2017_cytolytic", "Riaz_2017_RECIST_v33"]
    for cname in targets:
        # HGNC 口径用于 Riaz（与 benchmark_v33.json 一致）
        if cname in ("Riaz_2017_cytolytic", "Riaz_2017_RECIST_v33"):
            COHORTS[cname]["h5ad"] = (
                "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad"
                if cname == "Riaz_2017_cytolytic" else
                "data/cohorts/Riaz_2017_RECIST/processed/"
                "Riaz_2017_RECIST_v33_HGNC.h5ad")
        adata, folds, _ = load_cohort(cname)
        print(f"===== {cname} =====", flush=True)
        for m in SIGNATURE_METHODS:
            try:
                obs, p = prediction_perm_p(adata, folds, m)
            except Exception as e:
                print(f"  {m:10s} ERROR: {e}", flush=True)
                continue
            perms.setdefault(cname, {})[m] = {
                "obs_auroc": round(obs, 4), "p": round(p, 5),
                "n_shuffle": N_SHUFFLE}
            print(f"  {m:10s} obs={obs:.3f}  permutation p={p:.5f} "
                  f"({N_SHUFFLE} shuffles)", flush=True)
        try:
            with open(OUT / "permutation_v33.json") as fh:
                merged = json.load(fh)
            for m, v in perms.get(cname, {}).items():
                merged.setdefault("cohorts", {}).setdefault(cname, {})[m] = v
            with open(OUT / "permutation_v33.json", "w") as fh:
                json.dump(merged, fh, indent=1)
        except OSError:
            pass
    # merge with any already-recorded cells
    with open(OUT / "permutation_v33.json") as fh:
        merged = json.load(fh)
    for c, mv in perms.items():
        for m, v in mv.items():
            merged.setdefault("cohorts", {}).setdefault(c, {})[m] = v
    with open(OUT / "permutation_v33.json", "w") as fh:
        json.dump(merged, fh, indent=1)
    print("\nMerged into results/benchmark/v33/permutation_v33.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
