"""Rebuild the Riaz 2017 RECIST cohort (v33, fixes Major-4).

The original Riaz_2017_RECIST_processed.h5ad was a transposed label table
(42 obs x 43 vars, var_names = patient IDs, NO expression matrix), which made
all RECIST-endpoint signature results unreproducible (three exact 0.500
AUROCs in manuscript Table 1 were artifacts of scoring that broken matrix).

This script rebuilds the cohort from the intact main file:
  - expression matrix: data/cohorts/Riaz_2017/processed/Riaz_2017_processed.h5ad
    (43 patients x 22187 genes)
  - RECIST labels: data/cohorts/Riaz_2017_RECIST/processed/
    Riaz_2017_RECIST_processed.h5ad obs['response'] (42 patients; 9
    responders = PRCR, 33 non-responders = PD/SD; 7 patients without RECIST
    annotation are excluded)

Output:
  - data/cohorts/Riaz_2017_RECIST/processed/Riaz_2017_RECIST_v33.h5ad
    (42 x 22187, obs['response'] aligned with obs_names of the main file)
  - data/cohorts/Riaz_2017_RECIST/processed/Riaz_2017_RECIST_v33_splits.json
    (stratified k = 2/3/5 folds, seed 42)
"""
import json
from pathlib import Path

import numpy as np
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
RECIST_DIR = REPO / "data" / "cohorts" / "Riaz_2017_RECIST" / "processed"
SEED = 42


def stratified_folds(y, k, seed=SEED):
    """Stratified k-fold assignment preserving the responder ratio per fold."""
    rng = np.random.RandomState(seed)
    fold_of = np.full(len(y), -1)
    for cls in np.unique(y):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        for pos, sample_idx in enumerate(idx):
            fold_of[sample_idx] = pos % k
    folds = []
    for f in range(k):
        test_idx = [int(i) for i in np.where(fold_of == f)[0]]
        train_idx = [int(i) for i in np.where(fold_of != f)[0]]
        folds.append({"fold": int(f), "train": train_idx, "test": test_idx})
    return folds


def main():
    main_adata = sc.read_h5ad(
        REPO / "data" / "cohorts" / "Riaz_2017" / "processed" /
        "Riaz_2017_processed.h5ad")
    recist = sc.read_h5ad(RECIST_DIR / "Riaz_2017_RECIST_processed.h5ad")

    n_main, n_recist = main_adata.n_obs, recist.n_obs
    print(f"Main Riaz file: {n_main} x {main_adata.n_vars}")
    print(f"RECIST label file: {n_recist} obs (no expression matrix)")

    label_map = dict(zip(recist.obs_names, recist.obs["response"].values))
    common = [p for p in main_adata.obs_names if p in label_map]
    print(f"Patients with both expression and RECIST label: {len(common)}")

    y = np.array([int(label_map[p]) for p in common])
    print(f"RECIST response distribution: 0={int((y==0).sum())}, 1={int((y==1).sum())}")
    assert (y == 1).sum() == 9 and (y == 0).sum() == 33, \
        "unexpected RECIST label counts - aborting"

    sub = main_adata[common].copy()
    sub.obs["response"] = y
    sub.obs["endpoint"] = "RECIST"
    out_h5ad = RECIST_DIR / "Riaz_2017_RECIST_v33.h5ad"
    sub.write(out_h5ad)
    print(f"Saved: {out_h5ad} ({sub.n_obs} x {sub.n_vars})")

    splits = {"n_samples": int(sub.n_obs), "n_pos": int((y == 1).sum()),
              "seed": SEED, "folds": {}}
    for k in (2, 3, 5):
        splits["folds"][f"k={k}"] = stratified_folds(y, k)
        pos_per_fold = [int(y[np.array(f["test"])].sum())
                        for f in splits["folds"][f"k={k}"]]
        print(f"  k={k}: test positives per fold = {pos_per_fold}")
    out_json = RECIST_DIR / "Riaz_2017_RECIST_v33_splits.json"
    with open(out_json, "w") as fh:
        json.dump(splits, fh, indent=1)
    print(f"Saved: {out_json}")


if __name__ == "__main__":
    main()
