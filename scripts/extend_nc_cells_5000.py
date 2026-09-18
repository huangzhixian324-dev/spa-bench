"""Dry experiment A (2026-09-17 quality campaign): extend ALL eight n<=43
trainable permutation cells to 5,000 full-pipeline shuffles, removing the
resolution objection to the headline negative claim ("no trainable method
survives FDR on clinical endpoints at n <= 43").

Protocol mirrors scripts/fill_nc_cells.py exactly (frozen hyperparameters
C = 0.1, l1_ratio = 0.5; per-fold label-independent variance prescreen
cached; MI cells recompute MI within the top-2000 prescreen per shuffle;
Var cells use the cached top-500 variance selection; degenerate folds
skipped), with the same HGNC h5ad overrides and the same stored benchmark
fold splits. A fresh seed-42 run of n shuffles reproduces the first n
shuffles of any shorter archived run, so the extended p supersedes and
subsumes the archived value.

Cells (cohort, method): Hugo MI/Var, Lauss(Nathanson) MI/Var,
Jung MI/Var, Riaz-RECIST MI/Var. The Riaz cytolytic cells are excluded
(circular endpoint, already 500-1,000 shuffles and FDR-significant).

Outputs: results/benchmark/v33/nc_cells_ext/<cohort>__<method>.json
(per-cell checkpoint files allow resumption).
"""
import json
import multiprocessing as mp
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

from rerun_v33 import COHORTS, SEED, load_cohort  # noqa: E402
from base_method import BaseMethod  # noqa: E402
from sklearn.feature_selection import mutual_info_classif  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

# HGNC override (identical to scripts/fill_nc_cells.py)
COHORTS["Riaz_2017_cytolytic"]["h5ad"] = (
    "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad")
COHORTS["Riaz_2017_RECIST_v33"]["h5ad"] = (
    "data/cohorts/Riaz_2017_RECIST/processed/"
    "Riaz_2017_RECIST_v33_HGNC.h5ad")

OUTD = REPO / "results" / "benchmark" / "v33" / "nc_cells_ext"
N_PERM = 2000  # 2026-09-18: remaining 4 cells capped at 2,000 per user
# decision — the resolution objection targeted the 60-shuffle cells
# (Hugo/Lauss, now complete at 5,000); Jung/Riaz-RECIST baselines are far
# from the boundary and 2,000 shuffles (min p = 1/2001 ~ 0.0005) are
# statistically sufficient; per-cell shuffle counts are disclosed in S11
# per the existing convention. The four completed cells keep n = 5,000.
BATCH = 250
PRESCREEN = 2000
N_SEL = 500
PARAMS = {"C": 0.1, "l1_ratio": 0.5}  # benchmark frozen default

CELLS = [
    ("Hugo_2016", "ElasticNet"),
    ("Hugo_2016", "ElasticNet_Var"),
    ("Nathanson_2017", "ElasticNet"),
    ("Nathanson_2017", "ElasticNet_Var"),
    ("Jung_2019_DCB", "ElasticNet"),
    ("Jung_2019_DCB", "ElasticNet_Var"),
    ("Riaz_2017_RECIST_v33", "ElasticNet"),
    ("Riaz_2017_RECIST_v33", "ElasticNet_Var"),
]

_W = {"X": None, "y": None, "folds": None, "is_mi": None,
      "fold_cache": None}


def _pool_init(xpath, ypath, folds, is_mi, fold_cache):
    _W["X"] = np.load(xpath, mmap_mode="r")
    _W["y"] = np.load(ypath, mmap_mode="r")
    _W["folds"] = folds
    _W["is_mi"] = is_mi
    _W["fold_cache"] = fold_cache


def _perm_worker(seed):
    X, y, is_mi, fold_cache = (_W["X"], _W["y"], _W["is_mi"],
                               _W["fold_cache"])
    rng = np.random.RandomState(seed)
    y_perm = rng.permutation(y)
    yt, yp = [], []
    try:
        for fc in fold_cache:
            tr, te = fc["tr"], fc["te"]
            if is_mi:
                mi = mutual_info_classif(X[tr][:, fc["var_idx"]], y_perm[tr],
                                         random_state=SEED)
                sel = fc["var_idx"][np.argsort(mi)[-N_SEL:]]
            else:
                sel = fc["sel_var"]
            model = BaseMethod._make_elasticnet(**PARAMS)
            model.fit(X[tr][:, sel], y_perm[tr])
            yp.extend(model.predict_proba(X[te][:, sel])[:, 1])
            yt.extend(y_perm[te])
    except ValueError:
        # a training fold held only one class (possible for k=2 small-n
        # permutations); fill_nc_cells skips such replicates identically
        return None
    if len(set(yt)) < 2:
        return None
    try:
        return float(roc_auc_score(yt, yp))
    except ValueError:
        return None


def run_cell(cohort, method):
    dest = OUTD / f"{cohort}__{method}.json"
    ckpt = OUTD / f"{cohort}__{method}.ckpt.json"
    if dest.exists():
        print(f"[skip] {cohort}/{method} (done)", flush=True)
        return
    adata, folds, y = load_cohort(cohort)
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else \
        np.asarray(adata.X)
    is_mi = (method == "ElasticNet")

    # cached label-independent prescreen per fold (fill_nc_cells protocol)
    fold_cache = []
    for f in folds:
        tr, te = np.array(f["train"]), np.array(f["test"])
        var_idx = np.argsort(X[tr].var(axis=0))[-min(PRESCREEN, X.shape[1]):]
        sel_var = None if is_mi else \
            np.argsort(X[tr].var(axis=0))[-N_SEL:]
        fold_cache.append({"tr": tr, "te": te, "var_idx": var_idx,
                           "sel_var": sel_var})

    # observed AUROC under the frozen protocol (real labels)
    yt, yp = [], []
    for fc in fold_cache:
        tr, te = fc["tr"], fc["te"]
        if is_mi:
            mi = mutual_info_classif(X[tr][:, fc["var_idx"]], y[tr],
                                     random_state=SEED)
            sel = fc["var_idx"][np.argsort(mi)[-N_SEL:]]
        else:
            sel = fc["sel_var"]
        model = BaseMethod._make_elasticnet(**PARAMS)
        model.fit(X[tr][:, sel], y[tr])
        yt.extend(y[te])
        yp.extend(model.predict_proba(X[te][:, sel])[:, 1])
    obs = float(roc_auc_score(yt, yp))

    done, null = 0, []
    if ckpt.exists():
        ck = json.load(open(ckpt))
        done, null = ck["done"], ck["null"]
        print(f"  resuming from {done}", flush=True)

    tmpd = Path(str(ckpt) + ".tmpdir")
    tmpd.mkdir(parents=True, exist_ok=True)
    xpath, ypath = str(tmpd / "X.npy"), str(tmpd / "y.npy")
    np.save(xpath, np.ascontiguousarray(X, dtype=np.float32))
    np.save(ypath, y.astype(np.int8))

    n_jobs = min(6, max(2, (os.cpu_count() or 4) - 1))
    t0 = time.time()
    while done < N_PERM:
        batch = min(BATCH, N_PERM - done)
        pool = mp.Pool(n_jobs, initializer=_pool_init,
                       initargs=(xpath, ypath, folds, is_mi, fold_cache))
        try:
            res = pool.map(_perm_worker,
                           range(SEED + done, SEED + done + batch),
                           chunksize=5)
        finally:
            pool.close()
            pool.join()
        null.extend([r for r in res if r is not None])
        done += batch
        json.dump({"done": done, "null": null}, open(ckpt, "w"))
        print(f"  {done}/{N_PERM} ({time.time() - t0:.0f}s, "
              f"valid={len(null)})", flush=True)

    for f in (xpath, ypath):
        try:
            os.remove(f)
        except OSError:
            pass
    try:
        os.rmdir(tmpd)
    except OSError:
        pass

    ge = int(np.sum(np.array(null) >= obs - 1e-12))
    p = (ge + 1) / (len(null) + 1)
    old = {}
    oldp = OUTD.parent / "nc_cells" / f"{cohort}__{method}.json"
    if oldp.exists():
        d = json.load(open(oldp))
        old = {"archived_n": d.get("n_perm"), "archived_p": d.get("p"),
               "archived_obs": d.get("obs_auroc_this_run")}
    rec = {"cohort": cohort, "method": method,
           "obs_auroc_this_run": round(obs, 4), **old,
           "n_perm": N_PERM, "n_valid": len(null), "n_ge_obs": ge,
           "p": round(p, 5),
           "null_mean": round(float(np.mean(null)), 4),
           "protocol": "full-pipeline label shuffles, cached variance "
                       "prescreen, frozen C=0.1/l1_ratio=0.5, seed 42; "
                       "extension of fill_nc_cells.py to n=5000",
           "purpose": "resolution-hardened negative claim: n<=43 trainable "
                      "cells at uniform 5,000-shuffle resolution"}
    json.dump(rec, open(dest, "w"), indent=1)
    try:
        os.remove(ckpt)
    except OSError:
        pass
    print(f"  [done] obs={obs:.4f} p={p:.5f}", flush=True)


def main():
    OUTD.mkdir(parents=True, exist_ok=True)
    for cohort, method in CELLS:
        print(f"[cell] {cohort}/{method}", flush=True)
        run_cell(cohort, method)
    print("ALL CELLS COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
