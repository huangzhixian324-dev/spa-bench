"""XGBoost Gide 2019 RECIST permutation test (fills the one untested
trainable cell; 2026-09-17 quality campaign).

Protocol mirrors fill_nc_cells.py / rerun_v33.permutation_p: full-pipeline
label shuffles. Per fold, the XGBoostWrapper protocol is reproduced exactly:
full-dimension mutual-information selection of the top-300 genes against the
(permutation) training labels, then XGBClassifier with the frozen
hyperparameters of the archived exclusion evidence (n_estimators=100,
max_depth=3, learning_rate=0.1, subsample=0.8, random_state=42,
eval_metric='logloss'). Pooled OOF AUROC per shuffle;
p = (#{null >= obs} + 1) / (n_valid + 1).

BH q is reported under both the six-cell within-cohort family (primary
matrix) and a seven-cell family that adds this cell, so the manuscript can
state the family consequence explicitly.

Resumable: progress checkpointed to
results/benchmark/v33/nc_cells/Gide_2019_cBio__XGBoost.ckpt.json every 50
shuffles.
"""
import json
import multiprocessing as mp
import os
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

from rerun_v33 import SEED, load_cohort, cv_predict  # noqa: E402
from sklearn.feature_selection import mutual_info_classif  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

OUTD = REPO / "results" / "benchmark" / "v33" / "nc_cells"
DEST = OUTD / "Gide_2019_cBio__XGBoost.json"
CKPT = OUTD / "Gide_2019_cBio__XGBoost.ckpt.json"
N_PERM = 1000
N_SELECT = 300

_W = {"X": None, "y": None, "folds": None}


def _pool_init(xpath, ypath, folds):
    _W["X"] = np.load(xpath, mmap_mode="r")
    _W["y"] = np.load(ypath, mmap_mode="r")
    _W["folds"] = folds


def _xgb_fit_fold(X_tr, y_tr, X_te):
    from xgboost import XGBClassifier
    mi = mutual_info_classif(X_tr, y_tr, random_state=42)
    top = np.argsort(mi)[-min(N_SELECT, X_tr.shape[1]):]
    model = XGBClassifier(n_estimators=100, max_depth=3, learning_rate=0.1,
                          subsample=0.8, random_state=42,
                          eval_metric="logloss")
    model.fit(X_tr[:, top], y_tr)
    return model.predict_proba(X_te[:, top])[:, 1]


def _perm_worker(seed):
    X, y, folds = _W["X"], _W["y"], _W["folds"]
    rng = np.random.RandomState(seed)
    y_perm = rng.permutation(y)
    yt, yp = [], []
    for f in folds:
        tr, te = np.array(f["train"]), np.array(f["test"])
        try:
            yp.extend(_xgb_fit_fold(X[tr], y_perm[tr], X[te]))
        except ValueError:
            return None  # degenerate training fold
        yt.extend(y_perm[te])
    try:
        return float(roc_auc_score(np.array(yt), np.array(yp)))
    except ValueError:
        return None


def bh_q(pvals):
    """Benjamini-Hochberg q-values, standard monotone enforcement."""
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        r = m - rank
        val = min(prev, p[idx] * m / r)
        q[idx] = val
        prev = val
    return q


def main():
    adata, folds, y = load_cohort("Gide_2019_cBio")
    print(f"Gide_2019_cBio: n={len(y)} pos={int(y.sum())} "
          f"genes={adata.n_vars}", flush=True)

    # observed AUROC via the same cv_predict path as the exclusion evidence
    yt, yp = cv_predict(adata, folds, "XGBoost")
    obs = float(roc_auc_score(yt, yp))
    archived = 0.7182  # xgb_exclusion_evidence_v33.json
    print(f"obs AUROC (this run) = {obs:.4f}  (archived evidence {archived})",
          flush=True)

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    tmpd = tempfile.mkdtemp(prefix="xgb_perm_")
    xpath, ypath = os.path.join(tmpd, "X.npy"), os.path.join(tmpd, "y.npy")
    np.save(xpath, np.ascontiguousarray(X, dtype=np.float32))
    np.save(ypath, y.astype(np.int8))

    done, null = 0, []
    if CKPT.exists():
        ck = json.load(open(CKPT))
        done, null = ck["done"], ck["null"]
        print(f"resuming from {done} shuffles", flush=True)

    n_jobs = min(6, max(2, (os.cpu_count() or 4) - 1))
    t0 = time.time()
    while done < N_PERM:
        batch = min(100, N_PERM - done)
        pool = mp.Pool(n_jobs, initializer=_pool_init,
                       initargs=(xpath, ypath, folds))
        try:
            res = pool.map(_perm_worker,
                           range(SEED + done, SEED + done + batch),
                           chunksize=4)
        finally:
            pool.close()
            pool.join()
        null.extend([r for r in res if r is not None])
        done += batch
        json.dump({"done": done, "null": null}, open(CKPT, "w"))
        el = time.time() - t0
        print(f"  {done}/{N_PERM} shuffles ({el:.0f}s, "
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

    # BH q under the 6-cell within-cohort family and a 7-cell family
    fam6 = {"GEP": 2e-05, "TIDE": 0.00014, "PD_L1": 2e-05,
            "IMPRES": 0.03168, "ElasticNet_MI": 0.0232,
            "ElasticNet_Var": 0.02}
    q7 = dict(zip(list(fam6) + ["XGBoost"],
                  bh_q(list(fam6.values()) + [p]).round(4)))
    q6_other = dict(zip(list(fam6),
                        bh_q(list(fam6.values())).round(4)))

    rec = {
        "cohort": "Gide_2019_cBio", "method": "XGBoost",
        "endpoint": "RECIST (n=73, 40 responders)",
        "obs_auroc_this_run": round(obs, 4),
        "archived_exclusion_evidence_auroc": archived,
        "obs_deviation_vs_archived": round(abs(obs - archived), 4),
        "params_frozen": {"n_estimators": 100, "max_depth": 3,
                          "learning_rate": 0.1, "subsample": 0.8,
                          "random_state": 42, "feature_selection":
                          "full-dim MI top-300 per fold (XGBoostWrapper)"},
        "n_perm_requested": N_PERM, "n_valid": len(null),
        "n_ge_obs": ge, "p": round(p, 5),
        "bh_q_7cell_family": q7,
        "bh_q_6cell_family_others_unchanged_check": q6_other,
        "protocol": "full-pipeline label shuffles; per-fold MI top-300 "
                    "selection on shuffled labels; frozen XGBClassifier "
                    "hyperparameters; seed 42; multiprocessing pool",
        "xgboost_version": __import__("xgboost").__version__,
    }
    with open(DEST, "w") as fh:
        json.dump(rec, fh, indent=1)
    try:
        os.remove(CKPT)
    except OSError:
        pass
    print(f"p = {p:.5f}  q(7-cell) = {q7['XGBoost']}", flush=True)
    print("WROTE", DEST, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
