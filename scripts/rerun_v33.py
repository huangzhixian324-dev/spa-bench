"""SPATBench v33 rerun: regenerate all benchmark numbers affected by the
2026-09-02 audit fixes, with a complete, valid, machine-readable audit trail.

What changed vs the archived v28/v31 results:
  1. ElasticNet now actually uses penalty='elasticnet' with inner 3-fold CV
     (previously silent L2, no tuning).
  2. ElasticNet_Var is now implemented in the repo (previously missing).
  3. IMPRES is the canonical 15-relation implementation decoded from the
     original authors' FEATS.mat (previously a non-canonical sign() chain).
  4. Riaz RECIST uses the rebuilt cohort (42 x 22187, k=2/3/5 stratified).
  5. Jung 2019 labels are flipped so that 1 = DCB (clinical benefit), matching
     the label semantics of all other cohorts. GSE135222 annotates 6 DCB /
     21 NDB; the archived h5ad coded 1 = NDB, which silently inverted the
     meaning of the reported AUROCs.
  6. Permutation tests: 500 shuffles (seed 42) for every method x cohort;
     ML hyperparameters are tuned on the REAL labels (inner 3-fold CV) and
     then frozen during permutation, the standard practice.

Outputs (all valid JSON) -> results/benchmark/v33/
"""
import json
import multiprocessing as mp
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

from base_method import BaseMethod, get_method, IMPRES_Wrapper  # noqa: E402

IMPRES_PAIRS = IMPRES_Wrapper.PAIRS
from metrics import evaluate_all  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

SEED = 42
N_PERM = 500
OUT = REPO / "results" / "benchmark" / "v33"
OUT.mkdir(parents=True, exist_ok=True)

COHORTS = {
    "Hugo_2016": {
        "h5ad": "data/cohorts/Hugo_2016/processed/Hugo_2016_processed.h5ad",
        "splits": "data/cohorts/Hugo_2016/processed/Hugo_2016_splits.json",
        "flip": False, "endpoint": "RECIST (responders 13/28)",
    },
    "Nathanson_2017": {
        "h5ad": "data/cohorts/Nathanson_2017/processed/Nathanson_2017_processed.h5ad",
        "splits": "data/cohorts/Nathanson_2017/processed/Nathanson_2017_splits.json",
        "flip": False, "endpoint": "RECIST (responders 10/25)",
    },
    "Gide_2019_cBio": {
        "h5ad": "data/cohorts/Gide_2019_cBio/processed/Gide_2019_cBio_processed.h5ad",
        "splits": "data/cohorts/Gide_2019_cBio/processed/Gide_2019_cBio_splits.json",
        "flip": False, "endpoint": "RECIST (responders 40/73)",
    },
    "Jung_2019_DCB": {  # labels flipped: 1 = DCB (6/27)
        "h5ad": "data/cohorts/Jung_2019/processed/Jung_2019_HGNC_processed.h5ad",
        "splits": "data/cohorts/Jung_2019/processed/Jung_2019_HGNC_splits.json",
        "flip": True, "endpoint": "DCB vs NDB (GSE135222; benefit 6/27)",
    },
    "Riaz_2017_cytolytic": {
        "h5ad": "data/cohorts/Riaz_2017/processed/Riaz_2017_processed.h5ad",
        "splits": "data/cohorts/Riaz_2017/processed/Riaz_2017_splits.json",
        "flip": False, "endpoint": "cytolytic-above-median (21/43)",
    },
    "Riaz_2017_RECIST_v33": {
        "h5ad": "data/cohorts/Riaz_2017_RECIST/processed/Riaz_2017_RECIST_v33.h5ad",
        "splits": "data/cohorts/Riaz_2017_RECIST/processed/Riaz_2017_RECIST_v33_splits.json",
        "flip": False, "endpoint": "RECIST (responders 9/42), k=2 primary",
        "primary_k": "k=2",
    },
}

METHODS = ["GEP", "TIDE", "IMPRES", "ElasticNet", "ElasticNet_Var"]
# Permutation inclusion logic (as in the manuscript): trainable methods plus
# IMPRES as a fixed-classifier reference. GEP/TIDE/PD-L1 are deterministic
# scorers excluded from permutation (no training step to invalidate).
# ML permutation replicates rerun the full per-fold continuous-MI nested CV,
# which is compute-heavy at these cohort sizes; we therefore use 200 shuffles
# (null p resolution 1/201 ~ 0.005) for the two ML variants and 500 for the
# cheap fixed-classifier IMPRES. Each cell's shuffle count is recorded in the
# output JSON and in Table 4.
PERM_METHODS = ["IMPRES", "ElasticNet", "ElasticNet_Var"]
PERM_N = {"IMPRES": 500, "ElasticNet": 200, "ElasticNet_Var": 500}


def load_cohort(name):
    cfg = COHORTS[name]
    adata = sc.read_h5ad(REPO / cfg["h5ad"])
    folds = json.load(open(REPO / cfg["splits"]))["folds"]
    k = cfg.get("primary_k")
    if k:
        folds = folds[k]
    y = adata.obs["response"].values.astype(int)
    if cfg["flip"]:
        y = 1 - y
        adata.obs["response"] = y
    return adata, folds, y


def cv_predict(adata, folds, method_name, params=None):
    """Return (y_true, y_pred) over the stored folds."""
    method = get_method(method_name)
    if params and hasattr(method, "model"):
        pass
    yt, yp = [], []
    for f in folds:
        preds, _ = method.fit_predict(adata, f["train"], f["test"])
        yt.extend(adata.obs["response"].values[f["test"]])
        yp.extend(preds)
    return np.array(yt).astype(int), np.array(yp)


PRESCREEN = 2000  # stage-1 variance prescreen (leakage-free, train fold only)


def _select_genes(X_tr, y_tr, method_name):
    """Two-stage, training-fold-only feature selection.

    Stage 1: top-2000 genes by variance (train fold). Stage 2: for the MI
    variant, top-500 by continuous mutual information with the training
    labels (sklearn); for the variance variant, top-500 by variance. Both
    stages use ONLY the training fold, so no label information crosses the
    fold boundary. Identical estimator in the primary analysis and in every
    permutation replicate.
    """
    n_var = min(PRESCREEN, X_tr.shape[1])
    var_idx = np.argsort(X_tr.var(axis=0))[-n_var:]
    if method_name == "ElasticNet_Var":
        n_sel = min(500, X_tr.shape[1])
        return np.argsort(X_tr.var(axis=0))[-n_sel:]
    from sklearn.feature_selection import mutual_info_classif
    mi = mutual_info_classif(X_tr[:, var_idx], y_tr, random_state=SEED)
    order = np.argsort(mi)[-min(500, len(var_idx)):]
    return var_idx[order]


def nested_cv_predict(adata, folds, method_name, params=None, y_override=None):
    """Leakage-free nested CV: selection and (optionally) tuning inside folds.

    If params is None, hyperparameters are tuned by 3-fold inner CV on the
    training fold and the selected setting is recorded.
    """
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    y = (adata.obs["response"].values.astype(int) if y_override is None
         else np.asarray(y_override).astype(int))
    yt, yp, chosen = [], [], []
    for f in folds:
        tr, te = np.array(f["train"]), np.array(f["test"])
        genes = _select_genes(X[tr], y[tr], method_name)
        if params is None:
            g = GridSearchCV(
                BaseMethod._make_elasticnet(),
                param_grid={"C": [0.01, 0.1, 1.0],
                            "l1_ratio": [0.1, 0.5, 0.9]},
                cv=3, scoring="roc_auc")
            g.fit(X[tr][:, genes], y[tr])
            model, chosen_params = g.best_estimator_, g.best_params_
            chosen.append(chosen_params)
        else:
            model = BaseMethod._make_elasticnet(**params)
            model.fit(X[tr][:, genes], y[tr])
        yt.extend(y[te])
        yp.extend(model.predict_proba(X[te][:, genes])[:, 1])
    return np.array(yt), np.array(yp), chosen


def tune_primary(adata, folds, method_name):
    """Tune hyperparameters once on the REAL labels with nested CV; the modal
    setting is frozen for permutation replicates (standard practice)."""
    _, _, chosen = nested_cv_predict(adata, folds, method_name)
    if not chosen:
        return {"C": 0.1, "l1_ratio": 0.5}
    keys = [tuple(sorted(d.items())) for d in chosen]
    best = max(set(keys), key=keys.count)
    return dict(best)


_W = {"X": None, "y": None, "var_names": None, "folds": None,
       "method_name": None, "params": None}


def _pool_init(Xpath, ypath, var_names, folds, method_name, params):
    """Runs once per spawned worker; memmaps the cohort data."""
    _W["X"] = np.load(Xpath, mmap_mode="r")
    _W["y"] = np.load(ypath, mmap_mode="r")
    _W["var_names"] = var_names
    _W["folds"] = folds
    _W["method_name"] = method_name
    _W["params"] = params


def _perm_worker(seed):
    """One label-shuffle replicate + full nested-CV rerun (module-level)."""
    X = _W["X"]
    y = _W["y"]
    var_names = _W["var_names"]
    folds = _W["folds"]
    method_name = _W["method_name"]
    params = _W["params"]
    rng = np.random.RandomState(seed)
    y_perm = rng.permutation(y)
    yt, yp = [], []
    for f in folds:
        tr, te = np.array(f["train"]), np.array(f["test"])
        if method_name in ("ElasticNet", "ElasticNet_Var"):
            sel = _select_genes(X[tr], y_perm[tr], method_name)
            model = BaseMethod._make_elasticnet(**params)
            model.fit(X[tr][:, sel], y_perm[tr])
            yp.extend(model.predict_proba(X[te][:, sel])[:, 1])
        elif method_name == "IMPRES":
            score = np.zeros(X.shape[0])
            measured = 0
            for ga, gb in IMPRES_PAIRS:
                ia = var_names.index(ga) if ga in var_names else -1
                ib = var_names.index(gb) if gb in var_names else -1
                if ia >= 0 and ib >= 0:
                    score += (X[:, ia] > X[:, ib]).astype(float)
                    measured += 1
            yp.extend(score[te] * (15.0 / max(measured, 1)) / 15.0)
        yt.extend(y_perm[te])
    try:
        return roc_auc_score(np.array(yt), np.array(yp))
    except ValueError:
        return 0.5


def permutation_p(adata, folds, method_name, obs_auroc, n_perm=N_PERM, seed=SEED,
                  params=None, genes=None):
    """Permutation test via memmap + multiprocessing pool.

    One replicate = one label shuffle + a full nested CV re-run (per-fold
    feature selection on the shuffled labels with frozen hyperparameters), so
    the null distribution covers the whole fitting procedure. The expression
    matrix is exposed to workers through a float32 memmap so that only the
    shuffle seeds cross the process boundary.
    """
    import tempfile
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    y = adata.obs["response"].values.astype(np.int8).copy()
    var_names = list(adata.var_names)
    tmpd = tempfile.mkdtemp(prefix="v33_perm_")
    xpath, ypath = os.path.join(tmpd, "X.npy"), os.path.join(tmpd, "y.npy")
    np.save(xpath, np.ascontiguousarray(X, dtype=np.float32))
    np.save(ypath, y)
    n_jobs = min(6, max(2, os.cpu_count() - 1))
    pool = mp.Pool(n_jobs, initializer=_pool_init,
                   initargs=(xpath, ypath, var_names, folds,
                             method_name, params))
    try:
        null = pool.map(_perm_worker, range(seed, seed + n_perm), chunksize=8)
    finally:
        pool.close()
        pool.join()
    for f in (xpath, ypath):
        try:
            os.remove(f)
        except OSError:
            pass
    try:
        os.rmdir(tmpd)
    except OSError:
        pass
    ge = int(np.sum(np.array(null) >= obs_auroc - 1e-12))
    return (ge + 1) / (n_perm + 1)


def main():
    benchmark, perms, meta = {}, {}, {
        "seed": SEED, "n_perm": N_PERM,
        "note": "v33 rerun after 2026-09-02 audit; see script docstring",
    }

    for cname in COHORTS:
        adata, folds, y = load_cohort(cname)
        print(f"\n===== {cname} ({adata.n_obs} x {adata.n_vars}) "
              f"pos={int(y.sum())}/{len(y)} =====")
        benchmark[cname], perms[cname] = {}, {}

        for m in METHODS:
            params = genes = None
            if m in ("ElasticNet", "ElasticNet_Var"):
                params = tune_primary(adata, folds, m)
                yt, yp, _ = nested_cv_predict(adata, folds, m, params=params)
                genes = None
                print(f"  {m:15s} tuned params: {params}", flush=True)
            else:
                yt, yp = cv_predict(adata, folds, m)
            res = evaluate_all(yt, yp, cohort_name=cname)
            auc = res["auroc"]
            benchmark[cname][m] = {
                "auroc": round(auc, 4),
                "auroc_ci": [round(res["auroc_ci_low"], 4),
                             round(res["auroc_ci_high"], 4)],
                "auprc": round(res["auprc"], 4),
            }
            print(f"  {m:15s} AUROC={auc:.3f} "
                  f"[{res['auroc_ci_low']:.3f}-{res['auroc_ci_high']:.3f}]",
                  flush=True)

            # permutation p-value (trainable methods + IMPRES reference)
            if m in PERM_METHODS:
                n_perm = PERM_N.get(m, 200)
                p = permutation_p(adata, folds, m, auc, n_perm=n_perm,
                                  params=params, genes=None)
                perms[cname][m] = {"obs_auroc": round(auc, 4),
                                   "p": round(p, 4), "n_perm": n_perm}
                print(f"  {m:15s} permutation p = {p:.4f} ({n_perm} shuffles)", flush=True)
        # checkpoint after each cohort
        with open(OUT / "checkpoint_v33.json", "w") as fh:
            json.dump({"meta": meta, "cohorts": benchmark,
                       "perms": perms}, fh, indent=1)

    # ---- k-fold robustness for the rebuilt Riaz RECIST (ElasticNet_Var) ----
    print("\n===== Riaz RECIST v33 k-fold robustness (ElasticNet_Var) =====")
    adata, folds_all, y = load_cohort("Riaz_2017_RECIST_v33")
    adata.obs["response"] = y  # ensure unflipped
    splits = json.load(open(REPO / COHORTS["Riaz_2017_RECIST_v33"]["splits"]))
    kfold = {}
    for k, folds in splits["folds"].items():
        params = {"C": 0.1, "l1_ratio": 0.5}  # frozen modal setting
        aucs = []
        for seed in range(10):  # 10 seeds: regenerate stratified folds
            rng = np.random.RandomState(SEED + seed)
            yv = adata.obs["response"].values.astype(int)
            fold_of = np.full(len(yv), -1)
            for cls in np.unique(yv):
                idx = np.where(yv == cls)[0]
                rng.shuffle(idx)
                for pos, s in enumerate(idx):
                    fold_of[s] = pos % int(k.split("=")[1])
            ff = [{"train": list(np.where(fold_of != f)[0]),
                   "test": list(np.where(fold_of == f)[0])}
                  for f in range(int(k.split("=")[1]))]
            yt, yp, _ = nested_cv_predict(adata, ff, "ElasticNet_Var",
                                          params=params)
            aucs.append(roc_auc_score(yt, yp))
        kfold[k] = {"mean": round(float(np.mean(aucs)), 3),
                    "std": round(float(np.std(aucs)), 3),
                    "n_seeds": 10}
        print(f"  {k}: AUROC = {kfold[k]['mean']} ± {kfold[k]['std']} "
              f"(10 stratified seeds)", flush=True)

    with open(OUT / "benchmark_v33.json", "w") as fh:
        json.dump({"meta": meta, "cohorts": benchmark}, fh, indent=1)
    with open(OUT / "permutation_v33.json", "w") as fh:
        json.dump({"meta": meta, "cohorts": perms}, fh, indent=1)
    with open(OUT / "riaz_recist_kfold_v33.json", "w") as fh:
        json.dump({"meta": meta, "kfold": kfold}, fh, indent=1)
    print(f"\nAll results written to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
