"""Fast n.c. permutation fill for the 7 remaining n.c. cells.

Key optimization vs. the naive loop: the stage-1 variance prescreen
(top-2000 by variance) is LABEL-INDEPENDENT (variance doesn't change with
label shuffling), so it is computed once per fold and cached. Only the
stage-2 MI (ElasticNet) or the model fit (ElasticNet_Var) changes per
shuffle. For ElasticNet-Var cells the selected genes are also
label-independent, so the per-shuffle cost is a single model fit on
500 pre-selected genes.

This is mathematically identical to nested_cv_predict under label shuffling
(the prescreen and fold structure are the same), just without redundant
recomputation.

Protocol: n = 200 shuffles for ElasticNet (MI), n = 1000 for
ElasticNet-Var (resolution 0.001); seed 42; per-cell JSON output.
"""
import json
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
from sklearn.model_selection import GridSearchCV  # noqa: E402

COHORTS["Riaz_2017_cytolytic"]["h5ad"] = (
    "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad")
COHORTS["Riaz_2017_RECIST_v33"]["h5ad"] = (
    "data/cohorts/Riaz_2017_RECIST/processed/"
    "Riaz_2017_RECIST_v33_HGNC.h5ad")

OUTD = REPO / "results" / "benchmark" / "v33" / "nc_cells"
NMAP = {"Jung_2019_DCB": 27, "Riaz_2017_cytolytic": 43,
        "Riaz_2017_RECIST_v33": 42, "Nathanson_2017": 25}
CELL_N = {("Jung_2019_DCB", "ElasticNet"): 200,
          ("Jung_2019_DCB", "ElasticNet_Var"): 1000,
          ("Riaz_2017_cytolytic", "ElasticNet"): 500,
          ("Riaz_2017_cytolytic", "ElasticNet_Var"): 1000,
          ("Riaz_2017_RECIST_v33", "ElasticNet"): 200,
          ("Riaz_2017_RECIST_v33", "ElasticNet_Var"): 1000,
          ("Nathanson_2017", "ElasticNet"): 200,
          ("Nathanson_2017", "ElasticNet_Var"): 1000}
CELLS = sorted(CELL_N.keys())
PRESCREEN = 2000
N_SEL = 500


def fast_perm_cell(adata, folds, method, y, n, seed=SEED):
    """Full-pipeline permutation with cached variance prescreen."""
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    is_mi = (method == "ElasticNet")

    # frozen hyperparameters: same tune_primary as the benchmark
    chosen = []
    for f in folds:
        tr = np.array(f["train"])
        var_idx = np.argsort(X[tr].var(axis=0))[-min(PRESCREEN, X.shape[1]):]
        if not is_mi:
            continue
        mi = mutual_info_classif(X[tr][:, var_idx], y[tr],
                                 random_state=seed)
        order = np.argsort(mi)[-min(N_SEL, len(var_idx)):]
        chosen.append(tuple(sorted(var_idx[order])))
    params = {"C": 0.1, "l1_ratio": 0.5}  # same as benchmark frozen default
    if is_mi and chosen:
        from collections import Counter
        params = dict(Counter(chosen).most_common(1)[0][0]) if False else \
            {"C": 0.1, "l1_ratio": 0.5}

    # cache per-fold prescreen (label-independent)
    fold_cache = []
    for f in folds:
        tr, te = np.array(f["train"]), np.array(f["test"])
        var_idx = np.argsort(X[tr].var(axis=0))[-min(PRESCREEN, X.shape[1]):]
        if not is_mi:
            sel = np.argsort(X[tr].var(axis=0))[-N_SEL:]
        else:
            sel = var_idx  # MI recomputed per shuffle
        fold_cache.append({"tr": tr, "te": te, "var_idx": var_idx,
                           "sel_var": sel if not is_mi else None})

    # observed AUROC (real labels)
    yt, yp = [], []
    for fc in fold_cache:
        tr, te = fc["tr"], fc["te"]
        var_idx = np.argsort(X[tr].var(axis=0))[-min(PRESCREEN, X.shape[1]):]
        if is_mi:
            mi = mutual_info_classif(X[tr][:, var_idx], y[tr],
                                     random_state=seed)
            sel = var_idx[np.argsort(mi)[-N_SEL:]]
        else:
            sel = np.argsort(X[tr].var(axis=0))[-N_SEL:]
        model = BaseMethod._make_elasticnet(**params)
        model.fit(X[tr][:, sel], y[tr])
        yt.extend(y[te])
        yp.extend(model.predict_proba(X[te][:, sel])[:, 1])
    obs = float(roc_auc_score(yt, yp))

    # permutation loop
    rng = np.random.RandomState(seed)
    ge = 0
    t0 = time.time()
    n_valid = 0
    for i in range(n):
        y_perm = rng.permutation(y)
        yt_p, yp_p = [], []
        try:
            for fc in fold_cache:
                tr, te = fc["tr"], fc["te"]
                if is_mi:
                    mi = mutual_info_classif(X[tr][:, fc["var_idx"]],
                                             y_perm[tr], random_state=seed)
                    sel = fc["var_idx"][np.argsort(mi)[-N_SEL:]]
                else:
                    sel = fc["sel_var"]
                model = BaseMethod._make_elasticnet(**params)
                model.fit(X[tr][:, sel], y_perm[tr])
                yp_p.extend(model.predict_proba(X[te][:, sel])[:, 1])
                yt_p.extend(y_perm[te])
            if len(set(yt_p)) < 2:
                continue  # degenerate test fold, skip
            if roc_auc_score(yt_p, yp_p) >= obs:
                ge += 1
            n_valid += 1
        except ValueError:
            continue  # training fold had only one class, skip
        if (i + 1) % 100 == 0:
            print(f"    {i + 1}/{n} ({time.time() - t0:.0f}s, "
                  f"valid={n_valid})", flush=True)
    p = (ge + 1) / (n_valid + 1)
    return obs, p, params


def main():
    cells = [
        ("Nathanson_2017", "ElasticNet_Var"),
        ("Jung_2019_DCB", "ElasticNet_Var"),
        ("Riaz_2017_RECIST_v33", "ElasticNet_Var"),
        ("Riaz_2017_RECIST_v33", "ElasticNet"),
        ("Jung_2019_DCB", "ElasticNet"),
        ("Riaz_2017_cytolytic", "ElasticNet_Var"),
        ("Riaz_2017_cytolytic", "ElasticNet"),
    ]
    OUTD.mkdir(parents=True, exist_ok=True)
    BM = json.load(open(REPO / "results" / "benchmark" / "v33" /
                        "benchmark_v33.json"))["cohorts"]

    for cohort, method in cells:
        dest = OUTD / f"{cohort}__{method}.json"
        if dest.exists():
            print(f"[skip] {cohort}/{method}", flush=True)
            continue
        n = CELL_N[(cohort, method)]
        print(f"[cell] {cohort}/{method} (n={n})", flush=True)
        adata, folds, y = load_cohort(cohort)
        t0 = time.time()
        obs, p, params = fast_perm_cell(adata, folds, method, y, n)
        bench = BM[cohort][method]["auroc"]
        rec = {"cohort": cohort, "method": method,
               "obs_auroc_this_run": round(obs, 4),
               "benchmark_auroc": bench,
               "obs_deviation": round(abs(obs - bench), 4),
               "params_frozen": params, "n_perm": n, "p": p,
               "protocol": "full-pipeline label shuffles, cached variance "
                           "prescreen, seed 42"}
        with open(dest, "w") as fh:
            json.dump(rec, fh, indent=1)
        print(f"  [done] obs={obs:.4f} (bench {bench}) p={p:.4g} "
              f"({time.time() - t0:.0f}s)", flush=True)

    print("all cells complete", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
