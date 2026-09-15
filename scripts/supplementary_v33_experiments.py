"""v33 supplementary-material recompute experiments.

Recomputes, under the v33 leakage-free protocol, every supplementary-table
quantity that is not already stored in results/benchmark/v33/*.json:

  S5  Cohen's kappa concordance (Hugo, OOF predictions, median split)
  S6  IMPRES leave-one-pair-out (canonical 15 pairs, Hugo)
  S7  preprocessing sensitivity (Hugo EN-MI; quantile row recomputable,
      TPM not derivable from the archived log2(FPKM+1) matrix - recorded)
  S8  CV sensitivity (Hugo EN-MI; k = 3/10/LOO; k = 5 from benchmark JSON)
  S10 immune deconvolution signature AUROCs (all 6 cohort-endpoints)
  S12 seed stability (Hugo EN-MI vs EN-Var, 10 stratified seeds)
  S13 learning curve (Hugo EN-Var subsampling, 10 seeds per size)
  Figs S2/S3  full-data top-500 MI vs variance gene sets (overlap 27%)
  Fig S4      Hugo KM step curves (v33 event convention, see
              scripts/hugo_survival_v33.py)
  CDS replication on the 4 archived endpoints + first stored Gide value
  (cds v1.1.0, declared endpoint genes GZMA/PRF1, same convention as the
  archived scores)

Frozen hyperparameters come from tune_primary() on the real labels (the same
convention as the v33 benchmark). Output:
  results/benchmark/v33/supplementary_v33_recompute.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

from rerun_v33 import (COHORTS, SEED, cv_predict, load_cohort,  # noqa: E402
                       nested_cv_predict, tune_primary)
from base_method import IMPRES_Wrapper  # noqa: E402
from metrics import evaluate_all  # noqa: E402
from sklearn.metrics import cohen_kappa_score, roc_auc_score  # noqa: E402
from sklearn.feature_selection import mutual_info_classif  # noqa: E402

# v33 fix (same as final_benchmarks_v33.py): HGNC-symbol matrices for Riaz.
COHORTS["Riaz_2017_cytolytic"]["h5ad"] = (
    "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad")
COHORTS["Riaz_2017_RECIST_v33"]["h5ad"] = (
    "data/cohorts/Riaz_2017_RECIST/processed/"
    "Riaz_2017_RECIST_v33_HGNC.h5ad")

OUT = REPO / "results" / "benchmark" / "v33" / \
    "supplementary_v33_recompute.json"
BM = json.load(open(REPO / "results" / "benchmark" / "v33" /
                    "benchmark_v33.json"))["cohorts"]

ALL_COHORTS = ["Hugo_2016", "Nathanson_2017", "Gide_2019_cBio",
               "Jung_2019_DCB", "Riaz_2017_cytolytic",
               "Riaz_2017_RECIST_v33"]

SIGNATURES = {
    "CD8 T cells": ["CD8A", "CD8B"],
    "CD4 T cells": ["CD4"],
    "B cells": ["CD19", "CD79A", "MS4A1"],
    "NK cells": ["NKG7", "KLRD1", "KLRF1"],
    "M1 Macrophages": ["NOS2", "IL12A", "TNF"],
    "M2 Macrophages": ["CD163", "MRC1", "MSR1"],
    "Tregs": ["FOXP3", "IL2RA"],
    "Dendritic cells": ["CD1C", "CLEC10A", "ITGAX"],
    "Neutrophils": ["ELANE", "MPO", "CEACAM8"],
    "Cytolytic score": ["GZMA", "PRF1", "GNLY"],
    "IFNG response": ["IFNG", "CXCL10", "CXCL9", "IDO1", "STAT1"],
    "TCR signaling": ["CD3D", "CD3E", "CD3G", "TRAC"],
    "Checkpoint": ["PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"],
    "TLS": ["CXCL13", "CCL19", "CCL21"],
}


def dense(adata):
    return (adata.X.toarray() if hasattr(adata.X, "toarray")
            else np.asarray(adata.X))


def stratified_folds(y, k, seed):
    rng = np.random.RandomState(seed)
    fold_of = np.full(len(y), -1)
    for cls in np.unique(y):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        for pos, s in enumerate(idx):
            fold_of[s] = pos % k
    return [{"train": list(np.where(fold_of != f)[0]),
             "test": list(np.where(fold_of == f)[0])} for f in range(k)]


def oof_auroc(adata, folds, method, params):
    yt, yp, _ = nested_cv_predict(adata, folds, method, params=params)
    return float(roc_auc_score(yt, yp))


# ------------------------------------------------------------------
def impres_lopo(adata, y):
    """Full canonical IMPRES score + leave-one-pair-out AUROCs."""
    X = dense(adata)
    genes = list(adata.var_names)
    pairs = [p for p in IMPRES_Wrapper.PAIRS
             if p[0] in genes and p[1] in genes]

    def score(skip=None):
        s = np.zeros(adata.n_obs)
        m = 0
        for pr in pairs:
            if pr == skip:
                continue
            s += (X[:, genes.index(pr[0])] > X[:, genes.index(pr[1])]
                  ).astype(float)
            m += 1
        return s * (15.0 / max(m, 1)) / 15.0

    full = float(roc_auc_score(y, score()))
    rows = []
    for pr in pairs:
        a = float(roc_auc_score(y, score(skip=pr)))
        rows.append({"dropped_pair": f"{pr[0]}-{pr[1]}",
                     "auroc": round(a, 3),
                     "delta_from_full": round(a - full, 3)})
    rows.sort(key=lambda r: -r["auroc"])
    return {"n_pairs": len(pairs), "full_auroc": round(full, 3),
            "rows": rows}


# ------------------------------------------------------------------
def kappa_table(adata, folds, y, params_mi):
    """OOF median-binarized calls for 5 methods -> Cohen's kappa."""
    preds = {}
    for m in ["IMPRES", "GEP", "TIDE", "PD_L1_IHC"]:
        _, yp = cv_predict(adata, folds, m)
        preds[m] = np.asarray(yp, dtype=float)
    _, yp, _ = nested_cv_predict(adata, folds, "ElasticNet",
                                 params=params_mi)
    preds["ElasticNet_MI"] = np.asarray(yp, dtype=float)
    calls = {m: (v > np.median(v)).astype(int) for m, v in preds.items()}
    names = list(calls)
    mat = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            mat[f"{a}|{b}"] = {
                "agreement": round(float(np.mean(calls[a] == calls[b])), 3),
                "kappa": round(float(cohen_kappa_score(calls[a], calls[b])),
                               3)}
    return {"convention": "OOF score binarized at its median",
            "pairs": mat}


# ------------------------------------------------------------------
def cv_sensitivity(adata, folds, y, params):
    out = {}
    for k in (3, 10):
        ff = stratified_folds(y, k, SEED)
        yt, yp, _ = nested_cv_predict(adata, ff, "ElasticNet", params=params)
        out[f"k={k}"] = {"auroc": round(float(roc_auc_score(yt, yp)), 3)}
    # LOO
    loo = [{"train": [i for i in range(len(y)) if i != t], "test": [t]}
           for t in range(len(y))]
    yt, yp, _ = nested_cv_predict(adata, loo, "ElasticNet", params=params)
    out["LOO"] = {"auroc": round(float(roc_auc_score(yt, yp)), 3)}
    # k=5 default row comes from the v33 benchmark JSON (same protocol)
    out["k=5"] = {"auroc": BM["Hugo_2016"]["ElasticNet"]["auroc"],
                  "source": "benchmark_v33.json"}
    return out


# ------------------------------------------------------------------
def quantile_normalize(X):
    """Per-sample quantile normalization to the mean ordered profile."""
    sx = np.sort(X, axis=1)
    ref = sx.mean(axis=0)
    out = np.empty_like(X)
    for i in range(X.shape[0]):
        r = np.argsort(np.argsort(X[i]))
        out[i] = ref[r]
    return out


def preprocessing(adata, folds, y, params):
    from scanpy import read_h5ad
    aq = adata.copy()
    aq.X = quantile_normalize(dense(adata))
    yt, yp, _ = nested_cv_predict(aq, folds, "ElasticNet", params=params)
    return {
        "fpkm_log2": {"auroc": BM["Hugo_2016"]["ElasticNet"]["auroc"],
                      "source": "benchmark_v33.json"},
        "quantile": {"auroc": round(float(roc_auc_score(yt, yp)), 3)},
        "tpm_log2": None,  # not derivable: archived matrix is log2(FPKM+1),
        # no raw counts or gene lengths are stored in the repository
    }


# ------------------------------------------------------------------
def seed_stability(adata, y, params_mi, params_var, n_seeds=10):
    mi_a, var_a = [], []
    for s in range(n_seeds):
        ff = stratified_folds(y, 5, SEED + s)
        mi_a.append(oof_auroc(adata, ff, "ElasticNet", params_mi))
        var_a.append(oof_auroc(adata, ff, "ElasticNet_Var", params_var))
    d = np.array(mi_a) - np.array(var_a)
    return {
        "mi_mean": round(float(np.mean(mi_a)), 3),
        "mi_sd": round(float(np.std(mi_a)), 3),
        "var_mean": round(float(np.mean(var_a)), 3),
        "var_sd": round(float(np.std(var_a)), 3),
        "delta_mean": round(float(np.mean(d)), 3),
        "delta_sd": round(float(np.std(d)), 3),
        "mi_wins": int((d > 0).sum()),
        "n_seeds": n_seeds,
        "mi_aurocs": [round(a, 3) for a in mi_a],
        "var_aurocs": [round(a, 3) for a in var_a],
    }


# ------------------------------------------------------------------
def learning_curve(adata, y, params, sizes=(8, 14, 21, 28), n_seeds=10):
    rows = []
    for n in sizes:
        aucs, fold_sds = [], []
        for s in range(n_seeds):
            rng = np.random.RandomState(SEED + s)
            idx = np.arange(len(y))
            keep = []
            for cls in np.unique(y):
                ci = idx[y == cls]
                rng.shuffle(ci)
                keep.extend(ci[:max(2, int(round(n * len(ci) / len(y))))])
            keep = np.array(sorted(keep))
            sub = adata[keep].copy()
            # manuscript convention: k = min(5, n_pos, n_neg) — avoids empty
            # test folds at small subsample sizes
            ys = y[keep]
            k_sub = min(5, int((ys == 1).sum()), int((ys == 0).sum()))
            ff = stratified_folds(ys, k_sub, SEED + s)
            yt, yp, _ = nested_cv_predict(sub, ff, "ElasticNet_Var",
                                          params=params)
            try:
                aucs.append(float(roc_auc_score(yt, yp)))
            except ValueError:
                pass
            fa = []
            for f in ff:
                if len(np.unique(yt[f["test"]])) == 2:
                    fa.append(roc_auc_score(yt[f["test"]], yp[f["test"]]))
            if fa:
                fold_sds.append(float(np.std(fa)))
        rows.append({"n": int(len(keep)), "target_n": n,
                     "pct_full": round(100.0 * len(keep) / len(y)),
                     "auroc_mean": round(float(np.mean(aucs)), 3),
                     "auroc_sd_seeds": round(float(np.std(aucs)), 3),
                     "mean_crossfold_sd": round(float(np.mean(fold_sds)), 3)})
    return rows


# ------------------------------------------------------------------
def deconvolution():
    out = {}
    for cname in ALL_COHORTS:
        adata, _, y = load_cohort(cname)
        X = dense(adata)
        genes = list(adata.var_names)
        res = {}
        for sig, gl in SIGNATURES.items():
            found = [g for g in gl if g in genes]
            if not found:
                res[sig] = {"auroc": None, "n_genes_found": 0, "n_genes": len(gl)}
                continue
            score = X[:, [genes.index(g) for g in found]].mean(axis=1)
            try:
                auc = round(float(roc_auc_score(y, score)), 3)
            except ValueError:
                auc = None
            res[sig] = {"auroc": auc, "n_genes_found": len(found),
                        "n_genes": len(gl)}
        out[cname] = res
        print(f"  deconv done: {cname}", flush=True)
    return out


# ------------------------------------------------------------------
def cds_endpoints():
    sys.path.insert(0, str(REPO / "cds_tool"))
    from cds import circularity_detection_score

    jobs = [("Hugo_2016", "Hugo 2016 RECIST"),
            ("Nathanson_2017", "Nathanson 2017 RECIST"),
            ("Gide_2019_cBio", "Gide 2019 RECIST"),
            ("Riaz_2017_cytolytic", "Riaz 2017 cytolytic"),
            ("Riaz_2017_RECIST_v33", "Riaz 2017 RECIST (v33)")]
    out = {}
    for cname, label in jobs:
        adata, _, y = load_cohort(cname)
        X = dense(adata)
        r = circularity_detection_score(X, y, ["GZMA", "PRF1"],
                                        list(adata.var_names))
        r = json.loads(json.dumps(r, default=float))
        out[label] = r
        print(f"  CDS {label}: {json.dumps(r)[:120]}", flush=True)
    return out


# ------------------------------------------------------------------
def km_hugo():
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    from scanpy import read_h5ad

    adata = read_h5ad(REPO / "data" / "cohorts" / "Hugo_2016" / "processed" /
                      "Hugo_2016_processed.h5ad")
    os_days = adata.obs["OS_days"].values.astype(float)
    y = adata.obs["response"].values.astype(int)
    keep = os_days > 0
    dur, ev = os_days[keep], (os_days[keep] <= 400).astype(int)

    def curve(groups):
        res = {}
        for name, mask in groups.items():
            kmf = KaplanMeierFitter()
            kmf.fit(dur[mask], ev[mask])
            res[name] = {"timeline": [round(float(t), 1) for t in kmf.survival_function_.index],
                         "survival": [round(float(v), 4) for v in
                                      kmf.survival_function_.iloc[:, 0]]}
        return res

    X = dense(adata)
    genes = list(adata.var_names)
    ful = np.zeros(adata.n_obs)
    m = 0
    for ga, gb in IMPRES_Wrapper.PAIRS:
        if ga in genes and gb in genes:
            ful += (X[:, genes.index(ga)] > X[:, genes.index(gb)]).astype(float)
            m += 1
    impres = ful * (15.0 / max(m, 1)) / 15.0
    hi = impres[keep] > np.median(impres[keep])
    curves = curve({"IMPRES_high": hi, "IMPRES_low": ~hi,
                    "RECIST_responder": y[keep] == 1,
                    "RECIST_nonresponder": y[keep] != 1})
    lr_i = logrank_test(dur[hi], dur[~hi], ev[hi], ev[~hi])
    resp = y[keep] == 1
    lr_r = logrank_test(dur[resp], dur[~resp], ev[resp], ev[~resp])
    return {"event_convention": "0 < OS_days <= 400 (see hugo_survival_v33.py)",
            "n": int(keep.sum()), "n_events": int(ev.sum()),
            "curves": curves,
            "logrank_impres_p": float(lr_i.p_value),
            "logrank_recist_p": float(lr_r.p_value)}


# ------------------------------------------------------------------
def mi_var_top500(adata):
    X = dense(adata)
    y = adata.obs["response"].values.astype(int)
    genes = list(adata.var_names)
    var = X.var(axis=0)
    mi = mutual_info_classif(X, y, random_state=SEED)
    top_var = set(np.argsort(var)[-500:])
    top_mi = set(np.argsort(mi)[-500:])
    inter = sorted(top_var & top_mi)
    return {
        "overlap_count": len(inter),
        "overlap_pct": round(100.0 * len(inter) / 500, 1),
        "top500_var_genes": [genes[i] for i in sorted(top_var)],
        "top500_mi_genes": [genes[i] for i in sorted(top_mi)],
        "var_scores_all": [round(float(v), 5) for v in var],
        "mi_scores_all": [round(float(v), 5) for v in mi],
        "genes_all": genes,
        "note": "descriptive full-data computation (v33 narrative: "
                "top-500, full data); not a predictive estimate",
    }


# ------------------------------------------------------------------
def main():
    print("=== loading Hugo ===", flush=True)
    adata, folds, y = load_cohort("Hugo_2016")
    params_mi = tune_primary(adata, folds, "ElasticNet")
    params_var = tune_primary(adata, folds, "ElasticNet_Var")
    print(f"frozen params: MI={params_mi} Var={params_var}", flush=True)

    out = {"meta": {
        "seed": SEED,
        "params_mi": params_mi, "params_var": params_var,
        "protocol": "v33 leakage-free nested CV (per-fold feature "
                    "selection); frozen hyperparameters from tune_primary",
    }}

    print("=== S6 IMPRES LOPO ===", flush=True)
    out["impres_lopo_hugo"] = impres_lopo(adata, y)
    print(f"  full IMPRES = {out['impres_lopo_hugo']['full_auroc']}", flush=True)

    print("=== S5 kappa ===", flush=True)
    out["kappa_hugo"] = kappa_table(adata, folds, y, params_mi)

    print("=== S8 CV sensitivity ===", flush=True)
    out["cv_sensitivity_hugo"] = cv_sensitivity(adata, folds, y, params_mi)
    print(f"  {out['cv_sensitivity_hugo']}", flush=True)

    print("=== S7 preprocessing ===", flush=True)
    out["preprocessing_hugo"] = preprocessing(adata, folds, y, params_mi)
    print(f"  {out['preprocessing_hugo']}", flush=True)

    print("=== S12 seed stability ===", flush=True)
    out["seed_stability_hugo"] = seed_stability(adata, y, params_mi,
                                                params_var)
    print(f"  {out['seed_stability_hugo']}", flush=True)

    print("=== S13 learning curve ===", flush=True)
    out["learning_curve_hugo"] = learning_curve(adata, y, params_var)
    print(f"  {out['learning_curve_hugo']}", flush=True)

    print("=== S10 deconvolution ===", flush=True)
    out["deconvolution"] = deconvolution()

    print("=== CDS endpoints ===", flush=True)
    out["cds_endpoints"] = cds_endpoints()

    print("=== Fig S4 KM curves ===", flush=True)
    out["km_hugo"] = km_hugo()

    print("=== Figs S2/S3 MI vs Var ===", flush=True)
    out["mi_var_top500_hugo"] = mi_var_top500(adata)

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nWROTE {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
