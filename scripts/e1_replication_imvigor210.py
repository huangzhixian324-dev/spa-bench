"""E1 (executed, second replication): response-definition effect in a
NON-melanoma cancer type — IMvigor210 (metastatic urothelial carcinoma,
anti-PD-L1, n = 297 response-annotated samples).

The pre-declared NSCLC target (GSE274975) is unreachable from this machine
(all NCBI endpoints blocked; no iAtlas NSCLC ICI cohort exists), and Liu
2019 (the first executed replication) is melanoma. IMvigor210 — bladder
carcinoma — upgrades the replication from two-cohort to cross-cancer and
its data is already local (data/external/blca_iatlas_imvigor210_2017/,
cBioPortal iAtlas harmonization, z-score expression).

Design (identical to Table S19; seed 42):
  endpoint A (molecular surrogate): cytolytic-above-median (GZMA/PRF1 mean)
  endpoint B (clinical): RECIST, CR/PR = responder (68/297)
  leakage-free nested CV, frozen hyperparameters from tune_primary,
  6 methods; exact 50,000-shuffle permutations for fixed scorers;
  CDS v1.1.0 on both endpoints + 10-replicate random-label null.

Output: results/benchmark/v33/e1_imvigor210_replication.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))
sys.path.insert(0, str(REPO / "cds_tool"))

from rerun_v33 import SEED, nested_cv_predict, tune_primary, cv_predict  # noqa: E402
from metrics import evaluate_all  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

EXT = REPO / "data" / "external" / "blca_iatlas_imvigor210_2017"
OUT = REPO / "results" / "benchmark" / "v33" / \
    "e1_imvigor210_replication.json"
RESPONDERS = ["COMPLETE RESPONSE", "PARTIAL RESPONSE"]


def auroc_fixed(score, y):
    r = rankdata(score)
    n1 = int(y.sum())
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def stratified_folds(y, k, seed):
    rng = np.random.RandomState(seed)
    fold_of = np.full(len(y), -1)
    for cls in np.unique(y):
        idx = np.where(y == cls)[0]
        rng.shuffle(idx)
        for pos, i in enumerate(idx):
            fold_of[i] = pos % k
    return [{"train": list(np.where(fold_of != f)[0]),
             "test": list(np.where(fold_of == f)[0])} for f in range(k)]


def perm_fixed(score, y, n=50000, seed=SEED):
    rng = np.random.RandomState(seed)
    obs = auroc_fixed(score, y)
    ge = 0
    for _ in range(n):
        if auroc_fixed(score, rng.permutation(y)) >= obs:
            ge += 1
    return obs, (ge + 1) / (n + 1)


def main():
    import scanpy as sc
    mat = pd.read_csv(EXT / "expression.tsv.gz", index_col=0)
    mat = mat.dropna(axis=0, how="any")  # complete-case genes (mirrors E5; the restored study-export z-score matrix contains unmeasured-gene rows)
    clin = pd.read_csv(EXT / "clinical.tsv", sep="\t", index_col=0)
    clin.index = clin.index.astype(str)
    resp = clin.get("RESPONSE")
    m = mat.T.join(resp.rename("response"), how="inner")
    m = m[m["response"].notna()]
    y_recist = m["response"].astype(str).str.strip().str.upper() \
        .isin(RESPONDERS).astype(int).values
    genes = list(mat.index)
    gi = {g: i for i, g in enumerate(genes)}
    cyt_genes = [g for g in ("GZMA", "PRF1") if g in gi]
    X = np.asarray(m.drop(columns=["response"]).values, dtype=float)
    cyt_score = X[:, [gi[g] for g in cyt_genes]].mean(axis=1)
    y_cyt = (cyt_score > np.median(cyt_score)).astype(int)

    adata = sc.AnnData(X=X)
    adata.var_names = genes
    print(f"cohort: n={len(y_recist)} RECIST responders={int(y_recist.sum())} "
          f"cytolytic-high={int(y_cyt.sum())}", flush=True)

    out = {"meta": {
        "cohort": "IMvigor210 (blca_iatlas_imvigor210_2017, cBioPortal "
                  "iAtlas)",
        "n": int(len(y_recist)),
        "endpoints": {
            "cytolytic_surrogate": {
                "definition": "GZMA/PRF1 mean, median split (as Riaz)",
                "n_high": int(y_cyt.sum()), "declared_genes": cyt_genes},
            "RECIST": {"definition": "CR/PR vs SD/PD",
                       "n_responders": int(y_recist.sum())}},
        "seed": SEED,
        "expression": "cBioPortal z-scores (affine-invariant for the "
                      "benchmark metrics)",
        "note": "cross-cancer (bladder carcinoma) replication; non-melanoma",
    }}
    bench_methods = ["GEP", "TIDE", "IMPRES", "PD_L1_IHC", "ElasticNet",
                     "ElasticNet_Var"]
    label = {"IMPRES": "IMPRES", "GEP": "GEP", "TIDE": "TIDE",
             "PD_L1_IHC": "PD-L1 (CD274)", "ElasticNet": "ElasticNet (MI)",
             "ElasticNet_Var": "ElasticNet (Var)"}

    for ename, y in (("cytolytic_surrogate", y_cyt), ("RECIST", y_recist)):
        adata.obs["response"] = y
        k = min(5, int((y == 1).sum()), int((y == 0).sum()))
        folds = stratified_folds(y, k, SEED)
        res = {"k": k}
        for meth in bench_methods:
            res[meth] = {}
            if meth in ("ElasticNet", "ElasticNet_Var"):
                params = tune_primary(adata, folds, meth)
                yt, yp, _ = nested_cv_predict(adata, folds, meth,
                                              params=params)
                res[meth]["params"] = params
            else:
                yt, yp = cv_predict(adata, folds, meth)
            r = evaluate_all(np.asarray(yt).astype(int),
                             np.asarray(yp, dtype=float))
            res[meth].update({
                "auroc": round(r["auroc"], 4),
                "auroc_ci": [round(r["auroc_ci_low"], 4),
                             round(r["auroc_ci_high"], 4)],
                "auprc": round(r["auprc"], 4)})
            print(f"  {ename} {label[meth]}: {r['auroc']:.3f} "
                  f"[{r['auroc_ci_low']:.3f}-{r['auroc_ci_high']:.3f}]",
                  flush=True)
        for meth in ("IMPRES", "GEP", "TIDE", "PD_L1_IHC"):
            yt_f, yp = cv_predict(adata, folds, meth)
            yt_f = np.asarray(yt_f).astype(int)
            obs, p = perm_fixed(np.asarray(yp, dtype=float), yt_f)
            res[meth]["perm_p"] = p
            res[meth]["n_shuffle"] = 50000
            print(f"  {ename} {label[meth]} perm p={p:.5f}", flush=True)
        out[ename] = res

    out["endpoint_switch_delta"] = {
        meth: round(out["cytolytic_surrogate"][meth]["auroc"]
                    - out["RECIST"][meth]["auroc"], 3)
        for meth in bench_methods}

    from cds import circularity_detection_score
    rng = np.random.RandomState(SEED)
    out.setdefault("cds", {})
    for ename, y in (("cytolytic_surrogate", y_cyt), ("RECIST", y_recist)):
        r = circularity_detection_score(X, y, cyt_genes, genes)
        r = json.loads(json.dumps(r, default=float))
        nulls = [float(circularity_detection_score(
            X, rng.permutation(y), cyt_genes, genes)["CDS"])
            for _ in range(10)]
        out["cds"][ename] = {"CDS": r["CDS"], "risk": r["risk_level"],
                             "components": r["components"],
                             "null_mean": round(float(np.mean(nulls)), 3),
                             "null_max": round(float(np.max(nulls)), 3)}
        print(f"  CDS {ename}: {r['CDS']:.3f} "
              f"(null mean {np.mean(nulls):.3f})", flush=True)

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print("WROTE", OUT, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
