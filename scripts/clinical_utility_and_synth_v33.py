"""v33 micro-recompute: (1) Hugo EN-Var clinical utility metrics under the
v33 protocol (the archived S9 block and manuscript DCA/NRI/IDI numbers were
v28-era); (2) synthetic-data framework validation for Note S1.

Output: results/benchmark/v33/clinical_utility_synth_v33.json
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

from rerun_v33 import SEED, load_cohort, nested_cv_predict, tune_primary  # noqa: E402
from metrics import (decision_curve_analysis,  # noqa: E402
                     integrated_discrimination_improvement,
                     net_reclassification_improvement)
from sklearn.metrics import (average_precision_score,  # noqa: E402
                             brier_score_loss, confusion_matrix,
                             f1_score, roc_auc_score)

OUT = REPO / "results" / "benchmark" / "v33" / "clinical_utility_synth_v33.json"


def main():
    out = {"meta": {"seed": SEED,
                    "protocol": "v33 leakage-free nested CV, frozen "
                                "hyperparameters from tune_primary"}}

    # ---- 1. Hugo EN-Var clinical utility ----
    adata, folds, y = load_cohort("Hugo_2016")
    params = tune_primary(adata, folds, "ElasticNet_Var")
    yt, yp, _ = nested_cv_predict(adata, folds, "ElasticNet_Var", params=params)
    baseline = np.full_like(yp, float(np.mean(yt)))  # random/prevalence baseline
    tn, fp, fn, tp = confusion_matrix(yt, yp >= 0.5).ravel()
    dca = decision_curve_analysis(yt, yp, thresholds=[0.10, 0.20, 0.50])
    nri = net_reclassification_improvement(yt, baseline, yp)
    idi = integrated_discrimination_improvement(yt, baseline, yp)
    out["hugo_envar_clinical_utility"] = {
        "params": params,
        "auroc": round(float(roc_auc_score(yt, yp)), 3),
        "auprc": round(float(average_precision_score(yt, yp)), 3),
        "f1": round(float(f1_score(yt, yp >= 0.5, zero_division=0)), 3),
        "sensitivity": round(float(tp / (tp + fn)) if (tp + fn) else None, 3),
        "specificity": round(float(tn / (tn + fp)) if (tn + fp) else None, 3),
        "brier": round(float(brier_score_loss(yt, yp)), 3),
        "confusion_tp_fp_tn_fn": [int(tp), int(fp), int(tn), int(fn)],
        "dca": {f"net_benefit_pt={t}": round(float(nb), 3) for t, nb in
                zip(dca["thresholds"], dca["net_benefit"])},
        "dca_treat_all": {f"net_benefit_pt={t}": round(float(nb), 3) for t, nb in
                          zip(dca["thresholds"], dca["net_benefit_all"])},
        "nri_vs_prevalence_baseline": {
            "nri_total": round(float(nri["nri_total"]), 3),
            "nri_pval": round(float(nri["nri_pval"]), 4)},
        "idi_vs_prevalence_baseline": {
            "idi_total": round(float(idi["idi_total"]), 4)},
    }
    print(json.dumps(out["hugo_envar_clinical_utility"], indent=1), flush=True)

    # ---- 2. Synthetic framework validation (Note S1) ----
    import scanpy as sc
    synth = sc.read_h5ad(REPO / "data" / "cohorts" / "SYNTH_TEST" /
                         "processed" / "SYNTH_TEST_processed.h5ad")
    splits = json.load(open(REPO / "data" / "cohorts" / "SYNTH_TEST" /
                            "processed" / "SYNTH_TEST_splits.json"))["folds"]
    y = synth.obs["response"].values.astype(int)
    res = {"n": int(synth.n_obs), "n_genes": int(synth.n_vars)}
    from base_method import get_method
    yt, yp = [], []
    m = get_method("GEP")
    for f in splits:
        preds, _ = m.fit_predict(synth, f["train"], f["test"])
        yt.extend(y[f["test"]])
        yp.extend(preds)
    res["gep_auroc"] = round(float(roc_auc_score(yt, yp)), 3)
    p = tune_primary(synth, splits, "ElasticNet")
    yt, yp, _ = nested_cv_predict(synth, splits, "ElasticNet", params=p)
    res["en_mi_auroc"] = round(float(roc_auc_score(yt, yp)), 3)
    res["en_mi_params"] = p
    out["synthetic_validation"] = res
    print(json.dumps(res, indent=1), flush=True)

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"WROTE {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
