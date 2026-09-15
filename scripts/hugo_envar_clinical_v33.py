"""Hugo 2016 ElasticNet-Var clinical-utility metrics under the v33 protocol.

The manuscript's DCA/NRI/IDI sentence for ElasticNet_Var on Hugo 2016 carried
v28 pipeline values (NRI +0.985, sensitivity 0.385, specificity 0.600).
This script recomputes them under the v33 leakage-free protocol with frozen
hyperparameters from tune_primary, using workflows/evaluation/metrics.py:

  - sensitivity/specificity/F1/Brier at the 0.5 threshold
  - DCA net benefit at decision thresholds 0.10 / 0.20 / 0.50 (+ treat-all)
  - NRI and IDI vs. the random baseline (constant 0.5), the same
    convention as the archived v28 table

Output: results/benchmark/v33/hugo_envar_clinical_v33.json
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

from rerun_v33 import load_cohort, nested_cv_predict, tune_primary  # noqa: E402
from metrics import (classification_metrics, decision_curve_analysis,  # noqa: E402
                     integrated_discrimination_improvement,
                     net_reclassification_improvement)

OUT = REPO / "results" / "benchmark" / "v33" / "hugo_envar_clinical_v33.json"


def main():
    adata, folds, y = load_cohort("Hugo_2016")
    params = tune_primary(adata, folds, "ElasticNet_Var")
    yt, yp, _ = nested_cv_predict(adata, folds, "ElasticNet_Var", params=params)
    yt = np.asarray(yt).astype(int)
    yp = np.asarray(yp, dtype=float)

    cm = classification_metrics(yt, yp)
    dca = decision_curve_analysis(yt, yp, thresholds=[0.10, 0.20, 0.50])
    base = np.full_like(yp, 0.5)
    nri = net_reclassification_improvement(yt, yp, base)
    idi = integrated_discrimination_improvement(yt, yp, base)

    out = {"meta": {"params": params,
                    "protocol": "v33 leakage-free nested CV, frozen params; "
                                "NRI/IDI vs constant-0.5 random baseline"},
           "classification": cm,
           "dca": {f"tb_{t:.2f}": {"model": float(nb), "treat_all": float(nba)}
                   for t, nb, nba in zip(dca["thresholds"], dca["net_benefit"],
                                         dca["net_benefit_all"])},
           "nri_vs_random": {k: float(v) for k, v in nri.items()},
           "idi_vs_random": {k: float(v) for k, v in idi.items()}}
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
