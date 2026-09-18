"""Dry experiment B: IMPRES-on-Hugo clinical-utility metrics (DCA/NRI/IDI)
under the identical v33 conventions as hugo_envar_clinical_v33.py.

Rationale: the manuscript's clinical-utility analysis currently uses
ElasticNet-Var on Hugo (AUROC 0.528 — a near-chance method, whose utility
failure is trivial). IMPRES is the best AUROC method on Hugo (0.795, the
head positive cell); computing the same metrics for it makes the
"AUROC does not capture clinical utility" argument non-trivial.

Output: results/benchmark/v33/hugo_impres_clinical_v33.json
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

from rerun_v33 import load_cohort, cv_predict  # noqa: E402
from metrics import (classification_metrics, decision_curve_analysis,  # noqa: E402
                     integrated_discrimination_improvement,
                     net_reclassification_improvement)

OUT = REPO / "results" / "benchmark" / "v33" / "hugo_impres_clinical_v33.json"


def main():
    adata, folds, y = load_cohort("Hugo_2016")
    yt, yp = cv_predict(adata, folds, "IMPRES")
    yt = np.asarray(yt).astype(int)
    yp = np.asarray(yp, dtype=float)

    cm = classification_metrics(yt, yp)
    dca = decision_curve_analysis(yt, yp, thresholds=[0.10, 0.20, 0.50])
    base = np.full_like(yp, 0.5)
    nri = net_reclassification_improvement(yt, yp, base)
    idi = integrated_discrimination_improvement(yt, yp, base)

    out = {"meta": {"method": "IMPRES (canonical 15-pair, fixed scorer)",
                    "protocol": "v33 OOF scores via cv_predict; metrics "
                                "identical to hugo_envar_clinical_v33.py; "
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
