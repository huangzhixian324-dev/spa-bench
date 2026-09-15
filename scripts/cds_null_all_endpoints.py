"""Random-label null distribution of CDS for all five evaluable endpoints.

Each endpoint gets a 10-replicate random-label null of CDS v1.1.0 on its
own feature space (declared genes GZMA/PRF1, same convention as Table
S4b). This is the reference frame the reworked CDS claims are stated
against (continuous, null-referenced) rather than absolute thresholds.

Output: results/benchmark/v33/cds_nulls_v33.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "cds_tool"))
from rerun_v33 import COHORTS, SEED, load_cohort  # noqa: E402
from cds import circularity_detection_score  # noqa: E402

# v33 HGNC matrices for Riaz (same as final_benchmarks_v33.py)
COHORTS["Riaz_2017_cytolytic"]["h5ad"] = (
    "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad")
COHORTS["Riaz_2017_RECIST_v33"]["h5ad"] = (
    "data/cohorts/Riaz_2017_RECIST/processed/"
    "Riaz_2017_RECIST_v33_HGNC.h5ad")

COHORTS_ = [("Hugo_2016", "Hugo 2016 RECIST"),
            ("Nathanson_2017", "Nathanson 2017 RECIST"),
            ("Gide_2019_cBio", "Gide 2019 RECIST"),
            ("Riaz_2017_cytolytic", "Riaz 2017 cytolytic"),
            ("Riaz_2017_RECIST_v33", "Riaz 2017 RECIST (v33)")]


def main():
    out = {}
    for cname, label in COHORTS_:
        adata, _, y = load_cohort(cname)
        X = (adata.X.toarray() if hasattr(adata.X, "toarray")
             else np.asarray(adata.X))
        rng = np.random.RandomState(SEED)
        vals = []
        for _ in range(10):
            r = circularity_detection_score(X, rng.permutation(y),
                                            ["GZMA", "PRF1"],
                                            list(adata.var_names))
            vals.append(float(r["CDS"]))
        out[label] = {"n": 10, "mean": round(float(np.mean(vals)), 3),
                      "sd": round(float(np.std(vals)), 3),
                      "max": round(float(np.max(vals)), 3),
                      "q95": round(float(np.quantile(vals, 0.95)), 3),
                      "values": [round(v, 3) for v in vals]}
        print(f"{label}: mean {out[label]['mean']} max {out[label]['max']}",
              flush=True)
    with open(REPO / "results" / "benchmark" / "v33" / "cds_nulls_v33.json",
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
