"""CDS v2.0.0 validation on real endpoints (reviewer items R1-M4/R2-M4/R3-M2).

Circular-by-construction endpoint: IMmotion150 cytolytic surrogate
(median split of GZMA/PRF1 mean on the HGNC matrix; archived v1 CDS 0.945).
Genuine-biology endpoint: Gide RECIST (v1 CDS 0.852, genuine signal).
Null endpoints: Hugo RECIST (0.449), Lauss RECIST (0.781).
permutation=200, hypothesis_free=True, n_hf_features=2000.

Output: results/benchmark/v33/cds_v2_validation.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

REPO = Path(r"C:\Users\hzx\projects\spa-bench")
sys.path.insert(0, str(REPO / "cds_tool"))
import cds as cds_mod  # noqa: E402

OUT = REPO / "results" / "benchmark" / "v33" / "cds_v2_validation.json"
PERM = 200
EP = ["GZMA", "PRF1"]


def load(path):
    a = sc.read_h5ad(path)
    X = (a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X))
    return X, a


def main():
    out = {"meta": {"cds_version": cds_mod.__version__,
                    "permutation": PERM,
                    "hypothesis_free": True,
                    "n_hf_features": 2000,
                    "declared_genes": EP,
                    "note": "v1.x composite unchanged; v2 block is additive"},
           "endpoints": {}}

    # 1) circular-by-construction: IMmotion150 cytolytic surrogate
    p = REPO / "data/cohorts/IMmotion150/processed/IMmotion150_HGNC.h5ad"
    X, a = load(p)
    gidx = [list(map(str, a.var_names)).index(g) for g in EP]
    cyt = X[:, gidx].mean(axis=1)
    y_cyt = (cyt > np.median(cyt)).astype(int)
    genes = [str(g) for g in a.var_names]

    def run(name, X, y, genes):
        print(f"== {name}: {X.shape[0]} x {X.shape[1]} pos={int(y.sum())}",
              flush=True)
        r = cds_mod.circularity_detection_score(
            X, y.astype(int), EP, genes, permutation=PERM,
            hypothesis_free=True)
        v1 = {"CDS": r["CDS"], "risk": r["risk_level"],
              "C1": r["components"]["C1_MI_endpoint_genes"],
              "C2": r["components"]["C2_rho_endpoint_response"],
              "C3": r["components"]["C3_MI_skew_ratio"],
              "C3_saturated": r["components"]["C3_saturated"]}
        v2 = r.get("v2", {})
        out["endpoints"][name] = {"v1": v1, "v2": v2}
        p3 = v2.get("permutation_c3", {}).get("empirical_p")
        dsp = v2.get("declared_gene_specificity", {}).get("empirical_p")
        dp = v2.get("declared_gene_specificity", {}).get("mean_percentile")
        hfp = v2.get("hypothesis_free", {}).get("empirical_p")
        print(f"   v1 CDS={v1['CDS']} [{v1['risk']}] C3={v1['C3']} "
              f"saturated={v1['C3_saturated']}")
        print(f"   v2 perm-C3 p={p3} | declared-gene pct={dp} p={dsp} | "
              f"hyp-free p={hfp}", flush=True)

    run("IMmotion150_cytolytic", X, y_cyt, genes)

    # 2) genuine biology: Gide RECIST
    p = REPO / "data/cohorts/Gide_2019_cBio/processed/Gide_2019_cBio_processed.h5ad"
    X, a = load(p)
    run("Gide_RECIST", X, a.obs["response"].values.astype(int),
        [str(g) for g in a.var_names])

    # 3) null: Hugo RECIST
    p = REPO / "data/cohorts/Hugo_2016/processed/Hugo_2016_processed.h5ad"
    X, a = load(p)
    run("Hugo_RECIST", X, a.obs["response"].values.astype(int),
        [str(g) for g in a.var_names])

    # 4) above-null-mean genuine case: Lauss RECIST
    p = REPO / "data/cohorts/Nathanson_2017/processed/Nathanson_2017_processed.h5ad"
    X, a = load(p)
    run("Lauss_RECIST", X, a.obs["response"].values.astype(int),
        [str(g) for g in a.var_names])

    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written:", OUT, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
