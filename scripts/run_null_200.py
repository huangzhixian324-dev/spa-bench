"""Upgrade the cohort-matched CDS nulls from 10 to 200 replicates for all
five manuscript endpoints (Table 3), using the restored ORIGINAL h5ad
files and the exact v1.x CDS pipeline. Motivation (reviewer R2): the max
of 10 replicates is an extremely noisy order statistic; 200 replicates
give a stable null mean/max/q95 and eliminate that criticism.

Output: results/benchmark/v33/cds_nulls_200_v33.json with, per endpoint,
the 200-replicate distribution AND a consistency comparison against the
archived 10-replicate values.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
os_chdir = None
import os
os.chdir(REPO)
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "cds_tool"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

from rerun_v33 import COHORTS, load_cohort  # noqa: E402
import importlib.util  # noqa: E402
# 显式加载 cds_tool/cds.py（教训 13：workflows/evaluation/cds.py 会遮蔽它）
_spec = importlib.util.spec_from_file_location(
    "cds_tool_cds", REPO / "cds_tool" / "cds.py")
_cds_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_cds_mod)
circularity_detection_score = _cds_mod.circularity_detection_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402

# endpoint label -> COHORTS key (Table 3's five endpoints)
ENDPOINTS = [
    ("Hugo 2016 RECIST", "Hugo_2016"),
    ("Lauss 2017 RECIST (ACT)", "Nathanson_2017"),
    ("Gide 2019 RECIST", "Gide_2019_cBio"),
    ("Riaz 2017 RECIST (v33)", "Riaz_2017_RECIST_v33"),
    ("Riaz 2017 cytolytic", "Riaz_2017_cytolytic"),
]
N_REP = 200
N_FEATURES = 500

# Riaz 端点必须用 HGNC 映射矩阵（默认 processed 矩阵是 Entrez ID，
# GZMA/PRF1 找不到会返回 UNKNOWN -> null 全 0；与 fill_nc_cells.py 的
# 路径覆盖先例一致）
COHORTS["Riaz_2017_cytolytic"]["h5ad"] = (
    "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad")
COHORTS["Riaz_2017_RECIST_v33"]["h5ad"] = (
    "data/cohorts/Riaz_2017_RECIST/processed/"
    "Riaz_2017_RECIST_v33_HGNC.h5ad")


def main():
    out_path = REPO / "results/benchmark/v33/cds_nulls_200_v33.json"
    done = {}
    if out_path.exists():
        done = json.load(open(out_path))
    results = {}

    for label, ckey in ENDPOINTS:
        if label in done and done[label].get("n_replicates") >= N_REP:
            print(f"[skip] {label} already has {N_REP} replicates", flush=True)
            results[label] = done[label]
            continue
        adata, folds, y = load_cohort(ckey)
        import scanpy as sc
        X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
        gene_list = list(adata.var_names)
        y = np.asarray(y).astype(int)
        rng = np.random.RandomState(42)
        t0 = time.time()
        vals = []
        for i in range(N_REP):
            y_perm = rng.permutation(y)
            r = circularity_detection_score(X, y_perm, ["GZMA", "PRF1"],
                                            gene_list,
                                            n_features=N_FEATURES,
                                            random_state=42)
            vals.append(float(r["CDS"]))
            if (i + 1) % 25 == 0:
                print(f"  {label}: {i+1}/{N_REP} "
                      f"({time.time()-t0:.0f}s)", flush=True)
        vals = np.array(vals)
        results[label] = {
            "n_replicates": N_REP,
            "mean": round(float(vals.mean()), 4),
            "sd": round(float(vals.std()), 4),
            "q95": round(float(np.percentile(vals, 95)), 4),
            "max": round(float(vals.max()), 4),
            "values": [round(float(v), 4) for v in vals],
        }
        print(f"[done] {label}: mean={results[label]['mean']} "
              f"max={results[label]['max']} "
              f"({time.time()-t0:.0f}s)", flush=True)
        # incremental save
        json.dump(results, open(out_path, "w"), indent=1)

    # consistency check vs the archived 10-replicate nulls
    arch = json.load(open(REPO / "results/benchmark/v33/cds_nulls_v33.json"))
    print("\n=== 10-rep (archived) vs 200-rep (this run) ===", flush=True)
    mapping = {"Hugo 2016 RECIST": "Hugo 2016 RECIST",
               "Lauss 2017 RECIST (ACT)": "Nathanson 2017 RECIST",
               "Gide 2019 RECIST": "Gide 2019 RECIST",
               "Riaz 2017 RECIST (v33)": "Riaz 2017 RECIST (v33)",
               "Riaz 2017 cytolytic": "Riaz 2017 cytolytic"}
    for label, akey in mapping.items():
        a, n = arch[akey], results[label]
        print(f"  {label:28s} 10rep mean/max = {a['mean']:.3f}/{a['max']:.3f}"
              f"   200rep mean/max = {n['mean']:.3f}/{n['max']:.3f}", flush=True)
    print("written", out_path, flush=True)
    return 0


if __name__ == "__main__":
    main()
