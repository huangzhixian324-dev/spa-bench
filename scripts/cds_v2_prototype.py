"""Dry experiment 1: CDS v2.0 prototype — permutation-based C3.

v1.x C3 uses a ratio (top-500-variance MI p90 / random-500 MI median) that
saturates at 1.0 on real bulk RNA-seq because the random-feature MI is near
zero. v2.0 replaces the ratio with a label-permutation z-score:

  For each of N permutations, permute y, recompute the 90th-percentile MI
  of the SAME top-500-variance features, and collect the permutation null.
  C3_v2 = (MI_obs_p90 - mean(MI_perm_p90)) / max(std(MI_perm_p90), eps)

This naturally handles immune-dominant feature spaces because the
permutation preserves the inter-gene correlation structure.

Runs on all five benchmark cohorts + IMvigor210, comparing v1.x and v2.0
on both the circular surrogate and clinical endpoints.

Output: results/benchmark/v33/cds_v2_prototype.json
"""
import json, sys, time, warnings
import numpy as np
import time
from pathlib import Path
from sklearn.feature_selection import mutual_info_classif
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "scripts"), str(REPO / "workflows" / "methods"),
                str(REPO / "workflows" / "evaluation")]
from rerun_v33 import SEED, COHORTS, load_cohort  # noqa: E402

# HGNC override (critical — see fill_nc_cells.py)
COHORTS["Riaz_2017_cytolytic"]["h5ad"] = "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad"
COHORTS["Riaz_2017_RECIST_v33"]["h5ad"] = "data/cohorts/Riaz_2017_RECIST/processed/Riaz_2017_RECIST_v33_HGNC.h5ad"

N_PERM = 500
N_FEATURES = 500
DECLARED = ["GZMA", "PRF1"]
OUT = REPO / "results" / "benchmark" / "v33" / "cds_v2_prototype.json"


def mi_p90(X, y, top_idx, rs):
    mi = mutual_info_classif(X[:, top_idx], y, random_state=rs)
    return float(np.percentile(mi, 90))


def mi_median_random(X, y, random_idx, rs):
    mi = mutual_info_classif(X[:, random_idx], y, random_state=rs)
    return float(np.median(mi))


def cds_v1(X, y, ep_idx, top_idx, random_idx, rs):
    # C1
    mi_ep = [mutual_info_classif(X[:, idx:idx+1], y, random_state=rs)[0]
             for idx in ep_idx]
    c1 = np.mean(mi_ep)
    # C2
    rho = []
    for idx in ep_idx:
        r, _ = __import__('scipy.stats', fromlist=['spearmanr']).spearmanr(
            X[:, idx], y.astype(float))
        rho.append(abs(r))
    c2 = np.mean(rho)
    # C3 (ratio-based, saturates)
    mi_t = mutual_info_classif(X[:, top_idx], y, random_state=rs)
    mi_r = mutual_info_classif(X[:, random_idx], y, random_state=rs)
    ratio = np.percentile(mi_t, 90) / max(np.median(mi_r), 0.0001)
    c3 = np.clip(ratio / 10.0, 0, 1)
    cds = 0.3 * np.clip(c1 / 0.05, 0, 1) + 0.3 * c2 + 0.4 * c3
    return {"CDS": round(float(cds), 3), "C1": round(float(c1), 4),
            "C2": round(float(c2), 4), "C3_v1": round(float(c3), 4),
            "C3_saturated": bool(ratio >= 10.0)}


def cds_v2(X, y, ep_idx, top_idx, random_idx, rs, n_perm=200):
    # C1, C2 same as v1
    mi_ep = [mutual_info_classif(X[:, idx:idx+1], y, random_state=rs)[0]
             for idx in ep_idx]
    c1 = np.mean(mi_ep)
    rho = []
    for idx in ep_idx:
        r, _ = __import__('scipy.stats', fromlist=['spearmanr']).spearmanr(
            X[:, idx], y.astype(float))
        rho.append(abs(r))
    c2 = np.mean(rho)

    # C3 v2: permutation-based z-score of top-feature MI p90
    mi_obs_p90 = mi_p90(X, y, top_idx, rs)
    rng = np.random.RandomState(rs + 1000)
    perm_p90 = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y)
        perm_p90.append(mi_p90(X, y, top_idx, rs) if False else
                        mi_p90(X, y_perm, top_idx, rs))
    mu = np.mean(perm_p90)
    sd = max(np.std(perm_p90), 1e-10)
    z = (mi_obs_p90 - mu) / sd
    # convert z to a 0-1 scale via the normal CDF
    from scipy.stats import norm
    c3_v2 = float(norm.cdf(z))  # P(Z <= z) = enrichment beyond chance
    c3_v2 = np.clip(c3_v2, 0, 1)

    cds = 0.3 * np.clip(c1 / 0.05, 0, 1) + 0.3 * c2 + 0.4 * c3_v2
    return {"CDS": round(float(cds), 3), "C1": round(float(c1), 4),
            "C2": round(float(c2), 4), "C3_v2": round(float(c3_v2), 4),
            "C3_z": round(float(z), 2),
            "MI_obs_p90": round(mi_obs_p90, 5),
            "MI_perm_mean": round(mu, 5),
            "MI_perm_std": round(sd, 5)}


def main():
    from rerun_v33 import load_cohort
    OUT = REPO / "results" / "benchmark" / "v33" / "cds_v2_prototype.json"
    results = {}

    # cohorts to test: (cohort_key, endpoint_name, endpoint_builder)
    # the surrogate is always the median split of GZMA/PRF1 mean
    cohorts_to_test = [
        ("Riaz_2017_cytolytic", "cytolytic_surrogate", "surrogate"),
        ("Riaz_2017_cytolytic", "RECIST_clinical", "clinical"),
        ("Gide_2019_cBio", "RECIST_clinical", "clinical"),
        ("Hugo_2016", "RECIST_clinical", "clinical"),
        ("Nathanson_2017", "RECIST_clinical", "clinical"),
        ("Jung_2019_DCB", "DCB_clinical", "clinical"),
    ]

    for cohort_key, ep_name, ep_type in cohorts_to_test:
        print(f"\n===== {cohort_key} / {ep_name} =====", flush=True)
        adata, folds, y = load_cohort(cohort_key)
        X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
        gene_list = list(adata.var_names)
        gi = {g: i for i, g in enumerate(gene_list)}

        if ep_type == "surrogate":
            cyt_idx = [gi[g] for g in DECLARED if g in gi]
            cyt_score = X[:, cyt_idx].mean(axis=1)
            y_bin = (cyt_score > np.median(cyt_score)).astype(int)
        else:
            y_bin = np.asarray(y).astype(int)

        ep_idx = [gi[g] for g in DECLARED if g in gi]
        top_idx = np.argsort(X.var(axis=0))[-N_FEATURES:]
        rng = np.random.RandomState(SEED)
        random_idx = rng.choice(X.shape[1], size=N_FEATURES, replace=False)

        rs = SEED
        v1 = cds_v1(X, y_bin, ep_idx, top_idx, random_idx, rs)
        t0 = time.time()
        v2 = cds_v2(X, y_bin, ep_idx, top_idx, random_idx, rs, n_perm=200)
        t_v2 = time.time() - t0

        is_circular = (ep_type == "surrogate")
        key = f"{cohort_key}__{ep_name}"
        results[key] = {
            "cohort": cohort_key, "endpoint": ep_name,
            "is_circular": is_circular, "n": int(len(y_bin)),
            "n_pos": int(y_bin.sum()),
            "v1": v1, "v2": v2, "v2_time_s": round(t_v2, 1),
        }
        sat = "SATURATED" if v1["C3_saturated"] else "ok"
        print(f"  v1: CDS={v1['CDS']:.3f} C3={v1['C3_v1']:.3f} [{sat}]")
        print(f"  v2: CDS={v2['CDS']:.3f} C3_v2={v2['C3_v2']:.3f} "
              f"z={v2['C3_z']:.1f} ({t_v2:.0f}s)")
        print(f"  circular={is_circular} | v1 {'flags' if v1['CDS'] > 0.70 else 'misses'} | "
              f"v2 {'flags' if v2['CDS'] > 0.70 else 'misses'}")

    # summary: does v2 separate circular from clinical better than v1?
    circ_v1 = [r["v1"]["CDS"] for k, r in results.items() if r["is_circular"]]
    clin_v1 = [r["v1"]["CDS"] for k, r in results.items() if not r["is_circular"]]
    circ_v2 = [r["v2"]["CDS"] for k, r in results.items() if r["is_circular"]]
    clin_v2 = [r["v2"]["CDS"] for k, r in results.items() if not r["is_circular"]]

    summary = {
        "v1_circular_mean": round(float(np.mean(circ_v1)), 3),
        "v1_clinical_mean": round(float(np.mean(clin_v1)), 3),
        "v1_separation": round(float(np.mean(circ_v1) - np.mean(clin_v1)), 3),
        "v2_circular_mean": round(float(np.mean(circ_v2)), 3),
        "v2_clinical_mean": round(float(np.mean(clin_v2)), 3),
        "v2_separation": round(float(np.mean(circ_v2) - np.mean(clin_v2)), 3),
        "v1_C3_saturated_count": sum(1 for r in results.values() if r["v1"]["C3_saturated"]),
        "total_endpoints": len(results),
    }
    print(f"\n===== v1 vs v2 separation =====")
    print(f"  v1: circular {summary['v1_circular_mean']} vs clinical {summary['v1_clinical_mean']} (Δ={summary['v1_separation']})")
    print(f"  v2: circular {summary['v2_circular_mean']} vs clinical {summary['v2_clinical_mean']} (Δ={summary['v2_separation']})")
    print(f"  C3 saturation: v1 {summary['v1_C3_saturated_count']}/{len(results)} endpoints")

    out = {"description": "CDS v2.0 prototype: permutation-based C3 replacing "
                          "the ratio-based C3 that saturates on real data",
           "n_perm": 200, "summary": summary, "results": results}
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print(f"\nWROTE {OUT}")
    return 0


if __name__ == "__main__":
    import time as _t
    sys.exit(main())
