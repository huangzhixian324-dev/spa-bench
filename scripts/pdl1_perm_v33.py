"""Complete the v33 permutation matrix with PD-L1 cells.

PD-L1 (CD274 expression) is a fixed scorer: its score is label-independent,
so the full-pipeline label-shuffle null is mathematically identical to
shuffling labels against the fixed score (the same equivalence used for
IMPRES/GEP/TIDE, verified on Hugo IMPRES). The v33 protocol
(fixed scorers = IMPRES/GEP/TIDE/PD-L1, 50,000 prediction shuffles)
therefore extends to PD-L1 exactly; these cells were simply missing from
permutation_v33.json.

For each cohort the observed AUROC must reproduce the benchmark value
(identical fixed score), which is asserted before writing.

Output: updates results/benchmark/v33/permutation_v33.json in place
(6 new PD_L1 cells + meta note). Deterministic: seed 42, 50,000 shuffles.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from rerun_v33 import COHORTS, SEED, load_cohort  # noqa: E402

COHORTS["Riaz_2017_cytolytic"]["h5ad"] = (
    "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad")
COHORTS["Riaz_2017_RECIST_v33"]["h5ad"] = (
    "data/cohorts/Riaz_2017_RECIST/processed/"
    "Riaz_2017_RECIST_v33_HGNC.h5ad")

PF = REPO / "results" / "benchmark" / "v33" / "permutation_v33.json"
BF = REPO / "results" / "benchmark" / "v33" / "benchmark_v33.json"
N_SHUFFLE = 50000
GENE = "CD274"


def oof_pdl1_score(adata, folds):
    """Replicate PD_L1_Wrapper exactly: CD274 min-max normalized within
    each test fold, pooled over folds (label-independent fixed vector)."""
    X = (adata.X.toarray() if hasattr(adata.X, "toarray")
         else np.asarray(adata.X))
    genes = list(adata.var_names)
    raw = X[:, genes.index(GENE)].astype(float) if GENE in genes \
        else np.zeros(adata.n_obs)
    out = np.zeros(adata.n_obs)
    for f in folds:
        te = np.array(f["test"], dtype=int)
        s = raw[te]
        rng_ = s.max() - s.min()
        out[te] = (s - s.min()) / rng_ if rng_ > 0 else np.zeros(len(te))
    return out


def auroc_fixed(score, y):
    r = rankdata(score)
    n1 = int(y.sum())
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def main():
    perms = json.load(open(PF))
    bm = json.load(open(BF))["cohorts"]
    rng = np.random.RandomState(SEED)

    for cname in bm:
        adata, folds, y = load_cohort(cname)
        score = oof_pdl1_score(adata, folds)
        obs = auroc_fixed(score, y)
        ref = bm[cname]["PD_L1"]["auroc"]
        assert abs(obs - ref) < 5e-4, f"{cname}: obs {obs:.4f} != bench {ref}"
        ge = 0
        for _ in range(N_SHUFFLE):
            yp = rng.permutation(y)
            if auroc_fixed(score, yp) >= obs:
                ge += 1
        p = (ge + 1) / (N_SHUFFLE + 1)
        perms["cohorts"].setdefault(cname, {})["PD_L1"] = {
            "obs_auroc": round(float(obs), 4), "p": p,
            "n_shuffle": N_SHUFFLE}
        print(f"{cname}: PD_L1 obs={obs:.4f} p={p:.5f}", flush=True)

    note = (" PD_L1 cells added 2026-09-04 (scripts/pdl1_perm_v33.py): "
            "fixed-scorer protocol, 50,000 label shuffles against the fixed "
            "CD274 score; observed AUROC asserted equal to the benchmark "
            "value.")
    perms["meta"]["note"] = perms["meta"].get("note", "") + note
    with open(PF, "w") as fh:
        json.dump(perms, fh, indent=1)
    print("updated", PF)
    return 0


if __name__ == "__main__":
    sys.exit(main())
