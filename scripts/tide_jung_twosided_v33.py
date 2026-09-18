"""Jung 2019 DCB TIDE permutation — lower-tail and two-sided statistics
(B5): the archived one-sided p = 0.952 establishes only that 0.278 is not
above the null envelope; the far-below-chance reading needs the lower
tail. Fixed-scorer protocol (label-independent score, identical to
Table S11): 50,000 label shuffles, full null distribution retained.

Writes results/benchmark/v33/tide_jung_twosided_v33.json.
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "workflows" / "methods"))
sys.path.insert(0, str(REPO / "workflows" / "evaluation"))

from rerun_v33 import SEED, load_cohort  # noqa: E402
from base_method import TIDE_Wrapper  # noqa: E402

OUT = REPO / "results" / "benchmark" / "v33" / "tide_jung_twosided_v33.json"


def auroc_fixed(score, y):
    r = rankdata(score)
    n1 = int(y.sum())
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def main():
    adata, folds, y = load_cohort("Jung_2019_DCB")
    y = np.asarray(y).astype(int)
    n = len(y)

    w = TIDE_Wrapper()
    idx_all = list(range(n))
    w.fit(adata, idx_all)
    scores = np.asarray(w.predict(adata, idx_all), dtype=float)

    obs = auroc_fixed(scores, y)
    print(f"n={n}, pos={int(y.sum())}, TIDE obs AUROC = {obs:.4f} "
          f"(archived benchmark 0.2778)", flush=True)

    rng = np.random.RandomState(SEED)
    N = 50000
    null = np.empty(N)
    for i in range(N):
        null[i] = auroc_fixed(scores, rng.permutation(y))
        if (i + 1) % 10000 == 0:
            print(f"  {i + 1}/{N}", flush=True)

    ge = int(np.sum(null >= obs - 1e-12))
    lt = int(np.sum(null < obs + 1e-12))
    p_upper = (ge + 1) / (N + 1)
    p_lower = (lt + 1) / (N + 1)
    p_two = float(min(1.0, 2 * min(p_upper, p_lower)))

    out = {
        "cohort": "Jung_2019_DCB", "method": "TIDE", "n": n,
        "n_pos": int(y.sum()),
        "obs_auroc": round(obs, 4),
        "archived_benchmark_auroc": 0.2778,
        "n_shuffle": N,
        "p_upper_onesided": round(p_upper, 5),
        "p_lower": round(p_lower, 5),
        "p_two_sided": round(p_two, 5),
        "null_mean": round(float(null.mean()), 4),
        "null_sd": round(float(null.std()), 4),
        "null_percentiles": {str(q): round(float(np.percentile(null, q)), 4)
                             for q in (1, 5, 25, 50, 75, 95, 99)},
        "protocol": "fixed-scorer label shuffles (50,000), seed 42; score "
                    "label-independent (tidepy v1.3.9, default vthres)",
    }
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"p_upper={p_upper:.5f}  p_lower={p_lower:.5f}  "
          f"p_two={p_two:.5f}", flush=True)
    print("WROTE", OUT, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
