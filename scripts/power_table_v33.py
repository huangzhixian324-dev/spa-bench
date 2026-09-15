"""Power analysis documentation (2026-09-07 review round).

Computes the minimum detectable delta-AUROC at 80% power (alpha = 0.05,
two-sided) using the Hanley-McNeil asymptotic variance of the AUC
(Mann-Whitney U), for a PAIRED comparison of two methods on the same
patients with AUC correlation r.  Baseline AUC A1 = 0.5 (the relevant
null for the small cohorts studied here).

Output: results/benchmark/v33/power_table_v33.json
  - grid: manuscript Table 7 rows recomputed under r in {0, 0.5, 0.75}
  - n_for_delta_010: minimum balanced n for delta = 0.10 per r

Motivation: the archived Table 7 values were computed ad hoc and are not
reproduced by any single-assumption model of this form; this script makes
the model and its sensitivity to the paired-correlation assumption
explicit and machine-readable.
"""
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "benchmark" / "v33" / "power_table_v33.json"

Z = 1.959963984540054 + 0.841621233572914  # z(0.975) + z(0.80) = 2.8016


def var_auc(a, n1, n0):
    q1 = a / (2.0 - a)
    q2 = 2.0 * a * a / (1.0 + a)
    return (a * (1 - a) + (n1 - 1) * (q1 - a * a)
            + (n0 - 1) * (q2 - a * a)) / (n1 * n0)


def var_delta(d, n1, n0, r):
    v1 = var_auc(0.5, n1, n0)
    v2 = var_auc(0.5 + d, n1, n0)
    return v1 + v2 - 2.0 * r * (v1 * v2) ** 0.5


def min_detectable(n1, n0, r):
    lo, hi = 1e-4, 0.49
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if Z * var_delta(mid, n1, n0, r) ** 0.5 < mid:
            hi = mid
        else:
            lo = mid
    return hi


def n_for_delta(delta, r):
    for n in range(10, 2001, 2):
        if min_detectable(n // 2, n - n // 2, r) <= delta:
            return n
    return None


ROWS = [(25, 0.40, "Nathanson 2017"), (28, 0.46, "Hugo 2016"),
        (42, 0.21, "Riaz RECIST"), (73, 0.55, "Gide 2019"),
        (100, 0.50, "-"), (200, 0.50, "-")]
RS = [0.0, 0.5, 0.75]

out = {"model": ("Hanley-McNeil asymptotic AUC variance; paired "
                 "two-method comparison on the same patients, baseline "
                 "AUC 0.5, alpha 0.05 two-sided, power 0.80"),
       "z": Z, "rows": [], "n_for_delta_010": {}}
for n, prev, example in ROWS:
    n1 = max(1, round(n * prev))
    n0 = n - n1
    row = {"n": n, "responder_pct": prev, "n1": n1, "n0": n0,
           "example": example}
    for r in RS:
        row[f"min_delta_r{r}"] = round(min_detectable(n1, n0, r), 3)
    out["rows"].append(row)
for r in RS:
    out["n_for_delta_010"][f"r{r}"] = n_for_delta(0.10, r)

with open(OUT, "w") as fh:
    json.dump(out, fh, indent=1)

print(f"{'n':>4} {'resp%':>6} | " + " | ".join(f"r={r:<4}" for r in RS)
      + " | Table 7 (archived)")
archived = {(25, 0.40): 0.25, (28, 0.46): 0.22, (42, 0.21): 0.27,
            (73, 0.55): 0.15, (100, 0.50): 0.12, (200, 0.50): 0.08}
for row in out["rows"]:
    print(f"{row['n']:>4} {row['responder_pct']:>6} | "
          + " | ".join(f"{row[f'min_delta_r{r}']:<5}" for r in RS)
          + f" | {archived[(row['n'], row['responder_pct'])]}")
print("n for delta=0.10 @80% power:", out["n_for_delta_010"])
print("wrote", OUT)
