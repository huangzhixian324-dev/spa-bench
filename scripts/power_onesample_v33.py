"""Corrected one-sample AUROC power table + the n>=50 benchmarking anchor
(A-defect #8 fix, 2026-09-17 quality campaign).

Two defects in the archived power_onesample_v33.json are corrected here:
  1. Responder counts. Three rows used placeholder (n_pos, n_neg) values
     (Hugo (16,12) vs true (13,15); Lauss (16,9) vs true (10,15);
     Riaz RECIST (21,21) vs true (9,33)). Gide/Jung/Riaz-cytolytic were
     correct.
  2. Test convention. The archived computation applied a 1.96 critical
     value while the manuscript states one-sided alpha = 0.05 (z = 1.645).
     1.96 corresponds to one-sided alpha = 0.025 (conservative). This
     script computes BOTH conventions; the JSON marks the one-sided
     1.645 column as primary, consistent with the one-sided permutation
     framework used throughout the manuscript, and keeps the 1.96 column
     for the audit trail.

It also computes the anchor for the Table 9 recommendation "n >= 50
minimum floor": the smallest cohort size at which a true AUROC in the
realistic-effect band (0.65-0.70) becomes detectable at 80% power
(one-sided alpha = 0.05, Hanley-McNeil), across responder proportions
40-55%.

Output: results/benchmark/v33/power_onesample_v33.json (regenerated;
archived values retained under 'archived_v33_values' for the audit trail).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "benchmark" / "v33" / "power_onesample_v33.json"


def se_auroc(A, n_pos, n_neg):
    """Hanley-McNeil (1982) variance for a single AUROC estimate."""
    Q1 = A / (2 - A)
    Q2 = 2 * A ** 2 / (1 + A)
    var = (A * (1 - A) + (n_pos - 1) * (Q1 - A ** 2) +
           (n_neg - 1) * (Q2 - A ** 2)) / (n_pos * n_neg)
    return float(np.sqrt(var))


def power(A, n_pos, n_neg, zcrit):
    return float(norm.cdf((A - 0.5) / se_auroc(A, n_pos, n_neg) - zcrit))


def min_auroc(n_pos, n_neg, zcrit, target=0.80):
    lo, hi = 0.5, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if power(mid, n_pos, n_neg, zcrit) < target:
            lo = mid
        else:
            hi = mid
    return float(hi)


TRUE_COUNTS = {
    # (n_pos, n_neg) as reported in the manuscript cohort table
    'Hugo_2016': (13, 15),
    'Nathanson_2017': (10, 15),
    'Gide_2019_cBio': (40, 33),
    'Jung_2019_DCB': (6, 21),
    'Riaz_2017_cytolytic': (21, 22),
    'Riaz_2017_RECIST_v33': (9, 33),
}
ARCHIVED_COUNTS = {
    'Hugo_2016': (16, 12), 'Nathanson_2017': (16, 9),
    'Gide_2019_cBio': (40, 33), 'Jung_2019_DCB': (6, 21),
    'Riaz_2017_cytolytic': (21, 22), 'Riaz_2017_RECIST_v33': (21, 21),
}
ARCHIVED_MINS = {
    'Hugo_2016': 0.756, 'Nathanson_2017': 0.767, 'Gide_2019_cBio': 0.675,
    'Jung_2019_DCB': 0.817, 'Riaz_2017_cytolytic': 0.72,
    'Riaz_2017_RECIST_v33': 0.722,
}


def main():
    table = {}
    for cohort, (n1, n0) in TRUE_COUNTS.items():
        a1, a0 = ARCHIVED_COUNTS[cohort]
        row = {
            'n_pos_true': n1, 'n_neg_true': n0,
            'archived_counts_were': [a1, a0],
            'archived_min_auroc_80pct_power': ARCHIVED_MINS[cohort],
        }
        for A in [0.55, 0.60, 0.63, 0.66, 0.70, 0.75]:
            row[f'AUROC={A}'] = {
                'se': round(se_auroc(A, n1, n0), 4),
                'power_onesided_z1645': round(power(A, n1, n0, 1.645), 3),
                'power_z1964': round(power(A, n1, n0, 1.96), 3),
            }
        row['min_auroc_80pct_power_onesided_z1645'] = round(
            min_auroc(n1, n0, 1.645), 3)
        row['min_auroc_80pct_power_z1964'] = round(min_auroc(n1, n0, 1.96), 3)
        table[cohort] = row
        print(f"{cohort}: min AUROC @80% one-sided 1.645 = "
              f"{row['min_auroc_80pct_power_onesided_z1645']} "
              f"(1.96: {row['min_auroc_80pct_power_z1964']}; "
              f"archived: {row['archived_min_auroc_80pct_power']})",
              flush=True)

    # ---- Table 9 anchor: smallest n with a realistic effect detectable ----
    anchor = {}
    for resp_frac in [0.40, 0.46, 0.50, 0.55]:
        row = {}
        for A_true in [0.65, 0.68, 0.70]:
            lo, hi = 20, 200
            for _ in range(60):
                n = (lo + hi) // 2
                n1 = int(round(n * resp_frac))
                n0 = n - n1
                if n1 < 2 or n0 < 2:
                    lo = n + 1
                    continue
                if power(A_true, n1, n0, 1.645) < 0.80:
                    lo = n + 1
                else:
                    hi = n
            row[f'AUROC={A_true}'] = hi
        anchor[f'responder_frac_{resp_frac}'] = row
        print(f"anchor responder_frac={resp_frac}: {row}", flush=True)

    out = {
        'method': 'Hanley-McNeil single-sample AUROC variance, normal '
                  'approximation; one-sided alpha = 0.05 (z = 1.645), '
                  'H0: AUROC = 0.5; the 1.96 column (one-sided 0.025) is '
                  'retained for the audit trail',
        'v33_correction': 'regenerated 2026-09-17 with the true responder '
                          'counts (archived rows for Hugo/Lauss/Riaz-RECIST '
                          'used placeholder counts) and the one-sided '
                          'critical value consistent with the permutation '
                          'framework',
        'table': table,
        'table9_anchor': {
            'question': 'smallest n at which a true AUROC becomes '
                        'detectable at 80% power (one-sided 0.05)',
            'sweep': anchor,
        },
    }
    json.dump(out, open(OUT, 'w'), indent=1)
    print('WROTE', OUT, flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
