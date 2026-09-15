"""Two minor items in one run:
(A) XGBoost exclusion evidence: evaluate XGBoost on all six cohort-endpoints
    with the same CV protocol, record AUROC and prediction degeneracy
    (unique predicted probabilities), archive as JSON so the exclusion
    ("degenerate by construction") is auditable.
(B) One-sample power: closed-form power to detect AUROC > 0.5 (Hanley-McNeil
    single-sample variance, normal approximation) at the cohort sizes used,
    for a range of true AUROCs. Writes power_onesample_v33.json.
"""
import json
import sys
import numpy as np
from pathlib import Path
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / 'scripts'), str(REPO / 'workflows' / 'methods'),
                str(REPO / 'workflows' / 'evaluation')]
from rerun_v33 import load_cohort, cv_predict  # noqa: E402

OUT = REPO / 'results/benchmark/v33'

# ---------- (A) XGBoost degeneracy archive ----------
COHORTS = ['Hugo_2016', 'Nathanson_2017', 'Gide_2019_cBio', 'Jung_2019_DCB',
           'Riaz_2017_cytolytic', 'Riaz_2017_RECIST_v33']
xgb_records = {}
for cohort in COHORTS:
    try:
        adata, folds, y = load_cohort(cohort)
        yt, yp = cv_predict(adata, folds, 'XGBoost')
        uniq = len(np.unique(np.round(yp, 6)))
        rec = {
            'n': int(len(yt)),
            'n_pos': int(yt.sum()),
            'auroc': (round(float(roc_auc_score(yt, yp)), 4)
                      if len(set(yt)) >= 2 else None),
            'unique_predictions': int(uniq),
            'degenerate': bool(uniq <= 1),
            'obs_auroc_this_run': (round(float(roc_auc_score(yt, yp)), 4)
                                   if len(set(yt)) >= 2 else None),
        }
        xgb_records[cohort] = rec
        print(f'XGBoost {cohort}: {rec}', flush=True)
    except Exception as e:
        xgb_records[cohort] = {'error': str(e)}
        print(f'XGBoost {cohort}: ERROR {e}', flush=True)

json.dump({
    'purpose': 'Auditable exclusion evidence for XGBoost (review items R1-m11/R2-m6): '
               'the method is excluded from all reported matrices because it is '
               'degenerate by construction at these sample sizes. This archive '
               'records the degeneracy determination (unique predicted values) and '
               'the AUROC of each run so the exclusion is verifiable, even though '
               'the prediction tables themselves are not reported.',
    'protocol': 'same patient-stratified CV and preprocessing as all other methods '
                '(rerun_v33.cv_predict), default XGBoost hyperparameters',
    'records': xgb_records,
    'verdict': 'excluded: degenerate (constant or single-value predictions) on '
               'cohorts where unique_predictions <= 1',
}, open(OUT / 'xgb_exclusion_evidence_v33.json', 'w'), indent=1)
print('saved xgb_exclusion_evidence_v33.json', flush=True)

# ---------- (B) One-sample power ----------
def se_auroc(A, n_pos, n_neg):
    """Hanley-McNeil (1982) variance for a single AUROC estimate."""
    Q1 = A / (2 - A)
    Q2 = 2 * A ** 2 / (1 + A)
    var = (A * (1 - A) + (n_pos - 1) * (Q1 - A ** 2) +
           (n_neg - 1) * (Q2 - A ** 2)) / (n_pos * n_neg)
    return float(np.sqrt(var))


N_BY_COHORT = {'Hugo_2016': (16, 12), 'Lauss_2017': (16, 9), 'Gide_2019_cBio': (40, 33),
               'Jung_2019_DCB': (6, 21), 'Riaz_2017_cytolytic': (21, 22),
               'Riaz_2017_RECIST_v33': (21, 21)}
power_tab = {}
for cohort, (n1, n0) in N_BY_COHORT.items():
    row = {}
    for A in [0.55, 0.60, 0.63, 0.66, 0.70, 0.75]:
        se = se_auroc(A, n1, n0)
        z = (A - 0.5) / se
        row[f'AUROC={A}'] = {'se': round(se, 4),
                             'power_onesample': round(float(norm.cdf(z - 1.96)), 3)}
    # minimum AUROC detectable with 80% power
    lo, hi = 0.5, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        se = se_auroc(mid, n1, n0)
        if norm.cdf((mid - 0.5) / se - 1.96) < 0.80:
            lo = mid
        else:
            hi = mid
    row['min_auroc_80pct_power'] = round(hi, 3)
    power_tab[cohort] = row
    print(f'{cohort} (n1={n1},n0={n0}): min AUROC @80% power = {row["min_auroc_80pct_power"]}',
          flush=True)

json.dump({
    'method': 'Hanley-McNeil single-sample AUROC variance, normal approximation, '
              'one-sided alpha 0.05, H0: AUROC = 0.5',
    'note': 'Complements Table 8 (two-method comparison power); answers the '
            'one-sample detectability question used in the Results text.',
    'table': power_tab,
}, open(OUT / 'power_onesample_v33.json', 'w'), indent=1)
print('saved power_onesample_v33.json', flush=True)
