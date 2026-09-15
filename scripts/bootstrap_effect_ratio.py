"""Effect-size bootstrap CI for Design Choice 1 (old-R1-M2 / partA-M4).

Patient-level bootstrap of the Riaz 2017 endpoint-switch effect:
- delta_collapse = AUROC(cytolytic surrogate) - AUROC(RECIST) for ElasticNet-Var
- gap = max - min method AUROC on the RECIST endpoint (same patients)
- ratio = delta / gap
Resampling is over patients (out-of-fold prediction scores, fixed models,
seed 42), 1,000 bootstrap replicates, percentile CIs. Both endpoints are run
on the SAME fold split so patient indices correspond across endpoints.
"""
import sys
import json
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / 'scripts'), str(REPO / 'workflows' / 'methods'),
                str(REPO / 'workflows' / 'evaluation')]
from rerun_v33 import load_cohort, cv_predict  # noqa: E402

METHODS = ['IMPRES', 'GEP', 'TIDE', 'PD_L1', 'ElasticNet', 'ElasticNet_Var']
PARAMS = {'ElasticNet': {'C': 0.1, 'l1_ratio': 0.5},
          'ElasticNet_Var': {'C': 0.1, 'l1_ratio': 0.5}}
SURR = 'Riaz_2017_cytolytic'
CLIN = 'Riaz_2017_RECIST_v33'
B = 1000

print('loading cohorts...', flush=True)
adata_s, folds_s, y_s = load_cohort(SURR)
adata_r, folds_r, y_r = load_cohort(CLIN)
# use ONE fold split (the clinical one) for both endpoints so patient order matches
folds = folds_r
# align the surrogate adata to the clinical adata's patient order
if list(adata_s.obs_names) != list(adata_r.obs_names):
    adata_s = adata_s[adata_r.obs_names].copy()
    print('surrogate adata reordered to clinical patient order', flush=True)

scores = {}
for name, ad, methods in [(SURR, adata_s, ['ElasticNet_Var']),
                          (CLIN, adata_r, ['IMPRES', 'GEP', 'TIDE', 'PD_L1_IHC',
                                           'ElasticNet', 'ElasticNet_Var'])]:
    scores[name] = {}
    for m in methods:
        print(f'  cv_predict {name} / {m}', flush=True)
        yt, yp = cv_predict(ad, folds, m, params=PARAMS.get(m))
        scores[name][m] = {'y': yt.tolist(), 'p': yp.tolist()}

# patient-order identity check: same obs order across the two adata objects
try:
    assert list(adata_s.obs_names) == list(adata_r.obs_names), 'patient order differs'
    print('patient order identity: OK', flush=True)
except Exception as e:
    print('patient order WARNING:', e, flush=True)

json.dump(scores, open(REPO / 'results/benchmark/v33/riaz_oof_scores_v33.json', 'w'))

n = len(scores[CLIN]['GEP']['y'])
rng = np.random.RandomState(42)
deltas, gaps, ratios = [], [], []
yv_s = {m: np.array(scores[SURR][m]['y']) for m in METHODS}
pv_s = {m: np.array(scores[SURR][m]['p']) for m in METHODS}
yv_r = {m: np.array(scores[CLIN][m]['y']) for m in scores[CLIN]}
pv_r = {m: np.array(scores[CLIN][m]['p']) for m in scores[CLIN]}

for b in range(B):
    idx = rng.randint(0, n, n)

    def auc(yv, pv, m):
        y, s = yv[m][idx], pv[m][idx]
        if len(set(y)) < 2:
            return np.nan
        return roc_auc_score(y, s)

    a_s = {m: auc(yv_s, pv_s, m) for m in METHODS}
    a_r = {m: auc(yv_r, pv_r, m) for m in scores[CLIN]}
    delta = a_s['ElasticNet_Var'] - a_r['ElasticNet_Var']
    vals = [v for v in a_r.values() if not np.isnan(v)]
    gap = max(vals) - min(vals) if vals else np.nan
    deltas.append(delta)
    gaps.append(gap)
    if gap and not np.isnan(gap) and gap > 0 and not np.isnan(delta):
        ratios.append(delta / gap)

out = {
    'n_patients': n,
    'bootstrap_replicates': B,
    'delta_collapse_ENVar': {
        'point': float(yv_s['ElasticNet_Var'].size and
                       (roc_auc_score(yv_s['ElasticNet_Var'], pv_s['ElasticNet_Var']) -
                        roc_auc_score(yv_r['ElasticNet_Var'], pv_r['ElasticNet_Var']))),
        'ci95': [float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))]},
    'gap_max_method_RECIST': {
        'ci95': [float(np.nanpercentile(gaps, 2.5)), float(np.nanpercentile(gaps, 97.5))]},
    'ratio_delta_over_gap': {
        'point': 3.5,
        'ci95': [float(np.percentile(ratios, 2.5)), float(np.percentile(ratios, 97.5))],
        'n_valid': len(ratios)},
    'seed': 42,
}
json.dump(out, open(REPO / 'results/benchmark/v33/effect_size_bootstrap_v33.json', 'w'),
          indent=1)
print(json.dumps(out, indent=1))
