"""SKCM sample-type sensitivity (B7): is the TCGA immune-signature
prognostic gradient a cancer-type effect or a metastatic-vs-primary
composition effect?

The v33 SKCM policy is all-tumour samples (n = 441; metastatic-dominated).
COAD (GDC primary-only, n = 458) and BRCA/LUAD (pan-cancer RSEM, primary-
dominant) are primary-heavy, so the cross-cancer gradient (11/12 -> 3/12 ->
1/12 -> 0/12) partially confounds cancer type with sample composition. This
script decomposes SKCM by TCGA sample-type code (06 = metastatic,
01 = primary) and reruns the identical 12-signature Cox + BH pipeline per
arm, so the gradient's interpretation can be conditioned on composition.

Arms: all (v33 baseline reproduction), metastatic-only ('06'),
primary-only ('01'). Output:
  data/tcga/SKCM/SKCM_cox_sampletype_sensitivity.json
"""
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

REPO = Path(__file__).resolve().parents[1]
TCGA_BARCODE = re.compile(r"^TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")

SIGS = {
    'IMPRES-like': ['PDCD1','CD27','CTLA4','CD28','CD86','CD80','CD274','LAG3','HAVCR2','TIGIT','TNFRSF9','ICOS'],
    'CD8 T cells': ['CD8A','CD8B'], 'B cells': ['CD19','CD79A','MS4A1'],
    'NK cells': ['NKG7','KLRD1','KLRF1'], 'IFNG response': ['IFNG','CXCL10','CXCL9','IDO1','STAT1'],
    'Cytolytic': ['GZMA','PRF1','GNLY'], 'Checkpoint': ['PDCD1','CTLA4','LAG3','TIGIT','HAVCR2'],
    'TLS': ['CXCL13','CCL19','CCL21'], 'Treg': ['FOXP3','IL2RA'],
    'M1 Macrophage': ['NOS2','IL12A','TNF'], 'Myeloid': ['CD14','CD33','ITGAM'],
    'Stromal': ['COL1A1','COL1A2','COL3A1','ACTA2','FAP'],
}


def run_arm(expr, clin, label):
    """Identical to run_tcga_cox.process_cohort, on a pre-filtered matrix."""
    sp = {}
    for c in expr.columns:
        parts = c.split('-')
        sp[c] = '-'.join(parts[:3]) if len(parts) >= 3 else c
    seen, cols = {}, []
    for c in expr.columns:
        pid = sp[c]
        if pid not in seen:
            seen[pid] = c
            cols.append(c)
    expr = expr[cols]
    pids = [sp[c] for c in expr.columns]
    clin_m = clin[clin.index.isin(pids)]
    pid_to_col = {sp[c]: c for c in expr.columns}
    matched_cols = [pid_to_col.get(pid) for pid in clin_m.index
                    if pid in pid_to_col]
    expr_m = expr[matched_cols]
    n = len(matched_cols)
    events = int(clin_m['OS_STATUS'].str.upper().str.contains('DECEASED').sum())
    print(f'  [{label}] n={n} events={events}', flush=True)
    X = np.log2(np.clip(expr_m.values, 0, None) + 1)
    genes = [str(g).upper() for g in expr_m.index]
    gidx = {g: i for i, g in enumerate(genes)}
    results = {}
    for sname, sgenes in SIGS.items():
        idx = [gidx[g] for g in sgenes if g.upper() in gidx]
        if len(idx) < 2:
            continue
        sig_vals = np.mean(X[idx, :], axis=0)
        dur = pd.to_numeric(clin_m['OS_MONTHS'], errors='coerce').values
        event = (clin_m['OS_STATUS'].astype(str).str.upper()
                 .str.contains('DECEASED').astype(int).values)
        df = pd.DataFrame({'dur': dur, 'event': event, 'sig': sig_vals})
        df = df.replace([np.inf, -np.inf], np.nan).dropna()
        df = df[(df['dur'] > 0) & (df['event'].isin([0, 1]))].copy()
        if df['event'].sum() < 5:
            continue
        cph = CoxPHFitter()
        cph.fit(df[['dur', 'event', 'sig']], duration_col='dur',
                event_col='event')
        results[sname] = {'HR': round(float(np.exp(cph.params_['sig'])), 3),
                          'p': float(cph.summary.loc['sig', 'p'])}
    tested = list(results)
    ps = np.array([results[s]['p'] for s in tested])
    order = np.argsort(ps)
    m = len(ps)
    bh = np.empty(m)
    bh[order] = np.minimum.accumulate(
        (ps[order] * m / (np.arange(m) + 1))[::-1])[::-1]
    for s, q in zip(tested, bh):
        results[s]['p_bh'] = float(min(q, 1.0))
    sig_n = sum(1 for v in results.values() if v['p'] < 0.05)
    sig_bh = sum(1 for v in results.values() if v['p_bh'] < 0.05)
    hrs = [v['HR'] for v in results.values() if v['p_bh'] < 0.05]
    print(f'  [{label}] significant: {sig_n}/{len(results)} nominal, '
          f'{sig_bh}/{len(results)} BH; protective HR range of significant '
          f'sigs: {min(hrs):.3f}-{max(hrs):.3f}' if hrs else
          f'  [{label}] no BH-significant signature', flush=True)
    return {'n': n, 'events': events, 'n_signatures': len(results),
            'n_significant_nominal': sig_n, 'n_significant_bh': sig_bh,
            'results': results}


def main():
    raw = REPO / 'data' / 'tcga'
    refetch = raw / 'SKCM' / 'refetch'
    if (refetch / 'expression.tsv.gz').exists():
        print('loading refetched data (API; validated by the all-tumour '
              'arm against the archived tcga_cox_all.json SKCM values)',
              flush=True)
        expr = pd.read_csv(refetch / 'expression.tsv.gz', sep='\t',
                           index_col=0)
        clin = pd.read_csv(refetch / 'clinical.tsv', sep='\t', index_col=0)
    else:
        expr = pd.read_csv(raw / 'SKCM' / 'raw' / 'data_mrna_seq_v2_rsem.txt',
                           sep='\t', comment='#', index_col=0)
        if 'Entrez_Gene_Id' in expr.columns:
            expr = expr.drop(columns=['Entrez_Gene_Id'])
        clin = pd.read_csv(raw / 'SKCM' / 'raw' /
                           'data_clinical_patient.txt',
                           sep='\t', comment='#', index_col=0)
    non_barcode = [c for c in expr.columns
                   if not TCGA_BARCODE.match(str(c).strip())]
    assert not non_barcode, f'non-barcode columns: {non_barcode[:3]}'

    code = {c: c.split('-')[3][:2] for c in expr.columns}
    arms = {
        'all_tumour_v33': expr,
        'metastatic_only': expr[[c for c in expr.columns
                                 if code[c] == '06']],
        'primary_only': expr[[c for c in expr.columns
                              if code[c] == '01']],
    }
    out = {'cohort': 'TCGA-SKCM',
           'note': 'sample-type decomposition of the v33 all-tumour policy; '
                   'identical 12-signature Cox + BH pipeline '
                   '(scripts/run_tcga_cox.py); answers whether the cancer-'
                   'type gradient confounds metastatic/primary composition',
           'sample_counts_raw': {k: int(v.shape[1]) for k, v in arms.items()},
           'arms': {}}
    for label, mat in arms.items():
        print(f'===== {label} =====', flush=True)
        out['arms'][label] = run_arm(mat, clin, label)

    dst = raw / 'SKCM' / 'SKCM_cox_sampletype_sensitivity.json'
    json.dump(out, open(dst, 'w'), indent=1)
    print('WROTE', dst, flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
