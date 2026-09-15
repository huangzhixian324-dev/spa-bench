"""TCGA Cox regression: immune signatures vs overall survival (Table S15).

v33 (2026-09-02): fixed IndentationError at line 54; p-values now reported in
scientific notation (no more round-to-zero); Benjamini-Hochberg FDR added
across the 12 signatures within each cancer type; paths resolved relative to
the repository root. Reproduces results into data/tcga/<CT>/<CT>_cox.json.
"""
import pandas as pd, numpy as np, json, os, sys, re
from pathlib import Path
from lifelines import CoxPHFitter

REPO = Path(__file__).resolve().parents[1]

# v33 fix: sample-type policy per cancer type. TCGA-SKCM is dominated by
# metastatic biopsies (367 of 444 samples are '-06'), so restricting to
# primary tumours ('01') silently changed n from 426 to 76 and produced a
# different result from the archived analysis. Melanoma therefore uses all
# tumour samples (primary + metastatic, as in the archived run); the other
# cohorts use primary tumours only.
SAMPLE_TYPE_POLICY = {
    "SKCM": None,   # None = all tumour samples
    "BRCA": "01",
    "LUAD": "01",
    "COAD": "01",
}
TCGA_BARCODE = re.compile(r"^TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")

# Per-cancer input files. SKCM/BRCA/LUAD use the original cBioPortal
# pan-can-atlas files (unchanged -> identical v33 results). COAD is
# rebuilt from the GDC current index (STAR counts, primary tumours;
# scripts/build_coad_gdc.py) after the original non-TCGA raw files forced
# the v33 withdrawal.
COHORT_FILES = {
    "SKCM": ("data_mrna_seq_v2_rsem.txt", "data_clinical_patient.txt"),
    "BRCA": ("data_mrna_seq_v2_rsem.txt", "data_clinical_patient.txt"),
    "LUAD": ("data_mrna_seq_v2_rsem.txt", "data_clinical_patient.txt"),
    "COAD": ("data_mrna_gdc_star_fpkm.txt", "data_clinical_patient_gdc.txt"),
}

SIGS = {
  'IMPRES-like': ['PDCD1','CD27','CTLA4','CD28','CD86','CD80','CD274','LAG3','HAVCR2','TIGIT','TNFRSF9','ICOS'],
  'CD8 T cells': ['CD8A','CD8B'], 'B cells': ['CD19','CD79A','MS4A1'],
  'NK cells': ['NKG7','KLRD1','KLRF1'], 'IFNG response': ['IFNG','CXCL10','CXCL9','IDO1','STAT1'],
  'Cytolytic': ['GZMA','PRF1','GNLY'], 'Checkpoint': ['PDCD1','CTLA4','LAG3','TIGIT','HAVCR2'],
  'TLS': ['CXCL13','CCL19','CCL21'], 'Treg': ['FOXP3','IL2RA'],
  'M1 Macrophage': ['NOS2','IL12A','TNF'], 'Myeloid': ['CD14','CD33','ITGAM'],
  'Stromal': ['COL1A1','COL1A2','COL3A1','ACTA2','FAP'],
}

def process_cohort(name, raw_dir):
  """Load expression + clinical, compute signatures, run Cox."""
  expr_fn, clin_fn = COHORT_FILES.get(
      name, ('data_mrna_seq_v2_rsem.txt', 'data_clinical_patient.txt'))
  expr_path = os.path.join(raw_dir, name, 'raw', expr_fn)
  clin_path = os.path.join(raw_dir, name, 'raw', clin_fn)
  print(f'\n===== {name} =====')
  expr = pd.read_csv(expr_path, sep='\t', comment='#', index_col=0)
  if 'Entrez_Gene_Id' in expr.columns:
      expr = expr.drop(columns=['Entrez_Gene_Id'])
  clin = pd.read_csv(clin_path, sep='\t', comment='#', index_col=0)
  # Map samples to patients
  sp = {}
  for c in expr.columns:
      parts = c.split('-')
      sp[c] = '-'.join(parts[:3]) if len(parts) >= 3 else c
  # v33 fix: verify that expression columns are TCGA barcodes before
  # splitting them. The local COAD file uses non-TCGA identifiers
  # ("01CO001", 107 samples) and cannot support the archived COAD analysis
  # (n = 588); previously this silently produced degenerate fits (HR = 0).
  non_barcode = [c for c in expr.columns
                 if not TCGA_BARCODE.match(str(c).strip())]
  if non_barcode:
      print(f'  [SKIP] {name}: {len(non_barcode)}/{len(expr.columns)} '
            f'expression columns are not TCGA barcodes '
            f'(e.g. {non_barcode[:3]}). Source data is incompatible with '
            f'this pipeline; re-download from cBioPortal.')
      return {}
  sample_type = SAMPLE_TYPE_POLICY.get(name, "01")
  # Keep tumour samples according to the per-cancer sample-type policy
  if sample_type is not None:
      primary = [c for c in expr.columns
                 if len(c.split('-')) > 3 and c.split('-')[3].startswith(sample_type)]
      if primary: expr = expr[primary]
  # Deduplicate patients
  seen, cols = {}, []
  for c in expr.columns:
      pid = sp[c]
      if pid not in seen: seen[pid] = c; cols.append(c)
  expr = expr[cols]
  # Match clinical
  pids = [sp[c] for c in expr.columns]
  clin_m = clin[clin.index.isin(pids)]
  # Order expression columns to match clinical
  pid_to_col = {sp[c]: c for c in expr.columns}
  matched_cols = [pid_to_col.get(pid) for pid in clin_m.index if pid in pid_to_col]
  if not matched_cols:
      print(f'  No matched patients!')
      return {}
  expr_m = expr[matched_cols]
  n = len(matched_cols)
  events = int(clin_m['OS_STATUS'].str.upper().str.contains('DECEASED').sum())
  print(f'  n={n} events={events}')
  # Compute signatures
  X = np.log2(np.clip(expr_m.values, 0, None) + 1)
  # Filter to valid gene symbols (remove rows where Hugo_Symbol is NaN)
  genes = [str(g).upper() for g in expr_m.index]
  gidx = {g: i for i, g in enumerate(genes)}
  results = {}
  for sname, sgenes in SIGS.items():
      idx = [gidx[g] for g in sgenes if g.upper() in gidx]
      if len(idx) < 2: continue
      sig_vals = np.mean(X[idx, :], axis=0)
      # v33 fix: coerce clinical fields to numeric and drop non-finite rows.
      # Non-numeric OS_MONTHS values ('[Not Available]') previously reached
      # lifelines and raised "cannot convert float NaN to integer".
      dur = pd.to_numeric(clin_m['OS_MONTHS'], errors='coerce').values
      event = (clin_m['OS_STATUS'].astype(str).str.upper()
               .str.contains('DECEASED').astype(int).values)
      df = pd.DataFrame({'dur': dur, 'event': event, 'sig': sig_vals})
      df = df.replace([np.inf, -np.inf], np.nan).dropna()
      df = df[(df['dur'] > 0) & (df['event'].isin([0, 1]))].copy()
      df['event'] = df['event'].astype(int)
      if df['event'].sum() < 5: continue
      try:
          cph = CoxPHFitter()
          cph.fit(df[['dur','event','sig']], duration_col='dur', event_col='event')
          results[sname] = {'HR': round(float(np.exp(cph.params_['sig'])), 3),
                            'p': float(cph.summary.loc['sig','p']),
                            'p_sci': '%.3e' % float(cph.summary.loc['sig','p']),
                            'n': n, 'events': events}
      except Exception as e:
          results[sname] = {'HR': 0, 'p': 1.0, 'p_sci': '1.000e+00', 'n': n, 'events': events, 'error': str(e)}
  # Benjamini-Hochberg FDR across the signatures tested in this cancer type
  tested = [s for s, v in results.items() if 'error' not in v]
  if tested:
      ps = np.array([results[s]['p'] for s in tested])
      order = np.argsort(ps)
      m = len(ps)
      bh = np.empty(m)
      bh[order] = np.minimum.accumulate((ps[order] * m / (np.arange(m) + 1))[::-1])[::-1]
      for s, q in zip(tested, bh):
          q = float(min(q, 1.0))
          results[s]['p_bh'] = q
          results[s]['p_bh_sci'] = '%.3e' % q  # keep precision for tiny values
  return results

def main():
  raw = REPO / 'data' / 'tcga'
  all_res = {}
  for ct in ['SKCM','BRCA','LUAD','COAD']:
      r = process_cohort(ct, raw)
      if r:
          all_res[ct] = r
          json.dump(r, open(os.path.join(raw, ct, f'{ct}_cox.json'), 'w'), indent=2)
          # Print summary (nominal p<0.05 and BH-adjusted p<0.05)
          sig_n = sum(1 for v in r.values() if v.get('p',1) < 0.05)
          sig_bh = sum(1 for v in r.values() if v.get('p_bh',1) < 0.05)
          print(f'  Significant (p<0.05): {sig_n}/{len(r)}; BH-adjusted: {sig_bh}/{len(r)}')
          for sname, vals in sorted(r.items(), key=lambda x: x[1].get('p',1)):
              pval = vals.get('p', 1)
              hrval = vals['HR']
              if pval < 0.001: sig = '***'
              elif pval < 0.01: sig = '**'
              elif pval < 0.05: sig = '*'
              else: sig = ''
              print('    %-18s HR=%.3f p=%s BH=%s %s' % (sname, hrval, vals.get('p_sci','%.3e'%pval), vals.get('p_bh_sci','-'), sig))
  json.dump(all_res, open(os.path.join(raw, 'tcga_cox_all.json'), 'w'), indent=2)
  print('\n===== SUMMARY =====')
  for ct in ['SKCM','BRCA','LUAD','COAD']:
      if ct in all_res:
          n_sig = sum(1 for v in all_res[ct].values() if v.get('p',1) < 0.05)
          n_bh = sum(1 for v in all_res[ct].values() if v.get('p_bh',1) < 0.05)
          n_total = len(all_res[ct])
          print(f'{ct}: {n_sig}/{n_total} significant (nominal), {n_bh}/{n_total} (BH-adjusted)')

if __name__ == '__main__':
  main()
