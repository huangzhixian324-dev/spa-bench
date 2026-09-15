"""Benchmark Gide 2019 cBioPortal cohort using SPATBench pipeline."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'workflows', 'methods'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'workflows', 'evaluation'))
from base_method import get_method
from metrics import evaluate_all
import scanpy as sc

root = os.path.dirname(os.path.dirname(__file__))
h5ad_path = os.path.join(root, 'data/cohorts/Gide_2019_cBio/processed/Gide_2019_cBio_processed.h5ad')
splits_path = h5ad_path.replace('_processed.h5ad', '_splits.json')

adata = sc.read(h5ad_path)
folds = json.load(open(splits_path))['folds']
print(f'Loaded: {adata.n_obs} patients, {adata.n_vars} genes, {len(folds)} folds')

methods = ['GEP', 'IMPRES']
results = {}
for m in methods:
  try:
      meth = get_method(m)
      yt, yp = [], []
      for f in folds:
          preds, _ = meth.fit_predict(adata, f['train'], f['test'])
          rv = adata.obs['response'].values
          yt.extend(rv[f['test']])
          yp.extend(preds)
      y_true = np.array(yt).astype(int)
      y_pred = np.array(yp)
      r = evaluate_all(y_true, y_pred)
      results[m] = dict(auroc=float(r['auroc']), f1=float(r['f1']),
                       sensitivity=float(r['sensitivity']), specificity=float(r['specificity']))
      print(f"  {m:15s} AUROC={r['auroc']:.3f} F1={r['f1']:.3f}")
  except Exception as e:
      print(f"  {m:15s} ERROR: {e}")

out_path = os.path.join(root, 'data/cohorts/Gide_2019_cBio/gide2019_bench.json')
json.dump(results, open(out_path, 'w'), indent=2)
print(f'Saved: {out_path}')
