#!/usr/bin/env python3
"""
Quick validation test for SPATBench framework.

Creates synthetic data and runs all non-spatial baseline methods
to verify the pipeline works end-to-end before downloading real data.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

# Set random seed
np.random.seed(42)

# Create synthetic bulk RNA-seq data mimicking an immunotherapy cohort
N_SAMPLES = 200
N_GENES = 1000

print("=" * 60)
print("SPATBench Validation Test")
print("=" * 60)

# --- 1. Create synthetic data ---
print("\n[1/4] Creating synthetic data...")

# Response labels (balanced)
response = np.random.binomial(1, 0.4, N_SAMPLES)

# Gene expression with differential signal for 50 "immune genes"
X = np.random.lognormal(mean=4, sigma=1.5, size=(N_SAMPLES, N_GENES))
immune_genes = np.arange(50)
X[response == 1, :50] *= 1.5  # Responders: higher immune expression

# Create gene names
gene_names = [f"GENE_{i}" for i in range(N_GENES)]
# Make first 10 genes "known" immune genes
immune_gene_names = ["CD8A", "GZMA", "GZMB", "IFNG", "CXCL9",
                      "CXCL10", "PRF1", "TBX21", "CD274", "PDCD1"]
for i, name in enumerate(immune_gene_names):
    gene_names[i] = name

# Create AnnData
adata = sc.AnnData(X=X)
adata.var_names = gene_names
adata.obs["response"] = response
adata.obs["PDL1_CPS"] = X[:, 8] * 0.1 + np.random.normal(0, 2, N_SAMPLES)
adata.obs["TMB"] = np.random.exponential(10, N_SAMPLES)

# Add TIDE-relevant genes (v33 fix: deterministic placement at indices
# 100-108 instead of the previous meaningless random-overwrite logic that
# also consumed RNG state; removed "N_GENES - 1 - ..." dead branch)
tide_genes = ["HAVCR2", "TIGIT", "LAG3", "CTLA4", "VEGFA", "TGFB1",
              "CCL2", "CXCL12", "CXCL8"]
for k, g in enumerate(tide_genes):
    gene_names[100 + k] = g

adata.var_names = gene_names

# Save
synth_dir = Path("data/cohorts/SYNTH_TEST/processed")
synth_dir.mkdir(parents=True, exist_ok=True)
adata.write(synth_dir / "SYNTH_TEST_processed.h5ad")

# Splits
import json
n_folds = 5
folds = []
for fold in range(n_folds):
    n_test = N_SAMPLES // n_folds
    test_idx = list(range(fold * n_test, (fold + 1) * n_test))
    train_idx = [i for i in range(N_SAMPLES) if i not in test_idx]
    folds.append({"fold": fold, "train": train_idx, "test": test_idx})

with open(synth_dir / "SYNTH_TEST_splits.json", "w") as f:
    json.dump({"folds": folds, "n_samples": N_SAMPLES, "n_folds": n_folds}, f)

print(f"   Created: {N_SAMPLES} samples × {N_GENES} genes")
print(f"   Responders: {response.sum()} ({response.mean()*100:.1f}%)")

# --- 2. Import and run non-spatial methods ---
print("\n[2/4] Running non-spatial baseline methods...")

import sys
sys.path.insert(0, "workflows/methods")
sys.path.insert(0, "workflows/evaluation")

from base_method import get_method
from metrics import evaluate_all, delong_test

methods = ["PD_L1_IHC", "TMB", "GEP", "TIDE", "IMPRES", "ElasticNet", "XGBoost"]
results = {}

for method_name in methods:
    print(f"   {method_name}...", end=" ")
    try:
        method = get_method(method_name)
        # 5-fold CV
        all_y_true = []
        all_y_pred = []
        for fold_data in folds:
            train_idx = fold_data["train"]
            test_idx = fold_data["test"]
            preds, timing = method.fit_predict(adata, train_idx, test_idx)
            all_y_true.extend(response[test_idx])
            all_y_pred.extend(preds)

        y_true = np.array(all_y_true).astype(int)
        y_pred = np.array(all_y_pred)

        eval_results = evaluate_all(y_true, y_pred, cohort_name="SYNTH_TEST")
        results[method_name] = eval_results
        print(f"AUROC={eval_results['auroc']:.4f} [{timing['fit_time']:.1f}s+{timing['predict_time']:.1f}s]")
    except Exception as e:
        print(f"ERROR: {e}")

# --- 3. Summarize ---
print("\n[3/4] Summary:")

summary = pd.DataFrame([
    {"Method": m, "AUROC": r["auroc"], "AUPRC": r["auprc"],
     "F1": r["f1"], "Sensitivity": r["sensitivity"],
     "Specificity": r["specificity"], "Brier": r["brier_score"]}
    for m, r in results.items()
]).sort_values("AUROC", ascending=False)

print(summary.to_string(index=False))

# Save results
import json
output_dir = Path("results/benchmark")
output_dir.mkdir(parents=True, exist_ok=True)

def convert(obj):
    if isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, dict):
        return {k: convert(v) for k, v in obj.items()}
    return obj

with open(output_dir / "benchmark_results_VALIDATION.json", "w") as f:
    json.dump(convert(results), f, indent=2)

# DeLong test: best vs second-best
best = summary.iloc[0]["Method"]
second = summary.iloc[1]["Method"]
y_pred_best = np.array([])
y_pred_second = np.array([])

for method_name in [best, second]:
    method = get_method(method_name)
    preds_list = []
    for fold_data in folds:
        preds, _ = method.fit_predict(adata, fold_data["train"], fold_data["test"])
        preds_list.append(preds)
    all_preds = np.concatenate(preds_list)
    if method_name == best:
        y_pred_best = all_preds
    else:
        y_pred_second = all_preds

y_true_all = response

dl = delong_test(y_true_all, y_pred_best, y_pred_second)
print(f"\n   DeLong test: {best} vs {second}")
print(f"   ΔAUROC = {dl['delta_auroc']:.4f}, p = {dl['p_value']:.4f}")

# --- 4. TLS test ---
print("\n[4/4] Testing TLS-MAP...")

# Create synthetic spatial data
N_SPOTS = 500
coords = np.random.rand(N_SPOTS, 2) * 2000  # 2000μm × 2000μm
cell_types = [
    'B_cell', 'Plasma', 'GC_B', 'CD8_T', 'CD4_T', 'Tfh', 'Treg',
    'Tumor', 'CAF', 'Endothelial', 'Myeloid', 'DC', 'NK'
]
cell_props = np.random.dirichlet(np.ones(len(cell_types)), N_SPOTS)

# Create TLS-like clusters
for cluster_center in [(500, 500), (1500, 1500), (800, 1200)]:
    cx, cy = cluster_center
    for i in range(N_SPOTS):
        dist = np.sqrt((coords[i, 0] - cx)**2 + (coords[i, 1] - cy)**2)
        if dist < 150:
            cell_props[i, 0] += 0.3  # B cell
            cell_props[i, 3] += 0.2  # CD8 T
            cell_props[i, 5] += 0.1  # Tfh

# Normalize
cell_props = cell_props / cell_props.sum(axis=1, keepdims=True)

# Create AnnData
# v33 fix: previously var_names was assigned cell_types twice and the
# functional genes (AICDA/CXCL13/MKI67) were only stored in obsm, so
# TLSMAPScorer._score_functionality could never find them and silently
# fell back to default 0.3 — the functionality dimension was never tested.
# Now the three genes are appended as extra columns so they are actually
# scored (they add ~0 proportion mass and do not affect B/T detection).
func_genes = ['AICDA', 'CXCL13', 'MKI67']
func_expr = np.random.lognormal(mean=2, sigma=1, size=(N_SPOTS, len(func_genes)))
all_vars = cell_types + func_genes
X_spatial = np.hstack([cell_props, func_expr / 100.0])  # small added mass
spatial_adata = sc.AnnData(X=X_spatial)
spatial_adata.var_names = all_vars
spatial_adata.obsm['spatial'] = coords
for gi, gene in enumerate(func_genes):
    spatial_adata.obsm[f'{gene}_expr'] = func_expr[:, gi]

spatial_adata.write(synth_dir.parent / "SYNTH_SPATIAL_processed.h5ad")
# For graph building test
(Path("data/spatial/SYNTH/SYNTH_SPATIAL/processed")).mkdir(parents=True, exist_ok=True)

# Test TLS detection
sys.path.insert(0, "workflows/tls")
from tls_map import TLSDetector, TLSMAPScorer

detector = TLSDetector(eps=100, min_cells=20)
tls_list = detector.detect(spatial_adata)
print(f"   Detected {len(tls_list)} TLSs in synthetic data")

scorer = TLSMAPScorer(spatial_adata)
for tls in tls_list:
    scorer.score(tls)

print(f"   TLS maturity levels: {[(t.tls_id, t.maturity_level, f'{t.tls_map_composite:.3f}') for t in tls_list]}")

# --- Done ---
print("\n" + "=" * 60)
print("VALIDATION PASSED: All components working correctly")
print("=" * 60)
print(f"   Benchmarks: {len(methods)} methods evaluated")
print(f"   TLS detection: {len(tls_list)} TLSs found")
print(f"   Results saved to: results/benchmark/")
