#!/usr/bin/env python3
"""
Preprocess a bulk RNA-seq immunotherapy cohort.

Steps:
  1. Load raw counts
  2. TPM normalization
  3. log2(TPM + 1) transform
  4. Gene symbol standardization
  5. QC filtering
  6. Save processed h5ad

Usage:
    python scripts/preprocess_bulk.py --cohort IMvigor210
"""

import argparse
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--raw_dir", default=None)
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--min_genes_per_sample", type=int, default=200)
    parser.add_argument("--min_samples_per_gene", type=int, default=10)
    parser.add_argument("--gene_type", default="protein_coding")
    args = parser.parse_args()

    # Paths
    base_dir = Path(f"data/cohorts/{args.cohort}")
    raw_dir = Path(args.raw_dir) if args.raw_dir else base_dir / "raw"
    output_dir = Path(args.output_dir) if args.output_dir else base_dir / "processed"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{args.cohort}] Loading raw counts from {raw_dir}...")

    # TODO: Adapt loading to cohort-specific format
    # GEO datasets: read GEO matrix files
    # EGA datasets: read after download with pyega3
    # This is a template - actual loading depends on the download format

    # Placeholder: Load count matrix + clinical data
    # counts = pd.read_csv(raw_dir / "counts.csv", index_col=0)
    # clinical = pd.read_csv(raw_dir / "clinical.csv", index_col=0)

    # --- Step 1: TPM normalization ---
    # gene_lengths = load_gene_lengths()
    # tpm = counts_to_tpm(counts, gene_lengths)

    # --- Step 2: log2(TPM + 1) ---
    # log_tpm = np.log2(tpm + 1)

    # --- Step 3: Gene symbol standardization ---
    # log_tpm = standardize_gene_symbols(log_tpm)

    # --- Step 4: QC ---
    # Filter lowly expressed genes and samples
    # gene_mask = (log_tpm > 0).sum(axis=1) >= args.min_samples_per_gene
    # sample_mask = (log_tpm > 0).sum(axis=0) >= args.min_genes_per_sample

    # --- Step 5: Create AnnData ---
    # adata = sc.AnnData(X=log_tpm.T)
    # adata.obs = clinical.loc[adata.obs_names]

    # --- Step 6: Save ---
    # adata.write(output_dir / f"{args.cohort}_processed.h5ad")

    print(f"[{args.cohort}] Preprocessing complete. Output: {output_dir}")
    return 0


if __name__ == "__main__":
    exit(main())
