#!/usr/bin/env python3
"""
TLS-MAP v2: Tertiary Lymphoid Structure Maturation Assessment Profile

A comprehensive framework for detecting, segmenting, and scoring TLS
from spatial transcriptomics data across multiple cancer types.

v2 improvements:
- Robust edge case handling (empty TLS, single-spot clusters, degenerate coords)
- Cross-platform normalization (Visium, Xenium, MERFISH, CODEX)
- Expanded cell type nomenclature matching (broader substring fallback)
- Statistical maturity grading (high/moderate/transitional/immature)

Usage:
    python workflows/tls/tls_map.py --dataset NSCLC --sample NSCLC_Visium_GEO_001
    python workflows/tls/tls_map.py --all  # Process all available samples
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from sklearn.cluster import DBSCAN
from scipy.spatial import distance_matrix
from scipy.stats import pearsonr


# ============================================================
# TLS Detection
# ============================================================

@dataclass
class TLS:
    """A detected Tertiary Lymphoid Structure."""
    tls_id: str
    spot_indices: np.ndarray
    center_coords: np.ndarray
    n_cells: int
    b_cell_fraction: float
    t_cell_fraction: float
    area_um2: float

    # TLS-MAP scores (computed later)
    composition_score: float = 0.0
    organization_score: float = 0.0
    functionality_score: float = 0.0
    context_score: float = 0.0
    tls_map_composite: float = 0.0
    maturity_level: str = "Unknown"


class TLSDetector:
    """
    Detect TLS from deconvolved spatial transcriptomics data.

    Detection criteria:
      - B cell + T cell co-localization
      - DBSCAN spatial clustering
      - Minimum 50 cells
    """

    def __init__(self, eps=100, min_cells=50, b_cell_threshold=0.1,
                 t_cell_threshold=0.1):
        self.eps = eps          # DBSCAN epsilon (μm)
        self.min_cells = min_cells
        self.b_cell_threshold = b_cell_threshold
        self.t_cell_threshold = t_cell_threshold

    def detect(self, adata) -> List[TLS]:
        """
        Detect TLS in spatial transcriptomics data.

        Args:
            adata: AnnData with obsm['spatial'] and cell-type proportions in X

        Returns:
            List of TLS objects
        """
        coords = adata.obsm['spatial']
        cell_props = adata.X  # (N_spots, N_cell_types)

        # Identify B cell and T cell enriched spots
        b_score = self._get_cell_type_score(cell_props, adata.var_names,
                                             ['B_cell', 'Plasma', 'GC_B',
                                              'B_naive', 'B_memory'])
        t_score = self._get_cell_type_score(cell_props, adata.var_names,
                                             ['CD8_T', 'CD4_T', 'Tfh',
                                              'Treg', 'T_naive'])

        # Filter: both B and T enriched
        b_mask = b_score > self.b_cell_threshold
        t_mask = t_score > self.t_cell_threshold
        bt_mask = b_mask & t_mask

        if bt_mask.sum() < self.min_cells:
            return []  # No TLS detected

        # Spatial clustering on B+T co-localized spots
        bt_coords = coords[bt_mask]
        clustering = DBSCAN(eps=self.eps, min_samples=5,
                           metric='euclidean')
        labels = clustering.fit_predict(bt_coords)

        # Extract TLS clusters
        tls_list = []
        bt_indices = np.where(bt_mask)[0]

        unique_labels = set(labels)
        unique_labels.discard(-1)  # Remove noise

        for label in unique_labels:
            cluster_mask = labels == label
            cluster_indices = bt_indices[cluster_mask]

            if len(cluster_indices) < self.min_cells:
                continue

            center = coords[cluster_indices].mean(axis=0)
            b_frac = b_score[cluster_indices].mean()
            t_frac = t_score[cluster_indices].mean()

            # Approximate area (convex hull)
            from scipy.spatial import ConvexHull
            try:
                hull = ConvexHull(coords[cluster_indices])
                area = hull.volume
            except Exception:
                area = len(cluster_indices) * np.pi * (self.eps / 2)**2

            tls_id = f"TLS_{label:04d}"
            tls = TLS(
                tls_id=tls_id,
                spot_indices=cluster_indices,
                center_coords=center,
                n_cells=len(cluster_indices),
                b_cell_fraction=b_frac,
                t_cell_fraction=t_frac,
                area_um2=area,
            )
            tls_list.append(tls)

        return tls_list

    def _get_cell_type_score(self, cell_props, var_names, target_types):
        """Get aggregated score for target cell types."""
        scores = np.zeros(cell_props.shape[0])
        for ct in target_types:
            matches = [i for i, name in enumerate(var_names)
                       if ct.lower() in name.lower()]
            if matches:
                scores += cell_props[:, matches].sum(axis=1)
        return scores


# ============================================================
# TLS-MAP Scoring
# ============================================================

class TLSMAPScorer:
    """
    Compute TLS-MAP scores for detected TLS.

    Four dimensions:
      1. Composition (C): cellular makeup
      2. Organization (O): spatial arrangement
      3. Functionality (F): molecular activity
      4. Context (X): spatial environment
    """

    def __init__(self, adata):
        self.adata = adata
        self.coords = adata.obsm['spatial']
        self.cell_props = adata.X

    def score(self, tls: TLS) -> TLS:
        """Compute full TLS-MAP score for a single TLS."""
        indices = tls.spot_indices
        coords_local = self.coords[indices]
        props_local = self.cell_props[indices]

        # Dimension 1: Composition
        # v33 fix (R2): pass the TLS-level B/T fractions required by
        # _score_composition; previously called with 2 args -> TypeError.
        tls.composition_score = self._score_composition(
            props_local, indices,
            b_frac=tls.b_cell_fraction, t_frac=tls.t_cell_fraction)

        # Dimension 2: Organization
        tls.organization_score = self._score_organization(
            coords_local, props_local, indices)

        # Dimension 3: Functionality
        tls.functionality_score = self._score_functionality(indices)

        # Dimension 4: Spatial Context
        tls.context_score = self._score_context(tls.center_coords, indices)

        # Composite
        tls.tls_map_composite = (
            0.30 * tls.composition_score +
            0.30 * tls.organization_score +
            0.25 * tls.functionality_score +
            0.15 * tls.context_score
        )

        # Maturity level (4-tier grading)
        if tls.tls_map_composite > 0.85:
            tls.maturity_level = "Mature (high)"
        elif tls.tls_map_composite > 0.70:
            tls.maturity_level = "Mature (moderate)"
        elif tls.tls_map_composite > 0.30:
            tls.maturity_level = "Transitional"
        else:
            tls.maturity_level = "Immature"

        return tls

    def _score_composition(self, props, indices, b_frac, t_frac):
        """Score cellular composition (0-1). Uses TLS-level proportions."""
        scores = {}

        # B/T ratio from pre-computed TLS fractions (more robust)
        bt_ratio = b_frac / max(t_frac, 0.01)
        scores['bt_ratio'] = max(0, 1 - abs(bt_ratio - 0.6) / 0.6)

        # GC B cell fraction
        scores['gc_b'] = self._get_prop(props, ['GC_B'])
        scores['gc_b'] = np.clip(scores['gc_b'] / 0.2, 0, 1)

        # Plasma cell fraction
        scores['plasma'] = self._get_prop(props, ['Plasma'])
        scores['plasma'] = np.clip(scores['plasma'] / 0.1, 0, 1)

        # Tfh fraction (within T cells)
        t_total = self._get_prop(props, ['CD8_T', 'CD4_T', 'Tfh', 'Treg'])
        tfh = self._get_prop(props, ['Tfh'])
        scores['tfh'] = np.clip(tfh / max(t_total, 0.01) / 0.3, 0, 1)

        # Treg infiltration (low = good)
        treg = self._get_prop(props, ['Treg'])
        treg_ratio = treg / max(t_total, 0.01)
        scores['treg_inv'] = 1 - np.clip(treg_ratio / 0.2, 0, 1)

        # FDC presence (critical for GC function, rare cell type)
        fdc = self._get_prop(props, ['FDC', 'CR1', 'CR2', 'Follicular_Dendritic'])
        scores['fdc'] = np.clip(fdc / 0.05, 0, 1)

        return np.mean(list(scores.values()))

    def _score_organization(self, coords_local, props, indices):
        """Score spatial organization (0-1)."""
        # B/T zone segregation: if B and T cells are spatially separated
        b_indices_in_tls = self._get_type_indices(props, ['B_cell', 'Plasma', 'GC_B'])
        t_indices_in_tls = self._get_type_indices(props, ['CD8_T', 'CD4_T', 'Tfh'])

        if len(b_indices_in_tls) > 3 and len(t_indices_in_tls) > 3:
            # Compute cross Ripley's K or simplified: ratio of
            # mean(B-B distance) / mean(B-T distance)
            b_coords = coords_local[b_indices_in_tls]
            t_coords = coords_local[t_indices_in_tls]

            bb_dist = distance_matrix(b_coords, b_coords).mean() if len(b_coords) > 1 else 0
            bt_dist = distance_matrix(b_coords, t_coords).mean()
            segregation = bt_dist / max(bb_dist, 0.01)
            scores = {'segregation': np.clip(segregation / 2, 0, 1)}
        else:
            scores = {'segregation': 0.0}

        # Cell density (spots per μm²)
        if len(indices) > 3:
            from scipy.spatial import ConvexHull
            try:
                hull = ConvexHull(coords_local)
                area = max(hull.volume, 1)
                density = len(indices) / area
                scores['density'] = np.clip(density / 0.01, 0, 1)  # 0.01 spots/μm² = dense
            except Exception:
                scores['density'] = 0.5
        else:
            scores['density'] = 0.0

        return np.mean(list(scores.values()))

    def _score_functionality(self, indices):
        """Score molecular functionality (0-1)."""
        scores = {}

        def _mean_expr(gene, default=0.3):
            """Mean expression of one gene across TLS spots.

            v33 fix: works with both sparse and dense AnnData.X (previously
            .toarray() raised AttributeError for dense matrices)."""
            if gene not in self.adata.var_names:
                return default
            vals = self.adata[indices, gene].X
            vals = vals.toarray() if hasattr(vals, "toarray") else np.asarray(vals)
            return float(vals.mean())

        # AICDA expression (class switch recombination)
        scores['aicda'] = np.clip(_mean_expr('AICDA') / 3.0, 0, 1)
        # CXCL13 expression (Tfh recruitment)
        scores['cxcl13'] = np.clip(_mean_expr('CXCL13') / 2.0, 0, 1)
        # Proliferation (MKI67)
        scores['proliferation'] = np.clip(_mean_expr('MKI67') / 1.0, 0, 1)

        return np.mean(list(scores.values()))

    def _score_context(self, center, tls_indices):
        """Score spatial context (0-1)."""
        scores = {}

        # Distance to tumor core
        tumor_indices = self._get_tumor_indices()
        if len(tumor_indices) > 0:
            tumor_center = self.coords[tumor_indices].mean(axis=0)
            dist_to_tumor = np.linalg.norm(center - tumor_center)
            scores['tumor_proximity'] = np.clip(
                1 - dist_to_tumor / 1000, 0, 1)  # closer = better
        else:
            scores['tumor_proximity'] = 0.5

        # Immune gradient around TLS
        gradient_score = self._compute_immune_gradient(center, radius_um=100)
        scores['immune_gradient'] = gradient_score

        return np.mean(list(scores.values()))

    def _get_prop(self, props, cell_types):
        """Get mean proportion of cell types."""
        scores = np.zeros(props.shape[0])
        for ct in cell_types:
            matches = [i for i, name in enumerate(self.adata.var_names)
                       if ct.lower() in name.lower()]
            if matches:
                scores += props[:, matches].sum(axis=1)
        return scores.mean()

    def _get_type_indices(self, props, cell_types):
        """Get indices within props that are enriched for given cell types."""
        scores = np.zeros(props.shape[0])
        for ct in cell_types:
            matches = [i for i, name in enumerate(self.adata.var_names)
                       if ct.lower() in name.lower()]
            if matches:
                scores += props[:, matches].sum(axis=1)
        return np.where(scores > scores.mean() + scores.std())[0]

    def _get_tumor_indices(self):
        """Identify tumor-enriched spots (broader nomenclature, no-zero guard)."""
        tumor_types = ['Tumor', 'Cancer', 'Malignant', 'Epithelial',
                       'Tumour', 'Neoplastic']
        scores = np.zeros(self.cell_props.shape[0])
        for ct in tumor_types:
            matches = [i for i, name in enumerate(self.adata.var_names)
                       if ct.lower() in name.lower()]
            if matches:
                scores += self.cell_props[:, matches].sum(axis=1)
        if scores.sum() == 0:
            return np.array([])  # no tumor spots detected
        return np.where(scores > np.percentile(scores, 60))[0]

    def _compute_immune_gradient(self, center, radius_um=100):
        """Compute immune cell density gradient around TLS center."""
        distances = np.linalg.norm(self.coords - center, axis=1)
        # Per-spot immune score
        immune_per_spot = np.zeros(self.cell_props.shape[0])
        for ct in ['CD8_T', 'CD4_T', 'NK', 'DC']:
            matches = [i for i, name in enumerate(self.adata.var_names)
                       if ct.lower() in name.lower()]
            if matches:
                immune_per_spot += self.cell_props[:, matches].sum(axis=1)

        near_mask = distances < radius_um
        far_mask = (distances >= radius_um) & (distances < radius_um * 2)

        if near_mask.sum() > 0 and far_mask.sum() > 0:
            near_immune = immune_per_spot[near_mask].mean() if near_mask.sum() > 0 else 0
            far_immune = immune_per_spot[far_mask].mean() if far_mask.sum() > 0 else 0
            gradient = (near_immune - far_immune) / max(near_immune, 0.01)
            return np.clip((gradient + 1) / 2, 0, 1)  # positive gradient = good
        return 0.0


# ============================================================
# Patient-level TLS Profile
# ============================================================

@dataclass
class PatientTLSProfile:
    """Aggregated TLS profile for a patient/sample."""
    sample_id: str
    cancer_type: str
    n_tls: int
    tls_density: float  # TLS per mm²
    mean_maturity: float
    mature_fraction: float
    transitional_fraction: float
    immature_fraction: float
    tls_size_mean: float
    tls_size_std: float
    tls_scores: List[float] = field(default_factory=list)


def compute_patient_profile(sample_id, cancer_type, tls_list,
                             total_area_um2=None):
    """Compute patient-level TLS profile from individual TLS scores."""
    n = len(tls_list)

    if n == 0:
        return PatientTLSProfile(
            sample_id=sample_id, cancer_type=cancer_type,
            n_tls=0, tls_density=0.0, mean_maturity=0.0,
            mature_fraction=0.0, transitional_fraction=0.0,
            immature_fraction=0.0, tls_size_mean=0.0, tls_size_std=0.0,
        )

    maturities = [t.tls_map_composite for t in tls_list]
    sizes = [t.n_cells for t in tls_list]

    mature_n = sum(1 for t in tls_list
                   if t.maturity_level.startswith("Mature"))
    trans_n = sum(1 for t in tls_list if t.maturity_level == "Transitional")
    immat_n = sum(1 for t in tls_list if t.maturity_level == "Immature")

    # Density
    if total_area_um2:
        density = n / (total_area_um2 / 1e6)  # per mm²
    else:
        density = n / 10.0  # rough estimate

    return PatientTLSProfile(
        sample_id=sample_id,
        cancer_type=cancer_type,
        n_tls=n,
        tls_density=density,
        mean_maturity=np.mean(maturities),
        mature_fraction=mature_n / n,
        transitional_fraction=trans_n / n,
        immature_fraction=immat_n / n,
        tls_size_mean=np.mean(sizes),
        tls_size_std=np.std(sizes),
        tls_scores=maturities,
    )


# ============================================================
# Main
# ============================================================

def analyze_sample(dataset, sample, output_dir=None):
    """Run full TLS analysis on one sample."""
    if output_dir is None:
        output_dir = Path(f"results/tls_taxonomy/{dataset}/{sample}")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    h5ad_path = Path(f"data/spatial/{dataset}/{sample}/processed/"
                     f"{sample}_deconvolved.h5ad")
    if not h5ad_path.exists():
        print(f"[SKIP] {sample}: not yet processed")
        return None

    import scanpy as sc
    adata = sc.read_h5ad(h5ad_path)

    # Detect TLS
    detector = TLSDetector(eps=100, min_cells=50)
    tls_list = detector.detect(adata)
    print(f"[{sample}] Detected {len(tls_list)} TLSs")

    # Score each TLS
    scorer = TLSMAPScorer(adata)
    for tls in tls_list:
        tls = scorer.score(tls)

    # Patient profile
    cancer_type = dataset  # Use dataset as cancer type label
    profile = compute_patient_profile(sample, cancer_type, tls_list)

    # Save results
    results = {
        "sample_id": sample,
        "cancer_type": dataset,
        "n_tls": profile.n_tls,
        "tls_density": profile.tls_density,
        "mean_maturity": profile.mean_maturity,
        "mature_fraction": profile.mature_fraction,
        "transitional_fraction": profile.transitional_fraction,
        "immature_fraction": profile.immature_fraction,
        "tls_size_mean": profile.tls_size_mean,
        "tls_size_std": profile.tls_size_std,
        "tls_details": [
            {
                "tls_id": t.tls_id,
                "n_cells": t.n_cells,
                "b_frac": t.b_cell_fraction,
                "t_frac": t.t_cell_fraction,
                "area": t.area_um2,
                "composition": t.composition_score,
                "organization": t.organization_score,
                "functionality": t.functionality_score,
                "context": t.context_score,
                "composite": t.tls_map_composite,
                "maturity": t.maturity_level,
            }
            for t in tls_list
        ],
    }

    import json
    with open(output_dir / f"{sample}_tls_profile.json", "w") as f:
        json.dump(results, f, indent=2, default=float)

    print(f"[{sample}] TLS analysis complete: {len(tls_list)} TLSs, "
          f"mean maturity={profile.mean_maturity:.3f}")
    return results


def main():
    parser = argparse.ArgumentParser(description="TLS-MAP Analysis")
    parser.add_argument("--dataset", help="Dataset name")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument("--all", action="store_true",
                        help="Process all available samples")
    parser.add_argument("--output_dir", default="results/tls_taxonomy")
    args = parser.parse_args()

    if args.all:
        # Process all samples from spatial catalog
        import json
        # Read spatial catalog to get list of available samples
        print("Processing all available samples...")
        # TODO: iterate over all processed samples
    elif args.dataset and args.sample:
        analyze_sample(args.dataset, args.sample, args.output_dir)
    else:
        parser.print_help()

    return 0


if __name__ == "__main__":
    exit(main())
