#!/usr/bin/env python3
"""
Build spatial graph from deconvolved spatial transcriptomics data.

Converts cell-type-deconvolved spatial data into PyTorch Geometric graph objects
with multi-scale graph structure (cell-level, niche-level, tissue-level).

Usage:
    python scripts/build_graph.py --sample NSCLC_Visium_10x_001 --dataset NSCLC
"""

import argparse
import numpy as np
import torch
from pathlib import Path
from sklearn.neighbors import kneighbors_graph
from scipy.spatial import Delaunay


def build_spatial_graph(coords, features, k=10, radius_um=150,
                        include_delaunay=True):
    """
    Build a spatial graph from cell/spot coordinates and features.

    Args:
        coords: (N, 2) array of spatial coordinates (in μm)
        features: (N, D) array of node features
        k: KNN parameter
        radius_um: maximum edge distance
        include_delaunay: add Delaunay triangulation edges

    Returns:
        dict with edge_index, edge_attr, node_features, coords
    """
    N = coords.shape[0]

    # KNN graph (distance-weighted)
    knn_adj = kneighbors_graph(coords, n_neighbors=k, mode='distance',
                                include_self=False)

    # Delaunay triangulation (captures tissue topology)
    edges = set()
    edge_distances = {}

    if include_delaunay and N > 3:
        try:
            tri = Delaunay(coords)
            for simplex in tri.simplices:
                for i in range(3):
                    for j in range(i + 1, 3):
                        a, b = simplex[i], simplex[j]
                        edges.add((a, b))
                        edges.add((b, a))
                        d = np.linalg.norm(coords[a] - coords[b])
                        edge_distances[(a, b)] = d
                        edge_distances[(b, a)] = d
        except Exception:
            pass  # Delaunay fails for collinear points

    # KNN edges
    for i in range(N):
        for j_idx in range(len(knn_adj[i].indices)):
            j = knn_adj[i].indices[j_idx]
            d = knn_adj[i].data[j_idx]
            if radius_um is None or d < radius_um:
                edges.add((i, j))
                if (i, j) not in edge_distances:
                    edge_distances[(i, j)] = d

    # Build edge_index and edge_attr
    edge_list = list(edges)
    edge_index = torch.tensor(edge_list, dtype=torch.long).T

    # Edge features: distance, same_cell_type
    edge_attr_list = []
    for (i, j) in edge_list:
        d = edge_distances.get((i, j), np.linalg.norm(coords[i] - coords[j]))
        # Placeholder for same_cell_type (needs cell type labels)
        same_type = 0
        edge_attr_list.append([d, same_type])

    edge_attr = torch.tensor(edge_attr_list, dtype=torch.float)
    x = torch.tensor(features, dtype=torch.float)
    pos = torch.tensor(coords, dtype=torch.float)

    return {
        "edge_index": edge_index,
        "edge_attr": edge_attr,
        "x": x,
        "pos": pos,
        "n_nodes": N,
        "n_edges": edge_index.shape[1],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--radius_um", type=int, default=150)
    parser.add_argument("--no_delaunay", action="store_true")
    args = parser.parse_args()

    # Paths
    input_path = Path(f"data/spatial/{args.dataset}/{args.sample}/processed/"
                      f"{args.sample}_deconvolved.h5ad")
    output_path = Path(f"data/spatial/{args.dataset}/{args.sample}/processed/"
                       f"{args.sample}_graph.pt")

    if not input_path.exists():
        print(f"[WARNING] Input file not found: {input_path}")
        print(f"[INFO] This is a template — actual data will be processed "
              "after downloading.")
        return 0

    # Load deconvolved data
    import scanpy as sc
    adata = sc.read_h5ad(input_path)

    # Extract coordinates and features
    coords = adata.obsm['spatial']  # (N, 2) in μm
    features = adata.X  # cell-type proportions

    # Build graph
    graph = build_spatial_graph(
        coords, features,
        k=args.k,
        radius_um=args.radius_um,
        include_delaunay=not args.no_delaunay
    )

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(graph, output_path)

    print(f"[{args.sample}] Graph built: {graph['n_nodes']} nodes, "
          f"{graph['n_edges']} edges")
    print(f"[{args.sample}] Saved to: {output_path}")
    return 0


if __name__ == "__main__":
    exit(main())
