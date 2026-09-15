"""Map the Riaz 2017 matrices from Entrez Gene IDs to HGNC symbols.

Why this exists (v33 finding): the Riaz 2017 expression matrices store
var_names as Entrez Gene IDs ('1', '10', '100', ...), so every symbol-based
signature method (GEP, IMPRES, PD-L1/CD274) silently falls back to a constant
score and returns AUROC = 0.500 exactly. In the archived v28 results these
cells were reported as "non-learning methods are completely insensitive
(AUROC 0.500, Delta = 0.000)" - that claim is an artefact of unmapped gene
identifiers, not a biological finding. v33 fixes it by mapping the matrix.

Mapping source: HGNC complete set (Approved symbol <-> NCBI Gene ID),
downloaded to scripts/_hgnc_entrez.txt.

Outputs:
  data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad          (43 x ~20k)
  data/cohorts/Riaz_2017_RECIST/processed/Riaz_2017_RECIST_v33_HGNC.h5ad
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO = Path(__file__).resolve().parents[1]
MAPFILE = REPO / "scripts" / "_hgnc_entrez.txt"


def load_map():
    df = pd.read_csv(MAPFILE, sep="\t", dtype=str)
    m = {}
    for sym, eg in zip(df["Approved symbol"], df["NCBI Gene ID(supplied by NCBI)"]):
        if pd.isna(sym) or pd.isna(eg):
            continue
        for e in str(eg).split(","):
            e = e.strip()
            if e and e not in m:
                m[e] = sym
    return m


def remap(src: Path, dst: Path, m: dict):
    adata = sc.read_h5ad(src)
    names = [str(v) for v in adata.var_names]
    symbols, idx = [], []
    for i, n in enumerate(names):
        sym = m.get(n)
        if sym:
            symbols.append(sym)
            idx.append(i)
    sub = adata[:, idx].copy()
    sub.var_names = symbols
    sub.var_names_make_unique()
    sub.write(dst)
    print(f"{dst.name}: {adata.n_obs} x {len(idx)} genes "
          f"({len(idx)}/{adata.n_vars} Entrez IDs mapped to symbols)")
    return sub


def main():
    m = load_map()
    print(f"HGNC Entrez->symbol entries: {len(m)}")
    remap(REPO / "data/cohorts/Riaz_2017/processed/Riaz_2017_processed.h5ad",
          REPO / "data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad", m)
    remap(REPO / "data/cohorts/Riaz_2017_RECIST/processed/"
                "Riaz_2017_RECIST_v33.h5ad",
          REPO / "data/cohorts/Riaz_2017_RECIST/processed/"
                 "Riaz_2017_RECIST_v33_HGNC.h5ad", m)
    return 0


if __name__ == "__main__":
    sys.exit(main())
