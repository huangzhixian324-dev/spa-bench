"""Build STAD_PRJEB25780 (gastric ICI) and RCC_Braun_2020 (RCC ICI) cohorts
from tigeR.data SummarizedExperiment .rda files.

Both are pre-processed FPKM/TPM expression matrices with clinical ICI
response annotation, sourced from the TIGER/tigeR R data package.

STAD_PRJEB25780: gastric cancer, anti-PD-1, n=78 (CR 6/PR 15/SD 25/PD 32)
RCC_Braun_2020:  RCC ccRCC, nivolumab arm, n=172 evaluable (CR 1/CRPR 25/PR 13/SD 64/PD 69)

Outputs:
  data/cohorts/STAD_PRJEB25780/processed/STAD_PRJEB25780_HGNC.h5ad + splits
  data/cohorts/RCC_Braun_2020/processed/RCC_Braun_2020_HGNC.h5ad + splits
"""
import hashlib
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
import scipy.sparse as sp

warnings.filterwarnings('ignore')

REPO = Path(__file__).resolve().parents[1]
CLONE = REPO / 'data/cohorts/tigeR_clone/data'
PROC = REPO / 'data/cohorts'
SEED = 42
DECLARED = ('GZMA', 'PRF1')


def parse_se(rda_path):
    """Parse SummarizedExperiment from .rda and return (expr_df, colData_df)."""
    result = rdata.read_rda(str(rda_path))
    obj = list(result.values())[0]

    # colData
    cd = obj.colData
    rownames = [str(x) for x in cd.rownames]
    ld = cd.listData
    cols = {}
    for k in ld.keys():
        k_str = str(k)
        vals = ld[k]
        if hasattr(vals, '__len__') and len(vals) == len(rownames):
            cols[k_str] = [str(v) for v in vals]
    clin = pd.DataFrame(cols, index=rownames)

    # assays -> expression matrix
    assays = obj.assays
    expr = None
    gene_names = None
    if hasattr(assays, 'listData'):
        ld2 = assays.listData
        if isinstance(ld2, dict):
            for ak, av in ld2.items():
                # av should be a matrix-like with dim
                if hasattr(av, 'dim') and len(av.dim) == 2:
                    nrows, ncols = av.dim
                    # rdata may return it as a dict of columns or a flat structure
                    if hasattr(av, 'listData') and hasattr(av.listData, 'keys'):
                        ld3 = av.listData
                        if isinstance(ld3, dict) and len(ld3) == ncols:
                            # each key is a column (sample), each value is a gene vector
                            col_names = list(ld3.keys())
                            # get gene names from the first column's names
                            first = list(ld3.values())[0]
                            if hasattr(first, 'index'):
                                gene_names = [str(x) for x in first.index]
                            elif hasattr(first, 'names') and first.names:
                                gene_names = [str(x) for x in first.names]
                            expr = np.array([list(ld3[c]) for c in col_names]).T
                            print(f'  assay "{ak}": {expr.shape} (genes x samples -> {nrows}x{ncols})')
                            break
                    elif hasattr(av, 'data') and hasattr(av.data, 'dict'):
                        pass
                break
    return expr, gene_names, clin


def build_cohort(name, rda_path, therapy_filter, therapy_label,
                 endpoint_description):
    print(f'\n{"="*60}\n{name}\n{"="*60}')
    expr, gene_names, clin = parse_se(rda_path)

    if expr is None:
        print('  [WARN] expression matrix not parsed; trying direct matrix')
        return None

    # map clinical
    resp = clin['response'].values if 'response' in clin.columns else None
    therapy_col = clin['Therapy'].values if 'Therapy' in clin.columns else None
    resp_nr = clin['response_NR'].values if 'response_NR' in clin.columns else None

    # filter
    keep_mask = np.ones(len(clin), dtype=bool)
    if therapy_filter:
        keep_mask &= np.isin(therapy_col, therapy_filter)
    # exclude NE
    keep_mask &= np.isin(resp, ('CR', 'PR', 'SD', 'PD', 'CRPR'))

    idx = np.where(keep_mask)[0]
    expr_k = expr[:, idx]
    clin_k = clin.iloc[idx].reset_index(drop=True)
    resp_k = resp[idx]
    n = len(idx)
    print(f'  after filter: n={n}')

    # response binary: CR/PR/CRPR = 1, SD/PD = 0
    y_recist = np.array([1 if r in ('CR', 'PR', 'CRPR') else 0 for r in resp_k])
    n_res = int(y_recist.sum())
    print(f'  responders: {n_res}/{n} ({100*n_res/n:.0f}%)')

    # log2 transform (already FPKM/TPM, apply log2(x+1))
    X = np.log2(np.maximum(expr_k, 0) + 1.0)
    X = X.astype(np.float32)
    print(f'  log2(FPKM+1): {X.shape}')

    # deduplicate gene symbols (mean)
    genes = [str(g) for g in gene_names]
    df = pd.DataFrame(X.T, columns=pd.Index(genes))
    agg_df = df.T.groupby(level=0).mean().T
    uniq = agg_df.columns.tolist()
    G = len(uniq)
    Xp = agg_df.values.astype(np.float32)
    print(f'  unique symbols: {G}')

    # DCB6mo proxy: no PFS in these tigeR datasets -> skip
    # cytolytic surrogate
    gi = {g: i for i, g in enumerate(uniq)}
    cyt_genes = [g for g in DECLARED if g in gi]
    cyt_score = Xp[:, [gi[g] for g in cyt_genes]].mean(axis=1) if len(cyt_genes) == 2 else None
    y_cyt = (cyt_score > np.median(cyt_score)).astype(int) if cyt_score is not None else None
    cyt_high = int(y_cyt.sum()) if y_cyt is not None else 0
    print(f'  cytolytic-high: {cyt_high}; declared genes present: {cyt_genes}')

    # h5ad
    import anndata as ad
    obs = pd.DataFrame({
        'response': y_recist,
        'responder_label': resp_k,
        'therapy': [str(t) for t in (therapy_col[idx] if therapy_col is not None else [''] * n)],
        'cancer_type': [name] * n,
    }, index=[f'{name}_{i}' for i in range(n)])
    adata = ad.AnnData(X=Xp, obs=obs)
    adata.var_names = uniq
    adata.uns['provenance'] = f'tigeR.data {name}; log2(FPKM+1); duplicate symbols mean-aggregated'

    # splits
    k = min(5, n_res, n - n_res)
    rng = np.random.RandomState(SEED)
    fold_of = np.full(n, -1)
    for cls in np.unique(y_recist):
        idx_cls = np.where(y_recist == cls)[0]
        rng.shuffle(idx_cls)
        for pos, ii in enumerate(idx_cls):
            fold_of[ii] = pos % k
    folds = [{'train': [int(i) for i in np.where(fold_of != f)[0]],
              'test': [int(i) for i in np.where(fold_of == f)[0]]} for f in range(k)]
    payload = {'cohort': name, 'k': k, 'seed': SEED,
               'endpoint': f'RECIST (responders {n_res}/{n})',
               'folds': folds}
    payload['splits_sha1_16'] = hashlib.sha1(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]

    out_dir = PROC / name / 'processed'
    out_dir.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_dir / f'{name}_HGNC.h5ad')
    json.dump(payload, open(out_dir / f'{name}_splits.json', 'w'), indent=1)
    print(f'  WROTE {out_dir}')

    return {'n': n, 'n_responders': n_res, 'n_symbols': G, 'k': k,
            'splits_sha1_16': payload['splits_sha1_16'],
            'cyt_genes': cyt_genes, 'cyt_high': cyt_high}


def main():
    results = {}

    # 1) STAD gastric
    r1 = build_cohort(
        'STAD_PRJEB25780',
        CLONE / 'STAD_PRJEB25780.rda',
        therapy_filter=None,  # all anti-PD-1
        therapy_label='anti-PD-1',
        endpoint_description='RECIST (CR/PR/SD/PD from tigeR annotation)')
    results['STAD_PRJEB25780'] = r1

    # 2) RCC Braun (anti-PD-1 only, exclude everolimus + NE)
    r2 = build_cohort(
        'RCC_Braun_2020',
        CLONE / 'RCC_Braun_2020.rda',
        therapy_filter=['anti-PD-1'],
        therapy_label='nivolumab (anti-PD-1)',
        endpoint_description='RECIST (CR/CRPR/PR vs SD/PD from tigeR annotation)')
    results['RCC_Braun_2020'] = r2

    print('\n=== SUMMARY ===')
    for k, v in results.items():
        if v:
            print(f'  {k}: n={v["n"]} responders={v["n_responders"]} '
                  f'genes={v["n_symbols"]} k={v["k"]} sha={v["splits_sha1_16"]}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
