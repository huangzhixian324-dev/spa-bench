# SPATBench Data Inventory

**Date:** 2026-09-14 | **Manuscript:** v34 (source) / v35 (Cell Systems submission) | **Status:** Pre-submission

## Scope

SPATBench is a benchmark-audit study of transcriptomic immunotherapy-response
prediction. Six prediction methods (GEP, TIDE, IMPRES, PD-L1/CD274,
ElasticNet-Var, ElasticNet-MI) are evaluated across five immunotherapy-treated
cohorts (**196 patients**; melanoma, NSCLC, bladder), plus two external
replication cohorts and a TCGA pan-cancer survival validation.

This document records what is present in this repository, what must be
re-downloaded from public archives, and which reported numbers each file backs.

## Benchmark cohorts (5 cohorts, 196 patients)

| Cohort | Accession | Cancer | n | Endpoint | Responders | k |
|--------|-----------|--------|---|----------|------------|---|
| Hugo 2016 | GSE78220 | Melanoma | 28 | RECIST v1.1 | 13 | 5 |
| Riaz 2017 | GSE91061 | Melanoma | 43 (42 with RECIST) | RECIST v1.1 / cytolytic median split | 9 RECIST / 21 cytolytic | 5 / 2 |
| Lauss 2017 (ACT) | GSE100797 | Melanoma | 25 | RECIST v1.1 | 10 | 5 |
| Jung 2019 | GSE135222 | NSCLC | 27 | DCB vs NDB | 6 DCB | 5 |
| Gide 2019 | mel_iatlas_gide_2019 (cBioPortal) | Melanoma | 73 | RECIST v1.1 | 40 | 5 |

Processed expression matrices and split files (patient-stratified; seed 42):

| Cohort | Matrix in repo | Size | Splits | Benchmark JSON |
|--------|----------------|------|--------|----------------|
| Hugo 2016 | `data/cohorts/Hugo_2016/processed/Hugo_2016_processed.h5ad` | 6.4 MB | `Hugo_2016_splits.json` | `cohort_Hugo_2016.json` |
| Riaz 2017 (cytolytic) | `data/cohorts/Riaz_2017/processed/Riaz_2017_HGNC.h5ad` | 8.0 MB | `Riaz_2017_splits.json` | `cohort_Riaz_2017_cytolytic.json` |
| Riaz 2017 (RECIST) | `data/cohorts/Riaz_2017_RECIST/processed/Riaz_2017_RECIST_v33_HGNC.h5ad` | 7.8 MB | `Riaz_2017_RECIST_v33_splits.json` (k = 2/3/5) | `cohort_Riaz_2017_RECIST_v33.json` |
| Lauss 2017 | `data/cohorts/Nathanson_2017/processed/Nathanson_2017_processed.h5ad` | 4.2 MB | `Nathanson_2017_splits.json` | `cohort_Nathanson_2017.json` |
| Jung 2019 | `data/cohorts/Jung_2019/processed/Jung_2019_HGNC_processed.h5ad` | 14.6 MB | `Jung_2019_HGNC_splits.json` | `cohort_Jung_2019_DCB.json` |
| Gide 2019 | `data/cohorts/Gide_2019_cBio/processed/Gide_2019_cBio_processed.h5ad` | 24.2 MB | `Gide_2019_cBio_splits.json` | `cohort_Gide_2019_cBio.json` |

All `cohort_*.json` files live in `results/benchmark/v33/`.

### Naming and provenance notes

1. **`Nathanson_2017` is a legacy directory alias for Lauss 2017 (GSE100797).**
   "Nathanson" was an attribution error in earlier versions; the key is retained
   only so archived filenames and machine-readable table rows remain stable. The
   manuscript refers to this cohort as Lauss throughout (Revision Transparency,
   Supplementary Note S7 item 11). It is adoptive TIL therapy, not checkpoint
   blockade, and is included as a deliberate context-dependence stress test.
2. **HGNC mapping.** The analysis inputs are HGNC-mapped. Riaz 2017 was mapped
   Entrez→HGNC (21,769 / 22,187 genes; `scripts/map_riaz_to_hgnc.py`, using the
   bundled `scripts/_hgnc_entrez.txt`); Jung 2019 was mapped Ensembl→HGNC
   (`scripts/rebuild_jung.py`, cache `scripts/_jung_ensembl_map.json`). The
   un-mapped `Riaz_2017_processed.h5ad` / `Jung_2019_processed.h5ad` are retained
   as the pre-mapping intermediates.
3. **Riaz 2017 provides two response definitions on the same patients** (RECIST
   v1.1, n = 42; cytolytic median split, n = 43), enabling the within-patient
   endpoint switch that is the study's headline contrast. The RECIST endpoint
   uses k = 2 (only 9 responders); k = 3 and k = 5 are archived as robustness
   checks (`results/benchmark/v33/riaz_recist_kfold_v33.json`).
4. **Jung 2019 labels are corrected.** The original release coded 1 = NDB; all
   analyses here use 1 = DCB (Revision Transparency, item 2).

## External replication cohorts

| Cohort | Source | n | Result file |
|--------|--------|---|-------------|
| Liu 2019 | cBioPortal `mel_iatlas_liu_2019` | 122 | `results/benchmark/v33/e1_liu2019_replication.json` |
| IMvigor210 | cBioPortal `blca_iatlas_imvigor210_2017` | 297 | `results/benchmark/v33/e1_imvigor210_replication.json` |
| GSE274975 | GEO (disqualified) | — | `scripts/e1_gse274975_replication.py` |

Both replications are pre-declared and reproduce the cytolytic→RECIST collapse on
identical patients and expression matrices (Supplementary Tables S19, S20).
GSE274975 was the pre-declared NSCLC candidate; it was investigated over FTP and
disqualified (no RECIST, PFS or ICI-treatment annotation). Its raw SOFT and series
matrix files (`data/cohorts/GSE274975/raw/`) and its complete pre-declared
pipeline are retained for reuse on any suitable NSCLC ICI cohort; Liu 2019 is the
disclosed substitution.

## TCGA validation (12 immune signatures × 4 cancer types)

- `data/tcga/clinical_{skcm,brca,luad,coad}.txt` — GDC clinical/survival tables
- `data/tcga/sample_meta_{skcm,brca,luad,coad}.txt` — sample metadata
- `data/tcga/{SKCM,BRCA,LUAD,COAD,COADREAD}_cox.json` — per-cancer Cox models
- `data/tcga/tcga_cox_all.json` — consolidated BH-adjusted results
- Analysis: `scripts/run_tcga_cox.py` (lifelines CoxPHFitter)
- Reported in Supplementary Table S11 and Supplementary Figure S11

## External clinical annotation tables

`data/external/` holds the cBioPortal/iAtlas clinical tables used for the
replication and survival analyses: `mel_iatlas_gide_2019_api/` (expression,
clinical and sample JSON; `expression_compact.json` is 60.7 MB),
`mel_iatlas_hugo_ucla_2016/clinical.tsv`, `mel_iatlas_liu_2019/clinical.tsv`,
`blca_iatlas_imvigor210_2017/clinical.tsv`,
`skcm_vanderbilt_mskcc_2015/clinical.tsv`.

## Machine-readable results and figures

- `results/benchmark/v33/` — `benchmark_v33.json`, `permutation_v33.json` (the
  complete 36-cell method × cohort-endpoint matrix), per-cell checkpoint files in
  `nc_cells/`, `power_table_v33.json`, `bh_sensitivity_v33.json`,
  `mi_var_overlap_v33.json`, `hugo_survival_true.json`, `cds_*.json`,
  `effect_size_bootstrap_v33.json`, `perm_obs_auroc_v33.json`
- `results/run_registry.json` — provenance for each reported number (script, seed
  policy, git commit, split-file digests)
- `results/figures/v33/` — every main and supplementary figure at 300 dpi
  (plus TIFF for submission and the graphical abstract)
- `docs/` — manuscript source (v34; typeset v35 Cell Systems), supplementary
  material, submission PDF and Word file, highlights/eTOC

## Not included in this repository

- Raw FASTQ / per-sample GEO archives — re-download from the accessions above
  (`scripts/download_cohorts.py` and `scripts/fetch_cbioportal_cohorts.py`)
- TCGA raw RNA-seq — re-download from the GDC Data Portal
- Real PD-L1 IHC scores (public metadata provides none; CD274 expression is used
  as a proxy) and TMB mutation calls (unavailable; TMB is a constant 0.5)
- XGBoost benchmark cells — recorded as excluded, with evidence in
  `results/benchmark/v33/xgb_exclusion_evidence_v33.json`

## Reproducing the reported numbers

```bash
pip install -r requirements.txt
python scripts/run_all_checks.py     # 12 checks; exit code 0 = all passed
python scripts/rerun_v33.py          # regenerates the benchmark and permutation tables
```

Results are deterministic under seed 42 (k-fold stability over seeds 42–51 in
Supplementary Table S12). `scripts/validate.py` is a full end-to-end self-test on
a synthetic cohort and runs from a fresh checkout with no external data.
