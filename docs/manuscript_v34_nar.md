# SPATBench: Methodological Choices Can Outweigh Algorithmic Choices in Transcriptomic Immunotherapy Prediction Benchmarking at Current Sample Sizes

**Zhixian Huang** (黄志贤)<sup>1</sup>

<sup>1</sup> Putian University, Putian, Fujian, China

**Lead contact:** Zhixian Huang, huangzhixian324@gmail.com (correspondence)

**Author affiliations:** <sup>1</sup> Putian University, Putian, Fujian, China

---

## SUMMARY

Transcriptomic predictors of immune checkpoint inhibitor (ICI) response are benchmarked on small, incompatible cohorts with inconsistent endpoint definitions, leaving unclear whether conclusions reflect algorithms or methodology. SPATBench, an open-source framework, isolates three design choices across six methods and five ICI-treated cohorts (196 patients). Switching from a molecular surrogate (cytolytic activity) to a clinical endpoint (RECIST v1.1) on the same patients collapses every method whose features overlap the endpoint-defining genes (ElasticNet-Var 0.985 to 0.465; GEP 0.968 to 0.503; IMPRES unchanged), 3.5x the largest method gap. Var- and MI-based feature selection yield 27%-overlapping gene sets at indistinguishable performance. The Circularity Detection Score (CDS) places the circular endpoint far above its cohort-matched null while its categorical labels prove non-specific. No trainable method survives FDR correction on clinical endpoints at n <= 43, whereas at n = 73 (Gide 2019) ElasticNet-MI is significant (p = 0.0232). SPATBench and CDS are open-source and Docker-reproducible.

**Keywords:** cancer immunotherapy; biomarker benchmarking; response definition; feature selection; endpoint circularity; transcriptomics; machine learning; immune checkpoint inhibitors

**Uncertainty on the effect-size comparison.** The 3.5x ratio compares point estimates from the tuning track (Delta = 0.520 against a same-endpoint method gap of 0.148). Under the frozen-hyperparameter pipeline on the full 42-patient RECIST-evaluable set, the collapse is Delta = 0.343 (patient-level bootstrap 95% CI 0.097 to 0.584, 1,000 resamples, seed 42) against a frozen-track method gap of 0.029 (CI 0.023 to 0.300); the frozen-track ratio (11.98, CI 0.46 to 14.32) is larger in point estimate but its interval is wide because the denominator (the same-endpoint method gap) is itself small and volatile at this n; the ratio is therefore reported as descriptive only. The inferential anchor is the Delta interval, which excludes zero; the tuning-track 3.5x is the conservative reading of the effect. (`results/benchmark/v33/effect_size_bootstrap_v33.json`)

---

## INTRODUCTION

Immune checkpoint inhibitors (ICIs) targeting PD-1, PD-L1, and CTLA-4 have transformed oncology, achieving durable responses in melanoma, non-small cell lung cancer (NSCLC), and other malignancies [1,2]. However, only 20–30% of unselected patients achieve objective responses, motivating intensive efforts to develop transcriptomic biomarkers for patient stratification [3–6]. 
Multiple gene-expression-based prediction methods have been proposed: the T-effector gene expression profile (GEP) [3], TIDE [4], IMPRES [5], and a growing number of machine learning (ML) models [6–8]. Each was evaluated on different cohorts with different preprocessing, response definitions, and evaluation protocols. This fragmentation has produced a literature in which nearly every method appears superior — but only within its own developmental context [9,10]. 
Three structural features of the field make method rankings unreliable: (i) publicly available immunotherapy cohorts with pre-treatment RNA-seq are small (typically n < 50), limiting statistical power [11–13]. (ii) response can be defined via clinical endpoints (RECIST v1.1), molecular surrogates (cytolytic activity [14]), or progression-based endpoints (PFS) — each capturing different underlying biology And (iii) methodological choices such as feature selection strategy are rarely treated as experimental variables. 
Several recent efforts have benchmarked immunotherapy prediction methods across multiple cohorts [15–18]; Supplementary Table S2 compares their published AUROCs with the present benchmark. Liang et al. (2026) found limited cross-cohort generalizability for nine predictors [15]. COMPASS (2026) demonstrated cross-cancer improvements on 1,133 patients [16]. The ICLR 2026 benchmark evaluated 35 tasks and found domain-specific models outperforming general-purpose ones [17]. EXPRESSO benchmarked 20 signatures across 3,729 patients [18]. These efforts share a common design: they hold methodology fixed and compare methods. Our approach is orthogonal — we hold methods fixed and systematically vary methodological choices, asking which experimental design decisions, rather than which algorithms, determine benchmark outcomes. 
Here we present SPATBench, a computational framework and benchmark resource designed to answer: *given that rankings are unreliable at available sample sizes, which methodological choices contribute most to that unreliability, and by how much?* Our contributions are:

1. An open-source benchmark framework enforcing patient-stratified cross-validation and standardized preprocessing across six prediction methods and five cohorts.
2. Empirical isolation and quantification of three design choices — response definition, feature selection, and endpoint circularity — with effect-size comparisons against method-level variation.
3. A Circularity Detection Score (CDS) for pre-benchmarking endpoint assessment — a null-referenced continuous diagnostic whose categorical labels are non-specific among structured endpoints and are advisory only, the validated deliverable being a reporting convention whose binding rule is the declaration of endpoint-defining genes, with CDS as its quantitative companion rather than an independent gate — evaluated on 5 real endpoints, 5 synthetic scenarios, 60 gene-program pseudo-endpoints, one retrospective external consistency check (Auslander et al. 2018), and a retrospective impact assessment on published benchmark studies. 4. A complete permutation test matrix with high-resolution validation for borderline cases, demonstrating that no trainable method is distinguishable from random label shuffling on any clinical endpoint at n ≤ 43 after FDR correction.
5. A power-analysis-derived checklist of minimum standards for reliable benchmarking. All evidence is retrospective on public cohorts; the response-definition result is demonstrated in melanoma and bladder carcinoma, with NSCLC replication an open gap.

---

## RESULTS

### Design Choice 1: Response Definition Exerts a Larger Effect Than Method Choice

The Riaz 2017 cohort provides two response definitions on the same patients with the same expression matrix.

**Table 1. Riaz 2017: RECIST (n = 42) versus cytolytic (n = 43) response definitions (v33, HGNC-mapped).**

| Method | RECIST AUROC (n=42) | Cytolytic AUROC (n=43) | ΔAUROC | Type |
|--------|-------------|-----------------|--------|------|
| ElasticNet (MI) | 0.569 [0.370–0.755] | 0.968 [0.909–1.000] | −0.399 | ML |
| ElasticNet (Var) | 0.465 [0.208–0.704] | 0.985 [0.954–1.000] | **−0.520** | ML |
| TIDE | 0.466 [0.222–0.723] | 0.470 [0.282–0.653] | −0.004 | Signature |
| IMPRES | 0.613 [0.381–0.835] | 0.660 [0.493–0.811] | −0.047 | Signature |
| GEP | 0.503 [0.305–0.710] | 0.968 [0.915–1.000] | −0.465 | Signature |
| PD-L1 (CD274) | 0.545 [0.330–0.760] | 0.868 [0.750–0.960] | −0.323 | Biomarker |

Switching response definitions transforms all ML methods from apparently state-of-the-art (Var 0.985 on cytolytic) to near-chance on RECIST (0.465); GEP collapses with them (0.968 → 0.503) because the cytolytic endpoint is defined by GZMA/PRF1, both members of the GEP signature. IMPRES — whose fifteen checkpoint-gene relations share no genes with the endpoint definition — is nearly insensitive (0.660 → 0.613, Δ = −0.047; permutation p = 0.032 on cytolytic and 0.154 on RECIST, each against its own label-shuffle null). Within the Riaz RECIST endpoint the largest between-method ΔAUROC is 0.147 (IMPRES vs TIDE), so the largest between-endpoint effect (0.520) exceeds it by 3.5× (effect-size comparison framework in Supplementary Table S14). The 0.520 contrast uses the benchmark (seed-42) split; across 10 stratified seeds the k = 2 RECIST mean is 0.389 ± 0.105 (Supplementary Table S16), so the collapse magnitude carries ≈0.1 seed standard deviation while its direction and scale are stable. 
**k-fold robustness:** The Riaz RECIST endpoint uses only 2-fold CV (constrained by 9 responders). To confirm the collapse magnitude is not a k-fold artifact, we re-ran ElasticNet (Var) with k = 3 and k = 5 on the rebuilt Riaz RECIST dataset (n = 42; 10 stratified random seeds per k). All RECIST AUROC values remain far below 0.60 (k = 2: 0.389 ± 0.105; k = 3: 0.397 ± 0.064; k = 5: 0.398 ± 0.088; ElasticNet-Var, 10 stratified seeds per k), confirming that the cytolytic→RECIST collapse is not an artifact of the k = 2 constraint (Supplementary Table S16; machine-readable in results/benchmark/v33/riaz_recist_kfold_v33.json).

**External replication (pre-declared):** The response-definition effect replicates on an independent cohort — Liu 2019 [33] (n = 122 annotated samples; DFCI anti-PD-1 melanoma; independent of all five benchmark cohorts; cBioPortal iAtlas harmonization), the pre-specified substitution for GSE274975, whose pipeline (scripts/e1_gse274975_replication.py) is included in the repository. On identical patients and expression matrix, the cytolytic-surrogate endpoint (GZMA/PRF1 mean, median split — the Riaz construction) yields GEP 0.959 and ElasticNet-MI 0.956, while the clinical RECIST endpoint yields 0.535 and 0.538 (Δ = +0.425 and +0.418); IMPRES — whose features do not include the endpoint genes — retains the most RECIST signal (0.601) and is the only nominally significant fixed scorer on the clinical endpoint (50,000 shuffles, p = 0.027; BH q = 0.108 over the four fixed scorers, not FDR-significant), mirroring the Hugo headline result; no trainable method's bootstrap CI excludes 0.5 on RECIST. CDS flags the surrogate endpoint HIGH (0.949 vs null mean 0.566) and places the clinical endpoint above its null mean (0.721 vs 0.591) — the same genuine-biology signature as Gide (Table S4b). Full results in Supplementary Table S19.

**Cross-cancer replication (pre-declared):** The response-definition effect further replicates in a non-melanoma tumour type — IMvigor210 [19] (metastatic urothelial carcinoma, n = 297 response-annotated samples; anti-PD-L1; cBioPortal iAtlas harmonization). On identical patients, the cytolytic-surrogate endpoint yields GEP 0.967 and ElasticNet-MI 0.951, while the clinical RECIST endpoint yields 0.613 and 0.657 (Δ = +0.354 and +0.294). PD-L1 shows the same pattern (0.863 → 0.566, Δ = +0.297). On the clinical endpoint, GEP is the only fixed scorer surviving FDR (50,000 shuffles, p = 0.0025; BH q = 0.010 over the four fixed scorers) — a weak effect (AUROC 0.613) detectable only at this sample size, consistent with the power analysis. IMvigor210 is an application cohort of the GEP (which Ayers et al. developed on pembrolizumab KEYNOTE-series cohorts [3]), so this significance reflects cross-cohort signal rather than a cohort-of-origin effect. CDS cleanly separates the two endpoints here: the surrogate is HIGH (0.953 vs null mean 0.514) while RECIST falls below its null mean (0.434 vs 0.492) — the cleanest CDS separation to date, confirming that the Gide prognostication conflation is cohort-specific rather than universal. Full results in Supplementary Table S20.

**Mechanism:** The cytolytic score is operationally defined by GZMA and PRF1 expression [14]. Learning methods predict this score from co-expressed immune-infiltration genes; GEP predicts it even better (0.968) because GZMA/PRF1 are members of the GEP signature itself. When the endpoint switches to RECIST — not operationally defined by any gene — models whose signal comes from immune-infiltration co-expression collapse (0.503, 0.465), whereas IMPRES, whose checkpoint-gene relations do not overlap the endpoint genes, retains near-chance performance on both definitions. CDS flagged the cytolytic endpoint as HIGH risk (0.954) precisely because endpoint-defining genes sit inside the feature space of most benchmarked methods.

On the HGNC-mapped Riaz matrix, TIDE is near chance for both endpoints (RECIST 0.466, cytolytic 0.470), consistent with the tidepy implementation checked at source (classification = correlation['TIDE'] < vthres, default vthres = 0).

### Design Choice 2: Feature Selection Determines Biological Interpretation

Within identical ElasticNet models on Hugo 2016 (n = 28, RECIST), we compared variance-based and MI-based feature selection.

**Table 2. ElasticNet on Hugo 2016 (n = 28): MI versus variance feature selection.**

| Metric | MI-based | Variance-based | Interpretation |
|--------|----------|----------------|----------------|
| AUROC (5-fold nested CV) | 0.631 | 0.528 | Not distinguishable (cross-seed Δ −0.008 ± 0.050, Table S12) |
| Gene overlap (top-500) | — | 27.0% (135/500) prescreen convention; 2.4% (12/500) full-data MI | Overlap is convention-dependent (v33; two-stage prescreen is the pipeline's actual selection path — MI is computed within the top-2000 variance prescreen; the full-data MI variant is the alternative convention) |
| Spearman ρ (OOF scores) | — | 0.75 (p < 0.001) | Highly correlated predictions |
| Observed ΔAUROC (MI − Var) | +0.103 (benchmark split) | −0.008 ± 0.050 (10 seeds, Table S12) | Includes zero |
| Permutation p (n=60 full pipeline) | 0.131 | 0.344 | Neither distinguishable from chance |

Under the corrected, leakage-free implementation (v33) the two strategies achieve indistinguishable discriminative performance and highly correlated out-of-fold predictions (Spearman ρ = 0.75; Supplementary Figure S2); their top-500 gene sets overlap by 27% under the pipeline's two-stage prescreen convention (2.4% under the full-data MI convention — the dependence of the overlap figure on the selection convention is itself a design-choice effect; Supplementary Figure S3). Neither variant is distinguishable from chance by permutation (p = 0.13 and 0.34 at n = 28). Feature-selection choice at n = 28 therefore changes the selected gene sets — and with them the biological narrative — without measurably changing predictive performance; the difference remains statistically unobservable at this sample size (cross-seed stability: Supplementary Figure S7; Table S12).
Learning curve analysis (Supplementary Figure S6; Table S13) shows the expected small-sample behaviour under the leakage-free protocol: point estimates remain near chance at every subsample size while variance inflates steeply as n decreases (across-seed SD 0.211 at n = 8 vs. 0.038 at n = 28). Caveat: n = 8 with k = 5 yields ~6 training samples per fold; this point characterizes estimator variance rather than providing a precise performance estimate. 
### Design Choice 3: Endpoint Circularity Can Be Detected Before Benchmarking

**Relation to Design Choice 1.** The response-definition effect quantified above is, mechanistically, an instance of endpoint circularity: the cytolytic surrogate is operationally defined by genes inside the feature space. Design Choice 1 measured the *magnitude* of that effect within patients; Design Choice 3 asks whether such circularity can be recognised *before* any benchmark is run, and develops the diagnostic for it. The two sections therefore analyse one phenomenon from the benchmark-design and the pre-benchmarking-diagnostic perspectives respectively; no non-circular response-definition contrast (e.g., RECIST versus PFS on the same patients) was available in any public cohort, and we flag this as a boundary of the response-definition claim (Limitations).

**Table 3. CDS with cohort-matched random-label nulls (200 replicates; declared genes GZMA/PRF1; cds v1.2.1 — CDS scores identical under v1.1.0; the archived 10-replicate nulls give the same readings and remain in cds_nulls_v33.json). 'ML AUROC Range' = trainable methods only (ElasticNet MI/Var), v33 values. Riaz RECIST is reported on the HGNC matrix (0.409); the archived 0.685 was computed on the pre-v33 Entrez matrix.**

| Endpoint | CDS | null mean | null max | CDS vs null | ML AUROC Range | Inflated? |
|---|---|---|---|---|---|---|
| Riaz 2017 cytolytic | 0.954 | 0.566 | 0.842 | far above null max | 0.968–0.985 | **Yes** |
| Gide 2019 RECIST | 0.852 | 0.555 | 0.800 | above null max | 0.605–0.630 | No |
| Lauss 2017 RECIST (ACT) | 0.781 | 0.578 | 0.866 | above null mean | 0.627–0.700 | No |
| Riaz 2017 RECIST (v33) | 0.409 | 0.560 | 0.836 | within null | 0.465–0.569 | No |
| Hugo 2016 RECIST | 0.449 | 0.586 | 0.834 | within null | 0.528–0.631 | No |

CDS places the only circular endpoint far above its random-label null (Riaz cytolytic 0.954 vs null mean 0.566, max 0.842) — exactly where an endpoint operationally defined by GZMA/PRF1 must sit — and ML methods reach AUROC 0.968–0.985 on it. The clinical endpoints sit at or inside their nulls (Hugo 0.449 vs null mean 0.586; Riaz RECIST 0.409 vs 0.560; Lauss 0.781 vs 0.578, max 0.866) and no trainable-method AUROC exceeds 0.70 on them. Gide 2019 (first stored value, Table S4b) returns 0.852, above its null maximum (0.800): the declared genes genuinely track response in this large cohort (GEP 0.830, q < 0.001), and CDS v1.x cannot separate that genuine biology from circularity. The categorical thresholds therefore do not carry the validated claim. The validated claims are: (a) the null-referenced ordering above, and (b) the procedural rule that any operationally gene-defined surrogate — of any program type (Supplementary Table S18) — is a circular prediction task by construction. 
**Within-study consistency check:** CDS scores for the endpoints with stored values at analysis time (Hugo, Riaz, Lauss) were computed before the corresponding ML benchmarks were run; all sat at or within their random-label nulls (Table 3) and no trainable-method AUROC exceeded 0.70 on them. No CDS value was stored for Gide 2019 at analysis time; the v33 recomputation (0.852, above its null maximum) is reported in Table S4b and discussed as a quantified limitation. Because the scores were stored before the benchmarks, the benchmarks constitute a subsequent verification of the CDS readings — a limited within-study check, not prospective validation. External prospective validation — applying CDS before running benchmarks on cohorts not analyzed in this study — remains an open challenge (see Limitations).

**External consistency check:** The Auslander et al. (2018) study provides an independent validation scenario. IMPRES was trained on Auslander's data and validated externally on Hugo 2016 (RECIST), achieving AUROC = 0.83 in the original paper (0.795 in our benchmark). If CDS had been computed before that validation, it would have returned CDS = 0.449 (below its cohort-matched random-label null mean of 0.586), consistent with the Hugo RECIST endpoint being non-circular and with the observed IMPRES AUROC reflecting biological signal in the checkpoint gene pairs rather than endpoint-gene dependency.

### Complete Permutation Test Matrix

**Table 4. Permutation tests (fixed scorers: 50,000 prediction shuffles; ML variants: full-pipeline label shuffles, per-cell n).**

| Method | Cohort | n | Obs AUROC | perm-obs AUROC | p | n_shuff | BH q | Sig |
|--------|--------|----|-----------|----------------|-----|--------|-------|-----|
| ElasticNet | Hugo_2016 | 28 | 0.631 | 0.631 | 0.1311 | 60 | 0.393 |  |
| ElasticNet | Lauss_2017 | 25 | 0.700 | 0.700 | 0.082 | 60 | 0.123 |  |
| ElasticNet | Gide_2019_cBio | 73 | 0.630 | 0.662 | 0.0232 | 5,000 | 0.028 | **Yes** |
| ElasticNet | Jung_2019_DCB | 27 | 0.476 | 0.564 | 0.1791 | 200 | 0.358 |  |
| ElasticNet | Riaz_2017_cytolytic | 43 | 0.968 | 0.933 | 0.002 | 500 | 0.003 | **Yes** |
| ElasticNet | Riaz_2017_RECIST_v33 | 42 | 0.569 | 0.616 | 0.0945 | 200 | 0.461 |  |
| ElasticNet_Var | Hugo_2016 | 28 | 0.528 | 0.528 | 0.3443 | 60 | 0.642 |  |
| ElasticNet_Var | Lauss_2017 | 25 | 0.627 | 0.627 | 0.036 | 1000 | 0.108 |  |
| ElasticNet_Var | Gide_2019_cBio | 73 | 0.605 | 0.605 | 0.02 | 500 | 0.0278 | **Yes** |
| ElasticNet_Var | Jung_2019_DCB | 27 | 0.476 | 0.524 | 0.2877 | 1000 | 0.432 |  |
| ElasticNet_Var | Riaz_2017_cytolytic | 43 | 0.985 | 0.974 | 0.001 | 1000 | 0.002 | **Yes** |
| ElasticNet_Var | Riaz_2017_RECIST_v33 | 42 | 0.465 | 0.546 | 0.25 | 1000 | 0.500 |  |
| GEP | Hugo_2016 | 28 | 0.446 | 0.6916 | 50000 | 0.755 |  | 
| GEP | Lauss_2017 | 25 | 0.693 | 0.0567 | 50000 | 0.113 |  | 
| GEP | Gide_2019_cBio | 73 | 0.830 | 2e-05 | 50000 | 0.000 | **Yes** |
| GEP | Jung_2019_DCB | 27 | 0.786 | 0.01692 | 50000 | 0.102 |  | 
| GEP | Riaz_2017_cytolytic | 43 | 0.968 | 2e-05 | 50000 | 0.000 | **Yes** |
| GEP | Riaz_2017_RECIST_v33 | 42 | 0.503 | 0.4889 | 50000 | 0.587 |  | 
| IMPRES | Hugo_2016 | 28 | 0.795 | 0.0029 | 50000 | 0.017 | **Yes** |
| IMPRES | Lauss_2017 | 25 | 0.610 | 0.1893 | 50000 | 0.227 |  | 
| IMPRES | Gide_2019_cBio | 73 | 0.626 | 0.03168 | 50000 | 0.038 | **Yes** |
| IMPRES | Jung_2019_DCB | 27 | 0.583 | 0.4046 | 50000 | 0.486 |  | 
| IMPRES | Riaz_2017_cytolytic | 43 | 0.660 | 0.03224 | 50000 | 0.039 | **Yes** |
| IMPRES | Riaz_2017_RECIST_v33 | 42 | 0.613 | 0.1535 | 50000 | 0.461 |  | 
| PD_L1 | Hugo_2016 | 28 | 0.523 | 0.428 | 50000 | 0.642 |  | 
| PD_L1 | Lauss_2017 | 25 | 0.780 | 0.0096 | 50000 | 0.058 |  | 
| PD_L1 | Gide_2019_cBio | 73 | 0.791 | 2e-05 | 50000 | 0.000 | **Yes** |
| PD_L1 | Jung_2019_DCB | 27 | 0.659 | 0.1299 | 50000 | 0.358 |  | 
| PD_L1 | Riaz_2017_cytolytic | 43 | 0.868 | 2e-05 | 50000 | 0.000 | **Yes** |
| PD_L1 | Riaz_2017_RECIST_v33 | 42 | 0.545 | 0.3463 | 50000 | 0.519 |  | 
| TIDE | Hugo_2016 | 28 | 0.426 | 0.7552 | 50000 | 0.755 |  | 
| TIDE | Lauss_2017 | 25 | 0.500 | 0.5103 | 50000 | 0.510 |  | 
| TIDE | Gide_2019_cBio | 73 | 0.752 | 0.00014 | 50000 | 0.000 | **Yes** |
| TIDE | Jung_2019_DCB | 27 | 0.278 | 0.9521 | 50000 | 0.952 |  | 
| TIDE | Riaz_2017_cytolytic | 43 | 0.470 | 0.6365 | 50000 | 0.637 |  | 
| TIDE | Riaz_2017_RECIST_v33 | 42 | 0.466 | 0.6223 | 50000 | 0.622 |  | 

**Table 4 note (revised).** "Perm-obs AUROC" is the observed AUROC of the frozen-hyperparameter pipeline against which each permutation p was computed; "Obs AUROC" is the tuning-track benchmark estimate printed for reference. p-values are therefore conditional on the frozen pipeline (deviations of 0.011-0.087 disclosed in the preceding footnote).

PD-L1 cells use the fixed-scorer protocol (50,000 label shuffles against the fold-normalized CD274 score, which exactly reproduces the benchmark AUROC; scripts/pdl1_perm_v33.py). The seven previously n.c. ElasticNet/ElasticNet_Var cells were completed by scripts/fill_nc_cells.py (per-cell n; seed 42; frozen hyperparameters) and merged by scripts/merge_nc_cells.py; all 36 of 36 method × cohort-endpoint cells now carry a recorded permutation p, and within-cohort BH families (six cells each) were recomputed (complete per-cell results in Supplementary Table S11). Permutation-run observed AUROCs deviate from the benchmark AUROC on six trainable cells by 0.011–0.087 (frozen hyperparameters differ from the tune_primary selections; deviations are recorded in the per-cell JSON files under results/benchmark/v33/nc_cells/). 
Permutation p-values in this table come from two equivalent-in-intent protocols: for fixed scorers (IMPRES/GEP/TIDE/PD-L1) the score is label-independent, so the exact full-pipeline null equals 50,000 label shuffles against the fixed out-of-fold scores; for the trainable ML variants the full fitting pipeline (per-fold selection + fit) is re-run under label shuffles with per-cell n (60–1,000). After within-cohort Benjamini-Hochberg correction (families of six cells per cohort-endpoint), the following survive. On Hugo 2016, IMPRES alone (0.795, q = 0.017; per-method pairwise concordance on Hugo in Supplementary Table S5, and a leave-one-pair-out sensitivity analysis of IMPRES in Supplementary Table S6). On Gide 2019 (n = 73), GEP (0.830, q < 0.001), TIDE (0.752, q < 0.001), PD-L1 (0.791), ElasticNet-Var (0.605, q = 0.030) and IMPRES (q = 0.038) And on the circular Riaz cytolytic endpoint, ElasticNet (0.968, q = 0.003), ElasticNet-Var (0.985, q = 0.002), GEP and PD-L1 (0.868 Both q < 0.001) and IMPRES (q = 0.039) — on a HIGH-circularity endpoint even the methods whose features include the endpoint-defining genes survive, which is the circularity effect itself. The GSE100797 (Lauss) PD-L1 single-gene biomarker (0.780) is nominal-only after the completed six-test family correction (p = 0.0096, q = 0.058), and no other Lauss cell reaches q < 0.10. No trainable method survives FDR correction on any clinical endpoint at n ≤ 43. On Gide (n = 73), the fixed scorers GEP, TIDE and PD-L1 are family-robust, and ElasticNet-Var becomes detectable under the within-cohort family (q = 0.030) though not under the global family (q = 0.066) — the first sample size at which any method's signal becomes detectable, consistent with the power analysis. IMPRES on Hugo carries the Carter et al. (2019) training-set-bias caveat [30]. Family-definition sensitivity (Supplementary Table S21): the Gide ElasticNet-Var, Gide IMPRES and Riaz-cytolytic IMPRES cells lose significance under global or method-class-stratified families, and the Lauss PD-L1 cell gains significance under them; the primary conclusions — no trainable-method FDR survival on clinical endpoints at n ≤ 43, circular-endpoint inflation, IMPRES–Hugo and GEP/TIDE/PD-L1–Gide — hold under all three definitions. Cohort display names in this table match the Methods cohort table; internal dataset keys retain the legacy identifier `Nathanson_2017` for machine readability (see Methods, Cohorts and Revision Transparency). 
### Integrated Multi-Cohort Results

**Table 5. AUROC across four clinical-endpoint cohorts (v33, HGNC-mapped; 95% CIs from 1,000 bootstrap iterations where shown).**

| Method | Hugo 2016 (28) | Lauss 2017 (25) | Gide 2019 (73) | Riaz RECIST (42) |
|--------|----------------|---------------------|----------------|------------------|
| IMPRES | **0.795** [0.597-0.947] | 0.610 | 0.626 | **0.613** |
| GEP | 0.446 | 0.693 | **0.830** | 0.503 |
| TIDE | 0.426 | 0.500 | 0.752 | 0.466 |
| PD-L1 (CD274) | 0.523 | **0.780** | 0.791 | 0.545 |
| ElasticNet (Var) | 0.528 [0.318-0.761] | 0.627 | 0.605 | 0.465 [0.208-0.704] |
| ElasticNet (MI) | 0.631 [0.417-0.857] | 0.700 | 0.630 | 0.569 [0.370-0.755] |

†Permutation status and BH q-values in Table 4. Riaz RECIST n = 42 (rebuilt cohort); AUROC and 95% bootstrap CIs for every method–cohort cell are plotted in Supplementary Figure S8. Bold = best AUROC per column for context only; no method ranking is claimed (all bootstrap CIs of the trainable methods on clinical endpoints overlap 0.5 at n ≤ 43). 
**Table 6. DCB endpoint (Jung 2019, n = 27, NSCLC; corrected labels: 1 = durable clinical benefit, 6/27; GSE135222).**

| Method | AUROC | 95% CI |
|---|---:|---|
| IMPRES | 0.583 | 0.380-0.750 |
| GEP | 0.786 | 0.573-0.960 |
| TIDE | 0.278 | 0.083-0.527 |
| PD-L1 (CD274) | 0.659 | 0.464-0.861 |
| ElasticNet (MI) | 0.476 | 0.152-0.792 |
| ElasticNet (Var) | 0.476 | 0.100-0.818 |

TIDE scores far below chance on this cohort (0.278); with only 6 DCB responders this deviation sits inside the permutation null envelope (p = 0.952), and we report the value as computed under the tidepy default threshold convention rather than re-tuning it to the cohort. 
### Survival Analysis: AUROC Alone Does Not Capture Clinical Utility

**Table 7. Overall survival stratification on Hugo 2016 (true events).**

| Cohort | n | Events | Median OS | IMPRES KM p | RECIST KM p | GEP KM p |
|--------|----|--------|-----------|-------------|-------------|----------|
| Hugo 2016 (true events) | 26 | 12 | 439 days | 0.554 | **0.0003** | 0.475 |

Hugo values come from true OS events (iAtlas OS_STATUS, 26 patients / 12 deaths; matching validated by 24/27 local OS_days agree with iAtlas OS_MONTHS (±60 d); results/benchmark/v33/hugo_survival_true.json). The local GSE100797 (Lauss) matrix contains no OS metadata; the archived Lauss survival row from earlier versions could not be regenerated under the v33 protocol and is therefore excluded from this table (archived values retained in Supplementary Table S9 for continuity). 
RECIST response strongly stratifies OS (overall survival) on Hugo (log-rank p = 0.0003; Supplementary Figure S4). IMPRES — the best AUROC method on Hugo 2016 — does not significantly stratify OS (KM p = 0.554 with true events). Cox regression with z-standardized signatures (HR per 1 SD) confirms: IMPRES HR = 0.631 (95% CI [0.325–1.227], p = 0.175) on Hugo, GEP HR = 1.03 (95% CI [0.562–1.887], p = 0.923). AUROC for binary response prediction does not guarantee clinically meaningful survival stratification. 
Clinical utility metrics for ElasticNet_Var on Hugo 2016 (v33 recomputed; Supplementary Table S9): **DCA** net benefit 0.345 vs treat-all 0.405 at p_t = 0.10 and 0.259 vs 0.330 at 0.20 — positive but consistently below "treat all"; at 0.50 the model is at 0. **NRI vs. random baseline** = -0.113 (p = 0.658, not significant). **IDI** = -0.029 (sensitivity 0.538, specificity 0.533). 
### Power Analysis: Minimum Sample Size Requirements

**One-sample detectability.** Table 8 quantifies two-method comparisons; the Results also make one-sample statements (a method's AUROC exceeding 0.5), which we now calibrate directly using the Hanley-McNeil single-sample variance (normal approximation, one-sided alpha 0.05; `results/benchmark/v33/power_onesample_v33.json`). At the n <= 43 cohort sizes, 80% power requires a true AUROC of 0.72-0.82 (Riaz 0.72; Hugo 0.76; Lauss 0.77; Jung 0.82), whereas at n = 73 (Gide) the requirement falls to 0.68. This supplies the quantitative version of the detectability argument: effect sizes in the 0.60-0.70 range that are genuine at n = 73 are below the detection floor at n <= 43, which is what the permutation matrix shows.

**Table 8. Minimum detectable ΔAUROC at 80% power (α = 0.05, two-sided; Hanley–McNeil paired model, r = 0.75, the most favourable assumption; exact values from `scripts/power_table_v33.py`).**

| n | Responder % | Min Detectable ΔAUROC | Example Cohort |
|----|------------|----------------------|----------------|
| 25 | 40% | 0.23 | Lauss 2017 (ACT) |
| 28 | 46% | 0.21 | Hugo 2016 |
| 42 | 21% | 0.21 | Riaz RECIST |
| 73 | 55% | 0.13 | Gide 2019 |
| 100 | 50% | 0.11 | — |
| 200 | 50% | 0.08 | — |

At n = 25–43, the minimum detectable ΔAUROC is 0.21–0.23. For reliable pairwise comparison at ΔAUROC = 0.10 (80% power), the required sample size is ≈130 under the most favourable correlation assumption (paired, r = 0.75) and ≈260–510 under weaker or no correlation (power_table_v33.json); we therefore report n ≥ 50 below as a minimum practical floor, not as a sufficient condition for adequately powered pairwise comparison. 
**Table 9. Minimum requirements for reliable immunotherapy prediction benchmarking.**

| Requirement | Current Practice | Recommendation |
|-------------|-----------------|----------------|
| Patients per cohort | n = 25–43 | n ≥ 50 minimum floor; ≈130 for pairwise ΔAUROC = 0.10 at 80% power (paired, r = 0.75) |
| Independent cohorts | 1–2 | ≥3 with consistent endpoint definitions |
| Response definitions | Mixed (RECIST, PFS, surrogate) | ≥1 clinical endpoint; surrogates labeled with defining genes |
| Feature selection | Often unreported | ≥2 strategies pre-specified; report overlap and cross-seed stability |
| Circularity check | Not performed | Compute CDS; report continuous value; calibrate thresholds on synthetic data |
| Non-learning baselines | Often omitted | Include IMPRES + GEP as negative controls |
| Permutation test | Rarely performed | Full permutation matrix; FDR correction; ≥5000 permutations for borderline cases† |
| Statistical reporting | AUROC only; p < 0.05 | Exact p-values; 95% CIs; effect sizes; cross-fold SD at n < 50 |
| Reproducibility | Variable | Docker container + Zenodo archived data |

†Our own trainable-method cells use 60–1,000 full-pipeline shuffles because each shuffle re-fits the pipeline (≈60–90 s); the ≥5,000 recommendation applies to fixed scorers and to borderline cases wherever compute permits, and the resolution ceiling this creates at n ≤ 43 is disclosed in Methods. 
**Priority guidance:** Among these recommendations, the three most critical — because they address the largest sources of benchmark variance — are: (1) use at least one clinical endpoint alongside any molecular surrogates; (2) pre-specify at least two feature selection strategies and report their gene set overlap; (3) include non-learning baselines (IMPRES, GEP) as negative controls against circular prediction. The remaining recommendations provide additional rigor but address smaller-magnitude sources of variation.

---

## DISCUSSION

### This Is a Design-Choice Isolation Study, Not a Method Ranking

### A predictive account of when gene-defined surrogates mislead

The three modules of this study compose into one statement. A gene-defined surrogate misleads in proportion to (i) the immune-infiltration-driven co-expression supply that links endpoint genes to the rest of the feature space (high in melanoma, low in COAD-like low-infiltration contexts, where the TCGA gradient attenuates), (ii) the overlap between declared endpoint genes and the feature space (complete for the cytolytic surrogate, absent for IMPRES), and (iii) sample size, which sets both the detectability floor (nothing is detectable below n ≈ 43) and the power to separate genuine signal from construction (detectable from n ≈ 73. Two-method comparisons from n ≈ 130, Table 8). These anchors are stated as a testable claim: in a low-infiltration tumour type the endpoint-switch collapse should shrink toward the method gap even at adequate n. The 36-cell matrix and the CDS nulls are the instruments that make that claim checkable. 
SPATBench asks a different question from existing benchmarks. COMPASS [16], Liang et al. [15], EXPRESSO [18], and the ICLR benchmark [17] all ask "which method performs best?" SPATBench asks: "given that rankings are unreliable at available sample sizes, which design choices cause that unreliability, and can the answer be turned into actionable benchmark-design rules?" The framework is applicable beyond immunotherapy. Any transcriptomic prediction task that uses a gene-defined molecular surrogate endpoint is subject to the same circularity risk, and the CDS diagnostic generalizes to any setting where an expression matrix is paired with a binary label. The answer: response definition exerts a larger effect (largest within-cohort Δ = 0.520, ElasticNet-Var) than any between-method difference measured on the same RECIST endpoint (max Δ = 0.147); across different cohorts the between-method spread can reach 0.508 (Jung), which is why the endpoint effect must be demonstrated within-cohort. In practical terms, at n = 43, *how you define response matters more than which algorithm you choose* — and the mechanism ports across three further cohorts spanning two cancer types (each using the same author-constructed cytolytic surrogate, so these establish portability of the construction rather than independent prevalence; the published-study audit in Tables S17 and S20 establishes prevalence). 
We use the term "design choice" rather than "confounder" throughout, because the three variables we identify are not confounders in the strict causal inference sense (they are operational variables directly under experimental control). They are more precisely described as *design choices that act as variance-dominating factors in small-sample benchmarks*. This terminological precision matters because it clarifies the prescriptive implication: these factors can and should be controlled through experimental design, not merely adjusted for statistically. 
### The Response Definition Problem

The cytolytic → RECIST collapse (largest Δ = 0.520, ElasticNet-Var) has a general mechanism: any molecular surrogate endpoint operationally defined by genes in the prediction feature space creates a circular prediction task. We recommend: (1) always include at least one clinical endpoint; (2) label molecular surrogates with their defining genes; (3) compute CDS for gene-expression-derived endpoints; (4) include non-learning baselines as circularity detectors; and (5) for endpoints with CDS > 0.80 (recalibrated HIGH threshold), assume ML methods are learning the endpoint-defining genes until proven otherwise through independent clinical endpoint validation. 
### Feature Selection as a Narrative Design Choice

The substantial overlap of MI- and variance-selected gene sets (27% of the top-500 genes) with highly correlated out-of-fold predictions (ρ = 0.75) demonstrates that feature selection determines the biological story even when AUROC is statistically indistinguishable at n = 28. We recommend: (1) pre-specify ≥2 feature selection strategies; (2) report gene set overlap; (3) at n < 50, assess cross-seed stability; (4) accompany biological interpretation with sensitivity analysis. **Priority note:** Among feature-selection recommendations, gene set overlap reporting (item 2) is the single most impactful, as it directly reveals whether biological conclusions are strategy-dependent. 
### On Method Context and Generalizability

Our data are consistent with the principle that method performance is inherently context-dependent [15,17]. IMPRES's AUROC ranges from 0.795 (Hugo) to 0.583 (Jung DCB) across cohorts, treatments and tumour types. Across clinical endpoints, GEP ranges from 0.446 (Hugo) to 0.830 (Gide). The context-sensitivity of IMPRES is consistent with Carter et al. (2019), who showed IMPRES failed in an independent melanoma cohort (AUROC 0.52–0.55) [30] and identified potential training-set bias in the original IMPRES feature selection (34% non-random repeats) — an interpretative caveat for IMPRES's Hugo performance (0.795). We recommend that benchmark studies transparently report both performance and known limitations of each method's original development. 
### Independent Prognostic Evidence: TCGA Cancer-Type Gradient

The TCGA analysis below is **prognostic, not predictive**: univariate Cox models of overall survival without clinical covariates (stage, age), evaluating whether immune-signature associations are context-dependent across cancer types. It cannot and does not validate treatment-response prediction; it provides large-scale independent evidence that immune-signature behaviour emerges from method–context interaction. 
Signature significance showed a pronounced cancer-type gradient (Figure S11. BH within each cancer type). The B-cell signature was significant in SKCM, BRCA and LUAD but not in the re-acquired COAD (HR = 0.904 / 0.885 / 0.873 vs 1.153) The CD8 T-cell and TLS signatures were significant in SKCM and BRCA but not LUAD or COAD And the IFNG-response signature was significant only in SKCM. Across the four cancer types the immune signatures point protectively in SKCM/BRCA/LUAD while COAD HRs hover at or slightly above 1 (e.g. CD8 T cells 1.177): the gradient is in the strength and direction of association, and COAD reproduces the archived finding that immune infiltration is not prognostic there. 
**v33 correction (data-provenance audit).** The archived TCGA analysis reported 2,597 patients across four cancer types. Re-running the pipeline against the data currently in the repository reproduces three of the four cohorts and reveals two provenance defects that we disclose here rather than silently re-use the archived numbers:

- **SKCM sample-type policy.** Of 444 SKCM samples, 367 are metastatic biopsies (barcode suffix `-06`) and 76 are primary (`-01`). The archived analysis used all tumour samples (n = 426); the script's primary-only filter would have analysed n = 76. We restored the all-tumour policy and now report **n = 441 patients, 212 deaths**, with **11/12 immune signatures significant after Benjamini–Hochberg correction** (HR 0.795–0.928 for the 11 protective signatures; all BH-adjusted p ≤ 2.7e-3). Notably, the *stromal* signature is the only non-significant one and the only one with HR > 1 (1.086, BH-adjusted p = 0.090), an internal negative control consistent with immune-specific biology rather than a generic prognostic effect. - **COAD source data was not TCGA-formatted; re-acquired from the GDC.** The local COAD expression file used non-TCGA identifiers ("01CO001"; 107 columns) and could not support the archived COAD result (n = 588, 119 deaths), so v33 withdrew it. The cohort has since been re-acquired from the GDC current index (STAR counts, GENCODE v36, primary tumours; clinical from the GDC demographic API; scripts/build_coad_gdc.py), yielding **n = 458 patients, 102 deaths**. On these data **0/12 immune signatures reach significance after BH correction (all HR ≈ 1)**, independently reproducing the archived COAD finding. The expression source (GDC STAR FPKM) differs from the other three cancers' cBioPortal RSEM — a disclosed caveat that does not affect within-cancer univariate statistics. The TCGA validation now comprises **2,491 patients** (SKCM 441 + BRCA 1,082 + LUAD 510 + COAD 458). - **Multiple testing.** The original report used nominal p < 0.05 only. We now report BH-adjusted values within each cancer type: SKCM 11/12, BRCA 3/12, LUAD 1/12 and COAD 0/12 signatures remain significant at FDR < 0.05 — the SKCM conclusion is unchanged and in fact strengthened because it survives FDR correction at 48 tests. 
Current TCGA results (12 signatures × 4 cancer types, 2,491 patients): SKCM 11/12 significant (BH-adjusted), BRCA 3/12, LUAD 1/12, COAD 0/12 (Supplementary Table S15; machine-readable in data/tcga/tcga_cox_all.json).

A pronounced cancer-type gradient emerged (Supplementary Table S15. BH within each cancer type): 11/12 signatures significant in SKCM (per-SD HR 0.795–0.928 for the 11 protective signatures The stromal signature is the only non-significant one and the only one with HR > 1, 1.086 — an internal negative control), dropping to 3/12 in BRCA (CD8 T cells, B cells, TLS), 1/12 in LUAD (B cells) and 0/12 in the re-acquired COAD (all HR ≈ 1) — an ordered decrease across the four cancer types, with the caveat that the intermediate steps rest on one to three signatures each and should be read as descriptive rather than as a calibrated dose–response. No time-dependent AUC is reported under v33 (the archived 0.33–0.60 range was not regenerated). This pattern provides large-scale independent evidence that immune signature performance is not intrinsic to a method but emerges from method–context interaction, consistent with Usset et al. (2024, Nature Genetics), who showed that thousands of published biomarkers collapse into five latent factors [32]. Full methods in Supplementary Note S4. 
### CDS: Current Status and Validation Roadmap

We recognize the epistemological tension inherent in this work: CDS is itself subject to the same validation standards we recommend for benchmarking tools. Its current status can be summarized as follows:

- **What CDS can do:** Quantify endpoint–gene association and place endpoints on a null-referenced scale: the operationally defined cytolytic endpoint scores 0.954 versus null mean 0.566 / max 0.842 (Riaz; 200-replicate null), while Hugo RECIST (0.449) and Riaz RECIST (0.409) sit below their null means — no detectable gene–label dependency beyond chance. Provide a lightweight (<10 seconds) pre-benchmarking diagnostic requiring only an expression matrix, binary labels and declared endpoint genes, with no dependencies beyond NumPy/SciPy/scikit-learn and independent pip installability. Support the procedural circularity rule: any endpoint defined by genes in the feature space is a circular task by construction, whatever its score — a rule the negative-control stress test (Table S18) shows is program-agnostic (60/60 gene-program pseudo-endpoints HIGH). Retrospective application to the Auslander et al. (2018) external validation of IMPRES on Hugo 2016 constitutes a retrospective consistency check: CDS = 0.449 (below its null mean) is consistent with IMPRES AUROC = 0.795 being genuine signal, not circular inflation [5,30]. In addition, we identified five published immunotherapy prediction benchmark studies (Table S17) and retrospectively assessed what CDS would have told the authors. In all cohorts the cytolytic positive controls sit far above their nulls, while the clinical endpoints sit at or within theirs (E1 pre-registered runs: IMvigor210 CDS 0.434 vs null mean 0.492 / max 0.688, positive control 0.953; Liu 2019 CDS 0.721 vs null mean 0.591 / max 0.739, positive control 0.949; an independent 10-replicate E3 run gives the same ordering with null means 0.449/0.623 and maxima 0.496/0.728 — the two null runs are reconciled and both reported in Tables S17/S19/S20; the CTLA-4 cohort of Van Allen et al. [34] remains expression-inaccessible in the harmonization) — the internal finding extends to independent studies, and the IMvigor210 case supersedes the earlier 'not publicly accessible' statement. (The Gide 2019 CDS value of Supplementary Table S4b, 0.852, is a v33 recomputation of a score never stored at analysis time; see Results.)

- **Practical pre-benchmarking workflow:** (1) Declare which genes, if any, operationally define the endpoint; a gene-defined surrogate is a circular task by construction — switch to a clinical endpoint (RECIST, OS) or accept that ML AUROC will be inflated by endpoint–gene co-expression. (2) Compute CDS on the declared genes before training any model (<10 seconds). (3) Report the continuous CDS value against a 200-replicate random-label null (the archived 10-replicate nulls give the same readings) on the same cohort and feature space; a value far above the null maximum flags strong gene–label association, which is either the circular construction itself or genuine biology — distinguishable only by endpoint provenance. (4) Treat the shipped categorical labels (HIGH > 0.70, MODERATE 0.30–0.70, LOW ≤ 0.30) as advisory: the negative-control stress test (Table S18) shows they carry no specificity among structured endpoints. This workflow executes in under 10 seconds per cohort and requires no ML model training. 
- **What CDS cannot yet do:** Prospectively validate on independent multi-endpoint cohorts (no such cohort is publicly available with the required dual-endpoint design). Provide well-calibrated categorical thresholds (current thresholds are provisional, recalibrated from synthetic data with n = 30 null replicates). Distinguish between MODERATE-risk endpoints with meaningful AUROC implications (the C3 saturation problem limits discrimination in the 0.50–0.80 range). Be computed on cohorts where the raw expression data is not publicly accessible — a limitation that also serves as a transparency enforcement mechanism. 
**CDS: Known Failure Modes.** Three failure modes should be considered when interpreting CDS results. (1) **Undeclared endpoint genes (false negatives):** CDS requires the user to specify which genes define the endpoint. Unknown or undisclosed molecular surrogates evade CDS detection. Diagnosis: if CDS is LOW but ML AUROC > 0.80, suspect undeclared endpoint-gene dependency. Mitigation: CDS v2.0 will include a hypothesis-free mode based on permutation testing of the full feature-label MI distribution. (2) **Label noise (false negatives):** If response labels contain substantial misclassification, CDS may return LOW risk because label noise attenuates all gene-label correlations. Diagnosis: compare CDS scores between the suspect labels and a known gold-standard endpoint on the same cohort; a large discrepancy suggests label quality issues. (3) **C3 saturation in immune-dominant contexts (loss of discrimination):** In highly immune-infiltrated cancers (melanoma, NSCLC), the current C3 formulation saturates at 1.000 on every real cohort and on all 60 gene-program pseudo-endpoints of the negative-control stress test (Table S18). Diagnosis: C3 = 1.000 indicates that C1 and C2 are the sole discrimination sources. CDS v2.0 will use permutation-based C3 to estimate the null distribution of feature-space enrichment.

- **Validation roadmap:** (1) Re-calibrate C3 using pathway-stratified feature selection (planned for CDS v2.0). (2) Validate on external multi-endpoint cohorts as they become available (the pre-declared candidate GSE274975 was investigated and disqualified; see Revision Transparency, item 10). (3) Third-party replication by independent groups applying CDS to their own benchmark endpoints before running full evaluations. (4) Pre-registration of CDS scores before benchmark results are published, enabling prospective validation of CDS predictions against subsequently observed AUROC values. 
### Limitations

1. **Sample size:** All 95% CIs for the trainable methods on clinical endpoints at n ≤ 43 overlap 0.5 (the sole fixed-scorer exception is IMPRES on Hugo — the headline positive result — which does not bear on method ranking among trained models. The GSE100797 (Lauss) PD-L1 single-gene biomarker is nominal-only after the completed six-test family correction, p = 0.0096, q = 0.058). IMPRES CI width = 0.350 on Hugo 2016 with 1,000 bootstrap iterations — genuine uncertainty at n = 28. No ML method is statistically superior for n ≤ 43. 
2. **Response-definition evidence base:** The cytolytic→RECIST collapse was demonstrated on Riaz 2017, replicated on two independent cohorts — Liu 2019 (melanoma, n = 122; Table S19) and IMvigor210 (bladder carcinoma, n = 297; Table S20) — using the same author-constructed surrogate (the Riaz construction applied to each cohort), and is k-fold robust (Supplementary Table S16). The collapse magnitude (largest Δ = 0.520) may be cohort-specific: in tumour types with weaker immune infiltration, where fewer immune signatures are prognostic (in TCGA, 0-1 of 12 signatures reach FDR-adjusted significance in LUAD and none in COAD), the circular prediction effect would be attenuated because there is less immune co-expression for ML methods to leverage. NSCLC replication remains an open gap: the pre-declared candidate (GSE274975) was investigated and found to lack clinical endpoints (see Revision Transparency, item 10); a suitable NSCLC ICI cohort with dual endpoints and public expression data has not yet been identified. 
3. **CDS calibration, program-agnosticity and C3 formulation:** The CDS thresholds are provisional, and the negative-control stress test (Supplementary Table S18) shows the categorical labels have no specificity among structured endpoints: all 60 gene-program pseudo-endpoints score HIGH (0.79–0.95), and random-label replicates on the real cohorts reach 0.87 (Lauss; Table S4b). Under the shipped 0.70 HIGH boundary, 20–36% of 200 random-label replicates exceed it on every real cohort (`cds_nulls_200_v33.json`), versus ≈3% on the synthetic calibration null — so the shipped labels overstate specificity on real data by roughly an order of magnitude. CDS's validated content is the continuous value read against a cohort-matched random-label null; the circularity protection is procedural (declare the endpoint-defining genes). The Gide 2019 RECIST value (0.852, above its 200-replicate null maximum of 0.800) quantifies the core limitation: CDS v1.x cannot separate genuine biology (GEP 0.830 on the same cohort) from circularity. The C3 component saturates at 1.000 on all real cohorts because highly variable genes in bulk tumour RNA-seq are predominantly immune-related; the composite is therefore bounded to [0.40, 1.00] on real tumour RNA-seq, making the shipped LOW tier (≤ 0.30) mathematically unreachable there, and C1 and C2 carry the entire discrimination burden. A permutation-based C3 is planned for CDS v2.0. 
4. **Data access and cross-cohort validation:** In earlier versions IMPRES, GEP, and TIDE could not be evaluated on Jung 2019 (Ensembl IDs without HGNC mapping). As of v33 the HGNC-mapped matrix is used for all methods on this cohort, and we additionally corrected the Jung label orientation (the release codes 1 = no durable benefit, the inverse of every other cohort All analyses here use 1 = DCB). Leave-one-cohort-out validation remains precluded by gene identifier heterogeneity across cohorts — a barrier that is itself a finding, and CDS's data-access requirement (expression matrix + labels must be publicly available) serves as a de facto transparency enforcement mechanism. 
5. **Clinical biomarker proxies:** TMB and PD-L1 IHC results use gene expression proxies (constant baseline and CD274 expression, respectively).

6. **Cohort composition:** The immunotherapy cohorts are predominantly melanoma (169/196 patients), and one of the five is adoptive T-cell therapy rather than ICI (Methods). The TCGA analysis provides cross-cancer validation (2,491 patients, 4 cancer types), but CDS calibration and the response-definition collapse have only been demonstrated in melanoma. Gide 2019 includes both monotherapy and combination therapy arms, introducing treatment heterogeneity. 
7. **Learning curve caveat:** The n = 8 subsample analysis (k = 5, ~6 training samples per fold) is diagnostic rather than a precise estimate.

8. **IMPRES interpretative caveat:** IMPRES's Hugo performance (0.795, q = 0.017) may partially reflect inherited feature selection bias from the original IMPRES training protocol (34% non-random repeats) [30]. The extent to which this reflects genuine checkpoint biology versus inherited selection bias is unresolved. 


---

## REVISION TRANSPARENCY

This manuscript has undergone a multi-round audit (v28 → v35) that identified and corrected data-layer and implementation defects — including gene-identifier mapping (Riaz Entrez IDs), label orientation (Jung), matrix structure (Riaz RECIST rebuild), implementation fidelity (IMPRES canonical form), feature-selection pipeline corrections, survival-event approximation (Hugo), TCGA-COAD provenance, a disqualified pre-declared replication cohort (GSE274975), a cohort-citation error (GSE100797 is Lauss et al. 2017, not the CTLA-4 Nathanson study), and a Benjamini–Hochberg implementation error in the reporting scripts. The full itemized correction log is Supplementary Note S7; superseded numbers are not reproducible from the current repository. 
---

## STAR METHODS

### RESOURCE AVAILABILITY

**Lead contact**

Further information and requests for resources should be directed to and will be fulfilled by the lead contact, Zhixian Huang (huangzhixian324@gmail.com), the sole author of this study.

**Materials availability**

This study did not generate new materials.

**Data and code availability**

- Data: All raw cohort data are publicly available from GEO (GSE78220, GSE91061, GSE100797, GSE135222) and cBioPortal (mel_iatlas_gide_2019, blca_iatlas_imvigor210_2017, mel_iatlas_liu_2019).
- Code: All original code is publicly available. SPATBench pipeline (v33): https://github.com/huangzhixian324-dev/spa-bench (MIT). CDS package (v1.2.1): https://github.com/huangzhixian324-dev/cds (MIT); the packaged source is also included in the Zenodo archive (`cds_tool/`, installable with `pip install ./cds_tool`). A versioned archive with a DOI is deposited on Zenodo (DOI: [ZENODO DOI — to be inserted at acceptance]). The Docker image (spa-bench:v33) reproduces the synthetic validation test (`scripts/validate.py`, exit code 0 verified) and ships the complete pipeline code; reproducing the reported tables additionally requires the Zenodo-archived inputs and `scripts/rerun_v33.py`. - Any additional information required to reanalyze the data reported in this paper is available from the lead contact upon request.

### EXPERIMENTAL MODEL AND STUDY PARTICIPANT DETAILS

This study is purely computational and generated no new experimental data: no experimental models, human participants, or animals were used. All cohorts are publicly available, de-identified bulk RNA-seq datasets from published studies (GEO accessions GSE78220, GSE91061, GSE100797 and GSE135222; cBioPortal studies mel_iatlas_gide_2019, mel_iatlas_liu_2019 and blca_iatlas_imvigor210_2017; TCGA via the GDC), analyzed under the terms of their original publications.

### METHOD DETAILS


### Cohort Selection and Preprocessing

Five publicly available cohorts with pre-treatment bulk RNA-seq data and annotated clinical outcomes were analyzed (Supplementary Table S3):

| Cohort | Accession | Cancer | n | Treatment | Endpoint | Responders |
|--------|----------|--------|----|-----------|----------|------------|
| Hugo 2016 | GSE78220 [8] | Melanoma | 28 | Pembrolizumab | RECIST v1.1 | 13 (46%) |
| Riaz 2017 | GSE91061 [7] | Melanoma | 43 (42 with RECIST) | Nivolumab | RECIST + Cytolytic | 9 RECIST (of 42) / 21 Cytolytic (of 43) |
| Lauss 2017 (ACT) | GSE100797 [20] | Melanoma | 25 | Adoptive TIL | RECIST v1.1 | 10 (40%) |
| Jung 2019 | GSE135222 [24] | NSCLC | 27 | Anti-PD-1/PD-L1 | DCB vs NDB (RECIST v1.1, ≥6-mo follow-up) | 6 DCB (22%) |
| Gide 2019 | mel_iatlas_gide_2019 [21] | Melanoma | 73 | Anti-PD-1 ± anti-CTLA-4 | RECIST v1.1 | 40 (55%) |

**Total: 196 patients** (28 + 43 + 25 + 27 + 73). Cohorts were selected based on availability of pre-treatment bulk RNA-seq data with public clinical annotations. One cohort is adoptive T-cell therapy rather than checkpoint blockade: GSE100797 is the cohort of Lauss et al. (2017) [20] — referred to as "Nathanson" in parts of the benchmarking literature, an attribution error in earlier manuscript versions that is retained only as a legacy identifier in internal dataset filenames and machine-readable table rows (Revision Transparency, Supplementary Note S7 item 11). It is included deliberately as a stress test of method context-dependence beyond ICI; no headline claim depends on it, and we refer to the five cohorts as immunotherapy-treated throughout. The Riaz 2017 cohort provides two response definitions on the same patients — RECIST v1.1 from GEO SOFT metadata (of 49 pre-treatment biopsies with RECIST annotations, 9 had PRCR and 33 had PD/SD; 7 had unknown response and were excluded, leaving **n = 42** for the RECIST endpoint; the cytolytic endpoint uses all n = 43) and cytolytic activity (median-split of geometric mean of GZMA and PRF1; 21 above median [14]) — enabling a within-cohort, within-patient comparison of response definition effects. Cohorts were not pooled for meta-analysis due to endpoint heterogeneity (RECIST, DCB/PFS, and molecular surrogates capture different biology).

**Label semantics:** Throughout this work label 1 = clinical benefit (RECIST responder / DCB). For Jung 2019 (GSE135222), the original release coded 1 = NDB — the inverse of every other cohort; all Jung analyses use corrected labels (1 = DCB, 6/27; see Revision Transparency, item 2).

All expression data were standardized to log2(FPKM + 1). Patient-stratified k-fold cross-validation used k = 5 for all five cohorts except the Riaz RECIST endpoint, where the small responder count (9 of 42) constrained the design to k = 2. To assess k-fold robustness for the Riaz RECIST endpoint, we repeated the analysis with k = 3 and k = 5 (10 stratified random seeds each) and confirmed that the cytolytic→RECIST collapse magnitude is not an artifact of the k-fold choice (Supplementary Table S16; regenerated values in results/benchmark/v33/riaz_recist_kfold_v33.json). Gene symbols were mapped to HGNC where possible (Riaz 2017: Entrez-to-HGNC, 21,769/22,187 genes; Jung 2019: Ensembl-to-HGNC; see Revision Transparency, items 1–3, for data-layer corrections). 
### Methods Evaluated

**Gene signatures (no training required):** GEP (the 8-gene Ayers IFN-γ/T-effector signature: CD8A, GZMA, GZMB, IFNG, CXCL9, CXCL10, PRF1, TBX21 [3]; the circularity mechanism analysed in this paper is specific to gene sets containing GZMA/PRF1 — the 18-gene expanded GEP does not contain them and would not be expected to show the same collapse), TIDE (tidepy v1.3.9), IMPRES (the 15 canonical pairwise checkpoint-gene relations of Auslander et al. [5]; the relations were decoded from the authors' published feature set (github.com/noamaus/IMPRES-codes, FEATS.mat) and cross-checked against the 9 pairs named in the paper text (9/9 match); the remaining 6 relations are not independently verified beyond the authors' released file. For a sample, the IMPRES score is the fraction of the 15 logical relations (e.g., CD40 > PD-1) whose expression ordering holds, linearly rescaled when features are missing, exactly as described in the original Methods. No retraining is performed, mirroring the original design in which features were selected once on neuroblastoma spontaneous-regression data.)

**Gene expression biomarker:** PD-L1 (CD274 expression; proxy for IHC — real IHC scores unavailable in public metadata).

**Machine learning (training required):** ElasticNet logistic regression (sklearn saga solver, penalty = 'elasticnet') with variance-based feature selection (top-500 genes by expression variance in the training fold), ElasticNet with mutual information-based feature selection (top-500 genes by MI with the response label in the training fold), and XGBoost gradient boosting. ElasticNet hyperparameters (C, l1_ratio) were tuned via 3-fold inner cross-validation within each training fold; the selected hyperparameters were frozen during permutation testing, as is standard.

**Supplementary baselines:** TMB (constant 0.5; real mutation data unavailable) [22]. Immune deconvolution signatures (14 cell types, 43 unique genes as implemented in `run_tcga_cox.py`; Supplementary Table S10) [23].

The six benchmarked prediction methods span four categories. All feature selection used only training-fold data.

### Evaluation Protocol

Patient-stratified k-fold cross-validation with identical splits across all methods within each cohort. Ten metrics were computed. Here we report AUROC as the primary metric with 95% bootstrap confidence intervals (1,000 iterations) [28]. Expression-scale, preprocessing and cross-validation-strategy sensitivity analyses on Hugo 2016 (ElasticNet-MI) are provided in Supplementary Figure S5 and Tables S7 and S8. AUROC, AUPRC and bootstrap CIs for every method–cohort cell are provided in Supplementary Table S9, with clinical-utility metrics (Brier score, DCA, NRI, IDI) for a representative trainable cell; per-fold values are fully determined by the archived split files (Supplementary Table S1). For cohorts with n ≤ 28, we additionally report cross-fold standard deviation (SD) as a stability metric. 
### Permutation Testing

Full permutation tests were performed for all trainable method–cohort pairs. Exact p-values = (number of permuted AUROCs ≥ observed + 1) / (n_permutations + 1). Benjamini-Hochberg FDR correction was applied across all tests within each cohort [26].

**Inclusion logic:** Two protocols are used, matched to each method class. (a) Fixed scorers — IMPRES, GEP, TIDE, and PD-L1 (CD274 expression) — produce scores that are independent of the labels, so the full-pipeline label-shuffle null is mathematically identical to shuffling labels against the fixed out-of-fold scores; each cell therefore uses 50,000 prediction shuffles (seconds of compute, exact to 1/50,001). This equivalence was verified by running both protocols on Hugo IMPRES. (b) Trainable methods — ElasticNet (MI) and ElasticNet (Var) — are re-fit end-to-end under every label shuffle (per-fold feature selection with frozen hyperparameters), with a per-cell number of shuffles recorded in Table 4. XGBoost is excluded from benchmark reporting (degenerate on all clinical endpoints; its only non-degenerate AUROC sits on the circular cytolytic endpoint — development-phase evaluation only. Correction (this revision): under the current protocol XGBoost is not degenerate — re-runs produce non-degenerate predictions (Hugo 0.554, Lauss 0.493, Jung 0.611, Riaz-RECIST 0.434, Gide 0.718, and 0.991 on the circular cytolytic endpoint; all six cohorts re-run; `results/benchmark/v33/xgb_exclusion_evidence_v33.json`) and its inclusion would not change any conclusion, since it shows no signal at n <= 43 and its only higher value lies at n = 73, outside the scope of the negative claim. XGBoost is absent from the reported matrix because the development-phase runs were not retained under the current protocol, not because of degeneracy.) are retained for audit even though the prediction tables are not).

**Permutation resolution:** p-values are exact given the recorded per-cell shuffle count: 50,000 for fixed scorers (minimum detectable p = 2 × 10⁻⁵) and 60–5,000 full-pipeline shuffles for the trainable methods, where the compute cost of re-fitting under every shuffle (≈60–90 s per shuffle) is the binding constraint. We state the resolution ceiling explicitly: at 60 shuffles the smallest achievable exact p is 1/61 ≈ 0.016, so within a six-test cohort family the smallest attainable BH q is ≈ 0.10; the non-survival of trainable methods on clinical endpoints at n ≤ 43 is therefore established at this resolution and is corroborated independently by bootstrap 95% CIs that all overlap 0.5 (Supplementary Table S9). Trainable cells with p ≤ 0.08 at the initial budget were extended where compute permitted (up to 1,000 shuffles; the Gide ElasticNet-MI cell to 5,000). The one remaining borderline trainable cell, ElasticNet-MI on Gide 2019, was subsequently extended to 5,000 full-pipeline shuffles (checkpointed resumable protocol, frozen hyperparameters): the observed AUROC under this protocol is 0.662 (0.032 above the benchmark 0.630 — frozen-versus-tuned hyperparameters, the same deviation class disclosed in the Table 4 footnote), giving p = 0.0232 and within-cohort BH q = 0.028 — **FDR-significant**. ElasticNet-MI therefore joins ElasticNet-Var as a detectable trainable signature on Gide (n = 73), strengthening the reading that Gide is the smallest sample size at which trainable transcriptomic signatures become detectable on a clinical endpoint; our headline claim is unchanged, since it is scoped to n ≤ 43, where no trainable method survives (Gide sits above that threshold). The smallest non-significant within-cohort BH q-value on any clinical endpoint at n ≤ 43 is 0.058 (Lauss PD-L1).

### Multiple Comparison Framework

Three levels of inference: (i) Exploratory: AUROC point estimates with 95% bootstrap CIs reported without multiplicity correction (descriptive quantities). (ii) Permutation testing: all evaluated method–cohort pairs with within-cohort FDR correction. (iii) Pairwise comparisons of AUROCs via paired bootstrap tests of ΔAUROC (bootstrap implementation. The analytic DeLong covariance method was considered but the bootstrap version was used throughout — we refer to these as paired bootstrap tests, not DeLong tests). The Bonferroni threshold (p < 0.05/35 ≈ 0.0014) is unattainable at current sample sizes — we treat this inaccessibility as a finding. Because family definition is an analysis degree-of-freedom, Supplementary Table S21 reports the full 36-cell matrix under three family definitions (within-cohort 6×6; global 36; method-class-stratified 24+12): the primary conclusions are family-robust, and the four family-dependent cells are disclosed there. 
### Effect-Size Comparison Framework

To compare the relative influence of methodological choices against method choice, we computed the range of AUROC values attributable to each source on the Riaz 2017 cohort: (a) between-endpoint ΔAUROC for each method; (b) between-method ΔAUROC within each endpoint; (c) between-feature-selection ΔAUROC within ElasticNet. We report these as effect-size comparisons rather than variance-component percentages, which are unreliable with limited independent data points. 
Between-method spreads differ across cohorts, so the endpoint effect must be compared with method effects at the same level of control. Within the same patients and the same RECIST endpoint (Riaz), the largest between-method ΔAUROC is 0.147 (IMPRES 0.613 vs TIDE 0.466), so the largest between-endpoint effect (Var: 0.985 → 0.465, Δ = 0.520) exceeds it by 3.5×; the same comparison yields ≈2.2× on Liu 2019 (0.425/0.195) and ≈2.0× on IMvigor210 (0.354/0.180) — the direction replicates, with a smaller ratio. The effect is further supported by non-overlapping bootstrap confidence intervals. The cytolytic-endpoint AUROC CIs for ElasticNet-Var ([0.954, 1.000]) and ElasticNet-MI ([0.909, 1.000]) do not overlap the RECIST-endpoint CIs ([0.208, 0.704] and [0.370, 0.755] respectively), and the same non-overlap pattern holds for the external replications (Liu 2019: cytolytic GEP [0.928, 0.986] vs RECIST [0.427, 0.642]. IMvigor210: cytolytic [0.948, 0.981] vs RECIST [0.536, 0.686]). Across different cohorts the between-method spread can be as large as 0.508 (Jung: GEP 0.786 vs TIDE 0.278), so a cross-cohort comparison cannot support the claim; the controlled within-cohort comparison on Riaz is the appropriate design, because it fixes patients, treatment, and expression data and varies only the endpoint definition. 
### Power Analysis

Minimum detectable ΔAUROC at 80% power (α = 0.05, two-sided) as a function of sample size and responder proportion, using the Hanley–McNeil asymptotic variance of the AUC (Mann–Whitney U) for a paired two-method comparison on the same patients (baseline AUC 0.5) with AUC correlation r = 0.75 — the upper end of the correlation observed between method pairs on out-of-fold scores (Supplementary Table S12). The computation is archived as `scripts/power_table_v33.py` with a machine-readable sensitivity grid over r ∈ {0, 0.5, 0.75} (`power_table_v33.json`); Table 8 reports the r = 0.75 (most favourable) column, rounded to 0.01, and regenerates exactly from the archived script. Weaker correlation assumptions inflate the required n roughly two- to four-fold. Full power curves in Supplementary Figure S9. 
### Circularity Detection Score (CDS)

CDS quantifies endpoint-feature dependency as a weighted composite: **C1 (30%):** mean MI between endpoint-defining genes and the response label, normalized as min(C1_MI / 0.05, 1). **C2 (30%):** mean |Spearman ρ| between endpoint-defining genes and the response **C3 (40%):** min(1, R/10), where R is the ratio of the 90th-percentile MI of the top-500-variance features to the median MI of 500 random features. The composite is CDS = 0.3·C1n + 0.3·C2 + 0.4·C3. These weights and normalization constants are those of the released tool (v1.2.1; source and packaging at https://github.com/huangzhixian324-dev/cds and in the Zenodo archive); every CDS value in this manuscript is exactly reproducible from them (CDS scores are identical between v1.1.0 and v1.2.1 — v1.2.1 changes only the output semantics: it adds a null-referenced reading mode, a caveats block shipping the calibration limits below, and condition-specific guidance). The CDS concept builds on earlier work by Zhao et al. (2011, BMC Genomics), who introduced a "consistency degree" metric quantifying correlation between disease endpoints and gene expression profiles [29]. We note two properties of this formulation: (a) C3 saturates at 1.0 for all evaluated endpoints and random-label replicates — highly variable genes in bulk tumor RNA-seq are predominantly immune-related, so the top-to-random MI ratio consistently hits the clipping ceiling regardless of endpoint circularity. And (b) we evaluated a re-weighting (C1 40% / C2 40% / C3 20%) as a candidate refinement, but the released tool and all reported values use the 30/30/40 formulation Permuting C3's role (v2.0) is planned to improve MODERATE-range discrimination. 
**Synthetic data calibration:** To validate CDS independently of our three real endpoints, we generated synthetic expression data (n = 200, 5,000 genes) with controlled circularity levels (0–100% endpoint-gene signal, with 50 independent biological-signal genes as non-circular signal). CDS ranks the five scenarios monotonically (CDS 0.697 < 0.820 < 0.830 < 0.832 < 0.845; Supplementary Table S4d). Scenarios with CDS ≥ 0.80 show inflated AUROC (0.691–0.845), while the pure-random scenario scores 0.697 — above the level of an unstructured label, reflecting C3 saturation. 

**Retrospective external consistency check:** The Auslander et al. (2018) study provides an independent validation context for CDS. In that study, IMPRES was trained on Auslander's own data and validated externally on Hugo 2016 (RECIST endpoint), achieving AUROC = 0.83 in the original paper (0.795 in our benchmark). If CDS had been computed before that validation was performed, it would have returned CDS = 0.449 (MODERATE under the shipped thresholds, and below its cohort-matched random-label null mean of 0.586), indicating that the Hugo RECIST endpoint is not circular — any observed AUROC reflects genuine biological signal in the IMPRES gene pairs, not endpoint-gene dependency. This counterfactual application of CDS to a published external validation scenario constitutes a retrospective consistency check: CDS's reading (non-circular endpoint → AUROC not inflated by endpoint-gene dependency) is consistent with the observed outcome (IMPRES achieved significant AUROC = 0.795 on Hugo RECIST, permutation p = 0.0029). The same logic applies to Lauss 2017 (ACT) (CDS = 0.781, MODERATE under the recalibrated boundary; no trainable AUROC exceeds 0.70, no inflation). For Gide 2019 no CDS value was stored at analysis time; the v33 recomputation (Supplementary Table S4b) returns CDS = 0.852 — the C1/C2 components read the strong association between the declared genes (GZMA/PRF1) and the clinical label in this large melanoma cohort as endpoint-gene dependency, i.e. the known inability of CDS v1.x to separate genuine immune-biology signal (GEP AUROC = 0.830, q < 0.001) from circularity. CDS correctness claims are therefore restricted to the endpoints with stored scores.

**Null distribution, stress test, and reporting convention:** We generated 30 independent replicates of random binary labels to estimate the CDS null distribution. Mean null CDS = 0.50 (SD = 0.09), 95th percentile = 0.65, 99th percentile = 0.71. Approximately 3% of random endpoints were flagged as HIGH under the shipped threshold (CDS > 0.70) due to C3 saturation at 1.0 — a known limitation of the current C3 formulation; no random endpoint achieved CDS < 0.30. We disclose that these 30 null replicates were generated **after** the four real benchmark endpoints had been inspected; threshold calibration is therefore retrospective, not prospective, and categorical risk labels must be treated as advisory. We deliberately ship a single threshold set (HIGH > 0.70, MODERATE 0.30–0.70, LOW ≤ 0.30) — identical in the manuscript, in the released pip package, and in all archived result files — so that users replicating our analysis obtain exactly our labels. A more conservative HIGH threshold of 0.80 (false-positive rate 0/30 on the synthetic null) is reported for transparency and will be adopted in CDS v2.0 together with a permutation-based C3. We recommend reporting the continuous CDS value; categorical risk assignments are advisory.

 To characterise what CDS measures, we generated 30 pseudo-endpoints per cohort from real gene programs [25] — ECM/stromal (collagens, fibroblast, endothelial, pericyte and basement-membrane markers), tumour-purity/proliferation (cell cycle, E2F/MYC targets, glycolysis, oxidative phosphorylation) and lineage/signalling programs (melanocytic differentiation, keratinocyte basal, hypoxia, interferon, antigen presentation) — each dichotomised at the median of the program-mean expression, the same construction as the cytolytic surrogate, with the program's own genes declared as endpoint genes. All 60 pseudo-endpoints (Gide 2019, n = 73; Hugo 2016, n = 28) scored above the HIGH boundary (CDS 0.792–0.946): CDS flags any operationally gene-defined surrogate regardless of program identity. Random-label nulls computed per endpoint on the real feature spaces (10 replicates each; Table 3) show that the shipped 0.70 threshold also captures 3/10 null replicates at n = 28. CDS is therefore validated as a continuous, null-referenced measure of endpoint–gene association; the categorical HIGH label carries no specificity among structured endpoints, and the protection it enables is procedural — declaring the endpoint-defining genes identifies the circular construction itself. Full design and per-endpoint results in Supplementary Table S18 and Note S6. 
 CDS values are reported alongside the mean, maximum and 95th percentile of a 200-replicate random-label null computed on the same cohort and feature space (Table 3); the archived 10-replicate nulls give the same readings for every endpoint. A CDS far above the null maximum indicates that the declared genes track the label beyond chance — expected for operationally gene-defined endpoints (circularity), but also possible for genuine biology when the declared genes are true correlates of response (the Gide 2019 case, Table S4b). The two are distinguishable only by endpoint provenance, which is why the declaration of endpoint-defining genes is the primary control. 
**Epistemological note:** We recognize the inherent tension: CDS is itself subject to the same validation standards we recommend for benchmarking tools — prospective evaluation on independent multi-endpoint cohorts, transparent reporting of limitations, and third-party replication. Endpoint–feature overlap of the kind CDS detects is one of the leakage classes formalized in the machine-learning reproducibility literature [31]. To facilitate this, CDS is released as an open-source, independently installable package with a documented API and bundled synthetic calibration data. Its current status is experimental, and its provisional thresholds should be refined as more multi-endpoint cohorts become available.

### IMPRES Implementation

IMPRES uses the canonical implementation decoded from the original authors' released feature set: 15 fixed logical relations (list in Supplementary Note S5), score = fraction fulfilled, no retraining. This matters because IMPRES on Hugo 2016 is the headline positive result of this study — it is the only FDR-significant non-circular-endpoint result at n ≤ 43 in Table 4; that result must rest on the published algorithm, not on an ad hoc variant. On the canonical implementation, IMPRES on Hugo 2016 yields AUROC = 0.795 (permutation p = 0.0029 at 50,000 prediction shuffles under the fixed-scorer protocol; q = 0.017 after BH). We note that Carter et al. (2019) identified potential training-set bias in the original IMPRES feature selection (34% non-random repeats) [30]. The extent to which IMPRES's Hugo performance reflects genuine biology versus inherited feature selection bias is unresolved; we report the AUROC as observed and flag this interpretative caveat. 
---

### QUANTIFICATION AND STATISTICAL ANALYSIS

Permutation testing uses full-pipeline label shuffles with frozen hyperparameters; exact p-values are (ge + 1)/(n + 1). Within-cohort multiple testing uses Benjamini-Hochberg FDR (q < 0.05). AUROC 95% CIs use the paired bootstrap (1,000 iterations, the convention used throughout and in Table S9); correlated-AUROC comparisons likewise use the paired bootstrap rather than DeLong. Power analysis and minimum-sample guidance are derived from observed effect sizes (Supplementary Table S9). 
### ADDITIONAL RESOURCES

None.

### KEY RESOURCES TABLE

| Resource | Source | Identifier |
|---|---|---|
| SPATBench pipeline (v33) | This paper; GitHub (https://github.com/huangzhixian324-dev/spa-bench) | spa-bench:v33 |
| circularity-detection-score | This paper; GitHub (https://github.com/huangzhixian324-dev/cds) | v1.2.1 |
| tidepy | Jingxin Fu | v1.3.9 |
| scikit-learn | scikit-learn.org | >=1.0 |
| GSE78220 / GSE91061 / GSE100797 / GSE135222 | GEO | GSE78220; GSE91061; GSE100797; GSE135222 |
| mel/blca/liu iatlas cohorts | cBioPortal | mel_iatlas_gide_2019; blca_iatlas_imvigor210_2017; mel_iatlas_liu_2019 |

## SUPPLEMENTARY DATA

Supplementary Data are provided as a separate PDF document (SPATBench_CellSystems_supplementary_v35.pdf; Markdown source supplementary_material_v33.md, generated by `scripts/synthesize_v33_supplementary.py` from the machine-readable outputs in `results/benchmark/v33/`). Contents:

- **Table S1:** Per-fold cross-validation setup
- **Table S2:** Literature comparison — published vs. SPATBench AUROC
- **Table S3:** Cohort characteristics (detailed)
- **Table S4:** CDS — full component breakdown, synthetic calibration data, null distribution analysis, recalibrated thresholds
- **Table S5:** Method concordance — Cohen's κ, full 10-pair matrix (Hugo 2016)
- **Table S6:** IMPRES leave-one-pair-out analysis
- **Table S7:** Preprocessing sensitivity
- **Table S8:** Cross-validation sensitivity
- **Table S9:** Full five-cohort benchmark results (AUROC/CI/AUPRC, clinical-utility metrics, survival stratification)
- **Table S10:** Immune deconvolution marker genes (43 genes, 14 signatures)
- **Table S11:** Complete permutation test results
- **Table S12:** Random seed stability analysis (10 seeds)
- **Table S13:** Learning curve analysis (Hugo 2016, n = 8–28)
- **Table S14:** Effect-size comparison — Riaz 2017
- **Table S15:** TCGA Cox regression (12 signatures × 4 cancer types, 2,491 patients; COAD re-acquired from the GDC current index)
- **Table S16:** k-fold robustness analysis — Riaz RECIST (k = 2, 3, 5)
- **Table S17:** CDS retrospective assessment on five published benchmark studies + computed CDS for the iAtlas-accessible cohorts (E3)
- **Table S18:** CDS negative-control stress test — 60 real gene-program pseudo-endpoints (E2)
- **Table S19:** External replication — response-definition effect on Liu 2019 (independent cohort, E1)
- **Table S20:** Cross-cancer replication — IMvigor210 (bladder carcinoma, non-melanoma, E1)
- **Table S21:** Multiple-comparison family sensitivity — BH q under within-cohort / global / method-class-stratified family definitions (36 cells)
- **Figures S1–S12**
- **Notes S1–S7** (Note S7: itemized revision-transparency log)

## FIGURE LEGENDS

**Figure 1. Response-definition collapse on Riaz 2017.** AUROC per method on the same patients under the cytolytic-surrogate endpoint versus the clinical RECIST endpoint (HGNC-mapped matrix, v33; 95% bootstrap CIs). Every method whose feature set overlaps the endpoint-defining genes (GZMA/PRF1) collapses; IMPRES, with no overlap, is nearly insensitive. The single FDR-significant clinical-endpoint cell (IMPRES on Hugo, 0.795) carries the Carter et al. training-set-bias caveat stated in Methods and Limitations and should not be read as validated biology until that question is resolved. The identical file is provided as Supplementary Figure S1.

Alt text: Bar chart comparing AUROC of six prediction methods on the same Riaz 2017 patients under two endpoints: the gene-defined cytolytic surrogate (all methods with overlapping feature sets score high, 0.83-0.99) versus the clinical RECIST endpoint (all collapse to 0.40-0.57). IMPRES, whose feature set does not overlap the endpoint genes, is nearly unchanged between endpoints.

**Figure 2. CDS versus observed trainable-ML performance.** Circularity Detection Score (v1.2.1; scores identical to the v1.1.0 recomputation, Table S4b) against the best trainable-ML AUROC per endpoint, with cohort-matched random-label nulls; the circular cytolytic endpoint separates from clinical endpoints. The identical file is provided as Supplementary Figure S10.

Alt text: Scatter plot of Circularity Detection Score versus best trainable-ML AUROC across five endpoints, with cohort-matched random-label nulls. The circular cytolytic endpoint (CDS 0.954) separates from the clinical endpoints, which sit at or inside their nulls except Gide RECIST (CDS 0.852, above its null maximum, reflecting genuine biology the categorical label cannot separate).

**Figure 3. Complete permutation test matrix.** Observed AUROC and permutation p per method–cohort cell; asterisks mark within-cohort BH q < 0.05. No trainable method survives FDR on any clinical endpoint at n ≤ 43. The identical file is provided as Supplementary Figure S12.

Alt text: Heatmap of the complete 36-cell permutation test matrix showing observed AUROC and permutation p per method-cohort cell; asterisks mark within-cohort BH q below 0.05. No trainable method survives FDR correction on any clinical endpoint with n at most 43.

## ACKNOWLEDGEMENTS

The author thanks the maintainers of the public datasets used in this study.

## AUTHOR CONTRIBUTIONS

The sole author was responsible for all aspects of this work: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Resources, Data curation, Writing - original draft, Writing - review and editing, Visualization, Supervision, Project administration, and Funding acquisition.

## FUNDING

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## DECLARATION OF INTERESTS

The author declares no competing interests.

## DECLARATION OF GENERATIVE AI AND AI-ASSISTED TECHNOLOGIES IN THE WRITING PROCESS

During the preparation of this work the author used ChatGPT (OpenAI) in order to improve the clarity and language of the manuscript. After using this tool, the author reviewed and edited the content as needed and takes full responsibility for the content of the published article.

## REFERENCES

1. Topalian, S.L., Hodi, F.S., Brahmer, J.R., Gettinger, S.N., Smith, D.C., McDermott, D.F., Powderly, J.D., Carvajal, R.D., Sosman, J.A., Atkins, M.B., et al. (2012). Safety, activity, and immune correlates of anti-PD-1 antibody in cancer. *N Engl J Med*, 366, 2443–2454. doi:10.1056/NEJMoa1200690
2. Hodi, F.S., O'Day, S.J., McDermott, D.F., Weber, R.W., Sosman, J.A., Haanen, J.B., Gonzalez, R., Robert, C., Schadendorf, D., Hassel, J.C., et al. (2010). Improved survival with ipilimumab in patients with metastatic melanoma. *N Engl J Med*, 363, 711–723. doi:10.1056/NEJMoa1003466
3. Ayers, M., Lunceford, J., Nebozhyn, M., Murphy, E., Loboda, A., Kaufman, D.R., Albright, A., Cheng, J.D., Kang, S.P., Shankaran, V., et al. (2017). IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade. *J. Clin. Invest.* 127, 2930–2940. doi:10.1172/JCI91190
4. Jiang, P., Gu, S., Pan, D., Fu, J., Sahu, A., Hu, X., Li, Z., Traugh, N., Bu, X., Li, B., et al. (2018). Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response. *Nat. Med.* 24, 1550–1558. doi:10.1038/s41591-018-0136-1
5. Auslander, N., Zhang, G., Lee, J.S., Frederick, D.T., Miao, B., Moll, T., Tian, T., Wei, Z., Madan, S., Sullivan, R.J., et al. (2018). Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma. *Nat. Med.* 24, 1545–1549. doi:10.1038/s41591-018-0157-9
6. Chowell, D., Yoo, S.K., Valero, C., Pastore, A., Krishna, C., Lee, M., Hoen, D., Shi, H., Kelly, D.W., Patel, N., et al. (2022). Improved prediction of immune checkpoint blockade efficacy across multiple cancer types. *Nat Biotechnol*, 40, 499–506. doi:10.1038/s41587-021-01070-8
7. Riaz, N., Havel, J.J., Makarov, V., Desrichard, A., Urba, W.J., Sims, J.S., Hodi, F.S., Martín-Algarra, S., Mandal, R., Sharfman, W.H., et al. (2017). Tumor and microenvironment evolution during immunotherapy with nivolumab. *Cell* 171, 934–949.e16. doi:10.1016/j.cell.2017.09.028
8. Hugo, W., Zaretsky, J.M., Sun, L., Song, C., Moreno, B.H., Hu-Lieskovan, S., Berent-Maoz, B., Pang, J., Chmielowski, B., Cherry, G., et al. (2016). Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma. *Cell* 165, 35–44. doi:10.1016/j.cell.2016.02.065
9. Simon, R., Radmacher, M.D., Dobbin, K. and McShane, L.M. (2003). Pitfalls in the use of DNA microarray data for diagnostic and prognostic classification. *J Natl Cancer Inst*, 95, 14–18. doi:10.1093/jnci/95.1.14
10. Vabalas, A., Gowen, E., Poliakoff, E. and Casson, A.J. (2019). Machine learning algorithm validation with a limited sample size. *PLoS One* 14, e0224365. doi:10.1371/journal.pone.0224365
11. Button, K.S., Ioannidis, J.P.A., Mokrysz, C., Nosek, B.A., Flint, J., Robinson, E.S.J. and Munafò, M.R. (2013). Power failure: why small sample size undermines the reliability of neuroscience. *Nat. Rev. Neurosci.* 14, 365–376. doi:10.1038/nrn3475
12. Shi, L., Campbell, G., Jones, W.D., Campagne, F., Wen, Z., Walker, S.J., Su, Z., Chu, T.M., Goodsaid, F.M., Pusztai, L., et al. (2010). The MicroArray Quality Control (MAQC)-II study of common practices for the development and validation of microarray-based predictive models. *Nat. Biotechnol.* 28, 827–838. doi:10.1038/nbt.1665
13. Berrar, D. and Flach, P. (2012). Caveats and pitfalls of ROC analysis in clinical microarray research (and how to avoid them). *Brief. Bioinform.* 13, 83–97. doi:10.1093/bib/bbr008
14. Rooney, M., Shukla, S., Wu, C., Getz, G. and Hacohen, N. (2015). Molecular and genetic properties of tumors associated with local immune cytolytic activity. *Cell* 160, 48–61. doi:10.1016/j.cell.2014.12.033
15. Liang, Y., Chhuo, L., Argha, A., Farbehi, N., Chen, L., Alizadehsani, R., Hosseinzadeh, M., Beheshti, A., Porntaveetusm, T., Ye, Y., et al. (2026). Transcriptomic models for immunotherapy response prediction show limited cross-cohort generalisability. *arXiv* 2604.05478 [q-bio.GN]. doi:10.48550/arXiv.2604.05478
16. Shen, W., Moon, I., Nguyen, T.H., Li, M.M., Huang, Y., Nair, N., Marbach, D. and Zitnik, M. (2026). Generalizable AI predicts immunotherapy outcomes across cancers and treatments. *Nat. Med.* 32, 3010–3022. doi:10.1038/s41591-026-04502-7
17. El Kanbi, K., Cattan, Y., Marschall, P. et al. (2026). A transcriptomic benchmark for foundation models in immunology and inflammation drug development. *Learning Meaningful Representations of Life (LMRL) Workshop at ICLR 2026* (spotlight).
18. Pal, L.R., Gertz, E.M., Nair, N.U., Mukherjee, S., Patiyal, S., Cantore, T., Campagnolo, E.M., Chang, T.G., Dhruba, S.R., Kim, Y., et al. (2025). A machine learning framework for supervised treatment response prediction from tumor transcriptomics: a large-scale pan-cancer study. *bioRxiv*. doi:10.1101/2025.10.24.684491
19. Mariathasan, S., Turley, S.J., Nickles, D., Castiglioni, A., Yuen, K., Wang, Y., Kadel, E.E., III, Koeppen, H., Astarita, J.L., Cubas, R., et al. (2018). TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. *Nature* 554, 544–548. doi:10.1038/nature25501
20. Lauss, M., Donia, M., Harbst, K., Andersen, R., Mitra, S., Rosengren, F., Salim, M., Vallon-Christersson, J., Törngren, T., Kvist, A., et al. (2017). Mutational and putative neoantigen load predict clinical benefit of adoptive T cell therapy in melanoma. *Nat. Commun.* 8, 1738. doi:10.1038/s41467-017-01460-0
21. Gide, T.N., Quek, C., Menzies, A.M., Tasker, A.T., Shang, P., Holst, J., Madore, J., Lim, S.Y., Velickovic, R., Wongchenko, M., et al. (2019). Distinct immune cell populations define response to anti-PD-1 monotherapy and anti-PD-1/anti-CTLA-4 combined therapy. *Cancer Cell* 35, 238–255.e6. doi:10.1016/j.ccell.2019.01.003
22. Goodman, A.M., Kato, S., Bazhenova, L., Patel, S.P., Frampton, G.M., Miller, V., Stephens, P.J., Daniels, G.A. and Kurzrock, R. (2017). Tumor mutational burden as an independent predictor of response to immunotherapy in Diverse Cancers. *Mol Cancer Ther*, 16, 2598–2608. doi:10.1158/1535-7163.MCT-17-0386
23. Newman, A.M., Liu, C.L., Green, M.R., Gentles, A.J., Feng, W., Xu, Y., Hoang, C.D., Diehn, M. and Alizadeh, A.A. (2015). Robust enumeration of cell subsets from tissue expression profiles. *Nat. Methods* 12, 453–457. doi:10.1038/nmeth.3337
24. Jung, H., Kim, H.S., Kim, J.Y., Sun, J.M., Ahn, J.S., Ahn, M.J., Park, K., Esteller, M., Lee, S.H. and Choi, J.K. (2019). DNA methylation loss promotes immune evasion of tumours with high mutation and copy number load. *Nat. Commun.* 10, 4278. doi:10.1038/s41467-019-12159-9
25. Subramanian, A., Tamayo, P., Mootha, V.K., Mukherjee, S., Ebert, B.L., Gillette, M.A., Paulovich, A., Pomeroy, S.L., Golub, T.R., Lander, E.S., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl. Acad. Sci. U.S.A.* 102, 15545–15550. doi:10.1073/pnas.0506580102
26. Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B*, 57, 289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
27. DeLong, E.R., DeLong, D.M. and Clarke-Pearson, D.L. (1988). Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. *Biometrics*, 44, 837–845. doi:10.2307/2531595
28. Efron, B. and Tibshirani, R.J. (1993). *An Introduction to the Bootstrap*. Chapman & Hall, New York.
29. Zhao, C., Shi, L., Tong, W., Shaughnessy, J.D., Oberthuer, A., Pusztai, L., Deng, Y., Symmans, W.F. and Shi, T. (2011). Maximum predictive power of the microarray-based models for clinical outcomes is limited by correlation between endpoint and gene expression profile. *BMC Genomics* 12(Suppl 5), S3. doi:10.1186/1471-2164-12-S5-S3
30. Carter, J.A., Gilbo, P. and Atwal, G.S. (2019). IMPRES does not reproducibly predict response to immune checkpoint blockade therapy in metastatic melanoma. *Nat Med*, 25, 1833–1835. doi:10.1038/s41591-019-0671-4
31. Kapoor, S. and Narayanan, A. (2023). Leakage and the reproducibility crisis in machine-learning-based science. *Patterns* 4, 100804. doi:10.1016/j.patter.2023.100804
32. Usset, J., Rosendahl Huber, A., Andrianova, M.A., Batlle, E., Carles, J., Cuppen, E., Elez, E., Felip, E., Gómez-Rey, M., Lo Giacco, D., et al. (2024). Five latent factors underlie response to immunotherapy. *Nat. Genet.* 56, 2112–2120. doi:10.1038/s41588-024-01899-0
33. Liu, D., Schilling, B., Liu, D., Sucker, A., Livingstone, E., Jerby-Arnon, L., Zimmer, L., Gutzmer, R., Satzger, I., Loquai, C., et al. (2019). Integrative molecular and clinical modeling of clinical outcomes to PD1 blockade in patients with metastatic melanoma. *Nat. Med.* 25, 1916–1927. doi:10.1038/s41591-019-0654-5
34. Van Allen, E.M., Miao, D., Schilling, B., Shukla, S.A., Blank, C., Zimmer, L., Sucker, A., Hillen, U., Geukes Foppen, M.H., Goldinger, S.M., et al. (2015). Genomic correlates of response to CTLA-4 blockade in metastatic melanoma. *Science*, 350, 207–211. doi:10.1126/science.aad0095
