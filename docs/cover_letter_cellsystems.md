# Cover Letter — Cell Systems submission

Dear Editors,

We submit "SPATBench: Methodological Choices Can Outweigh Algorithmic Choices
in Transcriptomic Immunotherapy Prediction Benchmarking at Current Sample
Sizes" for consideration as an Article in Cell Systems.

The manuscript reports a quantitative audit of how benchmarking design
choices — not algorithms — drive conclusions in immunotherapy transcriptomic
prediction:

1. **Endpoint definition collapses rankings.** On the same patients, switching
   from the gene-defined cytolytic surrogate to clinical RECIST v1.1 drops
   ElasticNet-Var from 0.985 to 0.465 and GEP from 0.968 to 0.503, while IMPRES
   (no feature overlap) is unchanged; the between-endpoint effect is 3.5x the
   largest between-method gap.
2. **Feature selection rewrites the biological narrative.** Variance- and
   MI-based selection yield 27% overlapping gene sets (2.4% without the shared
   prescreen) with indistinguishable performance.
3. **A pre-benchmarking diagnostic.** The Circularity Detection Score (CDS,
   open-source, v1.2.1) places the circular endpoint far above its
   cohort-matched null while its categorical labels prove non-specific — the
   actionable protection is procedural (declare the endpoint genes).

High-resolution validation: ElasticNet (MI) on Gide 2019 (n = 73) becomes
FDR-significant at 5,000 full-pipeline shuffles (p = 0.0232), while no
trainable method survives FDR on any clinical endpoint at n <= 43.

We believe this fits Cell Systems' interest in quantitative reasoning about
scientific practice: the study treats benchmark methodology itself as a
measurable system, with a released diagnostic (CDS) and a Docker-reproducible
benchmark (SPATBench) as its instruments.

Suggested reviewers (no conflicts): to be completed by the authors
(at least six names, per Cell Systems requirements).

The authors declare no conflict of interest. All data are public (GEO,
cBioPortal, TCGA); code is open-source (MIT).
