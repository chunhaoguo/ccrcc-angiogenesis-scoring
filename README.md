# Angiogenesis score implementation in renal cancer trials

Research code accompanying the secondary analysis of angiogenesis gene-expression scoring and patient classification in renal cell carcinoma trials.

## Contents and reproducibility

- `scoring/operators.R`: the seven scoring operators used in the scoring-method comparison.
- `scoring/frozen_scorer.R`: functions for applying a qualified frozen reference.
- `analysis/`: analysis scripts for clinical outcomes, trial scoring, renal cohorts, biology and post hoc construct anchors. These historical source scripts document analyses and require adapting `/LOCAL_PROJECT_ROOT`, `/LOCAL_DATA_ROOT`, and `/LOCAL_EXTERNAL_DRIVE` to locally acquired, appropriately authorized source data and the original folder layout. This repository is not a turnkey reproduction of every manuscript result; not all upstream inputs, pipelines or parameter files are included.
- `test/test_synthetic.R`: a self-contained check using generated synthetic numbers. Run `Rscript test/test_synthetic.R` from the repository root. Additional operators require the R packages named in the source, including GSVA and BiocParallel.

No patient-level data, actual expression matrices, participant identifiers, credentials, or fitted reference parameters are provided. Redistribution rights for fitted references are being verified; absence of these parameters means the frozen scorer cannot be used as a complete frozen assay from this repository alone. Data must be obtained under their original access terms. High/low cutoffs and clinical use are not validated. The software is for research.

## Manuscript-to-code guide

| Manuscript item | Relevant code |
| --- | --- |
| Figure 1 and scoring comparisons | `scoring/operators.R`; `analysis/B1_operator_panel_trials/`; `analysis/B2_independent_cohorts/` |
| Figure 2 and clinical outcomes | `analysis/A_clinical_consequence/` |
| Figure 3 and biological context | `analysis/D_biology/`; `analysis/F_construct_anchors/` |
| Figure 4 and additional cohort analysis | `analysis/B2_independent_cohorts/`; `scoring/frozen_scorer.R` |

The original inputs and some upstream plotting and cellular workflows are not included; the mapping identifies relevant code rather than asserting a complete independent reproduction.

## Authorship

Chunhao Guo and coauthors. Source datasets and third-party libraries retain their original terms. No software license is specified pending an author decision.
