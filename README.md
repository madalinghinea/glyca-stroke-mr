# GlycA and ischemic stroke: two-sample Mendelian randomization

Analysis code and results for the study of circulating glycoprotein acetyls (GlycA) as a causal risk
factor for ischemic stroke and its subtypes. Everything reported in the article is produced by the
scripts in this repository and exported to a single table, `tables/manuscript_numbers.csv`.

## Main result

| Outcome | Instruments | IVW odds ratio per 1-SD GlycA (95% CI) | p |
|---|---|---|---|
| Ischemic stroke, MEGASTROKE (ebi-a-GCST006908) | 54 | 0.979 (0.901-1.064) | 0.623 |
| Ischemic stroke, GIGASTROKE (GCST90104540) | 45 | 0.970 (0.893-1.053) | 0.463 |
| Cardioembolic stroke, GIGASTROKE (GCST90104541), exploratory | 41 | 0.853 (0.754-0.966) | 0.0118 |

The two overall analyses are null. The cardioembolic estimate is exploratory, does not survive
correction for the five outcome analyses, and is not reproduced by a published analysis of the same
GIGASTROKE data; it is reported as a hypothesis, not a finding.

## Data

No GWAS summary statistics are stored here. Download them from the original sources:

| Dataset | Accession | Source |
|---|---|---|
| GlycA (exposure, 60 instruments) | `met-d-GlycA` | IEU OpenGWAS, https://gwas.mrcieu.ac.uk/ (access token required) |
| Ischemic stroke, MEGASTROKE | `ebi-a-GCST006908` | IEU OpenGWAS |
| Ischemic stroke, GIGASTROKE (overall and subtypes) | GCST90104540, GCST90104541, GCST90104542, GCST90104543 | GWAS Catalog, https://www.ebi.ac.uk/gwas/ |
| Cardioembolic stroke, MEGASTROKE | GCST006910 | GWAS Catalog |
| Atrial fibrillation | GCST006414 | GWAS Catalog |
| Body mass index; smoking initiation (multivariable MR) | `ieu-b-40`, `ieu-b-4877` | IEU OpenGWAS |

`00_setup.R` downloads the GWAS Catalog files automatically into `data/raw/` when they are missing
(about 2.8 GB). OpenGWAS queries need a personal token in the environment variable `OPENGWAS_JWT`
(put it in `.Renviron`, never in a script).

## How to run

```r
# from this folder, in R 4.6.0
source("run_all.R")                                                    # steps 01-10, about 45 min (MR-PRESSO dominates)
Sys.setenv(GLYCA_STEPS = "07_robustness.qmd,08_tables.qmd")            # or selected steps only
Sys.setenv(GLYCA_STEPS = "11_reverse_and_af_mvmr_token.qmd")           # reverse MR and GlycA + AF MVMR; needs OPENGWAS_JWT
```

Packages: TwoSampleMR 0.7.13, MRPRESSO 1.0, MVMR 0.4.9, mr.raps 0.4.3, ieugwasr 1.2.0.9000,
data.table, dplyr, tidyr, ggplot2, knitr. The exact versions of the run behind the article are in
`logs/sessionInfo.txt`. Quarto is not required: `run_all.R` converts each `.qmd` with `knitr::purl()`
and sources it. One seed (20261105) is set before every stochastic call.

## Layout

| Path | Content |
|---|---|
| `00_setup.R` | study metadata, helper functions, the analysis battery applied to every outcome |
| `01`-`10` | instruments, MEGASTROKE, GIGASTROKE, subtypes, multivariable MR, power, robustness and atrial fibrillation, tables, figures, number export |
| `11_reverse_and_af_mvmr_token.qmd` | reverse-direction MR and joint GlycA + atrial-fibrillation MVMR (OpenGWAS token) |
| `tables/` | every supplementary table (S1-S17) and `manuscript_numbers.csv` |
| `tables/PROVENANCE.md` | which table came from R and which from the Python tools |
| `figures/` | Figures 1-5 and S1-S5 (PDF and 300-dpi PNG) |
| `tools/` | independent Python re-implementations of the estimators and the OpenGWAS API analyses |
| `tools/CODE_REVIEW.md` | code review, known issues and the execution status of each script |
| `data/processed/` | harmonised data and saved model objects, so the analysis steps run offline |
| `logs/` | run log and `sessionInfo.txt` |
| `tables_baseline_AF130/` | superseded atrial-fibrillation tables from an earlier lead-SNP list, kept for provenance |

## Verification

`tools/independent_check.py` recomputes IVW, MR-Egger, Cochran Q, the Egger intercept and the
MR-PRESSO residual sum of squares from the SNP-level tables without R (40 comparisons, no
mismatches). `tools/independent_raw_check.py` repeats the harmonisation and the atrial-fibrillation
analyses from the raw GWAS Catalog files. `tables/TableS17c_MVMR_code_validation.csv` checks the
Python multivariable code against the R results.

## License and citation

Code: MIT (see `LICENSE`). Cite the article and, for the code, the archived release (see
`CITATION.cff`).
