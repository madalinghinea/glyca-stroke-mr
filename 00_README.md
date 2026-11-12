# GlycA and ischemic stroke - two-sample Mendelian randomization pipeline

Reproduces every number, table and figure of the manuscript. Run from the project root (this folder).

## Quick start

```r
source("run_all.R")          # offline stages from data/processed; ~45 min (MR-PRESSO, 10,000 simulations x 6 data sets)
Sys.setenv(GLYCA_STEPS = "06_power.qmd,10_manuscript_numbers.qmd"); source("run_all.R")   # selected steps only
```

Environment variables (all optional)

| Variable | Default | Effect |
|---|---|---|
| `GLYCA_RAW_DIR` | `data/raw` | folder for the downloaded GWAS Catalog files (fetched automatically when missing in this folder and in `../data/raw`, about 2.5 GB) |
| `GLYCA_PRESSO_N` | `10000` | MR-PRESSO simulations; lower values overwrite `data/processed/mr_*.rds` and change the PRESSO p-values |
| `GLYCA_RUN_EXTRACTION` | unset | `1` rebuilds the inputs that need the OpenGWAS token or the raw files |
| `GLYCA_STEPS` | all core steps | comma-separated list of scripts to run |

OpenGWAS token: set `OPENGWAS_JWT` in `.Renviron`; never write it into a script. It is needed only for stage A of 01, 02 and 05
and for script 11.

Required R packages: TwoSampleMR 0.7.13, MRPRESSO 1.0, MVMR 0.4.9, mr.raps 0.4.3, ieugwasr 1.2.0, data.table, dplyr, tidyr,
ggplot2, knitr. `00_setup.R` writes the exact versions to `logs/sessionInfo.txt` on every run.

## Pipeline

| Step | File | Content | Main outputs |
|---|---|---|---|
| 00 | `00_setup.R` | seed, directories, study metadata, outcome table, helpers, `run_battery()` | `logs/sessionInfo.txt` |
| 01 | `01_extract_exposure.qmd` | GlycA instruments (`met-d-GlycA`) | `exp_dat.rds` |
| 02 | `02_outcome_megastroke.qmd` | MR against MEGASTROKE | `mr_MEGASTROKE.rds` |
| 03 | `03_outcome_gigastroke.qmd` | MR against GIGASTROKE | `mr_GIGASTROKE.rds` |
| 04 | `04_subtype_gigastroke.qmd` | cardioembolic, large-artery and small-vessel stroke (exploratory) | `mr_CES/LAS/SVS.rds` |
| 05 | `05_mvmr.qmd` | MVMR with BMI and smoking initiation, conditional F, Q_A | `mvmr_summary.rds` |
| 06 | `06_power.qmd` | idealised and empirical minimum detectable effects | Table S8 |
| 07 | `07_robustness.qmd` | locus exclusion, common SNP set, Steiger filter, leave-one-locus-out, MR-RAPS, MEGASTROKE CES check, atrial fibrillation | Tables S11-S14, Figure S5, `robustness_summary.rds` |
| 08 | `08_tables.qmd` | supplementary tables S1-S7b, S9a-e, S10 | `tables/*.csv` |
| 09 | `09_figures.qmd` | Figures 1-5 and S1-S4 (PDF and 300-dpi PNG) | `figures/*` |
| 10 | `10_manuscript_numbers.qmd` | every number quoted in the manuscript | `tables/manuscript_numbers.csv` |
| 11 | `11_reverse_and_af_mvmr_token.qmd` | reverse MR (stroke to GlycA) and GlycA + AF MVMR; needs the token; runs only when named in `GLYCA_STEPS` | `reverse_mr.rds`, `mvmr_glyca_af_ces.rds` |

Each analysis script has Stage A (extraction or harmonisation, runs when its input file is missing or `GLYCA_RUN_EXTRACTION=1`)
and Stage B (analysis, offline). Scripts 06, 08, 09 and 10 only read saved results.

## Design rules

1. One seed (`SEED = 20261105`) is set immediately before every stochastic call (bootstraps inside `mr()`, MR-PRESSO, MR-RAPS).
2. All outcomes pass through `run_battery()`: IVW with random effects, MR-Egger, weighted median, simple and weighted mode,
   Cochran Q, Egger intercept, single-SNP and leave-one-out analyses, MR-PRESSO, Steiger filtering.
3. Steiger filtering uses explicit study metadata in `00_setup.R`; stale Steiger columns in an input file are dropped first,
   so the analysis stages need no OpenGWAS token.
4. MR-PRESSO outliers are identified by rsID, never by row index.
5. `check_harmonised()` stops on duplicated SNPs or non-finite effects and reports instruments whose allele frequency differs
   by more than 0.15 between exposure and outcome.

## Status (11 November 2026)

Version 10 is a refactoring of the validated v9 pipeline (renumbered scripts, one place for study metadata, stale-cache and stale-column fixes, computed figure labels). It has been run end to end in R 4.6.0 on the author's machine: steps 01-06 on 11 November 13:34-13:57, steps 07-10 at 15:12 (after the `pick()` / `dplyr::pick()` name clash of the first attempt was fixed) and step 11 at 15:19. All steps finished; the results reproduce the manuscript numbers (see `tools/CODE_REVIEW.md`, findings 11-13 and the execution status).

To re-run only part of the pipeline: `Sys.setenv(GLYCA_STEPS = "07_robustness.qmd,08_tables.qmd,09_figures.qmd,10_manuscript_numbers.qmd"); source("run_all.R")`; step 11 needs `OPENGWAS_JWT` and a valid token, and each call of `run_all.R` overwrites `logs/run_all.log`.

Independent checks (Python, no R): `tools/independent_check.py` (40 comparisons from the tables), `tools/independent_raw_check.py` (harmonisation and AF analyses from the raw GWAS Catalog files; 31 comparisons), `tools/af_update_python.py` (AF lead SNPs and AF to stroke estimates).

## Known limitations

* GIGASTROKE and subtype instruments are matched by exact rsID (no LD proxies); MEGASTROKE uses the OpenGWAS proxy search.
* AF lead SNPs are positionally pruned (largest chi-square per 1-Mb window), not LD-clumped.
* MR-PRESSO global-test p-values are Monte Carlo estimates (resolution about 1e-4 at 10,000 simulations).
* On Windows, MRPRESSO and mr.raps may need conda-installed `r-rootsolve` and `r-nortest` if R refuses to load DLLs from a user library.
* `data/processed` of older versions contained obsolete intermediates (`dat_megastroke.rds`, `dat_gigastroke.rds`, `het_*.rds`,
  `loo_*.rds`, `pleio_*.rds`, `presso_*.rds`, `res_*.rds`, `res_mv.rds`); the pipeline no longer reads them and they can be deleted.
