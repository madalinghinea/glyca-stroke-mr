# run_all.R -- runs the pipeline in order from the project root.
#   Rscript run_all.R                          offline stages from data/processed (~45 min, mostly MR-PRESSO)
#   GLYCA_RUN_EXTRACTION=1 Rscript run_all.R   also rebuild the inputs (needs OPENGWAS_JWT and downloads)
#   GLYCA_PRESSO_N=1000 Rscript run_all.R      smoke test; overwrites data/processed/mr_*.rds (manuscript = 10000)
#   GLYCA_STEPS="06_power.qmd,10_manuscript_numbers.qmd" Rscript run_all.R   selected steps only
#   GLYCA_RAW_DIR=/path/to/raw ...             folder with the downloaded GWAS Catalog files
# Each .qmd is converted with knitr::purl() and sourced, so Quarto is not required to run the pipeline.
# 11_reverse_and_af_mvmr_token.qmd needs an OpenGWAS token and runs only when named in GLYCA_STEPS.
steps <- c("01_extract_exposure.qmd", "02_outcome_megastroke.qmd", "03_outcome_gigastroke.qmd", "04_subtype_gigastroke.qmd",
           "05_mvmr.qmd", "06_power.qmd", "07_robustness.qmd", "08_tables.qmd", "09_figures.qmd", "10_manuscript_numbers.qmd")
sel <- Sys.getenv("GLYCA_STEPS"); if (nzchar(sel)) steps <- trimws(strsplit(sel, ",")[[1]])
stopifnot(file.exists("00_setup.R"), all(file.exists(steps)))
dir.create("logs", showWarnings = FALSE)

log_lines <- character()
say <- function(...) { l <- paste0(...); cat(l, "\n", sep = ""); log_lines <<- c(log_lines, l) }
say("Run started: ", format(Sys.time()))
for (f in steps) {
  say("\n######## ", f, " ########")
  tmp <- tempfile(fileext = ".R"); knitr::purl(f, output = tmp, documentation = 0, quiet = TRUE)
  t0 <- Sys.time(); ok <- TRUE
  out <- capture.output(ok <- tryCatch({ source(tmp, local = new.env(), echo = FALSE, print.eval = TRUE); TRUE },
                                       error = function(e) { cat("ERROR in", f, ":", conditionMessage(e), "\n"); FALSE }))
  if (length(out)) { cat(out, sep = "\n"); log_lines <- c(log_lines, out) }
  say(sprintf("-- %s %s (%.1f min)", f, if (ok) "finished" else "FAILED", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (!ok) break
}
say("\nRun ended: ", format(Sys.time()))
writeLines(log_lines, "logs/run_all.log")
if (!ok) { if (interactive()) stop("pipeline step failed - see logs/run_all.log") else quit(status = 1) }
