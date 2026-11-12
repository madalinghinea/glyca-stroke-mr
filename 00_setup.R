# 00_setup.R -- configuration, study metadata and helper functions shared by every script.
# Run all scripts from the project root; each script starts with source("00_setup.R").

SEED <- 20261105L                                                   # set before every stochastic call
DIR  <- list(raw = Sys.getenv("GLYCA_RAW_DIR", "data/raw"), proc = "data/processed",
             tab = "tables", fig = "figures", log = "logs")
for (d in DIR) dir.create(d, recursive = TRUE, showWarnings = FALSE)

RUN_EXTRACTION <- identical(Sys.getenv("GLYCA_RUN_EXTRACTION"), "1")  # 1 = rebuild inputs (OpenGWAS token / raw files)
PRESSO_N       <- as.integer(Sys.getenv("GLYCA_PRESSO_N", "10000"))   # MR-PRESSO simulations; manuscript = 10000

# ---- Study metadata ---------------------------------------------------------------------------
N_EXPOSURE <- 115078L                                               # met-d-GlycA (UK Biobank, SD units)
OUTCOMES <- data.frame(
  key      = c("MEGASTROKE", "GIGASTROKE", "CES", "LAS", "SVS", "MEGA_CES"),
  label    = c("Any ischemic stroke (MEGASTROKE)", "Any ischemic stroke (GIGASTROKE)", "Cardioembolic stroke (GIGASTROKE)",
               "Large-artery stroke (GIGASTROKE)", "Small-vessel stroke (GIGASTROKE)", "Cardioembolic stroke (MEGASTROKE)"),
  short    = c("MEGASTROKE AIS", "GIGASTROKE AIS", "GIGASTROKE CES (GCST90104541)", "GIGASTROKE LAS", "GIGASTROKE SVS",
               "MEGASTROKE CES (GCST006910)"),
  id       = c("ebi-a-GCST006908", "GCST90104540", "GCST90104541", "GCST90104542", "GCST90104543", "GCST006910"),
  ncase    = c(34217L, 62100L, 10804L, 6399L, 6811L, 7193L),
  ncontrol = c(406111L, 1234808L, 1234808L, 1234808L, 1234808L, 355468L),
  stringsAsFactors = FALSE)
KEYS <- c("MEGASTROKE", "GIGASTROKE", "CES", "LAS", "SVS")           # outcomes run through the full battery

GCST_FTP  <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
RAW_FILES <- list(
  GIGASTROKE = list(file = "GCST90104540.h.tsv.gz", path = "GCST90104001-GCST90105000/GCST90104540/harmonised/GCST90104540.h.tsv.gz"),
  CES        = list(file = "GCST90104541.h.tsv.gz", path = "GCST90104001-GCST90105000/GCST90104541/harmonised/GCST90104541.h.tsv.gz"),
  LAS        = list(file = "GCST90104542.h.tsv.gz", path = "GCST90104001-GCST90105000/GCST90104542/harmonised/GCST90104542.h.tsv.gz"),
  SVS        = list(file = "GCST90104543.h.tsv.gz", path = "GCST90104001-GCST90105000/GCST90104543/harmonised/GCST90104543.h.tsv.gz"),
  MEGA_CES   = list(file = "GCST006910.h.tsv.gz",   path = "GCST006001-GCST007000/GCST006910/harmonised/29531354-GCST006910-EFO_1001976.h.tsv.gz"),
  AF         = list(file = "GCST006414.h.tsv.gz",   path = "GCST006001-GCST007000/GCST006414/harmonised/30061737-GCST006414-EFO_0000275.h.tsv.gz"))

# Instrument sets used in the robustness analyses (labels appear in Table S11 and in manuscript_numbers.csv)
LOCUS_ANCHORS <- list(HP = "rs77303550", LPA = "rs10455872", MHC = c("rs1005599", "rs9270074"))
SET_LABELS <- c(primary = "Primary (all retained SNPs)", hp = "Excluding HP locus (chr16 71.6-72.9 Mb)",
                lpa = "Excluding LPA locus (chr6 160.5-161.0 Mb)", hlm = "Excluding HP + LPA + MHC loci",
                steiger = "Excluding Steiger-FALSE SNPs", common = "Common SNP set across 5 outcomes",
                af_p = "Excluding SNPs with AF p < 1e-3", af_near = "Excluding SNPs within 500 kb of an AF lead SNP",
                af_both = "Excluding both")

# ---- Packages and session log -------------------------------------------------------------------
suppressPackageStartupMessages({ library(TwoSampleMR); library(dplyr); library(tidyr); library(ggplot2); library(data.table) })
log_session <- function() {
  key <- c("TwoSampleMR", "MRPRESSO", "MVMR", "mr.raps", "ieugwasr", "data.table", "dplyr", "ggplot2")
  ver <- vapply(key, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED", "")
  ns  <- sort(loadedNamespaces())
  nsv <- vapply(ns, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "?"), "")
  writeLines(c(R.version.string, paste("Platform:", R.version$platform), "", "Key packages:", paste(names(ver), ver),
               "", "Loaded namespaces:", paste(ns, nsv)), file.path(DIR$log, "sessionInfo.txt"))
  invisible(ver)
}
log_session()

# ---- Formatting and file helpers ------------------------------------------------------------------
or_ci  <- function(b, se) data.frame(OR = exp(b), OR_lci = exp(b - 1.96 * se), OR_uci = exp(b + 1.96 * se))
fmt_or <- function(b, se, digits = 3) { x <- or_ci(b, se); sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), x$OR, x$OR_lci, x$OR_uci) }
fmt_p  <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))

M <- c(ivw = "Inverse variance weighted", egger = "MR Egger", wm = "Weighted median", simple = "Simple mode", wmode = "Weighted mode")
pick_method <- function(tab, m) tab[tab$method == M[[m]], ]               # one method's row of mr() / mr_heterogeneity() output

outcome_meta  <- function(key) OUTCOMES[OUTCOMES$key == key, ]
outcome_short <- function(key) OUTCOMES$short[match(key, OUTCOMES$key)]

need_file <- function(path, hint) {
  if (!file.exists(path)) stop(sprintf("Missing input: %s\n  -> %s", path, hint), call. = FALSE)
  invisible(path)
}
RAW_SEARCH <- unique(c(DIR$raw, file.path("..", "data", "raw")))      # folders searched for an already downloaded file
get_raw <- function(key) {                                          # harmonised GWAS Catalog file; downloaded to DIR$raw if absent everywhere
  f <- file.path(DIR$raw, RAW_FILES[[key]]$file)
  hit <- file.path(RAW_SEARCH, RAW_FILES[[key]]$file); hit <- hit[file.exists(hit)]
  if (length(hit)) return(hit[1])
  if (!file.exists(f)) {
    message("Downloading ", RAW_FILES[[key]]$file)
    old <- options(timeout = 7200); on.exit(options(old))
    download.file(file.path(GCST_FTP, RAW_FILES[[key]]$path), f, mode = "wb", quiet = TRUE)
  }
  f
}

# ---- Instrument loci: single-linkage clusters of instruments (gap <= 1 Mb) -------------------------------
define_loci <- function(exp_dat, gap = 1e6) {
  exp_dat$pos.exposure <- as.numeric(exp_dat$pos.exposure)
  exp_dat$F <- (exp_dat$beta.exposure / exp_dat$se.exposure)^2
  exp_dat %>% arrange(chr.exposure, pos.exposure) %>% group_by(chr.exposure) %>%
    mutate(new = (pos.exposure - lag(pos.exposure, default = -1e12)) > gap) %>% ungroup() %>%
    mutate(locus = cumsum(new)) %>% group_by(locus) %>%
    mutate(n_snps = n(), locus_label = paste0("chr", first(chr.exposure), ":", round(min(pos.exposure) / 1e6, 2), "-",
                                              round(max(pos.exposure) / 1e6, 2), " Mb")) %>% ungroup() %>%
    select(SNP, chr.exposure, pos.exposure, F, locus, n_snps, locus_label)
}

# ---- Quality control of harmonised data and Steiger directionality -----------------------------------------
check_harmonised <- function(dat, label) {
  d <- dat[dat$mr_keep, ]
  stopifnot(!anyDuplicated(d$SNP), all(is.finite(c(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome))), all(d$se.outcome > 0))
  de <- abs(d$eaf.exposure - d$eaf.outcome); flag <- which(de > 0.15)
  cat(sprintf("[%s] %d SNPs retained | |EAF exposure - EAF outcome| > 0.15: %s\n", label, nrow(d),
              if (length(flag)) paste0(d$SNP[flag], " (", round(de[flag], 2), ")", collapse = ", ") else "none"))
  d
}
STEIGER_COLS <- c("rsq.exposure", "rsq.outcome", "effective_n.exposure", "effective_n.outcome", "steiger_dir", "steiger_pval")
add_steiger <- function(d, key) {                                   # explicit metadata: no OpenGWAS call
  d[intersect(names(d), STEIGER_COLS)] <- NULL                      # never reuse Steiger columns from an earlier run
  m <- outcome_meta(key)
  d$samplesize.exposure <- N_EXPOSURE; d$units.exposure <- "SD"
  d$ncase.outcome <- m$ncase; d$ncontrol.outcome <- m$ncontrol
  d$samplesize.outcome <- m$ncase + m$ncontrol
  d$prevalence.outcome <- m$ncase / (m$ncase + m$ncontrol)
  d$units.outcome <- "log odds"
  suppressMessages(steiger_filtering(d))
}

# ---- MR-PRESSO summary (outliers identified by rsID) ------------------------------------------------------
presso_summary <- function(presso, d) {
  pr   <- presso[[1]]$`MR-PRESSO results`
  main <- presso[[1]]$`Main MR results`
  dist <- pr$`Distortion Test`
  idx  <- if (!is.null(dist) && is.numeric(dist$`Outliers Indices`)) dist$`Outliers Indices` else integer(0)
  data.frame(
    global_rss = pr$`Global Test`$RSSobs, global_p = pr$`Global Test`$Pvalue,
    n_outliers = length(idx), outlier_snps = paste(d$SNP[idx], collapse = ";"),
    distortion_p = if (!is.null(dist) && !is.null(dist$Pvalue)) dist$Pvalue else NA_real_,
    raw_b = main$`Causal Estimate`[1], raw_se = main$Sd[1], raw_p = main$`P-value`[1],
    corr_b = main$`Causal Estimate`[2], corr_se = main$Sd[2], corr_p = main$`P-value`[2], stringsAsFactors = FALSE)
}

# ---- The analysis battery applied identically to every outcome ----------------------------------------------
run_battery <- function(dat, key, presso_n = PRESSO_N) {
  d <- add_steiger(check_harmonised(dat, key), key)
  set.seed(SEED); res <- mr(d)                                      # IVW (random effects), Egger, weighted median, simple/weighted mode
  het    <- mr_heterogeneity(d)
  pleio  <- mr_pleiotropy_test(d)
  single <- mr_singlesnp(d)
  loo    <- mr_leaveoneout(d)
  set.seed(SEED)
  presso <- suppressMessages(run_mr_presso(d, NbDistribution = presso_n))
  out <- list(key = key, meta = outcome_meta(key), dat = d, res = res, het = het, pleio = pleio, single = single, loo = loo,
              presso = presso, presso_sum = presso_summary(presso, d))
  saveRDS(out, file.path(DIR$proc, paste0("mr_", key, ".rds")))
  out
}
print_battery <- function(b) {
  ivw <- pick_method(b$res, "ivw")
  cat(sprintf("\n[%s] %d SNPs | IVW OR %s, p = %s\n", b$key, ivw$nsnp, fmt_or(ivw$b, ivw$se), fmt_p(ivw$pval)))
  cat(sprintf("  Cochran Q (IVW) p = %s | Egger intercept %.4f, p = %s\n", fmt_p(pick_method(b$het, "ivw")$Q_pval), b$pleio$egger_intercept, fmt_p(b$pleio$pval)))
  ps <- b$presso_sum
  cat(sprintf("  MR-PRESSO global p = %s | outliers: %s\n", fmt_p(ps$global_p), ifelse(ps$n_outliers == 0, "none", ps$outlier_snps)))
  cat(sprintf("  Steiger FALSE: %s\n", ifelse(all(b$dat$steiger_dir), "none", paste(b$dat$SNP[!b$dat$steiger_dir], collapse = ", "))))
  invisible(b)
}
load_battery <- function(key) {
  f <- file.path(DIR$proc, paste0("mr_", key, ".rds"))
  need_file(f, "run the outcome script (02, 03 or 04) first")
  readRDS(f)
}
load_batteries <- function(keys = KEYS) setNames(lapply(keys, load_battery), keys)
