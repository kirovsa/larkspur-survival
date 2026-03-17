# tests/test_utils.R
# Unit tests for R/utils.R
# Run with: Rscript tests/test_utils.R
# Requires: survival, dplyr, tidyr, readr, yaml, tibble, purrr, stringr, ggplot2, survminer

# Resolve repo root from script path or working directory
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 1 && nchar(script_path) > 0) {
  repo_root <- normalizePath(file.path(dirname(script_path), ".."))
} else {
  repo_root <- getwd()  # fall back to working directory
}

cat("Repo root:", repo_root, "\n")
cat("Loading utils.R...\n")
source(file.path(repo_root, "R", "utils.R"))

# ── Helper: minimal pass/fail ─────────────────────────────────────────────────
pass <- 0L
fail <- 0L

expect <- function(desc, expr) {
  result <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("  ERROR in '%s': %s\n", desc, conditionMessage(e)))
    FALSE
  })
  if (result) {
    cat(sprintf("  PASS: %s\n", desc))
    pass <<- pass + 1L
  } else {
    cat(sprintf("  FAIL: %s\n", desc))
    fail <<- fail + 1L
  }
}

# ── Test 1: parse_cbio_status ─────────────────────────────────────────────────
cat("\n[1] parse_cbio_status\n")
expect("DECEASED -> 1",       parse_cbio_status("1:DECEASED")       == 1L)
expect("LIVING -> 0",         parse_cbio_status("0:LIVING")         == 0L)
expect("PROGRESSED -> 1",     parse_cbio_status("1:PROGRESSED")     == 1L)
expect("NOT PROGRESSED -> 0", parse_cbio_status("0:NOT PROGRESSED") == 0L)
expect("bare '1' -> 1",       parse_cbio_status("1")                == 1L)
expect("bare '0' -> 0",       parse_cbio_status("0")                == 0L)
expect("vectorised",          all(parse_cbio_status(c("1:DECEASED", "0:LIVING")) == c(1L, 0L)))

# ── Test 2: read_cbio_clinical ────────────────────────────────────────────────
cat("\n[2] read_cbio_clinical\n")
clin_path <- file.path(repo_root, "cbioportal", "brca_tcga_pub", "data_clinical_patient.txt")
expect("clinical file exists", file.exists(clin_path))
clin <- read_cbio_clinical(clin_path)
expect("returns tibble",         inherits(clin, "tbl_df"))
expect("PATIENT_ID present",     "PATIENT_ID"  %in% names(clin))
expect("OS_STATUS present",      "OS_STATUS"   %in% names(clin))
expect("OS_MONTHS present",      "OS_MONTHS"   %in% names(clin))
expect("PFS_STATUS present",     "PFS_STATUS"  %in% names(clin))
expect("PFS_MONTHS present",     "PFS_MONTHS"  %in% names(clin))
expect("150 patients",           nrow(clin) == 150L)
expect("no comment rows",        !any(startsWith(clin$PATIENT_ID, "#")))

# ── Test 3: read_cbio_expression ──────────────────────────────────────────────
cat("\n[3] read_cbio_expression\n")
expr_path <- file.path(repo_root, "cbioportal", "brca_tcga_pub", "data_mrna_seq_v2_rsem.txt")
expect("expression file exists", file.exists(expr_path))
expr <- read_cbio_expression(expr_path)
expect("returns tibble",         inherits(expr, "tbl_df"))
expect("Hugo_Symbol present",    "Hugo_Symbol"    %in% names(expr))
expect("Entrez_Gene_Id present", "Entrez_Gene_Id" %in% names(expr))
expect("15 genes",               nrow(expr) == 15L)
expect("152 columns (2+150)",    ncol(expr) == 152L)

# ── Test 4: compute_signature_score ───────────────────────────────────────────
cat("\n[4] compute_signature_score\n")
scores <- compute_signature_score(expr, c("BRCA1", "BRCA2", "TP53"), method = "mean")
expect("returns named numeric",  is.numeric(scores) && !is.null(names(scores)))
expect("150 scores",             length(scores) == 150L)
expect("names are sample IDs",   all(startsWith(names(scores), "TCGA-BX-")))
expect("finite values",          all(is.finite(scores)))

scores_med <- compute_signature_score(expr, c("MKI67", "CCND1"), method = "median")
expect("median method works",    is.numeric(scores_med) && length(scores_med) == 150L)

# Check that a warning is issued for a missing gene but still returns scores
w_caught <- FALSE
withCallingHandlers(
  {
    scores_w <- compute_signature_score(expr, c("BRCA1", "NOTREAL"))
  },
  warning = function(w) {
    w_caught <<- TRUE
    invokeRestart("muffleWarning")
  }
)
expect("warns on missing gene",          w_caught)
expect("still returns scores despite missing gene", is.numeric(scores_w))

# ── Test 5: prepare_survival_data ─────────────────────────────────────────────
cat("\n[5] prepare_survival_data\n")
sd_os  <- prepare_survival_data(clin, scores, endpoint = "OS")
expect("returns tibble",          inherits(sd_os, "tbl_df"))
expect("columns present",         all(c("patient_id","time","status","score","score_group") %in% names(sd_os)))
expect("status is 0/1 integer",   all(sd_os$status %in% c(0L, 1L)))
expect("time >= 0",               all(sd_os$time >= 0))
expect("score_group factor",      is.factor(sd_os$score_group))
expect("levels Low/High",         all(levels(sd_os$score_group) == c("Low", "High")))
expect("roughly balanced groups", {
  tbl <- table(sd_os$score_group)
  # With >= median split on 150 patients the groups should differ by at most 2
  abs(tbl["High"] - tbl["Low"]) <= 2L
})

sd_pfs <- prepare_survival_data(clin, scores, endpoint = "PFS")
expect("PFS endpoint works",      nrow(sd_pfs) > 0)

# ── Test 6: run_cox_analysis ──────────────────────────────────────────────────
cat("\n[6] run_cox_analysis\n")
cox_res <- run_cox_analysis(sd_os)
expect("returns list",            is.list(cox_res))
expect("cox_continuous present",  inherits(cox_res$cox_continuous, "coxph"))
expect("cox_binary present",      inherits(cox_res$cox_binary, "coxph"))
expect("ph_test present",         !is.null(cox_res$ph_test))
expect("summary_tbl tibble",      inherits(cox_res$summary_tbl, "tbl_df"))
expect("summary has HR column",   "hr" %in% names(cox_res$summary_tbl))
expect("HR is positive",          all(cox_res$summary_tbl$hr > 0))
expect("p_value in [0,1]",        all(cox_res$summary_tbl$p_value >= 0 &
                                        cox_res$summary_tbl$p_value <= 1))

# ── Test 7: fmt_pval ──────────────────────────────────────────────────────────
cat("\n[7] fmt_pval\n")
expect("p<0.001 -> '<0.001'",  fmt_pval(0.0001) == "<0.001")
expect("p=0.005 -> '0.005'",   fmt_pval(0.005)  == "0.005")
expect("p=0.05  -> '0.05'",    fmt_pval(0.05)   == "0.05")
expect("p=0.5   -> '0.50'",    fmt_pval(0.5)    == "0.50")

# ── Test 8: build_summary_table ───────────────────────────────────────────────
cat("\n[8] build_summary_table\n")
mock_results <- list(
  "Breast Cancer (TCGA)" = list(
    OS = list(
      proliferation = c(
        cox_res,
        list(surv_data = sd_os, sig_label = "Proliferation", endpoint = "OS")
      )
    )
  )
)
stbl <- build_summary_table(mock_results)
expect("returns tibble",      inherits(stbl, "tbl_df"))
expect("has Cohort column",   "Cohort"   %in% names(stbl))
expect("has HR column",       "HR"       %in% names(stbl))
expect("has P_value column",  "P_value"  %in% names(stbl))
expect("at least one row",    nrow(stbl) >= 1L)

# ── Summary ───────────────────────────────────────────────────────────────────
cat(sprintf("\n===========================\n%d passed, %d failed\n", pass, fail))
if (fail > 0L) quit(status = 1L)
