# =============================================================================
# utils.R — Helper functions for cBioPortal data loading and survival analysis
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(survival)
  library(survminer)
  library(yaml)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(stringr)
})

# Null-coalescing operator (available in rlang, defined here for portability)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ── Data Loading ──────────────────────────────────────────────────────────────

#' Read a cBioPortal clinical file (skips the 4 comment lines)
#'
#' @param path  Path to the data_clinical_patient.txt file.
#' @return A tibble with one row per patient.
read_cbio_clinical <- function(path) {
  # Count comment lines (lines starting with #)
  raw_lines <- readLines(path)
  n_comment <- sum(startsWith(raw_lines, "#"))

  read_tsv(
    path,
    skip      = n_comment,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
}

#' Read a cBioPortal mRNA expression file (Hugo_Symbol + Entrez_Gene_Id + samples)
#'
#' @param path  Path to the data_mrna_seq_v2_rsem.txt file.
#' @return A tibble with genes as rows and samples as columns
#'         (Hugo_Symbol and Entrez_Gene_Id retained).
read_cbio_expression <- function(path) {
  read_tsv(path, col_types = cols(.default = col_character()),
           show_col_types = FALSE)
}

# ── Signature Scoring ─────────────────────────────────────────────────────────

#' Compute per-sample signature scores from an expression matrix
#'
#' Expression values are z-scored per gene across all samples before
#' aggregation to make scores comparable across cohorts.
#'
#' @param expr_tbl    Tibble from read_cbio_expression().
#' @param genes       Character vector of Hugo gene symbols.
#' @param method      One of "mean", "median", or "pc1".
#' @return Named numeric vector of signature scores (names = sample IDs).
compute_signature_score <- function(expr_tbl, genes, method = "mean") {
  method <- match.arg(method, c("mean", "median", "pc1"))

  # Filter to requested genes (warn if any are missing)
  missing <- setdiff(genes, expr_tbl$Hugo_Symbol)
  if (length(missing) > 0) {
    warning(sprintf(
      "Signature gene(s) not found in expression data: %s",
      paste(missing, collapse = ", ")
    ))
  }

  sub <- expr_tbl %>%
    filter(Hugo_Symbol %in% genes) %>%
    select(-Hugo_Symbol, -any_of("Entrez_Gene_Id"))

  if (nrow(sub) == 0) {
    stop("None of the signature genes were found in the expression data.")
  }

  # Convert to numeric matrix (genes × samples)
  mat <- sub %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()

  # Z-score each gene across samples
  mat_z <- t(scale(t(mat)))  # scale works on columns, so transpose

  # Aggregate across genes
  scores <- switch(
    method,
    mean   = colMeans(mat_z, na.rm = TRUE),
    median = apply(mat_z, 2, median, na.rm = TRUE),
    pc1    = {
      pca <- prcomp(t(mat_z[!apply(mat_z, 1, anyNA), ]), scale. = FALSE)
      pca$x[, 1]
    }
  )

  scores
}

# ── Survival Data Preparation ─────────────────────────────────────────────────

#' Parse OS / PFS status from cBioPortal string coding to 0/1 integer
#'
#' Accepts both numeric ("0","1") and string ("0:LIVING","1:DECEASED") formats.
#'
#' @param x  Character vector of status values.
#' @return   Integer vector (0 = censored, 1 = event).
parse_cbio_status <- function(x) {
  as.integer(str_extract(x, "^[01]"))
}

#' Merge clinical data with signature scores and prepare survival data
#'
#' @param clinical    Tibble from read_cbio_clinical().
#' @param scores      Named numeric vector (sample IDs as names).
#' @param endpoint    "OS" or "PFS".
#' @return Tibble with columns: patient_id, time, status, score, score_group.
prepare_survival_data <- function(clinical, scores, endpoint = c("OS", "PFS")) {
  endpoint <- match.arg(endpoint)

  time_col   <- paste0(endpoint, "_MONTHS")
  status_col <- paste0(endpoint, "_STATUS")

  required <- c("PATIENT_ID", time_col, status_col)
  missing_cols <- setdiff(required, colnames(clinical))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Clinical file is missing required column(s): %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Build score tibble — match patients to samples
  # cBioPortal sample IDs are patient IDs with "-01" (or similar) suffix
  score_tbl <- tibble(
    sample_id = names(scores),
    score     = unname(scores)
  ) %>%
    mutate(patient_id = str_remove(sample_id, "-\\d+$"))

  clinical %>%
    select(
      patient_id = PATIENT_ID,
      time       = all_of(time_col),
      status     = all_of(status_col)
    ) %>%
    mutate(
      time   = as.numeric(time),
      status = parse_cbio_status(status)
    ) %>%
    filter(!is.na(time), !is.na(status), time >= 0) %>%
    inner_join(score_tbl, by = "patient_id") %>%
    mutate(
      # Patients at exactly the median are assigned to 'High'.
      # With tied median values this may slightly imbalance groups,
      # but preserves reproducibility (no random tie-breaking).
      score_group = if_else(score >= median(score, na.rm = TRUE), "High", "Low"),
      score_group = factor(score_group, levels = c("Low", "High"))
    )
}

# ── Cox PH Analysis ───────────────────────────────────────────────────────────

#' Run Cox proportional hazards model for a signature score
#'
#' Both a continuous score model and a dichotomised (High/Low) model are fit.
#' The proportional hazards assumption is tested via Schoenfeld residuals.
#'
#' @param surv_data   Tibble from prepare_survival_data().
#' @param covariates  Optional character vector of additional covariate column names.
#' @return A list with elements:
#'   \item{cox_continuous}{coxph fit for continuous score}
#'   \item{cox_binary}{coxph fit for High vs Low}
#'   \item{ph_test}{cox.zph result for cox_binary}
#'   \item{summary_tbl}{tidy summary tibble}
run_cox_analysis <- function(surv_data, covariates = NULL) {
  surv_obj <- with(surv_data, Surv(time, status))

  # Continuous score model
  cox_cont <- coxph(surv_obj ~ score, data = surv_data)

  # Binary High/Low model (± additional covariates)
  if (!is.null(covariates) && length(covariates) > 0) {
    valid_covs <- intersect(covariates, colnames(surv_data))
    if (length(valid_covs) < length(covariates)) {
      warning("Some covariates not found in data: ",
              paste(setdiff(covariates, colnames(surv_data)), collapse = ", "))
    }
    formula_bin <- as.formula(
      paste("surv_obj ~ score_group +", paste(valid_covs, collapse = " + "))
    )
  } else {
    formula_bin <- surv_obj ~ score_group
  }

  cox_bin <- coxph(formula_bin, data = surv_data)

  # Test PH assumption
  ph_test <- tryCatch(cox.zph(cox_bin), error = function(e) NULL)

  # Tidy summary
  s <- summary(cox_bin)
  coef_tbl <- as_tibble(s$conf.int, rownames = "term") %>%
    rename(hr = `exp(coef)`, ci_lower = `lower .95`, ci_upper = `upper .95`) %>%
    left_join(
      as_tibble(s$coefficients, rownames = "term") %>%
        select(term, p_value = `Pr(>|z|)`),
      by = "term"
    ) %>%
    select(term, hr, ci_lower, ci_upper, p_value)

  list(
    cox_continuous = cox_cont,
    cox_binary     = cox_bin,
    ph_test        = ph_test,
    summary_tbl    = coef_tbl
  )
}

# ── Plotting ──────────────────────────────────────────────────────────────────

#' Plot Kaplan–Meier curves for High vs Low signature groups
#'
#' @param surv_data   Tibble from prepare_survival_data().
#' @param title       Plot title string.
#' @param endpoint    "OS" or "PFS" (used for y-axis label).
#' @param palette     Two-colour character vector (Low, High).
#' @return A ggsurvplot object.
plot_km_curves <- function(surv_data,
                           title    = NULL,
                           endpoint = c("OS", "PFS"),
                           palette  = c("#2E86AB", "#E84855")) {
  endpoint <- match.arg(endpoint)
  ylab     <- if (endpoint == "OS") "Overall Survival Probability" else
    "Progression-Free Survival Probability"

  fit <- survfit(Surv(time, status) ~ score_group, data = surv_data)

  ggsurvplot(
    fit,
    data          = surv_data,
    palette       = palette,
    title         = title,
    xlab          = "Time (months)",
    ylab          = ylab,
    legend.title  = "Signature",
    legend.labs   = c("Low", "High"),
    pval          = TRUE,
    pval.method   = TRUE,
    conf.int      = TRUE,
    risk.table    = TRUE,
    risk.table.height = 0.25,
    ggtheme       = theme_bw(base_size = 12),
    surv.median.line = "hv"
  )
}

#' Build a forest-plot tibble from a list of Cox results
#'
#' @param results_list Named list of results from run_cox_analysis().
#'                     Names are used as row labels.
#' @return A ggplot object (forest plot).
plot_forest <- function(results_list) {
  tbl <- imap_dfr(results_list, function(res, label) {
    row <- res$summary_tbl %>%
      filter(str_detect(term, "score_groupHigh"))
    if (nrow(row) == 0) return(NULL)
    tibble(
      label    = label,
      hr       = row$hr,
      ci_lower = row$ci_lower,
      ci_upper = row$ci_upper,
      p_value  = row$p_value
    )
  })

  if (nrow(tbl) == 0) {
    return(ggplot() + labs(title = "No results to display"))
  }

  tbl <- tbl %>%
    mutate(
      label     = factor(label, levels = rev(unique(label))),
      p_label   = ifelse(p_value < 0.001, "<0.001",
                    ifelse(p_value < 0.01, sprintf("%.3f", p_value),
                           sprintf("%.2f", p_value))),
      sig       = case_when(
        p_value < 0.05  ~ "p < 0.05",
        TRUE            ~ "p ≥ 0.05"
      )
    )

  ggplot(tbl, aes(x = hr, y = label, colour = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25,
                   linewidth = 0.8) +
    geom_point(size = 3, shape = 15) +
    geom_text(aes(label = sprintf("HR=%.2f [%.2f–%.2f]\np=%s",
                                  hr, ci_lower, ci_upper, p_label)),
              hjust = -0.1, size = 3) +
    scale_colour_manual(values = c("p < 0.05" = "#E84855", "p ≥ 0.05" = "#2E86AB"),
                        name = NULL) +
    scale_x_log10(limits = c(
      min(tbl$ci_lower, na.rm = TRUE) * 0.7,
      max(tbl$ci_upper, na.rm = TRUE) * 2.5
    )) +
    labs(x = "Hazard Ratio (95% CI, log scale)", y = NULL,
         title = "Forest Plot — High vs Low Signature Score") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
}

# ── Reporting Helpers ─────────────────────────────────────────────────────────

#' Format a p-value for display in a table
fmt_pval <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "<0.001",
    p < 0.01  ~ sprintf("%.3f", p),
    TRUE      ~ sprintf("%.2f",  p)
  )
}

#' Build a concise results summary table across all cohort × signature × endpoint
#'
#' @param all_results  Named list (cohort → endpoint → signature → cox result list).
#'                     Each leaf must contain: summary_tbl, sig_label, endpoint,
#'                     and optionally sig_direction ("positive" or "negative").
#' @return A tibble suitable for kableExtra rendering.
#'
#' The `Direction_OK` column flags whether the observed HR direction matches the
#' expected direction defined in signatures.yaml:
#'   - "positive" signature → HR > 1 expected (High = worse)
#'   - "negative" signature → HR < 1 expected (High = better)
#'   - NA when no direction is specified.
build_summary_table <- function(all_results) {
  rows <- list()
  for (cohort in names(all_results)) {
    for (endpoint in names(all_results[[cohort]])) {
      for (sig in names(all_results[[cohort]][[endpoint]])) {
        res <- all_results[[cohort]][[endpoint]][[sig]]
        row <- res$summary_tbl %>%
          filter(str_detect(term, "score_groupHigh"))
        if (nrow(row) == 0) next

        # Evaluate whether HR direction matches the declared signature direction
        direction    <- res$sig_direction %||% NA_character_
        direction_ok <- dplyr::case_when(
          is.na(direction)             ~ NA,
          direction == "positive"      ~ row$hr >= 1,
          direction == "negative"      ~ row$hr <  1,
          TRUE                         ~ NA
        )

        rows[[length(rows) + 1]] <- tibble(
          Cohort        = cohort,
          Endpoint      = endpoint,
          Signature     = sig,
          HR            = round(row$hr, 2),
          CI_lower      = round(row$ci_lower, 2),
          CI_upper      = round(row$ci_upper, 2),
          P_value       = fmt_pval(row$p_value),
          Direction_OK  = direction_ok
        )
      }
    }
  }
  bind_rows(rows)
}
