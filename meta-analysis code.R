# ---- Setup ----
library(readxl)
library(metafor)
library(ggplot2)

# ---- Helper: publication-bias checks on one rma() fit ----
publication_bias_checks <- function(fit,
                                    label         = NULL,
                                    do_egger_kmin = 10,
                                    draw_funnel   = TRUE,
                                    funnel_levels = c(90, 95, 99),
                                    funnel_shade  = c("gray90", "gray80", "gray70")) {
  if (inherits(fit, "try-error") || is.null(fit)) return(invisible(NULL))
  k <- fit$k
  cat("\n====================\n",
      ifelse(is.null(label), "Publication-bias checks", label),
      "\n====================\n", sep = "")
  cat("Number of studies (k):", k, "\n")
  
  # 1) Egger’s regression test (SE as predictor)
  if (k >= do_egger_kmin) {
    cat("\nEgger's test (predictor = sei):\n")
    print(regtest(fit, model = "rma", predictor = "sei"))
  } else {
    cat("k <", do_egger_kmin, ": Skipping Egger's test (low power / inflated FP).\n")
  }
  
  # 2) Begg & Mazumdar rank correlation test
  if (k >= do_egger_kmin) {
    cat("\nBegg & Mazumdar rank correlation test:\n")
    print(ranktest(fit))
  } else {
    cat("k <", do_egger_kmin, ": Skipping Begg's test.\n")
  }
  
  # 3) Contour-enhanced funnel plot
  if (draw_funnel) {
    cat("\nDrawing contour-enhanced funnel plot...\n")
    funnel(fit,
           yaxis  = "sei",
           level  = funnel_levels,
           shade  = funnel_shade,
           transf = transf.ilogit,
           xlab   = "Effect (back-transformed proportion)",
           ylab   = "Standard error")
    legend("topright", inset = 0.02,
           legend = paste0("p < ", c(0.10, 0.05, 0.01)),
           fill   = funnel_shade, bty = "n")
  }
  
  # 4) Trim-and-fill pooled estimate
  tf <- try(trimfill(fit), silent = TRUE)
  if (!inherits(tf, "try-error")) {
    cat("\nTrim-and-fill pooled estimate (back-transformed):\n")
    print(predict(tf, transf = transf.ilogit))
  } else {
    cat("\nTrim-and-fill failed (algorithm did not converge or insufficient data).\n")
  }
  
  invisible(list(k = k,
                 egger_ok = k >= do_egger_kmin,
                 begg_ok  = k >= do_egger_kmin,
                 tf       = if (!inherits(tf, "try-error")) tf else NULL))
}

# ---- Helper: run one metric block inside a file/sheet ----
run_metric_block <- function(path, sheet, metric_col, color,
                             exclude_ids  = NULL,
                             unit_equals  = NULL,
                             show_summary = TRUE,
                             title_prefix = NULL,
                             show_pub_bias = TRUE) {
  
  # 1) Load once
  df <- read_excel(path, sheet = sheet)
  
  # 2) Optional filters
  if (!is.null(exclude_ids) && "ID" %in% names(df)) {
    df <- df[!(df$ID %in% exclude_ids), ]
  }
  if (!is.null(unit_equals) && "unit" %in% names(df)) {
    df <- df[df$unit == unit_equals, ]
  }
  
  # 3) Coerce types safely
  if (!("n" %in% names(df))) stop("Column 'n' not found in sheet: ", sheet)
  df$n <- suppressWarnings(as.numeric(df$n))
  
  if (!(metric_col %in% names(df))) {
    stop(sprintf("Column '%s' not found in sheet: %s", metric_col, sheet))
  }
  df[[metric_col]] <- suppressWarnings(as.numeric(df[[metric_col]]))
  
  # 4) Keep complete rows with 0<=metric<=1 and n>0
  keep <- is.finite(df$n) & df$n > 0 &
    is.finite(df[[metric_col]]) & df[[metric_col]] >= 0 & df[[metric_col]] <= 1
  d <- df[keep, , drop = FALSE]
  
  if (nrow(d) < 2) {
    warning(sprintf("Not enough rows after cleaning for %s | %s", sheet, metric_col))
    return(invisible(NULL))
  }
  
  # 5) Construct pseudo-events and fit model
  d$TP <- round(d[[metric_col]] * d$n)
  
  fit <- rma(measure = "PLO",
             xi = d$TP,
             ni = d$n,
             data = d,
             method = "REML")
  
  if (show_summary) print(summary(fit))
  
  # 6) Labels
  slab <- if ("Citation" %in% names(d)) d$Citation else paste0("Study ", seq_len(nrow(d)))
  xlab_txt <- sprintf("Pooled %s (proportion)", metric_col)
  if (!is.null(title_prefix)) xlab_txt <- sprintf("%s %s", title_prefix, xlab_txt)
  
  # 7) Forest
  forest(fit,
         slab   = slab,
         transf = transf.ilogit,
         xlab   = xlab_txt,
         colout = color,
         col    = color,
         refline = NA)
  
  # 8) Publication bias checks
  if (show_pub_bias) {
    publication_bias_checks(
      fit,
      label = sprintf("%s | %s | %s", basename(path), sheet, metric_col)
    )
  }
  
  invisible(fit)
}

# ---- Convenience wrapper for a set of metrics on one sheet ----
run_all_metrics <- function(path, sheet, color,
                            metrics      = c("Sensitivity", "F-1 score", "Precision",
                                             "Accuracy", "Specificity", "AUC"),
                            exclude_ids  = NULL,
                            unit_equals  = NULL,
                            title_prefix = NULL,
                            show_pub_bias = TRUE) {
  fits <- list()
  for (m in metrics) {
    message(sprintf("==> %s | %s | %s", basename(path), sheet, m))
    fits[[m]] <- try(
      run_metric_block(path, sheet, m, color,
                       exclude_ids  = exclude_ids,
                       unit_equals  = unit_equals,
                       title_prefix = title_prefix,
                       show_pub_bias = show_pub_bias),
      silent = TRUE
    )
  }
  invisible(fits)
}

## adjust path based on need
# 1) Screen patients
run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "screen patients",
  color = "blue3",
  metrics = c("Sensitivity", "F-1 score"),
  exclude_ids = NULL,
  unit_equals = NULL,
  title_prefix = NULL
)

# Then remaining metrics
run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "screen patients",
  color = "blue3",
  metrics = c("Precision", "Accuracy", "Specificity", "AUC"),
  exclude_ids = 162,
  unit_equals = NULL,
  title_prefix = NULL
)

# 2) Identify patients
run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "identify patient",
  color = "darkgreen",
  metrics = c("Sensitivity", "F-1 score", "Precision", "Accuracy", "Specificity", "AUC"),
  exclude_ids = NULL,
  unit_equals = "patients",
  title_prefix = NULL
)

# 3) Identify eligibility criteria from
run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "identify eligibility criteria",
  color = "darkorchid",
  metrics = c("Sensitivity", "F-1 score", "Precision", "Accuracy", "Specificity", "AUC"),
  exclude_ids = NULL,
  unit_equals = NULL,
  title_prefix = NULL
)

# 4) Classify eligibility criteria
run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "classify eligibility criteria",
  color = "darkred",
  metrics = c("Sensitivity", "F-1 score", "Precision", "Accuracy", "Specificity", "AUC"),
  exclude_ids = NULL,
  unit_equals = NULL,
  title_prefix = NULL
)

