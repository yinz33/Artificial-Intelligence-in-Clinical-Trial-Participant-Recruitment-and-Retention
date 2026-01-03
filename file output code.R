# ---- Setup ----
setwd("/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention")
library(readxl)
library(metafor)
library(ggplot2)

# small util for NULL-coalescing
`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---- Helper: run one metric block inside a file/sheet and DRAW the forest ----
run_metric_block <- function(path, sheet, metric_col, color,
                             exclude_ids = NULL,
                             unit_equals = NULL,
                             show_summary = TRUE,
                             title_prefix = NULL) {
  
  df <- read_excel(path, sheet = sheet)
  
  if (!is.null(exclude_ids) && "ID" %in% names(df)) {
    df <- df[!(df$ID %in% exclude_ids), ]
  }
  if (!is.null(unit_equals) && "unit" %in% names(df)) {
    df <- df[df$unit == unit_equals, ]
  }
  
  if (!("n" %in% names(df))) stop("Column 'n' not found in sheet: ", sheet)
  df$n <- suppressWarnings(as.numeric(df$n))
  
  if (!(metric_col %in% names(df))) {
    stop(sprintf("Column '%s' not found in sheet: %s", metric_col, sheet))
  }
  df[[metric_col]] <- suppressWarnings(as.numeric(df[[metric_col]]))
  
  keep <- is.finite(df$n) & df$n > 0 &
    is.finite(df[[metric_col]]) & df[[metric_col]] >= 0 & df[[metric_col]] <= 1
  d <- df[keep, , drop = FALSE]
  if (nrow(d) < 2) {
    warning(sprintf("Not enough rows after cleaning for %s | %s", sheet, metric_col))
    return(NULL)
  }
  
  d$TP <- round(d[[metric_col]] * d$n)
  
  fit <- rma(measure = "PLO",
             xi = d$TP,
             ni = d$n,
             data = d,
             method = "REML")
  
  if (show_summary) print(summary(fit))
  
  slab <- if ("Citation" %in% names(d)) d$Citation else paste0("Study ", seq_len(nrow(d)))
  xlab_txt <- sprintf("Pooled %s (proportion)", metric_col)
  if (!is.null(title_prefix)) xlab_txt <- sprintf("%s %s", title_prefix, xlab_txt)
  
  forest(fit,
         slab   = slab,
         transf = transf.ilogit,
         xlab   = xlab_txt,
         colout = color,
         col    = color,
         refline = NA)
  
  list(fit = fit, n_rows = nrow(d))
}

# ---- Run a set of metrics on one sheet and return fits ----
run_all_metrics <- function(path, sheet, color,
                            metrics = c("Sensitivity", "F-1 score", "Precision", "Accuracy", "Specificity", "AUC"),
                            exclude_ids = NULL,
                            unit_equals = NULL,
                            task_name = NULL,
                            title_prefix = NULL) {
  out <- list()
  for (m in metrics) {
    message(sprintf("==> %s | %s | %s", basename(path), sheet, m))
    res <- try(
      run_metric_block(path, sheet, m, color,
                       exclude_ids = exclude_ids,
                       unit_equals = unit_equals,
                       title_prefix = title_prefix),
      silent = TRUE
    )
    out[[m]] <- list(
      ok         = !inherits(res, "try-error") && !is.null(res),
      metric     = m,
      path       = path,
      sheet      = sheet,
      color      = color,
      task_name  = task_name %||% sheet,
      fit        = if (!inherits(res, "try-error") && !is.null(res)) res$fit else NULL
    )
  }
  out
}

# =========================
#     RUN & SAVE FORESTS
# =========================
pdf("all_forests.pdf", width = 8.5, height = 11)

fits_list <- list()

# 1) Screen patients
fits_list[["screen_patients_part1"]] <- run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "screen patients",
  color = "blue3",
  metrics = c("Sensitivity", "F-1 score"),
  task_name = "Screen patients"
)

fits_list[["screen_patients_part2"]] <- run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "screen patients",
  color = "blue3",
  metrics = c("Precision", "Accuracy", "Specificity", "AUC"),
  exclude_ids = 162,
  task_name = "Screen patients"
)

# 2) Identify patients (unit == "patients")
fits_list[["identify_patients"]] <- run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "identify patient",
  color = "darkgreen",
  metrics = c("Sensitivity", "F-1 score", "Precision", "Accuracy", "Specificity", "AUC"),
  unit_equals = "patients",
  task_name = "Identify patients"
)

# 3) Identify eligibility criteria
fits_list[["identify_ec"]] <- run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "identify eligibility criteria",
  color = "darkorchid",
  metrics = c("Sensitivity", "F-1 score", "Precision", "Accuracy", "Specificity", "AUC"),
  task_name = "Identify EC"
)

# 4) Classify eligibility criteria
fits_list[["classify_ec"]] <- run_all_metrics(
  path = "/Users/ayin/Desktop/Artificial Intelligence in Clinical Trial Participant Recruitment and Retention/workflow_category (6).xlsx",
  sheet = "classify eligibility criteria",
  color = "darkred",
  metrics = c("Sensitivity", "F-1 score", "Precision", "Accuracy", "Specificity", "AUC"),
  task_name = "Classify EC"
)

dev.off()  # writes "all_forests.pdf"

# =========================
#     SUMMARY TO CSV
# =========================
extract_row <- function(rec) {
  if (!isTRUE(rec$ok)) return(NULL)
  fit <- rec$fit
  pred <- predict(fit, transf = transf.ilogit)
  data.frame(
    Task         = rec$task_name,
    File         = basename(rec$path),
    Sheet        = rec$sheet,
    Metric       = rec$metric,
    K            = fit$k,
    Tau2         = fit$tau2,
    I2_percent   = fit$I2,
    Pooled       = as.numeric(pred$pred),
    CI_L         = as.numeric(pred$ci.lb),
    CI_U         = as.numeric(pred$ci.ub),
    stringsAsFactors = FALSE
  )
}

all_rows <- do.call(
  rbind,
  unlist(lapply(fits_list, function(group) lapply(group, extract_row)), recursive = FALSE)
)

if (!is.null(all_rows) && nrow(all_rows) > 0) {
  all_rows <- all_rows[order(all_rows$Task, all_rows$Metric), ]
  print(all_rows, row.names = FALSE, digits = 4)
  write.csv(all_rows, "meta_summary.csv", row.names = FALSE)
} else {
  message("No successful models to summarize.")
}

cat("\nSaved plots to: all_forests.pdf\n")
cat("Saved summary to: meta_summary.csv\n")

# ============================================================
#  SAVE ALL CONTOUR-ENHANCED FUNNELS TO PDF
# ============================================================
pdf("all_funnels.pdf", width = 8.5, height = 11)

draw_one_funnel <- function(rec) {
  if (!isTRUE(rec$ok)) return(invisible(NULL))
  fit <- rec$fit
  
  funnel(fit,
         yaxis  = "sei",
         level  = c(90, 95, 99),
         shade  = c("gray90", "gray80", "gray70"),
         transf = transf.ilogit,
         xlab   = "Effect (back-transformed proportion)",
         ylab   = "Standard Error")
  legend("topright", inset = 0.02,
         legend = c("p < 0.10", "p < 0.05", "p < 0.01"),
         fill   = c("gray90", "gray80", "gray70"), bty = "n")
  
  title(main = sprintf("%s | %s | %s",
                       rec$task_name, basename(rec$path), rec$metric))
  invisible(NULL)
}

invisible(lapply(fits_list, function(group) lapply(group, draw_one_funnel)))

dev.off()  # writes "all_funnels.pdf"
cat("Saved funnels to: all_funnels.pdf\n")

# ============================================================
#  PUBLICATION BIAS RESULTS TO CSV (robust)
#   - Adds Unadjusted pooled (Unadj_pred / CI)
#   - Egger (if k>=egger_kmin; default 8)
#   - Begg  (if k>=egger_kmin; default 8)
#   - Trim-and-Fill pooled (+ k0)
# ============================================================
pub_bias_row <- function(rec, egger_kmin = 8) {
  if (!isTRUE(rec$ok)) return(NULL)
  fit <- rec$fit
  k <- fit$k
  
  grab_num <- function(x, keys) {
    for (nm in keys) if (!is.null(x[[nm]])) return(as.numeric(x[[nm]]))
    if (!is.null(x$statistic)) return(as.numeric(unname(x$statistic)))
    if (!is.null(x$estimate))  return(as.numeric(unname(x$estimate)))
    if (!is.null(x$p.value))   return(as.numeric(x$p.value))
    NA_real_
  }
  
  # ---- Unadjusted pooled estimate (back-transformed) ----
  unadj_pred <- unadj_ci_lb <- unadj_ci_ub <- NA_real_
  unadj <- try(predict(fit, transf = transf.ilogit), silent = TRUE)
  if (!inherits(unadj, "try-error")) {
    unadj_pred  <- as.numeric(unadj$pred)
    unadj_ci_lb <- as.numeric(unadj$ci.lb)
    unadj_ci_ub <- as.numeric(unadj$ci.ub)
  }
  
  # ---- Egger (SE as predictor) ----
  egger_z <- egger_p <- egger_b <- egger_ci_lb <- egger_ci_ub <- NA_real_
  if (k >= egger_kmin) {
    eg <- try(metafor::regtest(fit, model = "rma", predictor = "sei"), silent = TRUE)
    if (!inherits(eg, "try-error")) {
      egger_z <- grab_num(eg, c("zval"))
      egger_p <- grab_num(eg, c("pval", "p.value"))
      egger_b <- grab_num(eg, c("b", "int", "estimate"))
      if (!is.null(eg$ci.lb)) egger_ci_lb <- as.numeric(eg$ci.lb)
      if (!is.null(eg$ci.ub)) egger_ci_ub <- as.numeric(eg$ci.ub)
      if (is.na(egger_ci_lb) && !is.null(eg$conf.int)) {
        egger_ci_lb <- as.numeric(eg$conf.int[1])
        egger_ci_ub <- as.numeric(eg$conf.int[2])
      }
    }
  }
  
  # ---- Begg & Mazumdar ----
  begg_tau <- begg_p <- NA_real_
  if (k >= egger_kmin) {
    bg <- try(metafor::ranktest(fit), silent = TRUE)
    if (!inherits(bg, "try-error")) {
      begg_tau <- grab_num(bg, c("tau", "estimate"))
      begg_p   <- grab_num(bg, c("pval", "p.value"))
    }
  }
  
  # ---- Trim-and-Fill ----
  tf_k0 <- NA_integer_; tf_pred <- tf_ci_lb <- tf_ci_ub <- NA_real_
  tf <- try(metafor::trimfill(fit), silent = TRUE)
  if (!inherits(tf, "try-error")) {
    tf_k0 <- as.integer(if (!is.null(tf$k0)) tf$k0 else 0L)
    tf_pred_obj <- try(predict(tf, transf = transf.ilogit), silent = TRUE)
    if (!inherits(tf_pred_obj, "try-error")) {
      tf_pred  <- as.numeric(tf_pred_obj$pred)
      tf_ci_lb <- as.numeric(tf_pred_obj$ci.lb)
      tf_ci_ub <- as.numeric(tf_pred_obj$ci.ub)
    }
  }
  
  data.frame(
    Task         = rec$task_name,
    File         = basename(rec$path),
    Sheet        = rec$sheet,
    Metric       = rec$metric,
    K            = k,
    
    # Unadjusted pooled estimate (original fit)
    Unadj_pred   = unadj_pred,
    Unadj_CI_L   = unadj_ci_lb,
    Unadj_CI_U   = unadj_ci_ub,
    
    # Egger
    Egger_kmin   = egger_kmin,
    Egger_z      = egger_z,
    Egger_p      = egger_p,
    Egger_b_int  = egger_b,
    Egger_CI_L   = egger_ci_lb,
    Egger_CI_U   = egger_ci_ub,
    
    # Begg
    Begg_tau     = begg_tau,
    Begg_p       = begg_p,
    
    # Trim & Fill (adjusted pooled)
    TF_k0_filled = tf_k0,
    TF_pred      = tf_pred,
    TF_CI_L      = tf_ci_lb,
    TF_CI_U      = tf_ci_ub,
    stringsAsFactors = FALSE
  )
}

pub_bias_rows <- do.call(
  rbind,
  unlist(lapply(fits_list, function(group) lapply(group, pub_bias_row)), recursive = FALSE)
)

if (!is.null(pub_bias_rows) && nrow(pub_bias_rows) > 0) {
  pub_bias_rows <- pub_bias_rows[order(pub_bias_rows$Task, pub_bias_rows$Metric), ]
  print(pub_bias_rows, row.names = FALSE, digits = 4)
  write.csv(pub_bias_rows, "publication_bias.csv", row.names = FALSE)
  cat("Saved publication-bias table to: publication_bias.csv\n")
} else {
  message("No publication-bias rows to save.")
}
