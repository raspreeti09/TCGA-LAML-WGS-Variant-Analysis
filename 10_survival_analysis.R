# =============================================================================
# 10_survival_analysis.R — Overall survival analysis
#
# Steps:
#   A. maftools::mafSurvival — KM curves for key driver genes
#   B. Custom KM curves with log-rank test (survminer)
#   C. TMB-stratified survival (high vs low)
#   D. Co-mutation survival (e.g. FLT3+/NPM1+)
#   E. Multivariate Cox proportional hazards model
#   F. Forest plot of Cox HR estimates
# =============================================================================

source("config.R")

library(maftools)
library(survival)
library(survminer)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

set.seed(CONFIG$seed)

laml_maf <- readRDS(
  file.path(CONFIG$data_proc, "maf", "laml_maf_filtered.rds")
)
maf_dt  <- as.data.table(laml_maf@data)
clin_dt <- as.data.table(getClinicalData(laml_maf))
colnames(clin_dt) <- tolower(colnames(clin_dt))

# ── Prepare survival data ─────────────────────────────────────────────────────
# Identify OS columns (may vary between TCGA releases)
os_time_col   <- grep("days_to_death|os_time|overall_surv",
                       colnames(clin_dt), value = TRUE, ignore.case = TRUE)[1]
os_status_col <- grep("vital_status|os_event|dead",
                       colnames(clin_dt), value = TRUE, ignore.case = TRUE)[1]
last_fup_col  <- grep("days_to_last_follow|last_follow",
                       colnames(clin_dt), value = TRUE, ignore.case = TRUE)[1]

if (is.na(os_time_col)) stop("Cannot find OS time column in clinical data.")
message("OS time column:   ", os_time_col)
message("OS status column: ", os_status_col)

# Merge survival + mutation info
clin_dt[[os_time_col]]   <- suppressWarnings(as.numeric(clin_dt[[os_time_col]]))
clin_dt[[os_status_col]] <- toupper(trimws(clin_dt[[os_status_col]]))

# Compute unified OS time
if (!is.na(last_fup_col)) {
  clin_dt[[last_fup_col]] <- suppressWarnings(as.numeric(clin_dt[[last_fup_col]]))
  clin_dt[, OS_time := ifelse(
    !is.na(get(os_time_col)),
    get(os_time_col),
    get(last_fup_col)
  )]
} else {
  clin_dt[, OS_time := get(os_time_col)]
}

clin_dt[, OS_event := as.integer(get(os_status_col) == "DEAD")]

# Convert days → years
clin_dt[, OS_years := OS_time / 365.25]

surv_data <- clin_dt[!is.na(OS_years) & OS_years >= 0,
                      .(barcode    = tumor_sample_barcode,
                        patient_id = substr(tumor_sample_barcode, 1, 12),
                        OS_years, OS_event)]

message("Patients with survival data: ", nrow(surv_data))

# ══════════════════════════════════════════════════════════════════════════════
# A. maftools::mafSurvival
# ══════════════════════════════════════════════════════════════════════════════
message("=== A. maftools mafSurvival ===")

driver_genes <- c("FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2",
                  "TET2", "RUNX1", "TP53", "CEBPA", "NRAS",
                  "ASXL1", "WT1", "SF3B1", "STAG2")
driver_genes <- intersect(driver_genes, unique(maf_dt$Hugo_Symbol))

for (gene in driver_genes) {
  tryCatch({
    save_base_fig(
      expr = mafSurvival(
        maf         = laml_maf,
        genes       = gene,
        time        = "OS_years",
        Status      = "OS_event",
        isTCGA      = TRUE,
        clinicalData = surv_data
      ),
      name   = paste0("mafSurv_", gene),
      subdir = "survival",
      w = 8, h = 6
    )
  }, error = function(e) message("mafSurvival failed for ", gene, ": ", e$message))
}

# ══════════════════════════════════════════════════════════════════════════════
# B. Custom KM curves with survminer
# ══════════════════════════════════════════════════════════════════════════════
message("=== B. Custom KM curves ===")

# Add mutation status for each gene
km_data <- copy(surv_data)

for (gene in driver_genes) {
  mut_pts <- unique(maf_dt[Hugo_Symbol == gene, Tumor_Sample_Barcode])
  set(km_data, j = gene,
      value = ifelse(km_data$barcode %in% mut_pts, "Mutated", "Wild-type"))
}

# Panel of KM plots
km_plots <- lapply(driver_genes, function(gene) {
  if (!gene %in% colnames(km_data)) return(NULL)
  df <- km_data[!is.na(OS_years)]
  df[[gene]] <- factor(df[[gene]], levels = c("Wild-type", "Mutated"))

  fit <- survfit(as.formula(paste0("Surv(OS_years, OS_event) ~ ", gene)),
                 data = df)

  n_mut   <- sum(df[[gene]] == "Mutated",   na.rm = TRUE)
  n_wt    <- sum(df[[gene]] == "Wild-type", na.rm = TRUE)

  ggsurvplot(
    fit,
    data           = df,
    palette        = c("#457B9D", "#E63946"),
    pval           = TRUE,
    pval.method    = TRUE,
    conf.int       = TRUE,
    risk.table     = TRUE,
    risk.table.y.text.col = TRUE,
    ncensor.plot   = FALSE,
    xlab           = "Time (years)",
    ylab           = "Overall Survival Probability",
    title          = paste0(gene,
                            " — Mutated (n=", n_mut,
                            ") vs WT (n=", n_wt, ")"),
    legend.labs    = c("Wild-type", "Mutated"),
    ggtheme        = theme_bw(base_size = 10)
  )
})
names(km_plots) <- driver_genes

# Save each plot
for (gene in driver_genes) {
  if (is.null(km_plots[[gene]])) next
  save_base_fig(
    expr   = print(km_plots[[gene]]),
    name   = paste0("KM_", gene),
    subdir = "survival",
    w = 9, h = 8
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# C. TMB-stratified survival
# ══════════════════════════════════════════════════════════════════════════════
message("=== C. TMB survival ===")

smry <- as.data.table(getSampleSummary(laml_maf))
smry[, TMB_perMb := total / 2800]

tmb_surv <- merge(km_data, smry[, .(Tumor_Sample_Barcode, TMB_perMb)],
                  by.x = "barcode", by.y = "Tumor_Sample_Barcode")

tmb_median <- median(tmb_surv$TMB_perMb, na.rm = TRUE)
tmb_surv[, TMB_group := ifelse(TMB_perMb > tmb_median, "TMB-High", "TMB-Low")]

fit_tmb <- survfit(Surv(OS_years, OS_event) ~ TMB_group, data = tmb_surv)

p_tmb_surv <- ggsurvplot(
  fit_tmb,
  data         = tmb_surv,
  palette      = c("#D73027", "#4393C3"),
  pval         = TRUE,
  pval.method  = TRUE,
  conf.int     = TRUE,
  risk.table   = TRUE,
  xlab         = "Time (years)",
  ylab         = "Overall Survival",
  title        = paste0("TMB-High vs TMB-Low (median = ",
                        round(tmb_median, 3), " mut/Mb)"),
  legend.labs  = c("TMB-High", "TMB-Low"),
  ggtheme      = theme_bw(base_size = 10)
)

save_base_fig(
  expr   = print(p_tmb_surv),
  name   = "TMB_high_vs_low_survival",
  subdir = "survival",
  w = 9, h = 8
)

# ══════════════════════════════════════════════════════════════════════════════
# D. Co-mutation survival (FLT3-ITD + NPM1)
# ══════════════════════════════════════════════════════════════════════════════
message("=== D. Co-mutation KM curves ===")

if (all(c("FLT3", "NPM1") %in% colnames(km_data))) {
  km_data[, FLT3_NPM1_group := paste0(
    "FLT3:", FLT3, " / NPM1:", NPM1
  )]

  group_counts <- km_data[, .N, by = FLT3_NPM1_group][order(-N)]
  message("FLT3/NPM1 co-mutation groups:\n",
          paste(capture.output(print(group_counts)), collapse = "\n"))

  fit_co <- survfit(Surv(OS_years, OS_event) ~ FLT3_NPM1_group,
                    data = km_data)

  p_co <- ggsurvplot(
    fit_co,
    data        = km_data,
    pval        = TRUE,
    conf.int    = FALSE,
    risk.table  = TRUE,
    xlab        = "Time (years)",
    ylab        = "Overall Survival",
    title       = "FLT3 / NPM1 Co-mutation Survival",
    ggtheme     = theme_bw(base_size = 10)
  )

  save_base_fig(
    expr   = print(p_co),
    name   = "FLT3_NPM1_comutation_survival",
    subdir = "survival",
    w = 10, h = 9
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# E. Multivariate Cox proportional hazards model
# ══════════════════════════════════════════════════════════════════════════════
message("=== E. Multivariate Cox model ===")

# Genes mutated in ≥ 10% of patients for stable Cox estimates
gene_smry    <- as.data.table(getGeneSummary(laml_maf))
n_patients   <- nrow(getSampleSummary(laml_maf))
cox_genes    <- gene_smry[MutatedSamples / n_patients >= 0.10, Hugo_Symbol]
cox_genes    <- intersect(cox_genes, colnames(km_data))

# Also include clinical covariates if available
age_col <- grep("^age", colnames(clin_dt), value = TRUE)[1]

cox_df <- copy(km_data)
if (!is.na(age_col)) {
  age_dt   <- clin_dt[, .(barcode = tumor_sample_barcode,
                            age    = get(age_col))]
  age_dt$age <- suppressWarnings(as.numeric(age_dt$age))
  cox_df   <- merge(cox_df, age_dt, by = "barcode", all.x = TRUE)
  cox_genes_final <- c(cox_genes, "age")
} else {
  cox_genes_final <- cox_genes
}

# Convert mutation status to 0/1
for (g in cox_genes) {
  cox_df[[g]] <- as.integer(cox_df[[g]] == "Mutated")
}

cox_formula <- as.formula(
  paste0("Surv(OS_years, OS_event) ~ ",
         paste(cox_genes_final, collapse = " + "))
)

cox_model <- tryCatch(
  coxph(cox_formula, data = cox_df),
  error = function(e) { message("Cox model error: ", e$message); NULL }
)

if (!is.null(cox_model)) {
  cox_summary <- as.data.table(
    broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
  )
  setnames(cox_summary,
           c("estimate", "conf.low", "conf.high", "p.value"),
           c("HR", "CI_low", "CI_high", "pval"))
  cox_summary[, padj := p.adjust(pval, method = "BH")]
  cox_summary <- cox_summary[order(pval)]

  fwrite(cox_summary,
         file.path(CONFIG$tables, "cox_multivariate_results.tsv"),
         sep = "\t")
  message("Cox model summary saved.")

  # ── F. Forest plot ────────────────────────────────────────────────────────
  message("=== F. Forest plot ===")

  cox_plot_dt <- cox_summary[!is.na(HR) & !is.infinite(HR)]
  cox_plot_dt[, sig := ifelse(padj < 0.05, "FDR < 0.05", "NS")]
  cox_plot_dt[, label := paste0(round(HR, 2),
                                 " (", round(CI_low, 2),
                                 "–",  round(CI_high, 2), ")")]

  p_forest <- ggplot(cox_plot_dt,
                     aes(x = HR, y = reorder(term, -HR),
                         xmin = CI_low, xmax = CI_high,
                         colour = sig)) +
    geom_errorbarh(height = 0.2, linewidth = 0.8) +
    geom_point(size = 3) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_text(aes(x = max(CI_high, na.rm = TRUE) * 1.05,
                  label = label), hjust = 0, size = 2.8) +
    scale_colour_manual(values = c("FDR < 0.05" = "firebrick",
                                   "NS"          = "grey50")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(title = "Multivariate Cox Model — TCGA-LAML WGS",
         x = "Hazard Ratio (95% CI)",
         y = NULL,
         colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  save_fig(p_forest, "cox_forest_plot", "survival", w = 12, h = max(6, nrow(cox_plot_dt) * 0.4 + 2))
}

message("\n=== Survival analysis complete ===")
