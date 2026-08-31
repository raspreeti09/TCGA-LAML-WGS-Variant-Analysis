# =============================================================================
# 09_clinical_correlation.R — Variant-clinical feature correlations
#
# Analyses:
#   A. Mutation frequency stratified by FAB classification
#   B. Age at diagnosis vs mutation burden (TMB)
#   C. Sex differences in mutation landscape
#   D. Cytogenetic risk group vs key mutations
#   E. Oncogenic mutations vs clinical variables (Fisher's / Wilcoxon)
#   F. Clinical oncoprint with annotation tracks
# =============================================================================

source("c:/Users/asus/OneDrive/Desktop/Project-LAML/config/config.R")

library(maftools)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

set.seed(CONFIG$seed)

laml_maf <- readRDS(
  file.path(CONFIG$data_proc, "maf", "laml_maf_filtered.rds")
)
maf_dt   <- as.data.table(laml_maf@data)
clin_dt  <- as.data.table(getClinicalData(laml_maf))

# Standardise clinical columns
colnames(clin_dt) <- tolower(colnames(clin_dt))

# Add patient_id
clin_dt[, patient_id := substr(tumor_sample_barcode, 1, 12)]

message("Clinical table: ", nrow(clin_dt), " patients, ", ncol(clin_dt), " vars")
message("Clinical columns: ", paste(colnames(clin_dt)[1:20], collapse = ", "))

# ── Helper: count of mutated patients per gene per group ──────────────────────
mut_freq_by_group <- function(gene, group_col) {
  # 1 = mutated, 0 = not
  gene_pts <- unique(maf_dt[Hugo_Symbol == gene, Tumor_Sample_Barcode])

  if (!group_col %in% colnames(clin_dt)) return(NULL)
  tmp <- clin_dt[!is.na(get(group_col)),
                  .(barcode = tumor_sample_barcode,
                    group   = get(group_col))]
  tmp[, mutated := as.integer(barcode %in% gene_pts)]
  tmp[, .(n_mutated = sum(mutated),
          n_total   = .N,
          pct       = mean(mutated) * 100),
      by = group]
}

# ══════════════════════════════════════════════════════════════════════════════
# A. FAB classification
# ══════════════════════════════════════════════════════════════════════════════
message("=== A. FAB vs mutation ===")

fab_col <- grep("fab", colnames(clin_dt), value = TRUE, ignore.case = TRUE)[1]

if (!is.na(fab_col)) {
  clin_dt[[fab_col]] <- toupper(trimws(clin_dt[[fab_col]]))
  fab_counts <- clin_dt[, .N, by = get(fab_col)][order(-N)]
  message("FAB groups: ", paste(fab_counts$get, collapse = ", "))

  key_genes <- c("NPM1", "FLT3", "DNMT3A", "IDH1", "IDH2",
                 "CEBPA", "RUNX1", "TP53", "TET2", "NRAS")

  fab_mut_list <- lapply(key_genes, function(g) {
    res <- mut_freq_by_group(g, fab_col)
    if (!is.null(res)) { res[, gene := g]; res }
  })
  fab_mut_dt <- rbindlist(fab_mut_list, fill = TRUE)
  setnames(fab_mut_dt, "group", "FAB")

  fwrite(fab_mut_dt,
         file.path(CONFIG$tables, "FAB_mutation_frequency.tsv"), sep = "\t")

  p_fab <- ggplot(fab_mut_dt[!is.na(FAB)],
                  aes(x = FAB, y = pct, fill = gene)) +
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = "Paired") +
    labs(title = "Mutation Frequency by FAB Classification",
         x = "FAB Class", y = "% Patients Mutated", fill = "Gene") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  save_fig(p_fab, "01_FAB_mutation_frequency", "landscape", w = 14, h = 7)
}

# ══════════════════════════════════════════════════════════════════════════════
# B. Age at diagnosis vs TMB
# ══════════════════════════════════════════════════════════════════════════════
message("=== B. Age vs TMB ===")

# Compute TMB
smry <- as.data.table(getSampleSummary(laml_maf))
smry[, TMB_perMb := total / 2800]

age_col <- grep("age_at_diag|age_diag|^age$",
                colnames(clin_dt), value = TRUE, ignore.case = TRUE)[1]

if (!is.na(age_col)) {
  clin_dt[[age_col]] <- suppressWarnings(as.numeric(clin_dt[[age_col]]))
  clin_dt[, age_yr := get(age_col) / 365.25]   # may be in days

  tmb_age <- merge(
    smry[, .(Tumor_Sample_Barcode, TMB_perMb)],
    clin_dt[, .(tumor_sample_barcode,
                age_yr,
                age_group = cut(age_yr,
                                breaks = c(0, 40, 60, 80, 200),
                                labels = c("<40", "40-60", "60-80", ">80"),
                                right  = FALSE))],
    by.x = "Tumor_Sample_Barcode",
    by.y = "tumor_sample_barcode"
  )

  p_age_tmb <- ggplot(tmb_age[!is.na(age_yr)],
                      aes(x = age_yr, y = log10(TMB_perMb + 0.001))) +
    geom_point(aes(colour = age_group), alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE,
                colour = "firebrick", linewidth = 1) +
    stat_cor(method = "pearson", label.x.npc = 0.6) +
    scale_colour_brewer(palette = "Set2") +
    labs(title = "Age at Diagnosis vs Tumour Mutation Burden",
         x = "Age (years)", y = "log10(TMB per Mb)",
         colour = "Age Group") +
    theme_bw(base_size = 12)

  save_fig(p_age_tmb, "02_age_vs_TMB", "landscape")

  fwrite(tmb_age,
         file.path(CONFIG$tables, "TMB_age_correlation.tsv"), sep = "\t")
}

# ══════════════════════════════════════════════════════════════════════════════
# C. Sex differences
# ══════════════════════════════════════════════════════════════════════════════
message("=== C. Sex differences ===")

sex_col <- grep("^sex$|gender", colnames(clin_dt), value = TRUE, ignore.case = TRUE)[1]

if (!is.na(sex_col)) {
  clin_dt[[sex_col]] <- toupper(trimws(clin_dt[[sex_col]]))

  sex_freq_list <- lapply(
    c("DNMT3A", "FLT3", "NPM1", "IDH1", "IDH2", "TP53"),
    function(g) {
      res <- mut_freq_by_group(g, sex_col)
      if (!is.null(res)) { res[, gene := g]; res }
    }
  )
  sex_freq_dt <- rbindlist(sex_freq_list, fill = TRUE)
  setnames(sex_freq_dt, "group", "sex")

  # Chi-squared test per gene
  sex_freq_dt[, p_chisq := mapply(function(g, s) {
    pts_g <- unique(maf_dt[Hugo_Symbol == g, Tumor_Sample_Barcode])
    tmp   <- clin_dt[!is.na(get(sex_col)),
                      .(sex = get(sex_col),
                        mutated = as.integer(
                          tumor_sample_barcode %in% pts_g))]
    tbl <- table(tmp$sex, tmp$mutated)
    if (prod(dim(tbl)) < 4 || any(tbl < 5)) return(1)
    chisq.test(tbl)$p.value
  }, gene, sex)]

  fwrite(sex_freq_dt,
         file.path(CONFIG$tables, "sex_mutation_frequency.tsv"), sep = "\t")
}

# ══════════════════════════════════════════════════════════════════════════════
# D. Cytogenetic risk group
# ══════════════════════════════════════════════════════════════════════════════
message("=== D. Cytogenetic risk groups ===")

cyto_col <- grep("cyto|risk|cytogenetics",
                 colnames(clin_dt), value = TRUE, ignore.case = TRUE)[1]

if (!is.na(cyto_col)) {
  clin_dt[[cyto_col]] <- trimws(clin_dt[[cyto_col]])

  cyto_freq_list <- lapply(
    c("FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2", "CEBPA",
      "TP53", "TET2", "RUNX1", "NRAS"),
    function(g) {
      res <- mut_freq_by_group(g, cyto_col)
      if (!is.null(res)) { res[, gene := g]; res }
    }
  )
  cyto_freq_dt <- rbindlist(cyto_freq_list, fill = TRUE)
  setnames(cyto_freq_dt, "group", "cyto_risk")
  cyto_freq_dt <- cyto_freq_dt[!is.na(cyto_risk)]

  p_cyto <- ggplot(cyto_freq_dt,
                   aes(x = cyto_risk, y = pct, fill = gene)) +
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "Mutation Frequency by Cytogenetic Risk Group",
         x = "Cytogenetic Risk", y = "% Patients Mutated",
         fill = "Gene") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  save_fig(p_cyto, "03_cytogenetic_risk_mutation", "landscape", w = 14, h = 7)

  fwrite(cyto_freq_dt,
         file.path(CONFIG$tables, "cytogenetic_risk_mutation.tsv"),
         sep = "\t")
}

# ══════════════════════════════════════════════════════════════════════════════
# E. Comprehensive mutation × clinical variable association table
# ══════════════════════════════════════════════════════════════════════════════
message("=== E. Mutation-clinical association ===")

clin_vars <- intersect(
  c(fab_col, sex_col, cyto_col, age_col,
    grep("wbc|blast|plt|hgb|ldh|creat", colnames(clin_dt),
         value = TRUE, ignore.case = TRUE)),
  colnames(clin_dt)
)
clin_vars <- clin_vars[!is.na(clin_vars)]

driver_genes <- c("FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2", "TET2",
                  "RUNX1", "TP53", "CEBPA", "NRAS", "KRAS", "WT1",
                  "ASXL1", "SF3B1", "SRSF2", "STAG2")
driver_genes <- intersect(driver_genes, unique(maf_dt$Hugo_Symbol))

assoc_results <- vector("list", length(driver_genes) * length(clin_vars))
k <- 1L

for (gene in driver_genes) {
  pts_mutated <- unique(maf_dt[Hugo_Symbol == gene, Tumor_Sample_Barcode])
  clin_dt[, mut_status := as.integer(tumor_sample_barcode %in% pts_mutated)]

  for (cvar in clin_vars) {
    x <- clin_dt[[cvar]]

    if (is.numeric(x)) {
      # Wilcoxon rank sum
      wt  <- tryCatch(
        wilcox.test(x[clin_dt$mut_status == 1],
                    x[clin_dt$mut_status == 0]),
        error = function(e) list(p.value = NA)
      )
      assoc_results[[k]] <- data.table(
        gene = gene, variable = cvar,
        test = "Wilcoxon",
        median_mut   = median(x[clin_dt$mut_status == 1], na.rm = TRUE),
        median_nomut = median(x[clin_dt$mut_status == 0], na.rm = TRUE),
        pval = wt$p.value
      )
    } else {
      # Fisher's exact
      tbl <- table(clin_dt[[cvar]], clin_dt$mut_status)
      if (prod(dim(tbl)) < 4) {
        assoc_results[[k]] <- NULL; k <- k + 1L; next
      }
      ft  <- tryCatch(fisher.test(tbl, simulate.p.value = TRUE),
                      error = function(e) list(p.value = NA))
      assoc_results[[k]] <- data.table(
        gene = gene, variable = cvar,
        test = "Fisher",
        median_mut = NA, median_nomut = NA,
        pval = ft$p.value
      )
    }
    k <- k + 1L
  }
}

assoc_dt <- rbindlist(assoc_results, fill = TRUE)
assoc_dt[, padj := p.adjust(pval, method = "BH")]
assoc_dt <- assoc_dt[order(padj)]
fwrite(assoc_dt,
       file.path(CONFIG$tables, "mutation_clinical_association.tsv"),
       sep = "\t")

message("Significant mutation-clinical associations (FDR < 0.1): ",
        sum(assoc_dt$padj < 0.1, na.rm = TRUE))

message("\n=== Clinical correlation complete ===")
