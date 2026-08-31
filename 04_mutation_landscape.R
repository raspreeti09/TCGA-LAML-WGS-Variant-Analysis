# =============================================================================
# 04_mutation_landscape.R — Somatic mutation landscape analysis
#
# Outputs:
#   • maftools summary statistics (TMB, variant types, classifications)
#   • Waterfall / oncoprint (top-N genes)
#   • Variant classification bar chart
#   • Transition / Transversion (TiTv) plot
#   • VAF raincloud plot per variant class
#   • Mutation burden distribution across samples
#   • Top-mutated genes barplot
# =============================================================================

source("c:/Users/asus/OneDrive/Desktop/Project-LAML/config/config.R")

library(maftools)
library(ggplot2)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(patchwork)
library(scales)

set.seed(CONFIG$seed)

# ── Load filtered MAF object ──────────────────────────────────────────────────
laml_maf <- readRDS(
  file.path(CONFIG$data_proc, "maf", "laml_maf_filtered.rds")
)
n_patients <- nrow(getSampleSummary(laml_maf))
message("Loaded MAF: ", n_patients,
        " patients, ", nrow(getGeneSummary(laml_maf)), " genes")

# ══════════════════════════════════════════════════════════════════════════════
# 1. MAF summary overview
# ══════════════════════════════════════════════════════════════════════════════
save_base_fig(
  expr    = plotmafSummary(maf = laml_maf, rmOutlier = TRUE,
                           addStat = "median", dashboard = TRUE,
                           titvRaw = FALSE),
  name    = "01_MAF_summary",
  subdir  = "landscape",
  w = 16, h = 10
)

# ══════════════════════════════════════════════════════════════════════════════
# 2. Oncoprint / Waterfall (top 30 genes)
# ══════════════════════════════════════════════════════════════════════════════
save_base_fig(
  expr = oncoplot(
    maf            = laml_maf,
    top             = 30,
    showTumorSampleBarcodes = FALSE,
    sortByAnnotation = TRUE,
    drawColBar      = TRUE,
    removeNonMutated = TRUE,
    fontSize        = 0.6,
    titleFontSize   = 1
  ),
  name   = "02_oncoprint_top30",
  subdir = "landscape",
  w = 18, h = 12
)

# ── Oncoprint with clinical annotation tracks ─────────────────────────────────
clin_dt <- getClinicalData(laml_maf)
annot_cols <- intersect(
  c("age_at_diagnosis", "gender", "race", "vital_status",
    "fab_classification", "cytogenetic_risk"),
  colnames(clin_dt)
)

if (length(annot_cols) > 0) {
  save_base_fig(
    expr = oncoplot(
      maf              = laml_maf,
      top              = 25,
      clinicalFeatures = annot_cols,
      sortByAnnotation = TRUE,
      removeNonMutated = TRUE,
      fontSize         = 0.55
    ),
    name   = "02b_oncoprint_clinical_annotation",
    subdir = "landscape",
    w = 20, h = 14
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Transition / Transversion (TiTv) ratios
# ══════════════════════════════════════════════════════════════════════════════
titv_obj <- titv(maf = laml_maf, plot = FALSE, useSyn = TRUE)

save_base_fig(
  expr   = plotTiTv(res = titv_obj),
  name   = "03_TiTv_ratio",
  subdir = "landscape",
  w = 14, h = 8
)

# Export TiTv summary
fwrite(titv_obj$TiTv.fractions,
       file.path(CONFIG$tables, "TiTv_per_sample.tsv"), sep = "\t")
fwrite(titv_obj$fraction.contribution,
       file.path(CONFIG$tables, "TiTv_fractions.tsv"), sep = "\t")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Tumour Mutation Burden (TMB)
# ══════════════════════════════════════════════════════════════════════════════
# WGS effective genome size ≈ 2,800 Mb (coding + non-coding)
# For TMB in WGS, divide by 2800 to get mut/Mb
smry <- getSampleSummary(laml_maf)
smry[, TMB_perMb := total / 2800]

p_tmb <- ggplot(smry, aes(x = reorder(Tumor_Sample_Barcode, -TMB_perMb),
                           y = TMB_perMb)) +
  geom_bar(stat = "identity", fill = "#2166AC", width = 0.85) +
  geom_hline(yintercept = median(smry$TMB_perMb, na.rm = TRUE),
             linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Tumour Mutation Burden — TCGA-LAML WGS",
       x = "Patient", y = "Mutations per Mb",
       caption = paste0("Dashed line = median (",
                        round(median(smry$TMB_perMb, na.rm = TRUE), 2),
                        " mut/Mb)")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

save_fig(p_tmb, "04_TMB_distribution", "landscape")

fwrite(smry[, .(Tumor_Sample_Barcode, total, TMB_perMb)],
       file.path(CONFIG$tables, "TMB_per_sample.tsv"), sep = "\t")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Variant Classification breakdown
# ══════════════════════════════════════════════════════════════════════════════
vc_summary <- as.data.table(smry)[, -(1:2)]  # remove barcode and total
vc_cols    <- setdiff(colnames(vc_summary), "TMB_perMb")

vc_long <- melt(
  as.data.table(smry)[, c("Tumor_Sample_Barcode", vc_cols), with = FALSE],
  id.vars       = "Tumor_Sample_Barcode",
  variable.name = "Variant_Class",
  value.name    = "Count"
)
vc_long <- vc_long[!is.na(Count) & Count > 0]

vc_totals <- vc_long[, .(Total = sum(Count)), by = Variant_Class][
  order(-Total)]

p_vc <- ggplot(vc_totals,
               aes(x = reorder(Variant_Class, Total), y = Total,
                   fill = Variant_Class)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = comma) +
  labs(title = "Variant Classification Summary",
       x = NULL, y = "Total Count") +
  theme_bw(base_size = 12)

save_fig(p_vc, "05_variant_classification", "landscape")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Top mutated genes — detailed barplot
# ══════════════════════════════════════════════════════════════════════════════
gene_smry <- as.data.table(getGeneSummary(laml_maf))
top30      <- gene_smry[1:min(30, nrow(gene_smry))]

top30[, pct := round(MutatedSamples / n_patients * 100, 1)]

p_genes <- ggplot(top30,
                  aes(x = reorder(Hugo_Symbol, MutatedSamples),
                      y = MutatedSamples)) +
  geom_col(fill = "#4477AA") +
  geom_text(aes(label = paste0(pct, "%")),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Top 30 Mutated Genes — TCGA-LAML WGS",
       x = NULL, y = "Number of Patients Mutated") +
  theme_bw(base_size = 12)

save_fig(p_genes, "06_top_mutated_genes", "landscape", w = 10, h = 10)

# ══════════════════════════════════════════════════════════════════════════════
# 7. VAF distribution per variant classification
# ══════════════════════════════════════════════════════════════════════════════
maf_dt <- as.data.table(laml_maf@data)

if ("t_VAF" %in% colnames(maf_dt)) {
  vaf_classes <- c("Missense_Mutation", "Nonsense_Mutation",
                   "Frame_Shift_Del", "Frame_Shift_Ins",
                   "Splice_Site", "In_Frame_Del", "In_Frame_Ins")
  vaf_plot_dt <- maf_dt[Variant_Classification %in% vaf_classes]

  p_vaf_box <- ggplot(vaf_plot_dt,
                      aes(x = reorder(Variant_Classification,
                                      -t_VAF, median),
                          y  = t_VAF,
                          fill = Variant_Classification)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),
                trim = TRUE, scale = "width", alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.2, alpha = 0.1) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "VAF Distribution by Variant Classification",
         x = NULL, y = "Variant Allele Frequency") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  save_fig(p_vaf_box, "07_VAF_by_class", "landscape", w = 12, h = 7)
}

# ══════════════════════════════════════════════════════════════════════════════
# 8. Rainfall plot (kataegis detection) for a representative sample
# ══════════════════════════════════════════════════════════════════════════════
# Pick the sample with the highest TMB
top_sample <- smry[which.max(total), Tumor_Sample_Barcode]

save_base_fig(
  expr   = rainfallPlot(maf = laml_maf,
                        tsb = top_sample,
                        detectChangePoints = TRUE,
                        pointSize = 0.4),
  name   = paste0("08_rainfall_", gsub("[^A-Za-z0-9]", "_", top_sample)),
  subdir = "landscape",
  w = 16, h = 6
)

message("\n=== Mutation landscape analysis complete ===")
