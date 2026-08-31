# =============================================================================
# 11_summary_visualization.R — Comprehensive summary figures
#
# Produces:
#   A. Publication-ready integrated oncoprint (mutations + CNV + clinical)
#   B. Circos plot — chromosomal mutation density
#   C. Signature attribution summary heatmap
#   D. Summary dashboard (patchwork multi-panel)
# =============================================================================

source("c:/Users/asus/OneDrive/Desktop/Project-LAML/config/config.R")

library(maftools)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(viridis)

set.seed(CONFIG$seed)

laml_maf <- readRDS(
  file.path(CONFIG$data_proc, "maf", "laml_maf_filtered.rds")
)
maf_dt   <- as.data.table(laml_maf@data)
clin_dt  <- as.data.table(getClinicalData(laml_maf))
colnames(clin_dt) <- tolower(colnames(clin_dt))

gene_smry   <- as.data.table(getGeneSummary(laml_maf))
n_patients  <- nrow(getSampleSummary(laml_maf))
smry        <- as.data.table(getSampleSummary(laml_maf))
smry[, TMB_perMb := total / 2800]

# ══════════════════════════════════════════════════════════════════════════════
# A. Integrated oncoprint with clinical annotation tracks
# ══════════════════════════════════════════════════════════════════════════════
message("=== A. Integrated oncoprint ===")

# Clinical annotation tracks
annot_cols <- intersect(
  c("age_at_diagnosis", "gender", "race", "vital_status",
    "fab_classification", "cytogenetic_risk",
    "days_to_death", "OS_event", "OS_years"),
  colnames(clin_dt)
)

save_base_fig(
  expr = oncoplot(
    maf              = laml_maf,
    top              = 30,
    clinicalFeatures = if (length(annot_cols) > 0) annot_cols else NULL,
    sortByAnnotation = TRUE,
    removeNonMutated = FALSE,
    showTumorSampleBarcodes = FALSE,
    drawColBar       = TRUE,
    fontSize         = 0.55,
    titleFontSize    = 1,
    legendFontSize   = 0.6,
    annotationFontSize = 0.6
  ),
  name   = "integrated_oncoprint_top30",
  subdir = "landscape",
  w = 22, h = 15
)

# ══════════════════════════════════════════════════════════════════════════════
# B. Circos plot — chromosomal mutation density
# ══════════════════════════════════════════════════════════════════════════════
message("=== B. Circos plot ===")

# Mutation density per 1 Mb bin
maf_dt[, chr_std := paste0("chr", sub("^chr", "", Chromosome))]
maf_dt[, bin_1mb := floor(Start_Position / 1e6)]

density_dt <- maf_dt[, .N, by = .(chr_std, bin_1mb)]
density_dt <- density_dt[chr_std %in% paste0("chr", c(1:22, "X", "Y"))]

# Set up circlize
save_base_fig(
  expr = {
    circos.clear()
    circos.par(
      "track.height"    = 0.15,
      cell.padding      = c(0.02, 0, 0.02, 0),
      gap.after         = c(rep(1, 22), 5, 5)
    )

    circos.initializeWithIdeogram(
      plotType  = c("ideogram", "labels"),
      species   = "hg38",
      chromosome.index = paste0("chr", c(1:22, "X", "Y"))
    )

    # Add mutation density track
    bed_density <- data.frame(
      chr   = density_dt$chr_std,
      start = density_dt$bin_1mb * 1e6,
      end   = (density_dt$bin_1mb + 1) * 1e6,
      value = log1p(density_dt$N)
    )

    circos.genomicTrack(
      bed_density,
      numeric.column = 4,
      panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "h",
                            col = "#E63946", ...)
      },
      track.height = 0.2
    )

    title("Somatic Mutation Density — TCGA-LAML WGS", cex.main = 1.2)
    circos.clear()
  },
  name   = "circos_mutation_density",
  subdir = "landscape",
  w = 12, h = 12
)

# ══════════════════════════════════════════════════════════════════════════════
# C. Signature attribution summary heatmap
# ══════════════════════════════════════════════════════════════════════════════
message("=== C. Signature heatmap ===")

contrib_file <- file.path(CONFIG$tables, "COSMIC_refit_contribution.tsv")
if (file.exists(contrib_file)) {
  contrib_dt <- fread(contrib_file)
  contrib_mat <- as.matrix(contrib_dt[, -1, with = FALSE])
  rownames(contrib_mat) <- contrib_dt$signature

  # Normalise to relative proportions per sample
  contrib_rel <- sweep(contrib_mat, 2, colSums(contrib_mat), "/") * 100

  # Filter signatures present in > 5% of samples
  present_sigs <- rownames(contrib_rel)[
    apply(contrib_rel, 1, function(x) sum(x > 1) / ncol(contrib_rel)) > 0.05
  ]
  if (length(present_sigs) == 0) present_sigs <- rownames(contrib_rel)

  col_sig <- colorRamp2(c(0, 25, 50, 100),
                         c("white", "#FEE08B", "#FC8D59", "#D73027"))

  ht_sig <- Heatmap(
    contrib_rel[present_sigs, , drop = FALSE],
    name             = "% contribution",
    col              = col_sig,
    cluster_rows     = TRUE,
    cluster_columns  = TRUE,
    show_column_names = FALSE,
    row_names_gp     = gpar(fontsize = 9),
    column_title     = "COSMIC Signature Contribution per Sample",
    heatmap_legend_param = list(title = "% Contribution")
  )

  save_base_fig(
    expr   = draw(ht_sig),
    name   = "signature_attribution_heatmap",
    subdir = "signatures",
    w = 18, h = 8
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# D. Summary dashboard (4-panel patchwork)
# ══════════════════════════════════════════════════════════════════════════════
message("=== D. Summary dashboard ===")

# Panel 1: Mutation type distribution
vc_dt <- maf_dt[, .N, by = Variant_Classification]
p1 <- ggplot(vc_dt, aes(x = reorder(Variant_Classification, -N),
                         y = N, fill = Variant_Classification)) +
  geom_col(show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "A) Variant Classification",
       x = NULL, y = "Count") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# Panel 2: TMB distribution
p2 <- ggplot(smry, aes(x = TMB_perMb)) +
  geom_histogram(bins = 35, fill = "#4477AA", colour = "white") +
  geom_vline(xintercept = median(smry$TMB_perMb, na.rm = TRUE),
             colour = "firebrick", linetype = "dashed") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "B) Tumour Mutation Burden",
       x = "TMB (mut/Mb, log10)", y = "# Patients") +
  theme_bw(base_size = 10)

# Panel 3: Top 15 genes
top15 <- gene_smry[1:min(15, nrow(gene_smry))]
p3 <- ggplot(top15, aes(x = reorder(Hugo_Symbol, MutatedSamples),
                         y = MutatedSamples / n_patients * 100)) +
  geom_col(fill = "#E63946") +
  coord_flip() +
  labs(title = "C) Top Mutated Genes",
       x = NULL, y = "% Patients") +
  theme_bw(base_size = 10)

# Panel 4: Ti/Tv summary
titv_obj <- titv(maf = laml_maf, plot = FALSE, useSyn = TRUE)
if (!is.null(titv_obj)) {
  titv_frac <- titv_obj$fraction.contribution
  titv_long <- melt(titv_frac, id.vars = character(0))
  titv_long <- as.data.table(titv_long)
  setnames(titv_long, c("type", "fraction"))

  p4 <- ggplot(titv_long, aes(x = type, y = fraction, fill = type)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "D) Ti/Tv Fractions",
         x = NULL, y = "Fraction") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
} else {
  p4 <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "TiTv data unavailable")
}

dashboard <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title    = "TCGA-LAML WGS — Variant Analysis Summary",
    subtitle = paste0("n = ", n_patients, " patients"),
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 10))
  )

save_fig(dashboard, "summary_dashboard", "", w = 14, h = 10)

# ══════════════════════════════════════════════════════════════════════════════
# E. Print final pipeline summary to console
# ══════════════════════════════════════════════════════════════════════════════
message("\n", paste(rep("=", 60), collapse = ""))
message("TCGA-LAML WGS Variant Analysis Pipeline — COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("  Patients analysed       : %d", n_patients))
message(sprintf("  Total somatic variants  : %d", nrow(maf_dt)))
message(sprintf("  Median TMB (mut/Mb)     : %.3f",
                median(smry$TMB_perMb, na.rm = TRUE)))
message(sprintf("  Unique genes mutated    : %d", nrow(gene_smry)))
message(sprintf("  Top mutated gene        : %s (%d%%)",
                gene_smry$Hugo_Symbol[1],
                round(gene_smry$AlteredIn[1] * 100)))

message("\n  Output directories:")
message("    Figures  → ", CONFIG$figs)
message("    Tables   → ", CONFIG$tables)
message("    Reports  → ", CONFIG$results)
message(paste(rep("=", 60), collapse = ""))
