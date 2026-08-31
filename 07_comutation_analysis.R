# =============================================================================
# 07_comutation_analysis.R — Co-mutation & mutual exclusivity analysis
#
# Methods:
#   A. maftools::somaticInteractions  (Fisher's exact, co-occurrence score)
#   B. Custom pairwise Fisher's test for AML-specific gene pairs
#   C. Co-mutation heatmap & network graph
#   D. Clonality / VAF-based clonal hierarchy inference
# =============================================================================

source("c:/Users/asus/OneDrive/Desktop/Project-LAML/config/config.R")

library(maftools)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(RColorBrewer)

set.seed(CONFIG$seed)

laml_maf <- readRDS(
  file.path(CONFIG$data_proc, "maf", "laml_maf_filtered.rds")
)

n_patients <- nrow(getSampleSummary(laml_maf))
gene_smry  <- as.data.table(getGeneSummary(laml_maf))

# ── Select genes mutated in ≥ 5% of patients ──────────────────────────────────
min_pct <- 0.05
freq_genes <- gene_smry[MutatedSamples / n_patients >= min_pct, Hugo_Symbol]
message("Genes mutated in ≥ ", min_pct * 100, "% of patients: ",
        length(freq_genes))

# ══════════════════════════════════════════════════════════════════════════════
# A. somaticInteractions (maftools)
# ══════════════════════════════════════════════════════════════════════════════
message("=== A. somaticInteractions ===")

si_result <- somaticInteractions(
  maf         = laml_maf,
  top         = min(length(freq_genes), 30),
  pvalue      = c(0.05, 0.1),
  returnAll   = TRUE,
  plotAll     = FALSE
)

save_base_fig(
  expr   = somaticInteractions(
    maf       = laml_maf,
    top       = min(length(freq_genes), 25),
    pvalue    = c(0.05, 0.1),
    fontSize  = 0.8
  ),
  name   = "01_somatic_interactions",
  subdir = "landscape",
  w = 14, h = 12
)

fwrite(as.data.table(si_result),
       file.path(CONFIG$tables, "somatic_interactions.tsv"), sep = "\t")

sig_pairs <- as.data.table(si_result)[
  pValue < 0.05 & !is.na(pValue)][order(pValue)]
message("Significant gene pairs (p < 0.05): ", nrow(sig_pairs))

# ══════════════════════════════════════════════════════════════════════════════
# B. Custom pairwise Fisher's exact test
# ══════════════════════════════════════════════════════════════════════════════
message("=== B. Custom Fisher's test ===")

# Build binary mutation matrix (patients × genes)
maf_dt  <- as.data.table(laml_maf@data)
all_pts <- getSampleSummary(laml_maf)$Tumor_Sample_Barcode

mut_binary <- dcast(
  maf_dt[Hugo_Symbol %in% freq_genes,
          .(Tumor_Sample_Barcode, Hugo_Symbol)][
    , mutated := 1L],
  Tumor_Sample_Barcode ~ Hugo_Symbol,
  value.var = "mutated",
  fill      = 0L
)

# Ensure all patients appear (add rows for unmutated patients)
missing_pts <- setdiff(all_pts, mut_binary$Tumor_Sample_Barcode)
if (length(missing_pts) > 0) {
  empty <- matrix(0L,
                  nrow = length(missing_pts),
                  ncol = ncol(mut_binary) - 1,
                  dimnames = list(NULL, setdiff(colnames(mut_binary),
                                                "Tumor_Sample_Barcode")))
  empty_dt <- cbind(
    data.table(Tumor_Sample_Barcode = missing_pts),
    as.data.table(empty)
  )
  mut_binary <- rbind(mut_binary, empty_dt)
}

mat <- as.matrix(mut_binary[, -1, with = FALSE])
rownames(mat) <- mut_binary$Tumor_Sample_Barcode

# Save binary matrix
fwrite(as.data.table(mat, keep.rownames = "patient"),
       file.path(CONFIG$tables, "binary_mutation_matrix.tsv"),
       sep = "\t")

# All gene pairs Fisher test
genes_test <- colnames(mat)
n_genes    <- length(genes_test)
results    <- vector("list", n_genes * (n_genes - 1) / 2)
idx        <- 1L

for (i in seq_len(n_genes - 1)) {
  for (j in (i + 1):n_genes) {
    g1 <- genes_test[i]; g2 <- genes_test[j]
    a  <- sum(mat[, g1] == 1 & mat[, g2] == 1)   # both mutated
    b  <- sum(mat[, g1] == 1 & mat[, g2] == 0)   # only g1
    cc <- sum(mat[, g1] == 0 & mat[, g2] == 1)   # only g2
    d  <- sum(mat[, g1] == 0 & mat[, g2] == 0)   # neither

    ft <- fisher.test(matrix(c(a, b, cc, d), 2, 2))
    results[[idx]] <- data.table(
      gene1 = g1, gene2 = g2,
      both = a, only1 = b, only2 = cc, neither = d,
      OR    = ft$estimate,
      pval  = ft$p.value,
      lower = ft$conf.int[1],
      upper = ft$conf.int[2]
    )
    idx <- idx + 1L
  }
}

pairs_dt <- rbindlist(results)
pairs_dt[, padj := p.adjust(pval, method = "BH")]
pairs_dt[, interaction := ifelse(OR > 1, "Co-occurrence", "Mutual exclusivity")]
pairs_dt <- pairs_dt[order(padj)]

fwrite(pairs_dt,
       file.path(CONFIG$tables, "fisher_pairwise_interactions.tsv"),
       sep = "\t")

sig_pairs_custom <- pairs_dt[padj < 0.1]
message("Significant pairs after BH correction (FDR < 0.1): ",
        nrow(sig_pairs_custom))

# ── Bubble chart of significant interactions ──────────────────────────────────
if (nrow(sig_pairs_custom) > 0) {
  sig_pairs_custom[, pair_label := paste(gene1, "vs", gene2)]

  p_bubble <- ggplot(sig_pairs_custom,
                     aes(x = log2(OR + 0.01),
                         y = -log10(padj),
                         size = abs(log2(OR + 0.01)),
                         colour = interaction,
                         label = pair_label)) +
    geom_point(alpha = 0.75) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed",
               colour = "grey40") +
    scale_colour_manual(values = c("Co-occurrence"        = "#E63946",
                                   "Mutual exclusivity"   = "#457B9D")) +
    scale_size_continuous(range = c(2, 8), guide = "none") +
    labs(title = "Pairwise Gene Co-mutation Analysis",
         subtitle = "TCGA-LAML WGS | Fisher's exact test, BH-adjusted",
         x = "log2(Odds Ratio)", y = "-log10(FDR)",
         colour = "Interaction Type") +
    theme_bw(base_size = 12)

  save_fig(p_bubble, "02_comutation_bubble", "landscape")
}

# ══════════════════════════════════════════════════════════════════════════════
# C. Co-mutation heatmap (ComplexHeatmap)
# ══════════════════════════════════════════════════════════════════════════════
message("=== C. Co-mutation heatmap ===")

# Build OR matrix for top genes
top_n    <- min(20, length(genes_test))
top_genes <- gene_smry[1:top_n, Hugo_Symbol]
top_genes <- intersect(top_genes, colnames(mat))

or_mat   <- matrix(NA, nrow = length(top_genes), ncol = length(top_genes),
                   dimnames = list(top_genes, top_genes))
pval_mat <- matrix(1,  nrow = length(top_genes), ncol = length(top_genes),
                   dimnames = list(top_genes, top_genes))

for (i in seq_along(top_genes)) {
  for (j in seq_along(top_genes)) {
    if (i == j) next
    g1 <- top_genes[i]; g2 <- top_genes[j]
    a  <- sum(mat[, g1] == 1 & mat[, g2] == 1)
    b  <- sum(mat[, g1] == 1 & mat[, g2] == 0)
    cc <- sum(mat[, g1] == 0 & mat[, g2] == 1)
    d  <- sum(mat[, g1] == 0 & mat[, g2] == 0)
    ft <- tryCatch(fisher.test(matrix(c(a, b, cc, d), 2, 2)),
                   error = function(e) list(estimate = NA, p.value = 1))
    or_mat[i, j]   <- as.numeric(ft$estimate)
    pval_mat[i, j] <- ft$p.value
  }
}

diag(or_mat) <- 1

log_or_mat      <- log2(or_mat + 0.01)
col_fun         <- colorRamp2(c(-3, 0, 3),
                               c("#457B9D", "white", "#E63946"))
pval_layer      <- ifelse(pval_mat < 0.05, "*", "")

ht <- Heatmap(
  matrix        = log_or_mat,
  name          = "log2(OR)",
  col           = col_fun,
  cell_fun      = function(j, i, x, y, w, h, fill) {
    if (pval_mat[i, j] < 0.05)
      grid.text("*", x, y, gp = gpar(fontsize = 10))
  },
  cluster_rows  = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_names_gp  = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  column_title  = "Co-mutation Heatmap (log2 OR) — TCGA-LAML WGS",
  heatmap_legend_param = list(title = "log2(OR)")
)

save_base_fig(
  expr   = draw(ht),
  name   = "03_comutation_heatmap",
  subdir = "landscape",
  w = 12, h = 12
)

# ══════════════════════════════════════════════════════════════════════════════
# D. Network graph of significant co-mutation pairs
# ══════════════════════════════════════════════════════════════════════════════
message("=== D. Co-mutation network ===")

network_pairs <- pairs_dt[padj < 0.1 &
                            gene1 %in% top_genes & gene2 %in% top_genes]

if (nrow(network_pairs) > 2) {
  g_net <- graph_from_data_frame(
    d        = network_pairs[, .(gene1, gene2, OR, padj, interaction)],
    directed = FALSE
  )

  V(g_net)$size        <- 10
  E(g_net)$color       <- ifelse(network_pairs$OR > 1, "#E63946", "#457B9D")
  E(g_net)$width       <- -log10(network_pairs$padj + 1e-5) / 2

  # Add mutation frequency as node attribute
  V(g_net)$freq <- gene_smry[
    match(V(g_net)$name, Hugo_Symbol), MutatedSamples / n_patients
  ] * 100

  save_base_fig(
    expr = {
      plot(g_net,
           vertex.label      = V(g_net)$name,
           vertex.label.cex  = 0.75,
           vertex.color      = colorRampPalette(
             c("lightblue", "firebrick"))(20)[
               cut(V(g_net)$freq, breaks = 20)],
           vertex.size       = sqrt(V(g_net)$freq) * 3,
           edge.color        = E(g_net)$color,
           edge.width        = E(g_net)$width,
           layout            = layout_with_fr(g_net),
           main              = "Co-mutation Network — TCGA-LAML WGS")
      legend("bottomleft",
             legend = c("Co-occurrence", "Mutual exclusivity"),
             col    = c("#E63946", "#457B9D"),
             lwd    = 3, cex = 0.8, bty = "n")
    },
    name   = "04_comutation_network",
    subdir = "landscape",
    w = 12, h = 12
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# E. Clonal hierarchy inference via VAF ordering
# ══════════════════════════════════════════════════════════════════════════════
message("=== E. Clonal hierarchy (VAF-based) ===")

maf_full <- as.data.table(laml_maf@data)

if ("t_VAF" %in% colnames(maf_full)) {
  # For key AML driver genes compute median VAF across patients
  key_genes <- intersect(
    c("DNMT3A", "TET2", "IDH1", "IDH2", "ASXL1",   # clonal / pre-leukaemic
      "NPM1", "FLT3", "CEBPA", "RUNX1", "NRAS",     # sub-clonal
      "TP53", "WT1"),
    unique(maf_full$Hugo_Symbol)
  )

  vaf_summary <- maf_full[Hugo_Symbol %in% key_genes,
                            .(median_VAF = median(t_VAF, na.rm = TRUE),
                              mean_VAF   = mean(t_VAF,   na.rm = TRUE),
                              n_muts     = .N),
                            by = Hugo_Symbol][order(-median_VAF)]

  fwrite(vaf_summary,
         file.path(CONFIG$tables, "VAF_clonal_hierarchy.tsv"), sep = "\t")

  p_clonal <- ggplot(vaf_summary,
                     aes(x = reorder(Hugo_Symbol, median_VAF),
                         y = median_VAF, colour = median_VAF)) +
    geom_segment(aes(xend = Hugo_Symbol, y = 0, yend = median_VAF),
                 linewidth = 1.2) +
    geom_point(size = 4) +
    geom_hline(yintercept = 0.25, linetype = "dashed",
               colour = "grey40", linewidth = 0.7) +
    annotate("text", x = 1.2, y = 0.27,
             label = "Putative clonal threshold (VAF ≈ 0.25)", size = 3) +
    scale_colour_viridis_c(option = "plasma", direction = -1) +
    coord_flip() +
    labs(title = "Median VAF per Gene — Clonal Hierarchy Inference",
         x = NULL, y = "Median VAF",
         colour = "Median VAF") +
    theme_bw(base_size = 12)

  save_fig(p_clonal, "05_VAF_clonal_hierarchy", "landscape", w = 10, h = 8)
}

message("\n=== Co-mutation analysis complete ===")
