# =============================================================================
# 05_driver_genes.R — Driver gene identification in TCGA-LAML WGS
#
# Methods:
#   A. maftools::oncodrive    — positional clustering (OncodriveFML-like)
#   B. dndscv                 — dN/dS ratio test (non-synonymous vs synonymous)
#   C. Known AML driver panel — frequency & lollipop plots
#   D. Hotspot visualisation  — protein domain plots
#   E. Pathway enrichment     — of significant driver genes
# =============================================================================

source("config/config.R")

library(maftools)
library(dndscv)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(CONFIG$seed)

laml_maf <- readRDS(
  file.path(CONFIG$data_proc, "maf", "laml_maf_filtered.rds")
)

# ══════════════════════════════════════════════════════════════════════════════
# A. OncodriveCLUST (maftools implementation)
# ══════════════════════════════════════════════════════════════════════════════
message("=== A. OncodriveCLUST ===")

sig_genes_oc <- oncodrive(
  maf            = laml_maf,
  AACol          = "HGVSp_Short",   # amino-acid change column
  minMut         = 5,               # minimum mutations per gene
  pvalMethod     = "zscore"
)

save_base_fig(
  expr   = plotOncodrive(res  = sig_genes_oc,
                         fdrCutOff  = 0.1,
                         useFraction = TRUE,
                         labelSize  = 1.0),
  name   = "01_oncodrive_clusters",
  subdir = "drivers",
  w = 13, h = 11
)

fwrite(sig_genes_oc,
       file.path(CONFIG$tables, "oncodrive_significant_genes.tsv"),
       sep = "\t")

# ══════════════════════════════════════════════════════════════════════════════
# B. dN/dS driver detection (dndscv)
# ══════════════════════════════════════════════════════════════════════════════
message("=== B. dN/dS analysis (dndscv) ===")

# ── 1. Prepare mutation table ─────────────────────────────────────────────────
snp_dt <- as.data.table(laml_maf@data)[
  Variant_Type == "SNP",
  .(sampleID = Tumor_Sample_Barcode,
    chr      = sub("^chr", "", Chromosome),
    pos      = Start_Position,
    ref      = Reference_Allele,
    mut      = Tumor_Seq_Allele2)
] |> unique()

message("SNPs passed to dndscv: ", nrow(snp_dt))

# ── 2. Run dndscv ─────────────────────────────────────────────────────────────
dnds_result <- dndscv(
  mutations                    = as.data.frame(snp_dt),
  refdb                        = "hg38",
  max_muts_per_gene_per_sample = 3,
  outp                         = 3
)

# ── 3. Save global dN/dS rates ───────────────────────────────────────────────
fwrite(
  as.data.table(dnds_result$globaldnds, keep.rownames = "rate"),
  file.path(CONFIG$tables, "dndscv_global_dnds.tsv"),
  sep = "\t"
)
message("Global dN/dS rates:")
print(dnds_result$globaldnds)

# ── 4. Extract per-gene results and significant drivers ───────────────────────
sel <- dnds_result$sel_cv                          # one row per gene

fwrite(as.data.table(sel),
       file.path(CONFIG$tables, "dndscv_all_genes.tsv"), sep = "\t")

sig_genes_dnds <- sel[sel$qallsubs_cv < 0.1, ]    # FDR < 10 %
message("Significant driver genes (q < 0.1): ", nrow(sig_genes_dnds))
print(sig_genes_dnds[, c("gene_name", "wmis_cv", "wnon_cv", "qallsubs_cv")])

fwrite(as.data.table(sig_genes_dnds),
       file.path(CONFIG$tables, "dndscv_significant_genes.tsv"), sep = "\t")

# ── 5. Volcano plot ───────────────────────────────────────────────────────────
plot_df <- as.data.table(sel) |>
  dplyr::mutate(
    log10q     = -log10(qallsubs_cv + 1e-10),
    significant = qallsubs_cv < 0.1
  )

p_dnds <- ggplot(plot_df,
                 aes(x = wmis_cv, y = log10q)) +

  # All genes — grey
  geom_point(data = dplyr::filter(plot_df, !significant),
             colour = "#B4B2A9", size = 1.8, alpha = 0.5) +

  # Significant genes — purple
  geom_point(data = dplyr::filter(plot_df, significant),
             colour = "#534AB7", size = 2.5, alpha = 0.9) +

  # Labels for significant genes only
  ggrepel::geom_label_repel(
    data          = dplyr::filter(plot_df, significant),
    aes(label     = gene_name),
    size          = 3.8,
    fontface      = "bold",
    colour        = "#26215C",
    fill          = "#EEEDFE",
    label.padding = unit(0.25, "lines"),
    label.r       = unit(0.15, "lines"),
    box.padding   = unit(0.5,  "lines"),
    point.padding = unit(0.4,  "lines"),
    segment.colour= "#534AB7",
    segment.size  = 0.4,
    max.overlaps  = Inf,
    seed          = 42
  ) +

  # dN/dS = 1 reference line (neutral evolution)
  geom_vline(xintercept = 1,
             linetype = "dashed", colour = "#E24B4A", linewidth = 0.5) +

  # FDR threshold line
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed", colour = "#888780", linewidth = 0.4) +

  annotate("text", x = max(plot_df$wmis_cv, na.rm = TRUE),
           y = -log10(0.1) + 0.15,
           label = "FDR = 0.1", size = 3,
           colour = "#888780", hjust = 1) +

  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +

  labs(
    title    = "dN/dS driver genes — TCGA-LAML WGS",
    subtitle = "Missense dN/dS (wmis) vs global FDR  ·  dashed red = neutral evolution (dN/dS = 1)",
    x        = "Missense dN/dS ratio (wmis_cv)",
    y        = expression(-log[10](q-value))
  ) +

  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, colour = "#26215C"),
    plot.subtitle    = element_text(size = 10, colour = "#5F5E5A"),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#D3D1C7", linewidth = 0.3)
  )

ggsave(
  filename = file.path(CONFIG$figs, "drivers", "02_dndscv_volcano.pdf"),
  plot     = p_dnds,
  width    = 12,
  height   = 10,
  device   = cairo_pdf
)
message("Saved: ", file.path(CONFIG$figs, "drivers", "02_dndscv_volcano.pdf"))

# ══════════════════════════════════════════════════════════════════════════════
# C. Known AML driver genes — frequency analysis
# ══════════════════════════════════════════════════════════════════════════════
message("=== C. Known AML driver genes ===")

# Curated list of well-established AML driver genes
aml_drivers <- c(
  "FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2", "TET2",
  "RUNX1", "TP53", "CEBPA", "WT1", "NRAS", "KRAS",
  "PTPN11", "CBL", "KIT", "ASXL1", "EZH2", "SF3B1",
  "SRSF2", "U2AF1", "ZRSR2", "STAG2", "RAD21", "SMC1A",
  "SMC3", "PHF6", "KDM6A", "NF1", "TP53", "MPL",
  "SETBP1", "MLL", "KMT2A", "BCOR", "BCORL1"
)

gene_smry  <- as.data.table(getGeneSummary(laml_maf))
driver_freq <- gene_smry[Hugo_Symbol %in% aml_drivers][
  order(-MutatedSamples)]

n_patients <- nrow(getSampleSummary(laml_maf))
driver_freq[, pct := round(MutatedSamples / n_patients * 100, 1)]

p_drivers <- ggplot(driver_freq,
                    aes(x = reorder(Hugo_Symbol, pct), y = pct,
                        fill = pct)) +
  geom_col() +
  geom_text(aes(label = paste0(pct, "%")), hjust = -0.15, size = 3) +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, 100)) +
  labs(title = "Known AML Driver Gene Mutation Frequency",
       subtitle = paste0("TCGA-LAML WGS (n = ", n_patients, " patients)"),
       x = NULL, y = "% Patients Mutated",
       fill = "% Mutated") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

save_fig(p_drivers, "03_AML_driver_frequency", "drivers", w = 10, h = 12)
fwrite(driver_freq,
       file.path(CONFIG$tables, "AML_known_driver_frequency.tsv"), sep = "\t")

# ══════════════════════════════════════════════════════════════════════════════
# D. Lollipop plots for key AML drivers
# ══════════════════════════════════════════════════════════════════════════════
message("=== D. Lollipop plots ===")

key_genes <- c("FLT3", "NPM1", "DNMT3A", "IDH1", "IDH2",
               "TET2", "RUNX1", "TP53", "CEBPA", "NRAS")

for (gene in key_genes) {
  if (!gene %in% gene_smry$Hugo_Symbol) next
  tryCatch({
    save_base_fig(
      expr   = lollipopPlot(
        maf              = laml_maf,
        gene             = gene,
        AACol            = "HGVSp_Short",
        showMutationRate = TRUE,
        labelPos         = "all",
        cBioPortal       = TRUE,
        repel            = TRUE
      ),
      name   = paste0("04_lollipop_", gene),
      subdir = "drivers",
      w = 14, h = 6
    )
  }, error = function(e) {
    message("Lollipop skipped for ", gene, ": ", e$message)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# E. Pathway enrichment of driver genes
# ══════════════════════════════════════════════════════════════════════════════
message("=== E. Pathway enrichment ===")

# Combine driver genes from both methods
all_drivers <- union(
  sig_genes_dnds$gene_name,
  driver_freq[pct >= 5, Hugo_Symbol]
)
all_drivers <- unique(all_drivers)
message("Driver genes for enrichment: ", paste(all_drivers, collapse = ", "))

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(all_drivers, fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)

if (nrow(entrez_ids) > 3) {
  # KEGG enrichment
  kegg_res <- enrichKEGG(
    gene         = entrez_ids$ENTREZID,
    organism     = "hsa",
    pvalueCutoff = 0.05
  )

  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    p_kegg <- dotplot(kegg_res, showCategory = 20, font.size = 10) +
      labs(title = "KEGG Pathway Enrichment — AML Driver Genes")
    save_fig(p_kegg, "05_KEGG_driver_enrichment", "drivers", w = 12, h = 9)
    fwrite(as.data.table(kegg_res@result),
           file.path(CONFIG$tables, "KEGG_driver_enrichment.tsv"), sep = "\t")
  }

  # WikiPathways enrichment (replaces ReactomePA — no extra package needed)
  wiki_res <- tryCatch(
    enrichWP(gene     = entrez_ids$ENTREZID,
             organism = "Homo sapiens",
             pvalueCutoff = 0.05),
    error = function(e) NULL
  )

  if (!is.null(wiki_res) && nrow(wiki_res) > 0) {
    p_wiki <- dotplot(wiki_res, showCategory = 20, font.size = 9) +
      labs(title = "WikiPathways Enrichment — AML Driver Genes")
    save_fig(p_wiki, "06_WikiPathways_driver_enrichment", "drivers",
             w = 13, h = 10)
    fwrite(as.data.table(wiki_res@result),
           file.path(CONFIG$tables, "WikiPathways_driver_enrichment.tsv"),
           sep = "\t")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# F. AML pathway / category grouping (PlotPathways)
# ══════════════════════════════════════════════════════════════════════════════
message("=== F. maftools pathway grouping ===")

# Define AML-relevant pathway gene sets
aml_pathways <- list(
  Signalling      = c("FLT3", "KIT", "NRAS", "KRAS", "PTPN11", "CBL", "NF1",
                      "MPL"),
  Transcription   = c("RUNX1", "CEBPA", "WT1", "ETV6", "ERG", "GATA2",
                      "KMT2A"),
  Epigenetics_DNA = c("DNMT3A", "TET2", "IDH1", "IDH2"),
  Epigenetics_Hist= c("EZH2", "ASXL1", "KDM6A", "SETD2", "KDM2B"),
  Splicing        = c("SF3B1", "SRSF2", "U2AF1", "ZRSR2"),
  Cohesion        = c("STAG2", "RAD21", "SMC1A", "SMC3"),
  TumourSuppressor= c("TP53", "WT1", "PHF6"),
  NPM1            = c("NPM1")
)

# Convert named list to data.frame (Pathway, Gene) — required by oncoplot
pathway_df <- data.frame(
  Pathway = rep(names(aml_pathways), sapply(aml_pathways, length)),
  Gene    = unlist(aml_pathways),
  row.names = NULL
)

# Keep only genes present in the MAF (oncoplot errors if < 2 genes found)
pathway_df_filtered <- pathway_df[
  pathway_df$Gene %in% gene_smry$Hugo_Symbol, ]
message("Pathway genes found in MAF: ",
        nrow(pathway_df_filtered), " / ", nrow(pathway_df))

if (nrow(pathway_df_filtered) >= 2) {
  genes_to_plot <- unique(pathway_df_filtered$Gene)

  save_base_fig(
    expr   = oncoplot(
      maf                     = laml_maf,
      genes                   = genes_to_plot,
      removeNonMutated        = FALSE,
      showTumorSampleBarcodes = FALSE,
      fontSize                = 0.6,
      titleText               = "AML Driver Genes by Pathway — TCGA-LAML WGS"
    ),
    name   = "07_AML_pathway_oncoplot",
    subdir = "drivers",
    w = 18, h = 10
  )
} else {
  message("Fewer than 2 pathway genes found in MAF — skipping pathway oncoplot.")
}

message("\n=== Driver gene analysis complete ===")
