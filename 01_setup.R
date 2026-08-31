# =============================================================================
# 01_setup.R — Package installation & environment verification
# Run ONCE before the rest of the pipeline.
# =============================================================================

source("config/config.R")

# ── 1. CRAN packages ──────────────────────────────────────────────────────────
cran_pkgs <- c(
  "tidyverse",       # dplyr, ggplot2, readr, tidyr, stringr, forcats, purrr
  "data.table",      # fast I/O and manipulation
  "survival",        # survival analysis
  "survminer",       # KM plots
  "ggrepel",         # non-overlapping labels
  "patchwork",       # compose ggplots
  "RColorBrewer",    # palettes
  "viridis",         # colour scales
  "pheatmap",        # heatmaps
  "corrplot",        # correlation matrix plots
  "ggpubr",          # publication-ready ggplots
  "scales",          # axis formatting
  "reshape2",
  "gridExtra",       # grid arrange
  "gt",              # publication tables
  "knitr",           # reporting
  "rmarkdown"        # HTML/PDF reports
)

# ── 2. Bioconductor packages ──────────────────────────────────────────────────
bioc_pkgs <- c(
  "TCGAbiolinks",              # TCGA GDC data download
  "maftools",                  # MAF parsing, landscape, driver genes, oncoprint
  "MutationalPatterns",        # de-novo & reference mutational signatures
  "BSgenome.Hsapiens.UCSC.hg38", # human reference genome (hg38)
  "GenomicRanges",             # genomic interval operations
  "GenomicFeatures",           # transcript/gene models
  "TxDb.Hsapiens.UCSC.hg38.knownGene",  # gene coordinates hg38
  "org.Hs.eg.db",              # gene ID mapping
  "VariantAnnotation",         # VCF / variant annotation
  "BiocParallel",              # parallel backend
  "ComplexHeatmap",            # advanced heatmaps / oncoprint
  "circlize",                  # circular plots (Circos-style)
  "clusterProfiler",           # pathway / GO enrichment
  "DOSE",                      # disease ontology enrichment
  "ReactomePA"                 # Reactome pathway enrichment
)

# ── 3. Install missing packages ───────────────────────────────────────────────

# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) == 0) {
    message("All ", ifelse(bioc, "Bioconductor", "CRAN"),
            " packages already installed.")
    return(invisible(NULL))
  }
  message("Installing: ", paste(missing, collapse = ", "))
  if (bioc) {
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing, repos = "https://cloud.r-project.org",
                     dependencies = TRUE)
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

# ── 3a. BSgenome.Hsapiens.UCSC.hg38 (large ~800 MB — explicit install) ────────
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  message("Installing BSgenome.Hsapiens.UCSC.hg38 (may take several minutes) ...")
  BiocManager::install(
    "BSgenome.Hsapiens.UCSC.hg38", ask = FALSE, update = FALSE
  )
}

# ── 3b. ComplexHeatmap (explicit, sometimes skipped in batch install) ─────────
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  message("Installing ComplexHeatmap ...")
  BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
}

# ── 3c. ReactomePA (depends on graphite and reactome.db — install first) ──────
if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  message("Installing ReactomePA and dependencies ...")
  BiocManager::install(c("graphite", "reactome.db", "ReactomePA"),
                       ask = FALSE, update = FALSE)
}
if (requireNamespace("ReactomePA", quietly = TRUE) &&
    !requireNamespace("reactome.db", quietly = TRUE)) {
  message("Installing missing ReactomePA dependency: reactome.db ...")
  BiocManager::install("reactome.db", ask = FALSE, update = FALSE)
}

# ── 3d. dndscv (CRAN; fall back to GitHub if CRAN build fails) ────────────────
if (!requireNamespace("dndscv", quietly = TRUE)) {
  message("Installing dndscv from CRAN ...")
  install.packages("dndscv", repos = "https://cloud.r-project.org",
                   dependencies = TRUE)
}
if (!requireNamespace("dndscv", quietly = TRUE)) {
  message("CRAN install failed — trying GitHub (im3sanger/dndscv) ...")
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
  remotes::install_github("im3sanger/dndscv")
}

# ── 4. Verify all packages load ───────────────────────────────────────────────
all_pkgs    <- c(cran_pkgs, bioc_pkgs,
                 "BSgenome.Hsapiens.UCSC.hg38", "dndscv")
load_ok     <- sapply(all_pkgs, requireNamespace, quietly = TRUE)
failed_pkgs <- all_pkgs[!load_ok]

if (length(failed_pkgs) > 0) {
  message("\n--- Packages that FAILED to install ---")
  for (pkg in failed_pkgs) {
    message("  FAILED: ", pkg)
  }
  message("\nTry installing them manually:")
  pkg_str <- paste(failed_pkgs, collapse = '", "')
  message('  BiocManager::install(c("', pkg_str, '"))')
  warning(length(failed_pkgs), " package(s) could not be loaded: ",
          paste(failed_pkgs, collapse = ", "))
} else {
  message("\nAll packages verified successfully.")
}

# ── 5. Ensure directory tree exists ───────────────────────────────────────────
dirs_to_check <- c(
  file.path(CONFIG$data_raw,  "maf"),
  file.path(CONFIG$data_raw,  "cnv"),
  file.path(CONFIG$data_raw,  "sv"),
  file.path(CONFIG$data_raw,  "clinical"),
  file.path(CONFIG$data_proc, "maf"),
  file.path(CONFIG$data_proc, "cnv"),
  file.path(CONFIG$data_proc, "filtered"),
  file.path(CONFIG$figs,      "landscape"),
  file.path(CONFIG$figs,      "signatures"),
  file.path(CONFIG$figs,      "survival"),
  file.path(CONFIG$figs,      "cnv"),
  file.path(CONFIG$figs,      "drivers"),
  CONFIG$tables,
  CONFIG$results
)
invisible(lapply(dirs_to_check, ensure_dir))
message("Directory tree verified.")

# ── 6. Session info ───────────────────────────────────────────────────────────
session_file <- file.path(CONFIG$results, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_file)
message("Session info saved to: ", session_file)
