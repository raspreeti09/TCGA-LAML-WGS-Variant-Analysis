# =============================================================================
# config.R — Central configuration for TCGA-LAML WGS Variant Analysis Pipeline
# =============================================================================

# ── Project root ─────────────────────────────────────────────────────────────
# Automatically determines the project root from the current working directory.
# This makes the pipeline portable across computers.

PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = FALSE)

CONFIG <- list(

  # Paths ----------------------------------------------------------------------
  root      = PROJECT_ROOT,
  data_raw  = file.path(PROJECT_ROOT, "data", "raw"),
  data_proc = file.path(PROJECT_ROOT, "data", "processed"),
  scripts   = file.path(PROJECT_ROOT, "scripts"),
  results   = file.path(PROJECT_ROOT, "results"),
  figs      = file.path(PROJECT_ROOT, "results", "figures"),
  tables    = file.path(PROJECT_ROOT, "results", "tables"),

  # TCGA project ---------------------------------------------------------------
  tcga_project = "TCGA-LAML",
  genome_build = "hg38",

  # Data types -----------------------------------------------------------------
  data_category_mut = "Simple Nucleotide Variation",
  data_type_maf     = "Masked Somatic Mutation",
  workflow_maf      = "Aliquot Ensemble Somatic Variant Merging and Masking",

  data_category_cnv = "Copy Number Variation",
  data_type_seg     = "Masked Copy Number Segment",
  workflow_cnv      = "DNAcopy",

  data_category_sv  = "Structural Variation",
  data_type_sv      = "Structural Rearrangement",

  # Filtering thresholds -------------------------------------------------------
  min_vaf       = 0.05,
  min_depth     = 10,
  min_alt_reads = 3,
  max_pop_af    = 0.001,

  # Signature analysis ---------------------------------------------------------
  n_signatures_min = 2,
  n_signatures_max = 10,
  n_bootstrap      = 100,
  sig_database     = "SBS",

  # Plotting -------------------------------------------------------------------
  fig_width  = 14,
  fig_height = 10,
  fig_dpi    = 300,
  fig_format = "pdf",

  # Survival analysis ----------------------------------------------------------
  os_time_col   = "days_to_death",
  os_status_col = "vital_status",
  alive_label   = "Alive",

  # Reproducibility ------------------------------------------------------------
  seed = 42
)

# ── Helper: create directory if absent ----------------------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  invisible(path)
}

# ── Helper: save a ggplot figure ----------------------------------------------
save_fig <- function(plot_obj, name, subdir = "",
                     w = CONFIG$fig_width,
                     h = CONFIG$fig_height,
                     dpi = CONFIG$fig_dpi,
                     fmt = CONFIG$fig_format) {

  dir_path <- if (nchar(subdir) > 0)
    file.path(CONFIG$figs, subdir)
  else
    CONFIG$figs

  ensure_dir(dir_path)

  filepath <- file.path(dir_path, paste0(name, ".", fmt))

  if (fmt == "pdf") {

    pdf(filepath, width = w, height = h)
    print(plot_obj)
    dev.off()

  } else {

    ggplot2::ggsave(
      filepath,
      plot = plot_obj,
      width = w,
      height = h,
      dpi = dpi
    )
  }

  message("Saved: ", filepath)
  invisible(filepath)
}

# ── Helper: save a base-R / maftools plot -------------------------------------
save_base_fig <- function(expr, name, subdir = "",
                          w = CONFIG$fig_width,
                          h = CONFIG$fig_height,
                          dpi = CONFIG$fig_dpi,
                          fmt = CONFIG$fig_format) {

  dir_path <- if (nchar(subdir) > 0)
    file.path(CONFIG$figs, subdir)
  else
    CONFIG$figs

  ensure_dir(dir_path)

  filepath <- file.path(dir_path, paste0(name, ".", fmt))

  if (fmt == "pdf") {

    pdf(filepath, width = w, height = h)

  } else {

    png(
      filepath,
      width = w * dpi,
      height = h * dpi,
      res = dpi,
      units = "px"
    )
  }

  force(expr)
  dev.off()

  message("Saved: ", filepath)
  invisible(filepath)
}

message("Config loaded.")
message("Project root: ", CONFIG$root)
