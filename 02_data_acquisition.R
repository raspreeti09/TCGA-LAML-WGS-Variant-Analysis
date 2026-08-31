# =============================================================================
# 02_data_acquisition.R — Download TCGA-LAML WGS data via TCGAbiolinks
#
# Downloads:
#   A. Somatic MAF  — consensus masked somatic mutations (WGS)
#   B. CNV segments — masked copy number segments (WGS)
#   C. SV           — structural rearrangements (WGS)
#   D. Clinical     — patient clinical / biospecimen metadata
# =============================================================================

source("config/config.R")

library(TCGAbiolinks)
library(dplyr)
library(data.table)

set.seed(CONFIG$seed)

# ── A. Somatic mutations (MAF) ─────────────────────────────────────────────────
message("\n=== A. Querying somatic mutation MAF ===")

query_maf <- GDCquery(
  project          = CONFIG$tcga_project,
  data.category    = CONFIG$data_category_mut,
  data.type        = CONFIG$data_type_maf,
  access           = "open",
  workflow.type    = CONFIG$workflow_maf
)

# Show available files
maf_results <- getResults(query_maf)
message("MAF files found: ", nrow(maf_results))
print(head(maf_results[, c("cases", "file_name", "data_type",
                            "experimental_strategy")]))

# Download
GDCdownload(
  query     = query_maf,
  method    = "api",
  directory = file.path(CONFIG$data_raw, "maf")
)

# Prepare (merge downloaded MAF files into a single object)
maf_data <- GDCprepare(
  query     = query_maf,
  directory = file.path(CONFIG$data_raw, "maf"),
  save      = TRUE,
  save.filename = file.path(CONFIG$data_proc, "maf", "TCGA_LAML_WGS_raw.RData")
)

message("MAF object class: ", class(maf_data))
message("MAF dimensions: ", nrow(maf_data), " variants x ",
        ncol(maf_data), " columns")

# Save as TSV for inspection
fwrite(as.data.table(maf_data),
       file.path(CONFIG$data_proc, "maf", "TCGA_LAML_WGS_raw.maf.gz"),
       sep = "\t", compress = "gzip")

# ── B. Copy Number Variation segments ─────────────────────────────────────────
message("\n=== B. Querying CNV segments ===")

query_cnv <- GDCquery(
  project       = CONFIG$tcga_project,
  data.category = CONFIG$data_category_cnv,
  data.type     = CONFIG$data_type_seg,
  access        = "open"
)

cnv_results <- getResults(query_cnv)
message("CNV segment files found: ", nrow(cnv_results))

GDCdownload(
  query     = query_cnv,
  method    = "api",
  directory = file.path(CONFIG$data_raw, "cnv")
)

cnv_data <- GDCprepare(
  query     = query_cnv,
  directory = file.path(CONFIG$data_raw, "cnv"),
  save      = TRUE,
  save.filename = file.path(CONFIG$data_proc, "cnv", "TCGA_LAML_CNV_raw.RData")
)

message("CNV data dimensions: ", nrow(cnv_data), " x ", ncol(cnv_data))

fwrite(as.data.table(cnv_data),
       file.path(CONFIG$data_proc, "cnv", "TCGA_LAML_CNV_segments.tsv"),
       sep = "\t")

# ── C. Structural Variants ─────────────────────────────────────────────────────
message("\n=== C. Querying Structural Variants ===")

query_sv <- tryCatch({
  GDCquery(
    project       = CONFIG$tcga_project,
    data.category = CONFIG$data_category_sv,
    data.type     = CONFIG$data_type_sv,
    access        = "open"
  )
}, error = function(e) {
  message("SV query failed (may not be available as open access): ", e$message)
  NULL
})

if (!is.null(query_sv)) {
  sv_results <- getResults(query_sv)
  message("SV files found: ", nrow(sv_results))

  GDCdownload(
    query     = query_sv,
    method    = "api",
    directory = file.path(CONFIG$data_raw, "sv")
  )

  sv_data <- GDCprepare(
    query     = query_sv,
    directory = file.path(CONFIG$data_raw, "sv"),
    save      = TRUE,
    save.filename = file.path(CONFIG$data_raw, "sv", "TCGA_LAML_SV_raw.RData")
  )

  fwrite(as.data.table(sv_data),
         file.path(CONFIG$data_proc, "maf", "TCGA_LAML_SV_raw.tsv"),
         sep = "\t")
} else {
  message("Skipping SV download.")
}

# ── D. Clinical data ──────────────────────────────────────────────────────────
message("\n=== D. Downloading clinical data ===")

# Primary clinical (GDC harmonised)
query_clin <- GDCquery(
  project       = CONFIG$tcga_project,
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"
)

GDCdownload(
  query     = query_clin,
  method    = "api",
  directory = file.path(CONFIG$data_raw, "clinical")
)

clinical_data <- GDCprepare(
  query     = query_clin,
  directory = file.path(CONFIG$data_raw, "clinical")
)

# clinical_data is a list of BCR Biotab tables
# Save each table
for (tbl_name in names(clinical_data)) {
  outfile <- file.path(CONFIG$data_raw, "clinical",
                       paste0("TCGA_LAML_", tbl_name, ".tsv"))
  fwrite(as.data.table(clinical_data[[tbl_name]]), outfile, sep = "\t")
  message("Saved clinical table: ", tbl_name)
}

# Save the patient clinical table as the simplified clinical file
# (GDCquery_clinic is deprecated in newer TCGAbiolinks — use prepared data directly)
clin_simple <- clinical_data[["clinical_patient_laml"]]
fwrite(as.data.table(clin_simple),
       file.path(CONFIG$data_raw, "clinical", "TCGA_LAML_clinical_simple.tsv"),
       sep = "\t")

message("\n=== Data acquisition complete ===")
message("MAF  → ", file.path(CONFIG$data_proc, "maf"))
message("CNV  → ", file.path(CONFIG$data_proc, "cnv"))
message("Clin → ", file.path(CONFIG$data_raw,  "clinical"))

### MD5 checksum verification

library(tools)   # for md5sum()

verify_md5 <- function(query, download_dir) {
  results  <- getResults(query)

  # GDCdownload saves files as: dir/project/category/.../file_id/filename
  downloaded <- list.files(download_dir, recursive = TRUE, full.names = TRUE)

  results$md5_status <- sapply(results$file_name, function(fname) {
    match <- downloaded[basename(downloaded) == fname]
    if (length(match) == 0) return("FILE NOT FOUND")

    computed <- md5sum(match[1])
    expected <- results$md5sum[results$file_name == fname]

    if (is.na(expected) || expected == "") return("NO REFERENCE MD5")
    if (computed == expected) "OK" else paste0("MISMATCH (got ", computed, ")")
  })

  results[, c("file_name", "md5sum", "md5_status")]
}

# Run on each query object (must still be in your R session)
md5_maf  <- verify_md5(query_maf,  file.path(CONFIG$data_raw, "maf"))
md5_cnv  <- verify_md5(query_cnv,  file.path(CONFIG$data_raw, "cnv"))
md5_clin <- verify_md5(query_clin, file.path(CONFIG$data_raw, "clinical"))

print(md5_maf)
print(md5_cnv)
print(md5_clin)
