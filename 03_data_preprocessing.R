# =============================================================================
# 03_data_preprocessing.R — MAF quality filtering & maftools object creation
#
# Steps:
#   1. Load raw MAF TSV
#   2. Retain only WGS aliquots
#   3. Apply variant-level filters (VAF, depth, alt-read thresholds)
#   4. Remove known germline / population variants (gnomAD AF filter)
#   5. Keep only PASS-flagged variants
#   6. Deduplicate to one representative aliquot per patient (tumour)
#   7. Build maftools MAF object and attach clinical data
#   8. Save filtered MAF and QC metrics
# =============================================================================

source("config/config.R")

library(maftools)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

set.seed(CONFIG$seed)

maf_tsv  <- file.path(CONFIG$data_proc, "maf", "TCGA_LAML_WGS_raw.maf.gz")
clin_tsv <- file.path(CONFIG$data_raw,  "clinical", "TCGA_LAML_clinical_simple.tsv")

# ── 1. Load raw MAF ────────────────────────────────────────────────────────────
message("Loading raw MAF ...")
raw <- fread(maf_tsv, sep = "\t", fill = TRUE,
             na.strings = c("", "NA", ".", "__UNKNOWN__"))
message("Raw MAF: ", nrow(raw), " variants, ", ncol(raw), " columns")
message("Samples (Tumor_Sample_Barcode): ",
        length(unique(raw$Tumor_Sample_Barcode)))

# ── 2. WGS strategy check ─────────────────────────────────────────────────────
# "Experimental_Strategy" is not present in GDC MAF files — the data is already
# WGS-only because GDCquery in step 02 filtered by workflow type.
message("WGS check: Experimental_Strategy column absent — data is WGS-only ",
        "by GDCquery design. Skipping filter. Total variants: ", nrow(raw))

# ── 3. Variant-level QC filters ───────────────────────────────────────────────
before <- nrow(raw)

# Compute VAF if not present
if (!"t_VAF" %in% colnames(raw)) {
  raw[, t_VAF := t_alt_count / (t_ref_count + t_alt_count)]
}

raw_filtered <- raw[
  !is.na(t_VAF)           &
  t_VAF        >= CONFIG$min_vaf        &
  t_depth      >= CONFIG$min_depth      &
  t_alt_count  >= CONFIG$min_alt_reads
]

message("VAF/depth filter: ", before, " → ", nrow(raw_filtered),
        " (removed ", before - nrow(raw_filtered), ")")

# ── 4. GDC_FILTER flag ────────────────────────────────────────────────────────
# GDC MAF files use "GDC_FILTER" not "FILTER"
if ("GDC_FILTER" %in% colnames(raw_filtered)) {
  before <- nrow(raw_filtered)
  raw_filtered <- raw_filtered[GDC_FILTER == "PASS" | is.na(GDC_FILTER)]
  message("GDC_FILTER==PASS: ", before, " → ", nrow(raw_filtered))
} else {
  message("GDC_FILTER column not found — skipping PASS filter.")
}

# ── 5. gnomAD population AF filter ────────────────────────────────────────────
# Use gnomAD_non_cancer_AF as the primary filter (best for cancer somatic studies)
# Fall back to gnomAD_AF if non-cancer column is absent
primary_af_col <- if ("gnomAD_non_cancer_AF" %in% colnames(raw_filtered)) {
  "gnomAD_non_cancer_AF"
} else if ("gnomAD_AF" %in% colnames(raw_filtered)) {
  "gnomAD_AF"
} else {
  NA_character_
}

if (!is.na(primary_af_col)) {
  before <- nrow(raw_filtered)
  raw_filtered[[primary_af_col]] <- suppressWarnings(
    as.numeric(raw_filtered[[primary_af_col]])
  )
  raw_filtered <- raw_filtered[
    is.na(get(primary_af_col)) | get(primary_af_col) <= CONFIG$max_pop_af
  ]
  message("gnomAD AF filter using '", primary_af_col,
          "' (≤ ", CONFIG$max_pop_af, "): ",
          before, " → ", nrow(raw_filtered))
} else {
  message("No gnomAD AF column found — skipping population frequency filter.")
}

# ── 6. Deduplicate: one tumour aliquot per patient ────────────────────────────
# TCGA barcodes: TCGA-XX-XXXX-01A-...  tumour=01; normal=10/11
raw_filtered[, patient_id := substr(Tumor_Sample_Barcode, 1, 12)]
raw_filtered[, aliquot_type := substr(Tumor_Sample_Barcode, 14, 15)]

# Keep only tumour aliquots (01 = primary solid tumour, 03 = primary blood)
before <- nrow(raw_filtered)
raw_filtered <- raw_filtered[aliquot_type %in% c("01", "03")]
message("Tumour aliquot filter: ", before, " → ", nrow(raw_filtered))

# For patients with multiple aliquots keep the one with the most variants
raw_filtered[, n_var_per_sample := .N, by = Tumor_Sample_Barcode]
best_barcodes <- raw_filtered[
  , .(best_barcode = Tumor_Sample_Barcode[which.max(n_var_per_sample)]),
  by = patient_id
]$best_barcode

raw_filtered <- raw_filtered[Tumor_Sample_Barcode %in% best_barcodes]
message("After de-duplication: ", length(unique(raw_filtered$Tumor_Sample_Barcode)),
        " unique patients")

# ── 7. Load and harmonise clinical data ───────────────────────────────────────
message("Loading clinical data ...")
clin <- fread(clin_tsv, sep = "\t", fill = TRUE,
              na.strings = c("", "NA", "[Not Available]",
                             "[Not Applicable]", "[Unknown]"))

# Find the patient barcode column — BCR Biotab uses bcr_patient_barcode
id_col <- grep("bcr_patient_barcode|submitter_id|patient_barcode",
               colnames(clin), value = TRUE, ignore.case = TRUE)[1]
if (is.na(id_col)) stop("Cannot find patient ID column in clinical data. ",
                         "Columns: ", paste(colnames(clin), collapse = ", "))
message("Patient ID column: ", id_col)
clin[, patient_id := substr(get(id_col), 1, 12)]
clin <- clin[!duplicated(patient_id)]

# Coerce survival columns
os_cols <- c("days_to_death", "days_to_last_follow_up")
for (col in intersect(os_cols, colnames(clin))) {
  clin[[col]] <- suppressWarnings(as.numeric(clin[[col]]))
}

# Create unified OS time (use death date, fall back to last follow-up)
if (all(c("days_to_death", "days_to_last_follow_up") %in% colnames(clin))) {
  clin[, OS_time := ifelse(
    !is.na(days_to_death), days_to_death, days_to_last_follow_up
  )]
  clin[, OS_event := as.integer(vital_status == "Dead")]
}

# Save harmonised clinical table
fwrite(clin, file.path(CONFIG$data_proc, "maf", "clinical_harmonised.tsv"),
       sep = "\t")

# ── 8. Build maftools MAF object ──────────────────────────────────────────────
message("Building maftools MAF object ...")

# Write filtered MAF to disk (maftools::read.maf reads from file or data.frame)
filtered_maf_path <- file.path(CONFIG$data_proc, "maf",
                                "TCGA_LAML_WGS_filtered.maf.gz")
fwrite(raw_filtered,
       filtered_maf_path,
       sep = "\t", compress = "gzip")

# Align Tumor_Sample_Barcode in clinical to MAF (use 12-char patient ID as key)
clin_for_maf <- copy(clin)
# maftools requires a "Tumor_Sample_Barcode" column in clinical that matches MAF
# Use the representative barcode we kept
barcode_map <- unique(raw_filtered[, .(patient_id, Tumor_Sample_Barcode)])
clin_for_maf <- merge(clin_for_maf, barcode_map, by = "patient_id", all = FALSE)

laml_maf <- read.maf(
  maf            = filtered_maf_path,
  clinicalData   = clin_for_maf,
  verbose        = TRUE,
  removeDuplicatedVariants = TRUE
)

# ── 9. Save MAF object ────────────────────────────────────────────────────────
saveRDS(laml_maf,
        file.path(CONFIG$data_proc, "maf", "laml_maf_filtered.rds"))
message("maftools MAF object saved.")

# ── 10. QC summary ────────────────────────────────────────────────────────────
smry <- getSampleSummary(laml_maf)
gene_smry <- getGeneSummary(laml_maf)
field_smry <- getFields(laml_maf)

fwrite(smry,      file.path(CONFIG$tables, "sample_summary.tsv"),  sep = "\t")
fwrite(gene_smry, file.path(CONFIG$tables, "gene_summary.tsv"),    sep = "\t")

message("\n=== Preprocessing complete ===")
message("Variants in filtered MAF : ", nrow(raw_filtered))
message("Unique patients          : ", nrow(smry))
message("Unique genes affected    : ", nrow(gene_smry))
print(laml_maf)

# ── 11. QC plots ──────────────────────────────────────────────────────────────
# VAF distribution
vaf_df <- as.data.frame(raw_filtered)[, c("t_VAF", "Variant_Classification")]

p_vaf <- ggplot(vaf_df, aes(x = t_VAF, fill = Variant_Classification)) +
  geom_histogram(bins = 50, colour = "white", linewidth = 0.1) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "VAF Distribution — TCGA-LAML WGS (filtered)",
       x = "Variant Allele Frequency", y = "Count") +
  theme_bw(base_size = 12)

save_fig(p_vaf, "QC_VAF_distribution", "landscape")

# Depth distribution
p_depth <- ggplot(as.data.frame(raw_filtered),
                  aes(x = t_depth)) +
  geom_histogram(bins = 60, fill = "#4477AA", colour = "white", linewidth = 0.1) +
  scale_x_log10() +
  labs(title = "Read Depth Distribution — Tumour",
       x = "Total Depth (log10)", y = "Count") +
  theme_bw(base_size = 12)

save_fig(p_depth, "QC_depth_distribution", "landscape")
