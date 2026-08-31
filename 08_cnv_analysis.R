# =============================================================================
# 08_cnv_analysis.R — Copy Number Variation analysis (WGS segments)
#
# Steps:
#   A. Load & QC copy number segment data
#   B. Segment-level analysis (chromosome arm gains/losses)
#   C. Focal CNA detection (GISTIC-like thresholds)
#   D. Recurrent CNA heatmap (ComplexHeatmap)
#   E. Known AML CNA regions (monosomy 7, del5q, +8, inv16, t(8;21))
#   F. CNV burden per sample
#   G. Integrate CNV with somatic mutations (comutation × CNA)
# =============================================================================

source("c:/Users/asus/OneDrive/Desktop/Project-LAML/config/config.R")

library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(maftools)
library(RColorBrewer)
library(patchwork)

set.seed(CONFIG$seed)

# ── Load CNV segment data ──────────────────────────────────────────────────────
seg_file <- file.path(CONFIG$data_proc, "cnv", "TCGA_LAML_CNV_segments.tsv")
message("Loading CNV segments from: ", seg_file)

seg_dt <- fread(seg_file, sep = "\t", fill = TRUE,
                na.strings = c("", "NA", ".", "__UNKNOWN__"))

# Standardise column names (GDC CNV segment files)
# Typical columns: GDC_Aliquot, Chromosome, Start, End, Num_Probes, Segment_Mean
# Harmonize column names robustly
colnames(seg_dt) <- tolower(gsub("[. ]", "_", colnames(seg_dt)))
message("CNV columns: ", paste(colnames(seg_dt), collapse = ", "))

# Map to standard names
name_map <- c(
  sample         = grep("sample|aliquot|case|barcode",
                        colnames(seg_dt), value = TRUE)[1],
  chr            = grep("^chr|chrom", colnames(seg_dt), value = TRUE)[1],
  start          = grep("^start", colnames(seg_dt), value = TRUE)[1],
  end            = grep("^end",   colnames(seg_dt), value = TRUE)[1],
  seg_mean       = grep("segment_mean|mean|log2", colnames(seg_dt),
                        value = TRUE)[1],
  num_probes     = grep("num_probes|probes|snp", colnames(seg_dt),
                        value = TRUE)[1]
)
name_map <- name_map[!is.na(name_map)]

setnames(seg_dt,
         old = unname(name_map),
         new = names(name_map),
         skip_absent = TRUE)

# Ensure numeric columns
for (col in c("start", "end", "seg_mean")) {
  if (col %in% colnames(seg_dt))
    seg_dt[[col]] <- suppressWarnings(as.numeric(seg_dt[[col]]))
}

seg_dt <- seg_dt[!is.na(seg_mean)]

# Clean patient ID (first 12 chars of TCGA barcode)
seg_dt[, patient_id := substr(sample, 1, 12)]

message("CNV segments: ", nrow(seg_dt), " | Patients: ",
        length(unique(seg_dt$patient_id)))

# ══════════════════════════════════════════════════════════════════════════════
# A. Genome-wide segment mean distribution
# ══════════════════════════════════════════════════════════════════════════════
p_seg_dist <- ggplot(seg_dt, aes(x = seg_mean)) +
  geom_histogram(bins = 80, fill = "#4477AA", colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed",
             colour = "firebrick", linewidth = 0.7) +
  labs(title = "CNV Segment Log2 Ratio Distribution — TCGA-LAML WGS",
       x = "Segment Mean (log2 ratio)", y = "Count",
       caption = "Dashed lines: ±0.3 threshold for amplification/deletion") +
  theme_bw(base_size = 12)

save_fig(p_seg_dist, "01_segment_mean_distribution", "cnv")

# ══════════════════════════════════════════════════════════════════════════════
# B. Chromosome arm-level analysis
# ══════════════════════════════════════════════════════════════════════════════
message("=== B. Chromosome arm analysis ===")

# Define chromosome arm boundaries (GRCh38 centromere positions, hg38)
centromeres_hg38 <- data.table(
  chr = c(paste0("chr", 1:22), "chrX", "chrY"),
  centromere = c(
    122503483, 92326171, 91218929, 49691432, 46280423, 58829970,
    58053104, 43838887, 42352429, 39686682, 51078348, 34856694,
    16000000, 16000000, 17000000, 35335801, 22263006, 15460898,
    24498980, 26436232, 10864560, 13000000, 60632012, 10104553
  ),
  chr_length = c(
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468, 156040895, 57227415
  )
)

# Assign each segment to chromosome arm
seg_dt[, chr_std := paste0("chr", sub("^chr", "", chr))]

seg_gr <- GRanges(
  seqnames = seg_dt$chr_std,
  ranges   = IRanges(seg_dt$start, seg_dt$end),
  seg_mean = seg_dt$seg_mean,
  patient  = seg_dt$patient_id
)

# Arm assignment function
assign_arm <- function(chr, start, end, centromeres) {
  cen_row <- centromeres[chr == chr]
  if (nrow(cen_row) == 0) return(NA_character_)
  mid <- (start + end) / 2
  ifelse(mid < cen_row$centromere, paste0(chr, "p"), paste0(chr, "q"))
}

seg_dt <- merge(seg_dt, centromeres_hg38[, .(chr = chr, centromere)],
                by.x = "chr_std", by.y = "chr", all.x = TRUE)
seg_dt[, arm := ifelse(
  (start + end) / 2 < centromere,
  paste0(chr_std, "p"),
  paste0(chr_std, "q")
)]

# Arm-level summary per patient (weighted mean by segment length)
arm_summary <- seg_dt[!is.na(arm),
  .(weighted_mean = weighted.mean(seg_mean, w = (end - start + 1)),
    total_bases   = sum(end - start + 1)),
  by = .(patient_id, arm)
]

arm_matrix <- dcast(arm_summary,
                    patient_id ~ arm,
                    value.var = "weighted_mean",
                    fill = 0)

arm_mat <- as.matrix(arm_matrix[, -1, with = FALSE])
rownames(arm_mat) <- arm_matrix$patient_id

# Sort arms
arm_order <- sort(colnames(arm_mat))
arm_mat   <- arm_mat[, arm_order, drop = FALSE]

# Threshold for display
arm_mat_thresh <- ifelse(arm_mat >  0.3,  1L,
                  ifelse(arm_mat < -0.3, -1L, 0L))

col_arm <- colorRamp2(c(-1, 0, 1),
                       c("#2166AC", "white", "#D73027"))

ht_arms <- Heatmap(
  arm_mat_thresh,
  name            = "CNA",
  col             = col_arm,
  cluster_rows    = TRUE,
  cluster_columns = FALSE,
  show_row_names  = FALSE,
  column_names_gp = gpar(fontsize = 7),
  column_title    = "Chromosome Arm CNA — TCGA-LAML WGS",
  heatmap_legend_param = list(
    title  = "CNA state",
    at     = c(-1, 0, 1),
    labels = c("Deletion", "Neutral", "Amplification")
  )
)

save_base_fig(
  expr   = draw(ht_arms),
  name   = "02_arm_CNA_heatmap",
  subdir = "cnv",
  w = 18, h = 10
)

# Arm frequency summary
arm_freq <- as.data.table(arm_mat_thresh, keep.rownames = "patient_id")
arm_freq_long <- melt(arm_freq, id.vars = "patient_id",
                       variable.name = "arm", value.name = "state")

arm_agg <- arm_freq_long[,
  .(pct_amp = mean(state ==  1) * 100,
    pct_del = mean(state == -1) * 100),
  by = arm
]

arm_agg_long <- rbind(
  arm_agg[, .(arm, type = "Amplification", pct = pct_amp)],
  arm_agg[, .(arm, type = "Deletion",      pct = -pct_del)]
)

p_arm_freq <- ggplot(arm_agg_long,
                     aes(x = arm, y = pct, fill = type)) +
  geom_col(position = "identity") +
  scale_fill_manual(values = c("Amplification" = "#D73027",
                               "Deletion"      = "#2166AC")) +
  geom_hline(yintercept = 0) +
  labs(title = "Chromosome Arm CNA Frequency — TCGA-LAML",
       x = "Chromosome Arm", y = "% Patients",
       fill = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))

save_fig(p_arm_freq, "03_arm_CNA_frequency", "cnv", w = 18, h = 7)

fwrite(arm_agg,
       file.path(CONFIG$tables, "arm_CNA_frequency.tsv"), sep = "\t")

# ══════════════════════════════════════════════════════════════════════════════
# C. Focal CNA detection
# ══════════════════════════════════════════════════════════════════════════════
message("=== C. Focal CNAs ===")

# Amplification: seg_mean > 0.3 (log2 ratio > 0.3 ≈ 3 copies if diploid)
# Deep deletion: seg_mean < -1 (log2 ratio < -1 ≈ hemizygous loss)
# Shallow del:   seg_mean < -0.3

seg_dt[, cna_class := ifelse(
  seg_mean >  1.0, "High-level Amplification",
  ifelse(seg_mean >  0.3, "Low-level Gain",
  ifelse(seg_mean < -1.0, "Homozygous Deletion",
  ifelse(seg_mean < -0.3, "Heterozygous Loss",
                           "Neutral"))))]

# Focal = segment < 5 Mb
seg_dt[, is_focal := (end - start) < 5e6]

focal_cna <- seg_dt[is_focal == TRUE &
                      cna_class %in% c("High-level Amplification",
                                        "Homozygous Deletion")]

message("Focal amplifications: ",
        nrow(focal_cna[cna_class == "High-level Amplification"]))
message("Focal homozygous deletions: ",
        nrow(focal_cna[cna_class == "Homozygous Deletion"]))

fwrite(focal_cna,
       file.path(CONFIG$tables, "focal_CNA_events.tsv"), sep = "\t")

# ══════════════════════════════════════════════════════════════════════════════
# D. Known AML CNA regions
# ══════════════════════════════════════════════════════════════════════════════
message("=== D. Known AML CNA regions ===")

# Chromosome counts for whole-chromosome events
chr_mean <- seg_dt[, .(mean_seg = weighted.mean(seg_mean,
                                                  w = (end - start + 1))),
                    by = .(patient_id, chr_std)]

# Monosomy 7 (chr7 deletion), trisomy 8 (+8)
chr_wide <- dcast(chr_mean, patient_id ~ chr_std,
                  value.var = "mean_seg", fill = 0)

count_chr_event <- function(chr_col, thresh_del = -0.5, thresh_amp = 0.5) {
  if (!chr_col %in% colnames(chr_wide)) return(NULL)
  x <- chr_wide[[chr_col]]
  data.table(
    event       = chr_col,
    n_amp       = sum(x >  thresh_amp, na.rm = TRUE),
    n_del       = sum(x <  thresh_del, na.rm = TRUE),
    pct_amp     = mean(x >  thresh_amp, na.rm = TRUE) * 100,
    pct_del     = mean(x <  thresh_del, na.rm = TRUE) * 100
  )
}

aml_chr_events <- rbindlist(lapply(
  paste0("chr", c(3, 5, 7, 8, 11, 12, 17, 21)),
  count_chr_event
))
fwrite(aml_chr_events,
       file.path(CONFIG$tables, "chr_level_CNA_AML.tsv"), sep = "\t")
print(aml_chr_events)

# ══════════════════════════════════════════════════════════════════════════════
# E. CNV burden per sample
# ══════════════════════════════════════════════════════════════════════════════
cnv_burden <- seg_dt[,
  .(n_segments     = .N,
    pct_amp_genome = sum((end - start + 1) *
                           (seg_mean >  0.3), na.rm = TRUE) /
                     sum(end - start + 1) * 100,
    pct_del_genome = sum((end - start + 1) *
                           (seg_mean < -0.3), na.rm = TRUE) /
                     sum(end - start + 1) * 100,
    mean_abs_log2  = mean(abs(seg_mean), na.rm = TRUE)),
  by = patient_id
]
cnv_burden[, total_altered := pct_amp_genome + pct_del_genome]

fwrite(cnv_burden,
       file.path(CONFIG$tables, "CNV_burden_per_sample.tsv"), sep = "\t")

p_burden <- ggplot(cnv_burden, aes(x = pct_amp_genome,
                                    y = pct_del_genome,
                                    colour = total_altered)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_colour_viridis_c(option = "inferno", direction = -1) +
  labs(title = "CNV Burden — TCGA-LAML WGS",
       x = "% Genome Amplified",
       y = "% Genome Deleted",
       colour = "% Total Altered") +
  theme_bw(base_size = 12)

save_fig(p_burden, "04_CNV_burden", "cnv")

# ══════════════════════════════════════════════════════════════════════════════
# F. maftools gistic plot (if GISTIC2 results are available)
# ══════════════════════════════════════════════════════════════════════════════
# If you run GISTIC2 externally and have all_lesions.conf_90.txt etc.,
# load results here:
gistic_all   <- file.path(CONFIG$data_raw, "cnv", "all_lesions.conf_90.txt")
gistic_amp   <- file.path(CONFIG$data_raw, "cnv", "amp_genes.conf_90.txt")
gistic_del   <- file.path(CONFIG$data_raw, "cnv", "del_genes.conf_90.txt")
gistic_score <- file.path(CONFIG$data_raw, "cnv", "scores.gistic")

if (all(file.exists(c(gistic_all, gistic_amp,
                       gistic_del, gistic_score)))) {
  message("GISTIC2 results found — loading ...")
  laml_gistic <- readGistic(
    gisticAllLesionsFile    = gistic_all,
    gisticAmpGenesFile      = gistic_amp,
    gisticDelGenesFile      = gistic_del,
    gisticScoresFile        = gistic_score,
    isTCGA                  = TRUE
  )

  save_base_fig(
    expr   = gisticChromPlot(gistic = laml_gistic, markBands = "all"),
    name   = "05_GISTIC_chromplot",
    subdir = "cnv",
    w = 18, h = 8
  )

  save_base_fig(
    expr   = gisticOncoPlot(
      gistic         = laml_gistic,
      maf            = readRDS(file.path(CONFIG$data_proc, "maf",
                                          "laml_maf_filtered.rds")),
      top            = 15
    ),
    name   = "06_GISTIC_oncoplot",
    subdir = "cnv",
    w = 16, h = 12
  )
} else {
  message("GISTIC2 result files not found — skipping GISTIC plots.")
  message("Place GISTIC2 output files in: ", file.path(CONFIG$data_raw, "cnv"))
}

message("\n=== CNV analysis complete ===")
