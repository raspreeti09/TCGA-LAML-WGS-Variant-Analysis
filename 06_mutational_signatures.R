# =============================================================================
# 06_mutational_signatures.R — Mutational signature analysis
#
# Steps:
#   A. Build trinucleotide mutation matrix (SBS 96-class)
#   B. De-novo NMF signature extraction (selecting optimal K)
#   C. Compare extracted signatures to COSMIC SBS v3.3
#   D. Decompose each sample into COSMIC reference signatures
#   E. Doublet-base (DBS) and indel (ID) signature matrices
#   F. Visualisations: profiles, cosine similarities, per-sample contributions
# =============================================================================

source("c:/Users/asus/OneDrive/Desktop/Project-LAML/config/config.R")

# ── Install BSgenome for hg38 (run once) ─────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")


library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
library(maftools)
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(NMF)            # NMF comes with MutationalPatterns if needed separately

set.seed(CONFIG$seed)

genome_ref <- BSgenome.Hsapiens.UCSC.hg38

# ── Load filtered MAF ─────────────────────────────────────────────────────────
maf_dt <- as.data.table(
  readRDS(file.path(CONFIG$data_proc, "maf", "laml_maf_filtered.rds"))@data
)

# ── A. Build GRanges from MAF (SNPs only for SBS) ─────────────────────────────
message("=== A. Building mutation catalogue ===")

snp_dt <- maf_dt[Variant_Type == "SNP",
                  .(sample    = Tumor_Sample_Barcode,
                    seqnames  = paste0("chr", sub("^chr", "", Chromosome)),
                    start     = Start_Position,
                    end       = End_Position,
                    REF       = Reference_Allele,
                    ALT       = Tumor_Seq_Allele2)]

# Build GRangesList: one GRanges per sample
gr_list <- split(
  GRanges(
    seqnames = snp_dt$seqnames,
    ranges   = IRanges(snp_dt$start, snp_dt$end),
    REF      = snp_dt$REF,
    ALT      = snp_dt$ALT,
    sample   = snp_dt$sample
  ),
  snp_dt$sample
)

message("Samples with SNPs: ", length(gr_list))

# Retrieve mutation context (trinucleotide)
mut_mat <- mut_matrix(
  vcf_list = gr_list,
  ref_genome = genome_ref
)

message("Mutation matrix dimensions: ",
        nrow(mut_mat), " trinucleotide contexts x ",
        ncol(mut_mat), " samples")

saveRDS(mut_mat,
        file.path(CONFIG$data_proc, "maf", "sbs96_mutation_matrix.rds"))

# Plot overall mutation profile
save_base_fig(
  expr   = plot_96_profile(mut_mat, condensed = TRUE, ymax = NULL),
  name   = "01_SBS96_profile_all_samples",
  subdir = "signatures",
  w = 16, h = 6
)

# ══════════════════════════════════════════════════════════════════════════════
# B. De-novo NMF signature extraction
# ══════════════════════════════════════════════════════════════════════════════
message("=== B. De-novo NMF extraction ===")

K_range <- CONFIG$n_signatures_min:CONFIG$n_signatures_max

estimate <- estimate_rank(
  mut_mat        = mut_mat,
  rank_range     = K_range,
  method         = "brunet",
  nrun           = CONFIG$n_bootstrap,
  verbose        = FALSE
)

# Plot NMF rank survey
save_base_fig(
  expr   = plot(estimate),
  name   = "02_NMF_rank_survey",
  subdir = "signatures",
  w = 12, h = 8
)

# Choose optimal K based on cophenetic coefficient (highest before sharp drop)
cophenetic_vals <- sapply(estimate$measures, function(m) m$cophenetic)
optimal_k       <- K_range[which.max(diff(cophenetic_vals) < -0.01) + 1]
if (is.na(optimal_k) || length(optimal_k) == 0) optimal_k <- 5

message("Optimal number of signatures (K) selected: ", optimal_k)

# Run NMF with optimal K
nmf_res <- extract_signatures(
  mut_mat   = mut_mat,
  rank      = optimal_k,
  nrun      = CONFIG$n_bootstrap,
  single_core = FALSE
)

# Rename signatures
colnames(nmf_res$signatures)    <- paste0("SBS_denovo_", seq_len(optimal_k))
rownames(nmf_res$contribution)  <- paste0("SBS_denovo_", seq_len(optimal_k))

# Plot de-novo signature profiles
save_base_fig(
  expr   = plot_96_profile(nmf_res$signatures, condensed = TRUE),
  name   = "03_denovo_signature_profiles",
  subdir = "signatures",
  w = 16, h = 4 * optimal_k
)

# Plot per-sample contribution
save_base_fig(
  expr   = plot_contribution(
    nmf_res$contribution,
    nmf_res$signatures,
    mode       = "relative",
    coord_flip = FALSE
  ),
  name   = "04_denovo_contribution_relative",
  subdir = "signatures",
  w = 18, h = 7
)

# ══════════════════════════════════════════════════════════════════════════════
# C. Match de-novo signatures to COSMIC SBS v3.3
# ══════════════════════════════════════════════════════════════════════════════
message("=== C. Matching to COSMIC SBS v3.3 ===")

# Load COSMIC SBS v3.3 signatures (built into MutationalPatterns)
cosmic_sbs <- get_known_signatures(
  source  = "COSMIC_v3.3",
  sig_type = "SBS"
)

cos_sim_matrix <- cos_sim_matrix(nmf_res$signatures, cosmic_sbs)

save_base_fig(
  expr   = plot_cosine_heatmap(
    cos_sim = cos_sim_matrix,
    col_order = NULL,
    cluster_rows = TRUE,
    cluster_cols = FALSE
  ),
  name   = "05_cosine_similarity_COSMIC",
  subdir = "signatures",
  w = 20, h = 6
)

# Write cosine similarity table
fwrite(as.data.table(cos_sim_matrix, keep.rownames = "denovo_sig"),
       file.path(CONFIG$tables, "cosine_similarity_COSMIC_SBS.tsv"),
       sep = "\t")

# Best COSMIC match per de-novo signature
best_match <- apply(cos_sim_matrix, 1, function(row) {
  names(row)[which.max(row)]
})
match_df <- data.frame(
  denovo_sig   = names(best_match),
  cosmic_match = best_match,
  cosine_sim   = apply(cos_sim_matrix, 1, max)
)
fwrite(as.data.table(match_df),
       file.path(CONFIG$tables, "denovo_to_COSMIC_match.tsv"),
       sep = "\t")
print(match_df)

# ══════════════════════════════════════════════════════════════════════════════
# D. Refitting samples to COSMIC reference signatures
# ══════════════════════════════════════════════════════════════════════════════
message("=== D. COSMIC refitting ===")

# Select biologically relevant COSMIC signatures for AML
# SBS1 (ageing), SBS2 (APOBEC), SBS13 (APOBEC), SBS5 (ageing),
# SBS6, SBS15, SBS20, SBS26 (MMR), SBS17 (oxidative damage),
# SBS18 (ROS), SBS22 (AA), SBS35 (platinum), SBS3 (HRD)
aml_relevant_sigs <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS6",
                        "SBS13", "SBS15", "SBS17a", "SBS17b",
                        "SBS18", "SBS22a", "SBS22b", "SBS26",
                        "SBS35", "SBS40a", "SBS40b", "SBS40c")

available_sigs <- intersect(aml_relevant_sigs, colnames(cosmic_sbs))
message("Using COSMIC signatures: ", paste(available_sigs, collapse = ", "))

cosmic_subset <- cosmic_sbs[, available_sigs, drop = FALSE]

refit <- fit_to_signatures(
  mut_matrix  = mut_mat,
  signatures  = cosmic_subset
)

# Absolute contributions
save_base_fig(
  expr   = plot_contribution(
    refit$contribution,
    cosmic_subset,
    mode       = "absolute",
    coord_flip = FALSE
  ),
  name   = "06_COSMIC_refit_absolute",
  subdir = "signatures",
  w = 18, h = 8
)

# Relative contributions
save_base_fig(
  expr   = plot_contribution(
    refit$contribution,
    cosmic_subset,
    mode       = "relative",
    coord_flip = FALSE
  ),
  name   = "07_COSMIC_refit_relative",
  subdir = "signatures",
  w = 18, h = 7
)

# Bootstrap refit for confidence intervals
refit_boot <- fit_to_signatures_bootstrapped(
  mut_matrix  = mut_mat,
  signatures  = cosmic_subset,
  n_boots     = 100,
  method      = "NNLS"
)

save_base_fig(
  expr   = plot_bootstrapped_contribution(
    refit_boot,
    mode       = "relative",
    dominant   = FALSE
  ),
  name   = "08_COSMIC_bootstrap_contribution",
  subdir = "signatures",
  w = 18, h = 8
)

# Save contribution tables
fwrite(as.data.table(refit$contribution, keep.rownames = "signature"),
       file.path(CONFIG$tables, "COSMIC_refit_contribution.tsv"),
       sep = "\t")

# ══════════════════════════════════════════════════════════════════════════════
# E. Doublet-base substitutions (DBS) matrix
# ══════════════════════════════════════════════════════════════════════════════
message("=== E. DBS signatures ===")

dbs_mat <- mut_matrix_stranded(
  vcf_list   = gr_list,
  ref_genome = genome_ref,
  ranges     = NULL,
  mode       = "DBS"
)

if (!is.null(dbs_mat) && nrow(dbs_mat) > 0) {
  cosmic_dbs <- get_known_signatures(source = "COSMIC_v3.3", sig_type = "DBS")

  refit_dbs <- fit_to_signatures(
    mut_matrix = dbs_mat,
    signatures = cosmic_dbs
  )

  save_base_fig(
    expr   = plot_contribution(refit_dbs$contribution, cosmic_dbs,
                               mode = "relative", coord_flip = FALSE),
    name   = "09_DBS_refit_contribution",
    subdir = "signatures",
    w = 18, h = 7
  )

  fwrite(as.data.table(refit_dbs$contribution, keep.rownames = "signature"),
         file.path(CONFIG$tables, "DBS_refit_contribution.tsv"),
         sep = "\t")
}

# ══════════════════════════════════════════════════════════════════════════════
# F. Indel (ID) signatures
# ══════════════════════════════════════════════════════════════════════════════
message("=== F. ID (indel) signatures ===")

# Extract indels from MAF
indel_dt <- maf_dt[Variant_Type %in% c("INS", "DEL"),
                    .(sample   = Tumor_Sample_Barcode,
                      seqnames = paste0("chr", sub("^chr", "", Chromosome)),
                      start    = Start_Position,
                      end      = End_Position,
                      REF      = Reference_Allele,
                      ALT      = Tumor_Seq_Allele2)]

if (nrow(indel_dt) > 0) {
  gr_indels <- split(
    GRanges(
      seqnames = indel_dt$seqnames,
      ranges   = IRanges(indel_dt$start, indel_dt$end),
      REF      = indel_dt$REF,
      ALT      = indel_dt$ALT
    ),
    indel_dt$sample
  )

  id_mat <- get_indel_context(gr_indels, ref_genome = genome_ref)
  id_counts <- count_indel_contexts(id_mat)

  save_base_fig(
    expr   = plot_indel_contexts(id_counts, condensed = TRUE),
    name   = "10_ID83_indel_profile",
    subdir = "signatures",
    w = 16, h = 6
  )

  cosmic_id <- get_known_signatures(source = "COSMIC_v3.3", sig_type = "ID")
  refit_id  <- fit_to_signatures(id_counts, cosmic_id)

  save_base_fig(
    expr   = plot_contribution(refit_id$contribution, cosmic_id,
                               mode = "relative", coord_flip = FALSE),
    name   = "11_ID_refit_contribution",
    subdir = "signatures",
    w = 18, h = 7
  )

  fwrite(as.data.table(refit_id$contribution, keep.rownames = "signature"),
         file.path(CONFIG$tables, "ID_refit_contribution.tsv"),
         sep = "\t")
}

message("\n=== Mutational signature analysis complete ===")
message("Optimal K = ", optimal_k)
print(match_df)
