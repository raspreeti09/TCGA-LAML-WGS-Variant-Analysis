# =============================================================================
# 00_run_pipeline.R — Master runner for the full TCGA-LAML WGS pipeline
#
# Usage:
#   Rscript scripts/00_run_pipeline.R             # run all steps
#   Rscript scripts/00_run_pipeline.R 3 4 5       # run only steps 3, 4, 5
#
# Steps:
#   1  — Package setup
#   2  — Data acquisition (TCGA GDC download)
#   3  — Data preprocessing & MAF object
#   4  — Mutation landscape
#   5  — Driver gene analysis
#   6  — Mutational signatures
#   7  — Co-mutation analysis
#   8  — CNV analysis
#   9  — Clinical correlations
#   10 — Survival analysis
#   11 — Summary visualisation
# =============================================================================

source("c:/Users/asus/OneDrive/Desktop/Project-LAML/config/config.R")

SCRIPT_DIR <- CONFIG$scripts

step_scripts <- c(
  "1"  = "01_setup.R",
  "2"  = "02_data_acquisition.R",
  "3"  = "03_data_preprocessing.R",
  "4"  = "04_mutation_landscape.R",
  "5"  = "05_driver_genes.R",
  "6"  = "06_mutational_signatures.R",
  "7"  = "07_comutation_analysis.R",
  "8"  = "08_cnv_analysis.R",
  "9"  = "09_clinical_correlation.R",
  "10" = "10_survival_analysis.R",
  "11" = "11_summary_visualization.R"
)

# Parse CLI args (which steps to run)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  steps_to_run <- names(step_scripts)
} else {
  steps_to_run <- args
}

# Validate step numbers
invalid <- setdiff(steps_to_run, names(step_scripts))
if (length(invalid) > 0)
  stop("Unknown steps: ", paste(invalid, collapse = ", "),
       "\nValid steps: ", paste(names(step_scripts), collapse = ", "))

# Run selected steps
log_df <- data.frame(
  step    = character(),
  script  = character(),
  status  = character(),
  elapsed = numeric(),
  stringsAsFactors = FALSE
)

for (step in steps_to_run) {
  script_path <- file.path(SCRIPT_DIR, step_scripts[step])
  message("\n", strrep("─", 60))
  message("STEP ", step, ": ", step_scripts[step])
  message(strrep("─", 60))

  t_start <- proc.time()

  result <- tryCatch({
    source(script_path, local = FALSE)
    "SUCCESS"
  }, error = function(e) {
    message("ERROR in step ", step, ": ", conditionMessage(e))
    paste0("ERROR: ", conditionMessage(e))
  }, warning = function(w) {
    message("WARNING in step ", step, ": ", conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  elapsed <- as.numeric((proc.time() - t_start)["elapsed"])

  log_df <- rbind(log_df, data.frame(
    step    = step,
    script  = step_scripts[step],
    status  = ifelse(is.null(result), "SUCCESS",
                     ifelse(grepl("ERROR", result), result, "SUCCESS")),
    elapsed = round(elapsed, 1),
    stringsAsFactors = FALSE
  ))
}

# Print run summary
message("\n", strrep("═", 60))
message("PIPELINE RUN SUMMARY")
message(strrep("═", 60))
print(log_df)

# Save log
log_path <- file.path(CONFIG$results, "pipeline_run_log.tsv")
write.table(log_df, log_path, sep = "\t", row.names = FALSE, quote = FALSE)
message("\nLog saved: ", log_path)

failed <- log_df[grepl("ERROR", log_df$status), "step"]
if (length(failed) > 0) {
  warning("The following steps FAILED: ", paste(failed, collapse = ", "))
} else {
  message("\nAll steps completed successfully.")
}
