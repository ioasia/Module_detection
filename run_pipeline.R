# Full module detection pipeline.
# Run from the repo root: Rscript run_pipeline.R
#                    or:  source("run_pipeline.R")  (from an R session)
#
# Steps:
#   1. network_reconstruction.R  — build co-expression network
#   2. basin_finder.py (x2)      — find basins of attraction
#   3. basin_analysis.R          — filter, NMF clustering, figures
#
# Pass --skip-step1 and/or --skip-step2 as command-line arguments to reuse
# previously computed network or basin files.

repo_dir <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

# ============================================================
# CONFIG — edit these before running
# ============================================================

# Abundance matrix: TSV, proteins x samples, first column = protein names
proteomics    <- "Data/my_proteomics.txt"

# Sample metadata: TSV, samples x annotation columns, first column = sample names
# Set to NULL to skip heatmap annotation
metadata      <- NULL

# CORUM complex pairs file (used to train the classifier)
corum         <- "Data/CORUM_complex_pairs.txt"

# Output directories
data_dir      <- "Data"
figures_dir   <- "Figures"

# Short name used as prefix for all output filenames
prefix        <- "my_dataset"

# Probability threshold for keeping a network edge (inspect _graph_quality.pdf first)
pred_threshold <- "0.8"

# Number of NMF components for the final decomposition.
# Inspect the _NMF_cophenetic.pdf from the first run, then set this to the optimal k.
nmf_k         <- "10"

# ============================================================
# (nothing below here should need editing)
# ============================================================

log <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

run <- function(cmd, args = character()) {
  status <- system2(cmd, args = args)
  if (status != 0) stop(sprintf("Command failed (exit %d): %s %s", status, cmd, paste(args, collapse = " ")))
}

cli_args   <- commandArgs(trailingOnly = TRUE)
skip_step1 <- "--skip-step1" %in% cli_args
skip_step2 <- "--skip-step2" %in% cli_args

network    <- file.path(data_dir, paste0(prefix, "_igraph_predictions_sig.txt"))
basins_pos <- file.path(data_dir, paste0(prefix, "_basins.pkl"))
basins_neg <- file.path(data_dir, paste0(prefix, "_basins_neg.pkl"))

# ---- step 1 -----------------------------------------------------------------
if (!skip_step1) {
  log("Step 1: Building protein co-expression network...")
  run("Rscript", c(
    "--vanilla",
    file.path(repo_dir, "Code", "network_reconstruction.R"),
    paste0("--proteomics=",     proteomics),
    paste0("--corum=",          corum),
    paste0("--data-dir=",       data_dir),
    paste0("--figures-dir=",    figures_dir),
    paste0("--prefix=",         prefix),
    paste0("--pred-threshold=", pred_threshold)
  ))
  log("Step 1 done.")
} else {
  log("Step 1 skipped.")
}

# ---- step 2 -----------------------------------------------------------------
if (!skip_step2) {
  log("Step 2a: Finding basins (positive abundances)...")
  run("python", c(
    file.path(repo_dir, "Code", "basin_finder.py"),
    paste0("--network=",     network),
    paste0("--proteomics=",  proteomics),
    paste0("--output-dir=",  data_dir),
    paste0("--prefix=",      prefix),
    "--no-negate"
  ))

  log("Step 2b: Finding basins (negated abundances)...")
  run("python", c(
    file.path(repo_dir, "Code", "basin_finder.py"),
    paste0("--network=",     network),
    paste0("--proteomics=",  proteomics),
    paste0("--output-dir=",  data_dir),
    paste0("--prefix=",      prefix)
  ))
  log("Step 2 done.")
} else {
  log("Step 2 skipped.")
}

# ---- step 3 -----------------------------------------------------------------
log("Step 3: Analysing basins and NMF clustering...")
step3_args <- c(
  "--vanilla",
  file.path(repo_dir, "Code", "basin_analysis.R"),
  paste0("--proteomics=",  proteomics),
  paste0("--network=",     network),
  paste0("--basins-pos=",  basins_pos),
  paste0("--basins-neg=",  basins_neg),
  paste0("--data-dir=",    data_dir),
  paste0("--figures-dir=", figures_dir),
  paste0("--prefix=",      prefix),
  paste0("--nmf-k=",       nmf_k)
)
if (!is.null(metadata)) step3_args <- c(step3_args, paste0("--metadata=", metadata))

run("Rscript", step3_args)
log("Step 3 done.")

log("Pipeline complete. Figures in ", figures_dir, "/, tables in ", data_dir, "/.")
