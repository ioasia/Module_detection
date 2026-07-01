# Full module detection pipeline.
# Run from the repo root: Rscript run_pipeline.R
#                    or:  source("run_pipeline.R")  (from an R session)
#
# Steps:
#   1. network_reconstruction.R  — build co-expression network
#   2. basin_finder.py (x2)      — find basins of attraction
#   3. basin_analysis.R          — filter, NMF clustering, figures
#
# Pass --from={basins|analysis} and/or --until={fdr|network|basins} to run
# a subset of the pipeline.

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
# Inspect _NMF_cophenetic.pdf from the first run, then set this to the optimal k.
nmf_k         <- "10"

# ============================================================
# (nothing below here should need editing)
# ============================================================

log <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

run <- function(cmd, args = character()) {
  status <- system2(cmd, args = args)
  if (status != 0) stop(sprintf("Command failed (exit %d): %s %s", status, cmd, paste(args, collapse = " ")))
}

# ---- parse --from / --until -------------------------------------------------
# Pipeline positions:  fdr(1) -> network(2) -> basins(3) -> analysis(4)
#
#   --from=basins    skip network reconstruction, start at basin finding
#   --from=analysis  skip steps 1 and 2, run only the analysis
#   --until=fdr      stop after FDR diagnostics (before applying threshold)
#   --until=network  stop after full network construction
#   --until=basins   stop after basin finding

cli_args  <- commandArgs(trailingOnly = TRUE)
from_arg  <- NA_character_
until_arg <- NA_character_
for (arg in cli_args) {
  if      (grepl("^--from=",  arg)) from_arg  <- sub("^--from=",  "", arg)
  else if (grepl("^--until=", arg)) until_arg <- sub("^--until=", "", arg)
  else stop("Unknown argument: ", arg,
            "\nValid: --from={basins|analysis}  --until={fdr|network|basins}")
}
if (!is.na(from_arg)  && !from_arg  %in% c("basins", "analysis"))
  stop("--from must be 'basins' or 'analysis'")
if (!is.na(until_arg) && !until_arg %in% c("fdr", "network", "basins"))
  stop("--until must be 'fdr', 'network', or 'basins'")

from_pos  <- if (is.na(from_arg))  1L else switch(from_arg,  basins = 3L, analysis = 4L)
until_pos <- if (is.na(until_arg)) 4L else switch(until_arg, fdr = 1L, network = 2L, basins = 3L)

if (from_pos > until_pos)
  stop("--from=", from_arg, " is after --until=", until_arg, " in the pipeline")

network    <- file.path(data_dir, paste0(prefix, "_igraph_predictions_sig.txt"))
basins_pos <- file.path(data_dir, paste0(prefix, "_basins.pkl"))
basins_neg <- file.path(data_dir, paste0(prefix, "_basins_neg.pkl"))

if (from_pos == 3L && !file.exists(network))
  stop("Network file not found: ", network,
       "\nRun step 1 first, or check your prefix / data-dir settings.")
if (from_pos >= 4L) {
  if (!file.exists(basins_pos))
    stop("Basins file not found: ", basins_pos,
         "\nRun step 2 first, or check your prefix / data-dir settings.")
  if (!file.exists(basins_neg))
    stop("Negated basins file not found: ", basins_neg,
         "\nRun step 2 first, or check your prefix / data-dir settings.")
}

run_step1  <- from_pos <= 2L
early_stop <- until_pos <= 1L
run_step2  <- from_pos <= 3L && until_pos >= 3L
run_step3  <- until_pos >= 4L

# ---- step 1 -----------------------------------------------------------------
if (run_step1) {
  log("Step 1: Building protein co-expression network...")
  step1_args <- c(
    "--vanilla",
    file.path(repo_dir, "Code", "network_reconstruction.R"),
    paste0("--proteomics=",     proteomics),
    paste0("--corum=",          corum),
    paste0("--data-dir=",       data_dir),
    paste0("--figures-dir=",    figures_dir),
    paste0("--prefix=",         prefix),
    paste0("--pred-threshold=", pred_threshold)
  )
  if (early_stop) step1_args <- c(step1_args, "--early-stop")
  run("Rscript", step1_args)
  log("Step 1 done.")
} else {
  log("Step 1 skipped.")
}

# ---- step 2 -----------------------------------------------------------------
if (run_step2) {
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
if (run_step3) {
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
} else {
  log("Step 3 skipped.")
}

log("Pipeline complete. Figures in ", figures_dir, "/, tables in ", data_dir, "/.")
