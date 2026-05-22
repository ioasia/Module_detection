#!/usr/bin/env bash
# Full module detection pipeline.
# Run from the repo root: bash run_pipeline.sh
#
# Steps:
#   1. network_reconstruction.R  — build co-expression network
#   2. basin_finder.py (x2)      — find basins of attraction
#   3. basin_analysis.R          — filter, NMF clustering, figures
#
# Pass --skip-step1 to reuse an existing network file.

set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================
# CONFIG — edit these before running
# ============================================================

# Abundance matrix: TSV, proteins x samples, first column = protein names
PROTEOMICS="Data/my_proteomics.txt"
#'~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Protein_landscape/ALL_cell_lines/ALL_cell_line_proteomics_median_centered.txt'


# Sample metadata: TSV, samples x annotation columns, first column = sample names
# Leave empty to skip heatmap annotation
METADATA=
#"Data/ALL_cell_line_meta_data_init.txt"

# CORUM complex pairs file (used to train the classifier)
CORUM="Data/CORUM_complex_pairs.txt"

# Output directories
DATA_DIR="Data"
FIGURES_DIR="Figures"

# Short name used as prefix for all output filenames
PREFIX="my_dataset"

# Probability threshold for keeping a network edge (inspect _graph_quality.pdf first)
PRED_THRESHOLD="0.8"

# Number of NMF components for the final decomposition.
# Inspect the cophenetic plot from the first run, then set this to the optimal k.
NMF_K="13"

# ============================================================
# (nothing below here should need editing)
# ============================================================

log() { echo "[$(date '+%H:%M:%S')] $*"; }

SKIP_STEP1=false
SKIP_STEP2=false
for arg in "$@"; do
  case $arg in
    --skip-step1) SKIP_STEP1=true ;;
    --skip-step2) SKIP_STEP2=true ;;
    *) echo "Unknown argument: $arg"; exit 1 ;;
  esac
done

# ---- step 1 -----------------------------------------------------------------
if [ "$SKIP_STEP1" = false ]; then
  log "Step 1: Building protein co-expression network..."
  Rscript --vanilla "$REPO_DIR/Code/network_reconstruction.R" \
    --proteomics="$PROTEOMICS" \
    --corum="$CORUM" \
    --data-dir="$DATA_DIR" \
    --figures-dir="$FIGURES_DIR" \
    --prefix="$PREFIX" \
    --pred-threshold="$PRED_THRESHOLD"
  log "Step 1 done."
else
  log "Step 1 skipped."
fi

NETWORK="$DATA_DIR/${PREFIX}_igraph_predictions_sig.txt"

# ---- step 2 -----------------------------------------------------------------
if [ "$SKIP_STEP2" = false ]; then
  log "Step 2a: Finding basins (positive abundances)..."
  python "$REPO_DIR/Code/basin_finder.py" \
    --network="$NETWORK" \
    --proteomics="$PROTEOMICS" \
    --output-dir="$DATA_DIR" \
    --prefix="$PREFIX" \
    --no-negate

  log "Step 2b: Finding basins (negated abundances)..."
  python "$REPO_DIR/Code/basin_finder.py" \
    --network="$NETWORK" \
    --proteomics="$PROTEOMICS" \
    --output-dir="$DATA_DIR" \
    --prefix="$PREFIX"

  log "Step 2 done."
else
  log "Step 2 skipped."
fi

BASINS_POS="$DATA_DIR/${PREFIX}_basins.pkl"
BASINS_NEG="$DATA_DIR/${PREFIX}_basins_neg.pkl"

# ---- step 3 -----------------------------------------------------------------
log "Step 3: Analysing basins and NMF clustering..."
Rscript --vanilla "$REPO_DIR/Code/basin_analysis.R" \
  --proteomics="$PROTEOMICS" \
  --network="$NETWORK" \
  --basins-pos="$BASINS_POS" \
  --basins-neg="$BASINS_NEG" \
  ${METADATA:+--metadata="$METADATA"} \
  --data-dir="$DATA_DIR" \
  --figures-dir="$FIGURES_DIR" \
  --prefix="$PREFIX" \
  --nmf-k="$NMF_K"
log "Step 3 done."

log "Pipeline complete. Figures in $FIGURES_DIR/, tables in $DATA_DIR/."
