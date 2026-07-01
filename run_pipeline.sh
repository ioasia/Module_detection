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
# Inspect _NMF_cophenetic.pdf from the first run, then set this to the optimal k.
NMF_K="13"

# ============================================================
# (nothing below here should need editing)
# ============================================================

log() { echo "[$(date '+%H:%M:%S')] $*"; }

NETWORK="$DATA_DIR/${PREFIX}_igraph_predictions_sig.txt"
BASINS_POS="$DATA_DIR/${PREFIX}_basins.pkl"
BASINS_NEG="$DATA_DIR/${PREFIX}_basins_neg.pkl"

# ---- parse --from / --until -------------------------------------------------
# Pipeline positions:  fdr(1) -> network(2) -> basins(3) -> analysis(4)
#
#   --from=basins    skip network reconstruction, start at basin finding
#   --from=analysis  skip steps 1 and 2, run only the analysis
#   --until=fdr      stop after FDR diagnostics (before applying threshold)
#   --until=network  stop after full network construction
#   --until=basins   stop after basin finding

FROM_POS=1
UNTIL_POS=4
for arg in "$@"; do
  case $arg in
    --from=basins)    FROM_POS=3 ;;
    --from=analysis)  FROM_POS=4 ;;
    --from=*)         echo "Unknown --from value: ${arg#--from=}"; exit 1 ;;
    --until=fdr)      UNTIL_POS=1 ;;
    --until=network)  UNTIL_POS=2 ;;
    --until=basins)   UNTIL_POS=3 ;;
    --until=*)        echo "Unknown --until value: ${arg#--until=}"; exit 1 ;;
    *)                echo "Unknown argument: $arg"; exit 1 ;;
  esac
done

if [ "$FROM_POS" -gt "$UNTIL_POS" ]; then
  echo "Error: --from is after --until in the pipeline."; exit 1
fi

if [ "$FROM_POS" -eq 3 ] && [ ! -f "$NETWORK" ]; then
  echo "Error: network file not found: $NETWORK"
  echo "Run step 1 first, or check your PREFIX / DATA_DIR settings."; exit 1
fi
if [ "$FROM_POS" -ge 4 ]; then
  if [ ! -f "$BASINS_POS" ]; then
    echo "Error: basins file not found: $BASINS_POS"
    echo "Run step 2 first, or check your PREFIX / DATA_DIR settings."; exit 1
  fi
  if [ ! -f "$BASINS_NEG" ]; then
    echo "Error: negated basins file not found: $BASINS_NEG"
    echo "Run step 2 first, or check your PREFIX / DATA_DIR settings."; exit 1
  fi
fi

if [ "$FROM_POS" -le 2 ]; then RUN_STEP1=true; else RUN_STEP1=false; fi
if [ "$UNTIL_POS" -le 1 ]; then EARLY_STOP=true; else EARLY_STOP=false; fi
if [ "$FROM_POS" -le 3 ] && [ "$UNTIL_POS" -ge 3 ]; then RUN_STEP2=true; else RUN_STEP2=false; fi
if [ "$UNTIL_POS" -ge 4 ]; then RUN_STEP3=true; else RUN_STEP3=false; fi

# ---- step 1 -----------------------------------------------------------------
if [ "$RUN_STEP1" = true ]; then
  log "Step 1: Building protein co-expression network..."
  STEP1_EXTRA=""
  if [ "$EARLY_STOP" = true ]; then STEP1_EXTRA="--early-stop"; fi
  Rscript --vanilla "$REPO_DIR/Code/network_reconstruction.R" \
    --proteomics="$PROTEOMICS" \
    --corum="$CORUM" \
    --data-dir="$DATA_DIR" \
    --figures-dir="$FIGURES_DIR" \
    --prefix="$PREFIX" \
    --pred-threshold="$PRED_THRESHOLD" \
    $STEP1_EXTRA
  log "Step 1 done."
else
  log "Step 1 skipped."
fi

# ---- step 2 -----------------------------------------------------------------
if [ "$RUN_STEP2" = true ]; then
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

# ---- step 3 -----------------------------------------------------------------
if [ "$RUN_STEP3" = true ]; then
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
else
  log "Step 3 skipped."
fi

log "Pipeline complete. Figures in $FIGURES_DIR/, tables in $DATA_DIR/."
