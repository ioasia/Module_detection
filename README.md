# Module Detection Pipeline

A pipeline for detecting protein modules from quantitative proteomics data using a basin-of-attraction approach on a protein co-abundance network.

## What it does

1. **Network reconstruction** — Computes pairwise Pearson correlations between proteins across samples and trains a logistic regression classifier (using CORUM protein complex pairs as ground truth) to assign each protein pair a probability of functional interaction. Edges above a chosen threshold form the network.

2. **Basin finding** — For each sample, protein abundances are used as a scalar function on the network. Proteins "flow" toward local minima; all proteins that converge to the same minimum form one module (basin of attraction). Run on both raw and negated abundances to capture modules defined by either high or low abundance.

3. **Basin analysis and NMF clustering** — Modules are filtered by size and statistical significance (permutation test). Per-module abundance profiles are computed across all samples and clustered using Non-negative Matrix Factorization (NMF) to identify recurrent sample archetypes.

## Dependencies

**R:** `tidyverse`, `igraph`, `Hmisc`, `ComplexHeatmap`, `RcppML`, `caret`, `ROCR`, `reticulate`, `rgexf`, `ggrepel`, `patchwork`, `circlize`, `reshape2`

**Python:** `numpy`, `pandas`, `networkx`, `tqdm`, `pickle` (stdlib)

## Input data

| File | Format |
|---|---|
| Proteomics abundance matrix | TSV, proteins × samples, first column = protein names |
| Sample metadata (optional) | TSV, samples × annotation columns, first column = sample names |
| CORUM complex pairs | TSV (provided in `Data/`) |
| Reference databases (optional) | STRING, Reactome, BioPlex etc. in `Data/` — used for network overlap figures |

## Usage

Edit the **CONFIG block** at the top of either wrapper and run from the repo root:

```bash
# Shell
bash run_pipeline.sh

# R
Rscript run_pipeline.R
```

### Controlling which steps run

Both wrappers accept `--from` and `--until` flags to run a subset of the pipeline. The stages in order are: **fdr → network → basins → analysis**.

| Flag | Effect |
|---|---|
| `--from=basins` | Skip network reconstruction; start at basin finding (requires existing network file) |
| `--from=analysis` | Skip steps 1 and 2; run only the analysis (requires existing network and basin files) |
| `--until=fdr` | Stop after writing FDR/AUC diagnostic plots; do **not** apply the threshold or write the network |
| `--until=network` | Stop after full network construction |
| `--until=basins` | Stop after basin finding; skip the analysis |

The flags can be combined, e.g. `--from=basins --until=basins` runs only the basin-finding step.

When using `--from`, the required input files must already exist in `DATA_DIR` with names matching the configured `PREFIX`:

- `--from=basins` needs `{prefix}_igraph_predictions_sig.txt`
- `--from=analysis` additionally needs `{prefix}_basins.pkl` and `{prefix}_basins_neg.pkl`

The wrapper will exit with a clear error message if any required file is missing.

### Recommended two-pass workflow

**Pass 1 — inspect the threshold:**
```bash
bash run_pipeline.sh --until=fdr
```
Examine `{prefix}_graph_quality.pdf` in `Figures/`, pick a probability threshold, and update `PRED_THRESHOLD` in the config.

**Pass 2 — run the full pipeline with the chosen threshold:**
```bash
bash run_pipeline.sh
```

### Choosing the NMF rank

On the first full run, Step 3 estimates NMF stability across ranks 2–20 and saves a cophenetic correlation plot (`{prefix}_NMF_cophenetic.pdf`). Inspect it to find the best rank, set `NMF_K` in the config, then re-run the analysis only:

```bash
bash run_pipeline.sh --from=analysis
```

## Output files

All outputs are prefixed with the `PREFIX` you set in the config.

| File | Description |
|---|---|
| `{prefix}_igraph_predictions_sig.txt` | Significant network edge list |
| `{prefix}_basins.pkl` / `{prefix}_basins_neg.pkl` | Basins of attraction per sample |
| `{prefix}_basins_combined.txt` | All filtered basins in long format |
| `{prefix}_module_median_abund.txt` | Per-module median abundance across samples |
| `{prefix}_NMF.RDS` | Cached cophenetic correlation table (re-used if present) |
| `{prefix}_nmf_scores.txt` | Module scores per NMF component |
| `{prefix}_membership.txt` | Per-sample NMF membership |
| `Figures/{prefix}_*.pdf` | All figures |

## Repository structure

```
Module_detection/
├── Code/
│   ├── classifier_functions.R       helper functions for network_reconstruction.R
│   ├── network_reconstruction.R     step 1 — build co-expression network
│   ├── basin_finder.py              step 2 — find basins of attraction
│   └── basin_analysis.R             step 3 — filter, NMF clustering, figures
├── Data/                            input and output data files
├── Figures/                         output figures
├── run_pipeline.sh                  shell wrapper
└── run_pipeline.R                   R wrapper
```
