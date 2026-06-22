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

Both wrappers support `--skip-step1` and `--skip-step2` to reuse previously computed network or basin files (useful when re-running the NMF step with a different rank `k`).

### Choosing the NMF rank

On the first run, Step 3 estimates stability across ranks 2–20 and saves a cophenetic correlation plot (`{prefix}_NMF_cophenetic.pdf`). Inspect this plot to identify the best rank, set `NMF_K` in the config to that value, then re-run with `--skip-step1 --skip-step2` to skip the slow network and basin steps and regenerate the final figures quickly.

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
