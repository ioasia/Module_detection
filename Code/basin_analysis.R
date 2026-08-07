### Basin analysis and NMF clustering
#
# Takes the basins-of-attraction identified by basin_finder.py, filters them
# by size and statistical significance, computes per-module abundance profiles,
# and clusters modules across samples using NMF.
#
# Usage (run from repo root):
#   Rscript Code/basin_analysis.R \
#     --proteomics=Data/my_proteomics.txt \
#     --network=Data/my_dataset_igraph_predictions_sig.txt \
#     --basins-pos=Data/my_dataset_basins.pkl \
#     --basins-neg=Data/my_dataset_basins_neg.pkl \
#     --metadata=Data/my_metadata.txt \
#     --nmf-k=10 \
#     --prefix=my_dataset
#
# Required arguments:
#   --proteomics    TSV, proteins x samples, first column = protein names
#   --network       TSV network edge list (Protein1 / Protein2 columns)
#   --basins-pos    pickle file of positive-abundance basins (from basin_finder.py)
#   --basins-neg    pickle file of negative-abundance basins (from basin_finder.py)
#
# Optional arguments:
#   --metadata      TSV, samples x annotation columns, first column = sample names
#                   used for heatmap annotation; omit to skip annotation
#   --data-dir      directory for output tables (default: Data)
#   --figures-dir   directory for output figures (default: Figures)
#   --prefix        prefix for all output filenames (default: dataset)
#   --nmf-k         NMF rank to use for final decomposition (default: 10)
#                   Run once first to inspect the cophenetic correlation plot,
#                   then re-run with the optimal k. The NMF estimation RDS is
#                   cached so re-runs with different k values are fast.
#
# Outputs (in data-dir):
#   {prefix}_basins_combined.txt
#   {prefix}_module_median_abund.txt
#   {prefix}_NMF.RDS                  (cached; re-used if already present)
#   {prefix}_nmf_scores.txt
#   {prefix}_membership.txt

library(reticulate)
if (!reticulate::py_module_available("pandas")) {
  reticulate::py_install("pandas")
}
library(tidyverse)
library(data.table)
library(igraph)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RcppML)

# ---- argument parsing -------------------------------------------------------
args <- (function() {
  a <- commandArgs(trailingOnly = TRUE)
  d <- list(
    "data-dir"    = "Data",
    "figures-dir" = "Figures",
    "prefix"      = "dataset",
    "nmf-k"       = "10",
    "metadata"    = ""
  )
  for (x in a) {
    m <- regmatches(x, regexec("^--([^=]+)=(.+)$", x))[[1]]
    if (length(m) == 3) d[[m[2]]] <- m[3]
  }
  d
})()

required <- c("proteomics", "basins-pos", "basins-neg")
missing  <- required[!required %in% names(args)]
if (length(missing) > 0) stop("Missing required arguments: --", paste(missing, collapse=", --"))

has_network <- !is.null(args$network) && nchar(args$network) > 0 && file.exists(args$network)

data_dir    <- args[["data-dir"]]
figures_dir <- args[["figures-dir"]]
prefix      <- args[["prefix"]]
K           <- as.integer(args[["nmf-k"]])

# ---- data -------------------------------------------------------------------
proteomics_quant <- read.delim(args$proteomics, sep = '\t', row.names = 1)
if (has_network) net <- read.delim(args$network, sep = '\t')

pd              <- import("pandas")
pickle_data_pos <- pd$read_pickle(args[["basins-pos"]])
pickle_data_neg <- pd$read_pickle(args[["basins-neg"]])

has_metadata <- nchar(args[["metadata"]]) > 0 && file.exists(args[["metadata"]])
if (has_metadata) {
  proteomics_meta <- read.delim(args[["metadata"]], sep = '\t', row.names = 1)
}

# ---- basin extraction -------------------------------------------------------
get.basins <- function(sampleIDs, pickle_dat, basin_length_thres_low,
                       basin_length_thres_high, sign_level) {
  do.call(rbind, lapply(sampleIDs, function(sampleID) {
    pickle_data_sample <- pickle_dat[[sampleID]]
    n_indices <- length(pickle_data_sample)

    all_basins <- do.call(rbind, lapply(1:n_indices, function(i) {
      basins        <- pickle_data_sample[[i]]
      basin_length  <- sapply(basins, function(b) length(b))
      tokeep        <- names(basin_length)[basin_length >= basin_length_thres_low &
                                           basin_length <= basin_length_thres_high]
      basins_filtered <- basins[tokeep]
      basins_filtered <- reshape2::melt(basins_filtered)
      basins_filtered$index <- i - 1  # 0-based, matching basin_finder.py's pass counter (index-0, index-1, ...)
      basins_filtered
    }))
    all_basins$sample <- sampleID
    colnames(all_basins)[1:2] <- c('basin', 'critical')
    all_basins$sign <- sign_level
    all_basins
  }))
}

# ---- significance filter ----------------------------------------------------
basins.thres <- function(quant_dat, thres_quant, basins_res, sign_level, iter) {
  basins_res_list <- split(basins_res, paste(basins_res$critical, basins_res$sample, sep = '_'))

  do.call(rbind, lapply(1:ncol(quant_dat), function(i) {
    print(i)
    sampleID          <- colnames(quant_dat)[i]
    basin_idx         <- grep(sampleID, names(basins_res_list))
    sample_basins_list <- basins_res_list[basin_idx]

    eff_size_basins <- do.call(rbind, lapply(names(sample_basins_list), function(k) {
      totest      <- sample_basins_list[[k]]
      basin_nodes <- totest$basin
      val         <- quant_dat[basin_nodes, sampleID]
      data.frame(mod = k, eff = abs(median(val, na.rm = TRUE)), n = length(basin_nodes))
    }))

    sample_basins_genes <- lapply(sample_basins_list, function(b) b$basin)
    gene_all <- if (sign_level == 'negative')
      rownames(quant_dat)[quant_dat[, sampleID] < 0] else
      rownames(quant_dat)[quant_dat[, sampleID] > 0]

    eff_size_basins_null <- do.call(rbind, lapply(1:iter, function(j) {
      set.seed(j)
      sample_basins_relist <- relist(sample(gene_all), skeleton = sample_basins_genes)
      do.call(rbind, lapply(names(sample_basins_relist), function(k) {
        basin_nodes <- sample_basins_relist[[k]]
        val         <- quant_dat[basin_nodes, sampleID]
        data.frame(eff = abs(median(val, na.rm = TRUE)), n = length(basin_nodes), iter = j)
      }))
    }))

    thres_group <- eff_size_basins_null %>% group_by(n) %>%
      summarize(thres = quantile(eff, thres_quant, na.rm = TRUE))
    eff_size_basins$thres <- thres_group$thres[match(eff_size_basins$n, thres_group$n)]

    eff_size_basins$pval <- sapply(1:nrow(eff_size_basins), function(j) {
      eff_tocheck <- eff_size_basins[j, ]
      eff_back    <- eff_size_basins_null[eff_size_basins_null$n == eff_tocheck$n, ]
      1 - sum(eff_tocheck$eff > eff_back$eff, na.rm = TRUE) / nrow(eff_back)
    })
    eff_size_basins$adjpval <- p.adjust(eff_size_basins$pval, method = 'BH')
    eff_size_basins
  }))
}

# ---- extract and filter basins ----------------------------------------------
all_samples_basins_neg <- get.basins(colnames(proteomics_quant), pickle_data_neg,
                                     basin_length_thres_low = 2,
                                     basin_length_thres_high = nrow(proteomics_quant),
                                     sign_level = 'negative')
all_samples_basins_neg <- all_samples_basins_neg[all_samples_basins_neg$index %in% c(0, 1), ]
all_samples_basins_neg <- all_samples_basins_neg[all_samples_basins_neg$critical != 'boundary', ]

basins_thres_neg     <- basins.thres(proteomics_quant, 0.95, all_samples_basins_neg, 'negative', 100)
basins_thres_neg_sig <- basins_thres_neg[basins_thres_neg$adjpval < 1.1, ]
all_samples_basins_neg <- all_samples_basins_neg[
  paste0(all_samples_basins_neg$critical, '_', all_samples_basins_neg$sample) %in% basins_thres_neg_sig$mod, ]

all_samples_basins_pos <- get.basins(colnames(proteomics_quant), pickle_data_pos,
                                     basin_length_thres_low = 2,
                                     basin_length_thres_high = nrow(proteomics_quant),
                                     sign_level = 'positive')
all_samples_basins_pos <- all_samples_basins_pos[all_samples_basins_pos$index %in% c(0, 1), ]
all_samples_basins_pos <- all_samples_basins_pos[all_samples_basins_pos$critical != 'boundary', ]

basins_thres_pos     <- basins.thres(proteomics_quant, 0.95, all_samples_basins_pos, 'positive', 100)
basins_thres_pos_sig <- basins_thres_pos[basins_thres_pos$adjpval < 1.1, ]
all_samples_basins_pos <- all_samples_basins_pos[
  paste0(all_samples_basins_pos$critical, '_', all_samples_basins_pos$sample) %in% basins_thres_pos_sig$mod, ]

all_samples_basins <- rbind(all_samples_basins_neg, all_samples_basins_pos)
write.table(all_samples_basins,
            file.path(data_dir, paste0(prefix, '_basins_combined.txt')), sep = '\t')

# ---- basin size and count plots ---------------------------------------------
basins_per_sample <- all_samples_basins %>%
  group_by(sample, index, sign) %>% summarise(n_unique_critical = n_distinct(critical))
critical_nodes_sum <- basins_per_sample %>% group_by(sample) %>% summarise(s = sum(n_unique_critical))

basins_per_sample$sample <- factor(basins_per_sample$sample,
                                   levels = critical_nodes_sum$sample[order(-critical_nodes_sum$s)])
basins_per_sample$col    <- factor(paste0(basins_per_sample$index, '_', basins_per_sample$sign),
                                   levels = c("0_positive", "1_positive", "0_negative", "1_negative"))
basins_per_sample$n_unique_critical[basins_per_sample$sign == 'negative'] <-
  -basins_per_sample$n_unique_critical[basins_per_sample$sign == 'negative']

p_number <- ggplot(basins_per_sample, aes(x = sample, y = n_unique_critical, fill = as.factor(col))) +
  geom_bar(stat = 'identity', col = 'black') +
  scale_y_continuous(expand = c(0, 0, 0.1, 0.1), labels = function(i) abs(i)) +
  scale_fill_manual(name = 'Index_abundance', values = c("#FFAAAA", "#FF0000", "#AAAAFF", "#0000FF")) +
  annotate(geom = 'text', x = Inf, y = round(median(critical_nodes_sum$s), 0), col = 'red',
           label = paste0('n = ', round(median(critical_nodes_sum$s), 0)), hjust = 1, vjust = 1.2) +
  labs(x = 'sample ID', y = '# modules') +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

pdf(file.path(figures_dir, paste0(prefix, '_modules_per_sample.pdf')), width = 10, height = 5)
plot(p_number); dev.off()

size_per_node <- all_samples_basins %>%
  group_by(sample, critical, sign) %>% summarise(n_unique_basin = n_distinct(basin))
median_size   <- size_per_node %>% group_by(sample) %>% summarise(m = median(n_unique_basin))
size_per_node$sample <- factor(size_per_node$sample,
                               levels = unique(median_size$sample[order(-median_size$m)]))

p_size <- ggplot(size_per_node, aes(x = sample, y = n_unique_basin)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_y_log10() +
  annotate(geom = 'text', x = Inf, y = Inf, col = 'red',
           label = paste0('median = ', round(median(size_per_node$n_unique_basin), 0)), hjust = 1, vjust = 1.2) +
  labs(x = 'sample ID', y = 'Module size') +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

pdf(file.path(figures_dir, paste0(prefix, '_module_size.pdf')), width = 10, height = 5)
plot(p_size); dev.off()

# use macOS native graphics device for any rasterized heatmap rendering
options(bitmapType = 'quartz')

# ---- correlation heatmap of module profiles ---------------------------------
all_samples_basins_list <- split(all_samples_basins,
                                 paste0(all_samples_basins$critical, '_', all_samples_basins$sample))
basin_sets <- lapply(all_samples_basins_list, function(df) unique(df$basin))

basins_profiles <- do.call(rbind, lapply(basin_sets, function(i) {
  as.numeric(apply(proteomics_quant[i, ], 2, median, na.rm = TRUE))
}))
colnames(basins_profiles) <- colnames(proteomics_quant)

basins_profiles_cor <- cor(t(basins_profiles), method = 'pearson', use = 'pairwise.complete.obs')
row_clust           <- hclust(as.dist(1 - basins_profiles_cor), method = 'ward.D2')

ht_module_cor <- Heatmap(basins_profiles_cor,
                         col = colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red')),
                         show_column_names = FALSE, show_row_names = FALSE,
                         cluster_columns = row_clust, cluster_rows = row_clust,
                         use_raster = TRUE,
                         raster_resize_mat = TRUE,
                         border = TRUE,
                         column_title = paste0(prefix, ' (n_modules = ', nrow(basins_profiles_cor), ')'),
                         heatmap_legend_param = list(title = 'Pearson cor.'))

pdf(file.path(figures_dir, paste0(prefix, '_module_cor_heatmap.pdf')), width = 10, height = 8)
draw(ht_module_cor, padding = unit(c(2, 20, 2, 2), "mm")); dev.off()

# ---- module abundance profiles ----------------------------------------------
merged_basins_list <- all_samples_basins_list
unique_basins_list <- lapply(merged_basins_list, function(i) {
  res <- i[!duplicated(i$basin), ]
  res$basin
})

unique_basins_profiles <- do.call(rbind, lapply(unique_basins_list, function(i) {
  as.numeric(apply(proteomics_quant[i, ], 2, median, na.rm = TRUE))
}))
colnames(unique_basins_profiles) <- colnames(proteomics_quant)

write.table(unique_basins_profiles,
            file.path(data_dir, paste0(prefix, '_module_median_abund.txt')), sep = '\t')

# ---- NMF clustering ---------------------------------------------------------
nmf_input <- as.matrix(unique_basins_profiles)
nmf_input <- nmf_input[complete.cases(nmf_input), ]

posneg <- function(x) rbind(pmax(x, 0), pmax(-x, 0))
nmf_input_pos <- posneg(nmf_input)

nmf_rds <- file.path(data_dir, paste0(prefix, '_NMF.RDS'))
if (file.exists(nmf_rds)) {
  message("Loading cached NMF results from ", nmf_rds)
  coph_cor_dt <- readRDS(nmf_rds)
} else {
  message("Running NMF rank estimation (ranks 2-20, 10 runs each)...")
  nmf_ranks <- 2:20
  coph_cor_dt <- data.frame(
    k        = nmf_ranks,
    coph_cor = sapply(nmf_ranks, function(k) {
      message("  rank k = ", k)
      n_samp    <- ncol(nmf_input_pos)
      consensus <- matrix(0, n_samp, n_samp)
      for (i in seq_len(10)) {
        fit_i <- RcppML::nmf(nmf_input_pos, k = k, seed = i, verbose = FALSE)
        cls   <- apply(fit_i$h, 2, which.max)
        consensus <- consensus + outer(cls, cls, "==") * 1.0
      }
      consensus <- consensus / 10
      d <- as.dist(1 - consensus)
      cor(d, cophenetic(hclust(d, method = "average")))
    })
  )
  saveRDS(coph_cor_dt, file = nmf_rds)
}

p_coph <- ggplot(coph_cor_dt, aes(x = k, y = coph_cor)) +
  geom_line() +
  geom_point(pch = 21, size = 3, fill = 'white') +
  geom_point(data = coph_cor_dt[coph_cor_dt$k == K, ], pch = 21, size = 6, color = 'red', stroke = 2) +
  labs(x = 'k', y = 'cophenetic correlation coefficient',
       title = paste0('NMF rank estimation (selected k = ', K, ')')) +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5))

pdf(file.path(figures_dir, paste0(prefix, '_NMF_cophenetic.pdf')), width = 6, height = 5)
plot(p_coph); dev.off()

# ---- NMF decomposition at chosen K ------------------------------------------
message("Fitting NMF at k = ", K, "...")
fit_final <- RcppML::nmf(nmf_input_pos, k = K, seed = 1234, verbose = FALSE)
H          <- as.data.frame(t(fit_final$h))
rownames(H) <- colnames(nmf_input_pos)
h_membership <- as.data.frame(t(apply(H, 1, function(i) i / sum(i))))
colnames(h_membership) <- paste0('NMF', 1:K)

nmf_clusters <- apply(h_membership, 1, which.max)
h_membership$class   <- factor(apply(h_membership, 1, function(i) colnames(h_membership)[which.max(i)]),
                                levels = paste0('NMF', 1:K))
h_membership$max_val <- apply(h_membership[, paste0('NMF', 1:K)], 1, max)

sample_order <- rownames(h_membership)[order(h_membership$class, -h_membership$max_val)]

write.table(h_membership[sample_order, ],
            file.path(data_dir, paste0(prefix, '_membership.txt')), sep = '\t')

# ---- NMF module scores ------------------------------------------------------
top_scores <- do.call(rbind, lapply(1:nrow(unique_basins_profiles), function(i) {
  val    <- as.numeric(scale(t(unique_basins_profiles))[, i])
  dt     <- data.frame(scale(H), val = val)
  sapply(1:K, function(j) {
    res <- lm(val ~ dt[[j]])
    summary(res)$coefficients[, 't value'][-1]
  })
}))
rownames(top_scores) <- rownames(unique_basins_profiles)
colnames(top_scores) <- paste0('NMF', 1:K)

write.table(top_scores,
            file.path(data_dir, paste0(prefix, '_nmf_scores.txt')), sep = '\t')

top_score_mods <- apply(top_scores, 2, function(k) {
  k_idx       <- which(abs(k) > 2)
  modules_sub <- rownames(top_scores)[k_idx]
  modules_sub[order(k[k_idx], decreasing = TRUE)]
})

top_score_mods_order <- reshape2::melt(top_score_mods)
top_score_mods_order$L1 <- factor(top_score_mods_order$L1, levels = paste0('NMF', 1:K))

# ---- combined heatmap -------------------------------------------------------
ht_memb <- Heatmap(t(h_membership[sample_order, paste0('NMF', 1:K)]),
                   col = colorRamp2(c(0, 0.5, 1), c('blue', 'white', 'red')),
                   column_split = h_membership[sample_order, 'class'],
                   column_title = NULL, show_column_names = FALSE,
                   cluster_columns = FALSE, cluster_rows = FALSE,
                   use_raster = FALSE,
                   row_names_side = 'left', border = TRUE,
                   heatmap_legend_param = list(title = 'membership'),
                   height = unit(4, "cm"))

ht_prof <- Heatmap(unique_basins_profiles[top_score_mods_order$value, sample_order],
                   row_split   = top_score_mods_order$Var2,
                   column_split = h_membership[sample_order, 'class'],
                   column_title = NULL, row_title = NULL,
                   show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
                   show_row_names = FALSE, row_gap = unit(0, "mm"), column_gap = unit(3, "mm"),
                   use_raster = FALSE,
                   border = TRUE, heatmap_legend_param = list(title = 'log2 ratio'),
                   height = unit(22, "cm"))

if (has_metadata) {
  meta_ann <- proteomics_meta[sample_order, , drop = FALSE]
  ha_meta  <- HeatmapAnnotation(df = meta_ann, show_annotation_name = TRUE,
                                na_col = 'white', border = TRUE, annotation_name_side = 'left')
  ht_all <- ha_meta %v% ht_memb %v% ht_prof
} else {
  ht_all <- ht_memb %v% ht_prof
}

pdf(file.path(figures_dir, paste0(prefix, '_NMF_heatmap.pdf')), width = 12, height = 12)
draw(ht_all, padding = unit(c(2, 20, 2, 2), "mm")); dev.off()

# ---- protein frequency vs network degree (only if network file was provided) -
if (!has_network) {
  message("--network not provided; skipping degree-vs-modules plot.")
} else {

freq_proteins <- as.data.frame(table(unlist(unique_basins_list)))
freq_proteins$Var1 <- factor(freq_proteins$Var1,
                             levels = freq_proteins$Var1[order(freq_proteins$Freq, decreasing = TRUE)])
freq_proteins <- freq_proteins[order(freq_proteins$Var1), ]

net_graph   <- graph_from_data_frame(net[, 1:2])
node_degree <- igraph::degree(net_graph)

freq_proteins$node_degree <- node_degree[as.character(freq_proteins$Var1)]
freq_proteins$score       <- freq_proteins$Freq * freq_proteins$node_degree
freq_proteins <- freq_proteins[order(freq_proteins$score, decreasing = TRUE), ]

p_degree_modules <- ggplot(freq_proteins, aes(x = node_degree, y = Freq)) +
  geom_point() +
  geom_point(data = freq_proteins[1:50, ], color = 'red') +
  geom_text_repel(data = freq_proteins[1:50, ], mapping = aes(label = Var1),
                  force = 20, min.segment.length = 0.1, size = 2, color = 'red', fontface = 'bold') +
  labs(x = 'node degree', y = '# modules') +
  theme_classic()

ggsave(file.path(figures_dir, paste0(prefix, '_degree_vs_modules.pdf')),
       plot = p_degree_modules, width = 10, height = 7)

} # end has_network

message("basin_analysis.R done. Results written to ", data_dir)
