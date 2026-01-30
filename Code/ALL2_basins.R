### ALL cell line basins


## Libraries
library(reticulate)
reticulate::py_install("pandas")
library(tidyverse)
library(data.table)


## Libraries
library(openxlsx)
library(igraph)
library(rgexf)
library(ggrepel)
library(patchwork)
library(umap)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(fgsea)
library(NMF)
library(GSA)
library(clusterProfiler)
library(tidyverse)
library(igraph)
library(ggraph)
library(ggrepel)
library(rgexf)
library(patchwork)


path_file <- '~/Downloads/ALL2_cell_line_proj/data/'
path_save <- '~/Downloads/ALL2_cell_line_proj/code/figures/'

# Data
net <- read.delim(paste0(path_file, 'ALL_cell_line_merged_igraph_predictions_sig.txt'), sep = '\t')
proteomics_quant <- read.delim(paste0(path_file, 'ALL_cell_line_proteomics_median_centered.txt'), sep = '\t', row.names = 1)
proteomics_meta <- read.delim(paste0(path_file, 'ALL_cell_line_meta_data_init.txt'),  sep = '\t')

# Basin of attraction
pd <- import("pandas")
pickle_data_neg <- pd$read_pickle(paste0(path_file, "ALL2_sign_change_condition.pkl"))
pickle_data_pos <- pd$read_pickle(paste0(path_file, "ALL2_sign_change_condition_negated.pkl"))

get.basins <- function(sampleIDs, pickle_dat, basin_length_thres_low, basin_length_thres_high, sign_level) {
  
  res <- do.call(rbind, lapply(sampleIDs, function(sampleID) {
    
    # sampleID <- "FORALL_158"
    
    pickle_data_sample <- pickle_dat[[sampleID]]
    
    # lapply(pickle_data_sample, function(i) grep('BCL2', i))
    n_indices <- length(pickle_data_sample)
    
    ## Focus on 0s and 1s
    all_basins <- do.call(rbind, lapply(c(1:n_indices) , function(i) {
      # print(i)
      # i = 1
      basins <- pickle_data_sample[[i]]
      basin_length <-  sapply(basins, function(i) length(i))
      
      # hist(log2(basin_length))
      
      ## Filters
      # tokeep <- basin_length >= 10 & basin_length <= 100
      tokeep <- basin_length >= basin_length_thres_low & basin_length <= basin_length_thres_high
      # tokeep <- basin_length >= 0
      tokeep <- names(basin_length)[tokeep]
      
      basins_filtered <- basins[tokeep]
      # basins_filtered <- basins['NUDT9']
      
      basins_filtered <- reshape2::melt(basins_filtered)
      basins_filtered$index <- i
      return(basins_filtered)
    }))
    all_basins$sample <- sampleID
    colnames(all_basins)[1:2] <- c('basin', 'critical')
    
  
    all_basins$sign <- sign_level
    return(all_basins)
    
  }))
}

basins.thres <- function(quant_dat, thres_quant, basins_res, sign_level, iter,n_cores) {
  
  basins_res_list <-  split(basins_res, paste(basins_res$critical, basins_res$sample, sep = '_')) 
  
  res <- do.call(rbind, lapply(1:ncol(quant_dat), function(i) {
    
    # i = 1
    print(i)
    
    sampleID <- colnames(quant_dat)[i]
    basin_idx <- grep(sampleID, names(basins_res_list))
    sample_basins_list <- basins_res_list[basin_idx]
    
    eff_size_basins <- do.call(rbind, lapply(names(sample_basins_list), function(k) {
      totest <- sample_basins_list[[k]]
      critical_node <- totest$critical[1]
      basin_nodes <- totest$basin
      val <- quant_dat[basin_nodes, sampleID]
      eff_size <- abs(median(val, na.rm = TRUE))
      data.frame(mod = k, eff = eff_size, n = length(basin_nodes))
    }))
    
    # plot(density(eff_size_basins$eff))
    # plot(log2(eff_size_basins$n), eff_size_basins$eff)
    
    sample_basins_genes <- lapply(sample_basins_list, function(i) i$basin)
    
    # gene_all <- unlist(sample_basins_genes)
    if(sign_level == 'negative') {
      
      gene_all <- rownames(quant_dat)[quant_dat[, sampleID] < 0]
      
    } else {
      gene_all <- rownames(quant_dat)[quant_dat[, sampleID] > 0]
    }
    
    
    eff_size_basins_null <- do.call(rbind, lapply(1:iter, function(j) {
      
      set.seed(j)
      sample_basins_relist <- relist(sample(gene_all), skeleton = sample_basins_genes)
      
      res <- do.call(rbind, lapply(names(sample_basins_relist), function(k) {
        basin_nodes <- sample_basins_relist[[k]]
        val <- quant_dat[basin_nodes, sampleID]
        eff_size <- abs(median(val, na.rm = TRUE))
        data.frame(eff = eff_size, n = length(basin_nodes), iter =j)
      }))
    }))
    
    
    thres_group <- eff_size_basins_null %>% group_by(n) %>%
      summarize(thres = quantile(eff, thres_quant, na.rm = TRUE))
    eff_size_basins$thres <- thres_group$thres[match(eff_size_basins$n, thres_group$n)]

    
    # P-value threshold
    eff_size_basins$pval <- sapply(1:nrow(eff_size_basins), function(j) {
      
      eff_tocheck <- eff_size_basins[j, ]
      eff_back <- eff_size_basins_null[eff_size_basins_null$n == eff_tocheck$n, ]
      1 - sum(eff_tocheck$eff > eff_back$eff, na.rm = TRUE)/nrow(eff_back)
      
    })
    eff_size_basins$adjpval <- p.adjust(eff_size_basins$pval, method = 'BH')
    
    eff_size_basins
  }))
  
  return(res)
}


all_samples_basins_neg <- get.basins(sampleIDs = colnames(proteomics_quant), pickle_dat = pickle_data_neg, basin_length_thres_low = 3, basin_length_thres_high = nrow(proteomics_quant), sign_level = 'negative')
all_samples_basins_neg <- all_samples_basins_neg[all_samples_basins_neg$index %in% c(1,2), ]
basins_thres_neg <- basins.thres(quant_dat = proteomics_quant, thres_quant = 0.95, basins_res = all_samples_basins_neg, sign_level = 'negative', iter = 100)

thres_pval <- 0.1
basins_thres_neg_sig <- basins_thres_neg[basins_thres_neg$adjpval < thres_pval, ]
all_samples_basins_neg <- all_samples_basins_neg[paste0(all_samples_basins_neg$critical, '_',  all_samples_basins_neg$sample) %in% basins_thres_neg_sig$mod, ]


all_samples_basins_pos <- get.basins(sampleIDs = colnames(proteomics_quant), pickle_dat = pickle_data_pos, basin_length_thres = 3,  basin_length_thres_high = nrow(proteomics_quant), sign_level = 'positive')
all_samples_basins_pos <- all_samples_basins_pos[all_samples_basins_pos$index %in% c(1,2), ]
basins_thres_pos <- basins.thres(quant_dat = proteomics_quant, thres_quant = 0.95, basins_res = all_samples_basins_pos, sign_level = 'positive', iter = 100)

basins_thres_pos_sig <- basins_thres_pos[basins_thres_pos$adjpval < thres_pval, ]
all_samples_basins_pos <- all_samples_basins_pos[paste0(all_samples_basins_pos$critical, '_',  all_samples_basins_pos$sample) %in% basins_thres_pos_sig$mod, ]


all_samples_basins <- rbind(all_samples_basins_neg, all_samples_basins_pos)

write.table(all_samples_basins, paste0(path_file, 'ALL_cell_line_basins_comb_non_outliers.txt'), sep = '\t')
################################################################################################################################
### Merge 
all_samples_basins <- read.delim(  paste0(path_file, 'ALL_cell_line_basins_comb_non_outliers.txt'), sep = '\t')


# Basins per sample
# Number of critical nodes
basins_per_sample <- all_samples_basins %>% group_by(sample, index, sign) %>%
  summarise(n_unique_critical = n_distinct(critical))

critical_nodes_sum <- basins_per_sample %>% group_by(sample) %>%
  summarise(s = sum(n_unique_critical))

basins_per_sample$sample <- factor(basins_per_sample$sample,
                                   levels = critical_nodes_sum$sample[order(-critical_nodes_sum$s)])
basins_per_sample$col <- paste0(basins_per_sample$index-1, '_', basins_per_sample$sign)

basins_per_sample$col <- factor(basins_per_sample$col, levels = c("0_positive", "1_positive", "0_negative", "1_negative"))


median_basins <- round(median(critical_nodes_sum$s), 0)
basins_per_sample$n_unique_critical[basins_per_sample$sign == 'negative'] <- -basins_per_sample$n_unique_critical[basins_per_sample$sign == 'negative']
p_number <- ggplot(basins_per_sample, aes(x = sample, y = n_unique_critical, fill = as.factor(col))) + 
  geom_bar(stat = 'identity', col = 'black') + 
  # geom_hline(yintercept = median_basins, col = 'red', lty = 1) + 
  scale_y_continuous(expand = c(0, 0, 0.1, 0.1), labels = function(i) abs(i)) + 
  scale_fill_manual(name = 'Index_abundance', values = c( "#FFAAAA", "#FF0000",  "#AAAAFF", "#0000FF")) + 
  annotate(geom = 'text', x = Inf, y = median_basins, col = 'red', label = paste0('n = ', median_basins), hjust = 1, vjust = 1.2) +
  labs(x = 'sample ID', y = '# modules') + 
  theme_classic() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

pdf(paste0(path_save, 'ALL_cell_line_critical_nodes_per_sample_non_outliers.pdf'), width = 10, height = 5)
plot(p_number)
dev.off()


# Size
size_per_node <- all_samples_basins %>% group_by(sample, critical, sign) %>%
  summarise(n_unique_basin = n_distinct(basin))

median_size <- size_per_node %>% group_by(sample) %>%
  summarise(m = median(n_unique_basin))

size_per_node$sample <- factor(size_per_node$sample, levels = unique(median_size$sample[order(-median_size$m)]))
total_median <- round(median(size_per_node$n_unique_basin), 0)

p_size <- ggplot(size_per_node, aes(x = sample, y = n_unique_basin), col = sign) + 
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) + 
  scale_y_log10() + 
  # geom_hline(yintercept = total_median, col = 'red') + 
  annotate(geom = 'text', x = Inf, y = Inf, col = 'red', label = paste0('median = ', total_median), hjust = 1, vjust = 1.2) +
  labs(x = 'sample ID', y = 'Module size') + 
  theme_classic() + 
  theme(
    # axis.text.x = element_text(hjust = 1, angle = 45),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

pdf(paste0(path_save, 'ALL_cell_line_basin_size_non_outliers', '.pdf'), width = 10, height = 5)
plot(p_size) 
dev.off()
################################################################################################################################
## Jaccard index
all_samples_basins_list <- split(all_samples_basins, paste0(all_samples_basins$critical, '_', all_samples_basins$sample))
# all_samples_basins_list <- split(all_samples_basins, paste0(all_samples_basins$critical))
length(all_samples_basins_list)


# Step 1: Precompute the unique basins for each sample
basin_sets <- lapply(all_samples_basins_list, function(df) unique(df$basin))
sample_names <- names(basin_sets)
n <- length(sample_names)

# Step 2: Initialize an empty matrix to store Jaccard similarities
all_samples_overlap <- matrix(NA, nrow = n, ncol = n,
                              dimnames = list(sample_names, sample_names))

# Step 3: Compute pairwise Jaccard similarities
for (i in seq_len(n)) {
  for (j in seq(i, n)) {  # Only compute upper triangle (symmetry)
    set_i <- basin_sets[[i]]
    set_j <- basin_sets[[j]]
    
    intersection_size <- length(intersect(set_i, set_j))
    union_size <- length(union(set_i, set_j))
    
    jaccard <- if (union_size == 0) 0 else intersection_size / union_size
    all_samples_overlap[i, j] <- jaccard
    all_samples_overlap[j, i] <- jaccard  # Fill symmetric entry
  }
}
write.table(all_samples_overlap, paste0(path_file, 'ALL_cell_line_basin_overlap_non_outliers.txt'), sep = '\t')

ht_overlap <- Heatmap(all_samples_overlap, col = colorRamp2(breaks = c(0, 1), colors = c('white', 'red')),
                      show_column_names = FALSE,
                      show_row_names = FALSE,
                      cluster_columns = TRUE,
                      cluster_rows = TRUE,
                      clustering_distance_rows = 'euclidean',
                      clustering_method_rows = 'ward.D2',
                      clustering_distance_columns = 'euclidean',
                      clustering_method_columns = 'ward.D2',
                      show_row_dend = TRUE,
                      column_gap = unit(3, "mm"),
                      border = TRUE, 
                      column_title = paste0('ALL cell lines cohort (n_modules = ', nrow(all_samples_overlap), ')'),
                      heatmap_legend_param = list(title = 'Jaccard index'),
                      column_title_side = "top",
                      column_title_gp = gpar(cex = 2, fontface = "bold"),
                      # row_title_rot = 0,
                      column_names_side = "top",
                      row_dend_width = unit(2, "cm"),
                      column_dend_height = unit(3, "cm"),
                      row_names_gp = gpar(fontsize = 12, col = 'black'),
                      column_names_gp = gpar(fontsize = 12, col = 'black'))

pdf(file = paste0(path_save, 'ALL_cell_line_heatmap_critical_node_comb_overlap_non_outliers.pdf'), width = 10, 
    height = 8)
draw(ht_overlap, padding = unit(c(2, 20, 2, 2), "mm"))
dev.off()
##################################################################################################################################################################
all_samples_overlap <- read.delim(paste0(path_file, 'ALL_cell_line_basin_overlap_non_outliers.txt'), sep = '\t', check.names = FALSE)

## Pairwise 
merge.basins <- function(basins_list, overlap_matrix, jac_thres, level = 0) {
  
  
  repeat {
    
    # Get edges above threshold
    sig_pairs <- which(overlap_matrix >= jac_thres, arr.ind = TRUE)
    sig_pairs <- sig_pairs[sig_pairs[, 1] < sig_pairs[, 2], , drop = FALSE]  # avoid duplicates
    
    if (nrow(sig_pairs) == 0) break
    
    level <- level + 1
    cat("Iteration:", level, "- merging basins\n")
    
    edges_df <- data.frame(
      from = rownames(overlap_matrix)[sig_pairs[, 1]],
      to = rownames(overlap_matrix)[sig_pairs[, 2]]
    )
    
    ## Add size 
    # edges_df$basin_size1 <- sapply(edges_df$from, function(i) nrow(basins_list[[i]]))
    # edges_df$basin_size2 <- sapply(edges_df$to, function(i) nrow(basins_list[[i]]))
    edges_df$jacc <- sapply(1:nrow(edges_df), function(k) overlap_matrix[ edges_df$from[k], edges_df$to[k] ])
    edges_df <- edges_df[order(edges_df$jacc, decreasing = TRUE), ]
    
    edge_unique <- unique(c(edges_df$from, edges_df$to))
    
    edge_unique_size <- sapply(edge_unique, function(i) nrow(basins_list[[i]]))
    idx_size <- order(edge_unique_size, decreasing = TRUE)
    edge_unique_iter <- edge_unique[idx_size]
    
    merged_df_all <- list()
    toremove <- c()
    
    while(length(edge_unique_iter) > 0 ) {
      
      # print(length(edge_unique_iter))
      edges_df_sub <- head(edges_df[edges_df$from %in% edge_unique_iter[1] | edges_df$to %in% edge_unique_iter[1], ],1)
      
      if(nrow(edges_df_sub) == 0) {
        edge_unique_iter <- setdiff(edge_unique_iter, edge_unique_iter[1])
        next
      }
      
      merged_df <- do.call(rbind, basins_list[c(edges_df_sub$from, edges_df_sub$to)])
      unique_critical <- unlist(lapply(strsplit(unique(merged_df$critical), ';'), function(k) gsub('_.*', '', k)))
      unique_samples <- unique(merged_df$sample)
      
      merged_df <- merged_df[!duplicated(merged_df$basin), ]
      # merged_name <- paste0("clique_", length(edge_unique_iter), "_", level)
      merged_name <- paste(paste(unique_critical, unique_samples, sep = '_'), collapse = '_')
      merged_df$critical <- merged_name
      merged_df <- list(merged_df)
      names(merged_df) <- merged_name
      merged_df_all <- c(merged_df_all, merged_df)
      
      edges_df <-  edges_df[!(edges_df$from[1] %in% edges_df_sub[, 1:2] | edges_df$to %in%edges_df_sub[, 1:2]), ]
      edge_unique_iter <- setdiff(edge_unique_iter, edge_unique_iter[[1]])
      toremove <- c(toremove, edges_df_sub$from, edges_df_sub$to)
    }
    
    
    # Update basins list
    basins_list <- c(merged_df_all, basins_list[!names(basins_list) %in% toremove])
    basin_sets <- lapply(basins_list, function(df) unique(df$basin))
    
    # Efficiently update overlap matrix
    new_names <- names(basins_list)
    new_n <- length(new_names)
    new_matrix <- matrix(0, nrow = new_n, ncol = new_n, dimnames = list(new_names, new_names))
    unchanged <- setdiff(new_names, names(merged_df_all))
    
    # Reuse existing overlaps
    if (exists("overlap_matrix")) {
      common <- intersect(unchanged, rownames(overlap_matrix))
      # new_matrix[common, common] <- overlap_matrix[common, common]
      new_matrix[common, common] <- 0
    }
    
    # Recompute only necessary overlaps
    for (i in names(merged_df_all)) {
      # i = 'clique_10041_1'
      # print(i)
      set_i <- basin_sets[[i]]
      # for (j in new_names) {
      for( j in c(names(merged_df_all), setdiff(edge_unique, toremove))) {
        # j = 'clique_10043_1'
        if (i == j) next
        set_j <- basin_sets[[j]]
        jac <- if (length(union(set_i, set_j)) == 0) 0 else length(intersect(set_i, set_j)) / length(union(set_i, set_j))
        new_matrix[i, j] <- jac
        new_matrix[j, i] <- jac
      }
    }
    overlap_matrix <- new_matrix
  }
  return(basins_list)
}



merged_basins_list <- merge.basins(basins_list = all_samples_basins_list, overlap_matrix = all_samples_overlap, jac_thres = 0.7, level = 0)
save(merged_basins_list, file = paste0(path_file, 'ALL_cell_line_merged_basin_list_non_outliers'))
saveRDS(merged_basins_list, file = paste0(path_file, 'ALL_cell_line_merged_basin_list_non_outliers.RDS'))
##################################################################################################################################################################
load(paste0(path_file, 'ALL_cell_line_merged_basin_list_non_outliers'))
# Size of modules
basin_length <- sapply(merged_basins_list, function(i) length(unique(i$basin)))
hist(log2(basin_length), breaks = 50)
sort(basin_length, decreasing = TRUE)[1:300]

unique_basins_list <- lapply(merged_basins_list, function(i) {
  
  res <- i[!duplicated(i$basin), ]
  res$basin
})
sort(table(gsub('_.*', '' , names(basin_length))), decreasing = TRUE)
length(unique_basins_list)
##################################################################################################################################################################
# Calculate profiles
unique_basins_profiles <- do.call(rbind, lapply(unique_basins_list, function(i) {
  
  res <- proteomics_quant[i, ]
  
  as.numeric(apply(res, 2, median, na.rm = TRUE))
}))
colnames(unique_basins_profiles) <- colnames(proteomics_quant)

basins_profiles_cor <- cor(t(unique_basins_profiles), method = 'pearson', use = 'pairwise.complete.obs')

row_clust <- hclust(as.dist(1- basins_profiles_cor), method = 'ward.D2')
ht_module_cor <- Heatmap(basins_profiles_cor, col = colorRamp2(breaks = c(-1,0, 1), colors = c('blue', 'white', 'red')),
                         show_column_names = FALSE,
                         show_row_names = FALSE,
                         cluster_columns = row_clust,
                         cluster_rows = row_clust,
                         show_row_dend = TRUE,
                         column_gap = unit(3, "mm"),
                         border = TRUE, 
                         column_title = paste0('ALL cell lines cohort (n_modules = ', nrow(basins_profiles_cor), ')'),
                         heatmap_legend_param = list(title = 'Peason cor.'),
                         column_title_side = "top",
                         column_title_gp = gpar(cex = 2, fontface = "bold"),
                         # row_title_rot = 0,
                         column_names_side = "top",
                         row_dend_width = unit(2, "cm"),
                         column_dend_height = unit(3, "cm"),
                         row_names_gp = gpar(fontsize = 12, col = 'black'),
                         column_names_gp = gpar(fontsize = 12, col = 'black'))
pdf(file = paste0(path_save, 'ALL_cell_line_heatmap_critical_node_comb_non_outliers.pdf'), width = 10, 
    height = 8)
draw(ht_module_cor, padding = unit(c(2, 20, 2, 2), "mm"))
dev.off()

write.table(unique_basins_profiles, paste0(path_file, 'ALL_cell_line_module_median_abund_non_outliers.txt'), sep = '\t')
##################################################################################################################################################################
## NMF clustering
unique_basins_profiles <- read.delim(paste0(path_file, 'ALL_cell_line_module_median_abund_non_outliers.txt'), sep = '\t')


### NMF
nmf_input <- as.matrix(unique_basins_profiles)

# Complete cases
nmf_input <- nmf_input[complete.cases(nmf_input), ]


enrich_heat <- Heatmap(nmf_input, clustering_distance_columns = 'pearson', 
                                       clustering_distance_rows = 'pearson', 
                                       clustering_method_columns = 'ward.D2', 
                                       clustering_method_rows = 'ward.D2',
                                       show_column_names = FALSE,
                                       show_row_names = FALSE,
                                       show_row_dend = TRUE,
                                       column_gap = unit(3, "mm"),
                                       border = TRUE, 
                                       column_title = paste0('ALL cell lines cohort (n_modules = ', nrow(nmf_input), ')'),
                                       heatmap_legend_param = list(title = 'log2 ratio'),
                                       column_title_side = "top",
                                       column_title_gp = gpar(cex = 2, fontface = "bold"),
                                       # row_title_rot = 0,
                                       column_names_side = "top",
                                       row_dend_width = unit(2, "cm"),
                                       column_dend_height = unit(3, "cm"),
                                       row_names_gp = gpar(fontsize = 12, col = 'black'),
                                       column_names_gp = gpar(fontsize = 12, col = 'black'))

pdf(file = paste0(path_save, 'ALL_cell_line_heatmap_basins_NMF_non_outliers.pdf'), width = 10, 
    height = 10)
draw(enrich_heat, padding = unit(c(2, 20, 2, 2), "mm"))
dev.off()

nmf_input_pos <- nmf_input %>%  posneg

estim.r <- nmf(x = nmf_input_pos,
               rank = 2:20,
               method = 'brunet',
               nrun = 30,
               seed = 1234,
               .options = 'p10')

saveRDS(estim.r, file = paste0(path_file, "ALL_cell_line_NMF_non_outliers.RDS"))
estim.r <- readRDS(file = paste0(path_file, "ALL_cell_line_NMF_non_outliers.RDS"))

pdf(paste0(path_save, 'ALL_cell_line_proteomics_NMF_res_non_outliers.pdf'), width = 12, height = 10)
plot(estim.r)
dev.off()

## Cophenetic correlation
K <- 6


## Membership index
# H matrix 
H <- as.data.frame(t(coef(estim.r$fit[[K-1]])))
h_membership <- as.data.frame(t(apply(H, 1, function(i) {
  i/sum(i)
})))
nmf_clusters <- apply(h_membership, 1, which.max)

proteomics_meta$nmf_clusters <- nmf_clusters[rownames(proteomics_meta)]

table(proteomics_meta$nmf_clusters, proteomics_meta$Grouped_Subtype)

# Cophonetic correlation
coph_cor <- sapply(1:19,function(i) {
  cophcor(estim.r$consensus[[i]])
})
coph_cor_dt <- data.frame(k = 2:20, coph_cor = coph_cor)

p_coph_cor <- ggplot(coph_cor_dt ,aes(x = k, y = coph_cor)) + 
  geom_line() + 
  geom_point(pch = 21, size = 3, fill = 'white') +
  geom_point(data = coph_cor_dt[coph_cor_dt$k == K, ], pch = 21, size = 6, color = 'red', stroke = 2) + 
  labs(x = 'k', y = 'cophenetic correlation coefficient', title = 'NMF') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))

pdf(paste0(path_save, 'ALL_cell_line_proteomics_NMF_basin_cophonetic_non_outliers.pdf'), width = 6, height = 5)
plot(p_coph_cor)
dev.off()



# B matrix 
B <- as.data.frame(basis(estim.r$fit[[K-1]]))


scaled_profiles <- t(scale(t(unique_basins_profiles)))
# scaled_profiles <- unique_basins_profiles

top_scores <- do.call(rbind, lapply(1:nrow(scaled_profiles), function(i) {

  val <- as.numeric(scaled_profiles[i, ])
  
  dt <- data.frame(scale(H), val = val)

  res <- lm(val ~ ., data = dt)
  sum_dt <- summary(res)
  sum_dt$coefficients[, 't value'][-1]
}))
rownames(top_scores) <- rownames(unique_basins_profiles)
colnames(top_scores) <- paste0('NMF', 1:K)


write.table(top_scores, paste0(path_file, 'ALL_cell_line_nmf_modules_scores_non_outliers.txt'), sep = '\t')


top_score_mods <- apply(top_scores, 2, function(k) {
  sorted_mods <- rownames(top_scores)[order(k, decreasing = TRUE)]
  c(head(sorted_mods, 5), tail(sorted_mods, 5))
})



top.scores.cor <- function(top_score_mods) {
  
  top_scores_pair <-unlist(apply(top_score_mods, 2, function(i) apply(combn(sort(i, decreasing = TRUE), m = 2), 2, paste, collapse = '_')))
  top_scores_pair_freq <- sort(table(top_scores_pair), decreasing = TRUE)
  
  res_all <- do.call(rbind, lapply(colnames(top_score_mods), function(i) {
    
    res <- cor(t(unique_basins_profiles[top_score_mods[,i], ]), use = 'pairwise.complete.obs', method = 'pearson')
    res <- reshape2::melt(res)
    res <- res[as.character(res$Var1) > as.character(res$Var2), ]
    res$spec <- 1/top_scores_pair_freq[paste(res$Var1, res$Var2, sep = '_')]
    # res$value[abs(res$value) < 0.5] <- NA
    
    unique_genes1 <- unique_basins_list[as.character(res$Var1)]
    unique_genes2 <- unique_basins_list[as.character(res$Var2)]
    
    nmf_samples <- names(nmf_clusters)[nmf_clusters == gsub('NMF', '', i)]
    res$abnd1 <- apply(unique_basins_profiles[as.character(res$Var1), nmf_samples], 1, median, na.rm =TRUE)
    res$abnd2 <- apply(unique_basins_profiles[as.character(res$Var2), nmf_samples], 1, median, na.rm =TRUE)
    
    res$n_genes1 <- sapply(unique_genes1, length)
    res$n_genes2 <- sapply(unique_genes2, length)
    res$genes1 <- sapply(unique_genes1, paste, collapse = ', ')
    res$genes2 <- sapply(unique_genes2, paste, collapse = ', ')
    
    res$cluster <- i
    
    res
  }))
  
  res_all
}

top_scores_cor <- top.scores.cor(top_score_mods = top_score_mods)
top_score_mods_order <- reshape2::melt(top_score_mods)
top_score_mods_high <- as.character(apply(top_score_mods, 2, function(i) i[1:5]))

# Enrichment analysis
msig_cp <- GSA.read.gmt(paste0(path_file, 'c2.cp.reactome.v2025.1.Hs.symbols.gmt.txt'))
msig_cp_list <- msig_cp$genesets
names(msig_cp_list) <- msig_cp$geneset.names
msig_cp_list <- reshape2::melt(msig_cp_list)
colnames(msig_cp_list) <- c('gene', 'name')
msig_cp_list <- msig_cp_list[, c('name', 'gene')]


module_enrichment <- do.call(rbind, lapply(as.character(top_score_mods), function(i) {
  
  # i = 'DNTT_FORALL_063'
  # print(i)
  
  basin_genes <- unique_basins_list[[i]]
  
  enrich_res <- enricher(gene =  basin_genes,
                  pvalueCutoff = 0.01,
                  pAdjustMethod = "BH",
                  universe = rownames(proteomics_quant),
                  minGSSize = 3,
                  maxGSSize = max(basin_length),
                  qvalueCutoff = 0.01,
                  TERM2GENE = msig_cp_list,
                  TERM2NAME = NA)
  
  enrich_res <- as.data.frame(enrich_res)
  enrich_res <- enrich_res[which(enrich_res$Count>=3), ]
  
  if(nrow(enrich_res) == 0) {
    enrich_res[1, ] = NA
    enrich_res$cluster <- 1
    enrich_res$Module <- i
    enrich_res$NMF <- top_score_mods_order$Var2[match(i, top_score_mods_order$value)]
    enrich_res$dir <- ifelse(i %in% top_score_mods_high, 'high', 'low')
    
    # write.table(enrich_res,  paste0(path_file,'ALL_cell_line_module_enrich_', i, '.txt'), sep = '\t')
    return(NULL)
  }
  
  if(nrow(enrich_res) == 1) {
    enrich_res$cluster <- 1
    enrich_res$Module <- i
    enrich_res$NMF <- top_score_mods_order$Var2[match(i, top_score_mods_order$value)]
    enrich_res$dir <- ifelse(i %in% top_score_mods_high, 'high', 'low')
    
    # write.table(enrich_res,  paste0(path_file,'ALL_cell_line_module_enrich_', i, '.txt'), sep = '\t')
    return(enrich_res)
  }

  
  # Ontology terms similarity matrix
  enrich_res_sig <- enrich_res
  enrich_res_sig_genes <- strsplit(enrich_res_sig$geneID, '\\/')
  names(enrich_res_sig_genes) <- enrich_res_sig$ID
  
  pairwise_overlap_terms <- sapply(enrich_res_sig_genes, function(i) {
    sapply(enrich_res_sig_genes, function(j) {
      # length(intersect(i,j))/length(union(i,j))
      length(intersect(i,j))/min(length(i), length(j))
    })
  })

  
  # Convert similarity matrix to a long-format data frame
  pairwise_overlap_terms <- reshape2::melt(pairwise_overlap_terms)
  
  # Remove self-comparisons (term vs itself)
  pairwise_overlap_terms <- pairwise_overlap_terms[pairwise_overlap_terms$Var1 != pairwise_overlap_terms$Var2, ]
  
  # Keep only unique term pairs (avoid duplicates like A–B and B–A)
  pairwise_overlap_terms <- pairwise_overlap_terms[as.character(pairwise_overlap_terms$Var1) > 
                                                     as.character(pairwise_overlap_terms$Var2), ]
  
  # Filter to retain only edges with strong similarity (> 0.5)
  pairwise_overlap_terms_sig <- pairwise_overlap_terms[pairwise_overlap_terms$value > 0.5, ]
  
  
  # Build an undirected graph where nodes = terms, edges = overlap relationships
  enrich_g <- graph_from_data_frame(pairwise_overlap_terms_sig, directed = FALSE)
  
  # Identify terms not connected to others (no edge > 0.5)
  isolated_terms <- setdiff(names(enrich_res_sig_genes), names(V(enrich_g)))
  
  # Add these isolated terms as standalone nodes
  enrich_g <- add_vertices(enrich_g, 
                           nv = length(isolated_terms), 
                           name = isolated_terms)
  
  # Assign node size based on the number of genes associated with each term
  V(enrich_g)$size <- enrich_res$Count[match(V(enrich_g)$name, enrich_res$ID)]
  
  
  
  # Detect functional clusters using the Louvain community detection algorithm
  term_cluster <- cluster_louvain(enrich_g)
  
  # Assign cluster labels (C1, C2, etc.) to nodes for color-coding in the plot
  V(enrich_g)$Cluster <- paste0('C', membership(term_cluster)[names(V(enrich_g))])
  
  # Visualize
  group_map <- data.frame('name' = names(V(enrich_g)), 
                          'cluster' = V(enrich_g)$Cluster)
  group_map$FoldEnrichment <- enrich_res_sig$FoldEnrichment[match(group_map$name, enrich_res_sig$ID)]
  
  top_terms <- group_map %>%
    group_by(cluster) %>%
    arrange(desc(abs(FoldEnrichment))) %>%
    slice(1) %>%
    ungroup()
  
  top_terms$padj <-  enrich_res_sig$p.adjust[match(top_terms$name, enrich_res_sig$ID)]
  top_terms$genes <- enrich_res_sig$geneID[match(top_terms$name, enrich_res_sig$ID)]
  top_terms$size <- enrich_res_sig$Count[match(top_terms$name, enrich_res_sig$ID)]
  top_terms$name <- factor(top_terms$name, levels = top_terms$name[order(top_terms$FoldEnrichment)])
  top_terms <- top_terms[order(top_terms$name), ]

  
  enrich_res_sig$cluster <- membership(term_cluster)[enrich_res_sig$ID]
  enrich_res_sig <- enrich_res_sig[order(enrich_res_sig$cluster, -enrich_res_sig$FoldEnrichment), ]
  enrich_res_sig$Module <- i
  enrich_res_sig$NMF <- top_score_mods_order$Var2[match(i, top_score_mods_order$value)]
  enrich_res_sig$dir <- ifelse(i %in% top_score_mods_high, 'high', 'low')
  
  
  return(enrich_res_sig)

  
  }))

# Choose one representative
module_enrichment_grouped <- module_enrichment %>% group_by(NMF, Module, cluster, dir) %>%
  slice(1) 
module_enrichment_grouped$Module <- factor(module_enrichment_grouped$Module, levels = unique(as.character(top_score_mods)))
module_enrichment_grouped$ID <- factor(module_enrichment_grouped$ID, levels = unique(module_enrichment_grouped$ID[order(module_enrichment_grouped$Module)]))

plot_enrich <- ggplot(module_enrichment_grouped, mapping = aes(x = Module, y = ID)) + 
         geom_tile(aes(fill = dir)) + 
  facet_grid(.~NMF, scale = 'free', space = 'free') + 
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45)) + 
  guides(fill = guide_legend(title = ''))


pdf(paste0(path_save, 'ALL_cell_line_top_modules_enrichment.pdf'), width = 18, height = 15)
plot(plot_enrich)
dev.off()


write.table(module_enrichment, paste0(path_file, 'ALL_cell_line_top_score_module_enrichment.txt'), sep = '\t', row.names = FALSE)
write.table(top_scores_cor, paste0(path_file, 'ALL_cell_line_top_score_cor_non_outliers.txt'), sep = '\t', row.names = FALSE)
################################################################################################################################
# Build graph with median abundance
lapply(colnames(top_score_mods), function(nmf) {
  
  top_scores_cor_sub <- top_scores_cor[top_scores_cor$cluster == nmf, ]
  
  # Build edge list with aesthetics included
  edges <- top_scores_cor_sub %>%
    select(Var1, Var2, value, spec) %>%
    mutate(
      edge_color = ifelse(value > 0.5, "positive", ifelse(value < -0.5, "negative", 'else')),
      edge_width =  as.numeric(spec))
  
  # Build node list
  nodes1 <- top_scores_cor_sub %>%
    select(name = Var1, med_abnd = abnd1, n_genes = n_genes1) %>%
    mutate(across(c(med_abnd, n_genes), as.numeric))
  nodes2 <- top_scores_cor_sub %>%
    select(name = Var2, med_abnd = abnd2, n_genes = n_genes2) %>%
    mutate(across(c(med_abnd, n_genes), as.numeric))
  
  nodes <- bind_rows(nodes1, nodes2) %>%
    group_by(name) %>%
    summarise(
      med_abnd = mean(med_abnd, na.rm = TRUE),
      n_genes = mean(n_genes, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Build igraph
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  E(g)$weight <- ifelse(E(g)$edge_color == "positive", 2, 0.5)
  g_filtered <- delete_edges(g, which(E(g)$edge_color == "else"))
  
  set.seed(123)  # ensures reproducibility
  layout_coords <- layout_with_fr(g_filtered, weights = E(g_filtered)$weight)
  
  # Plot
  graph_sub <- ggraph(g_filtered, layout =layout_coords) +
    geom_edge_link(aes(width = edge_width, color = edge_color)) +
    geom_node_point(aes(color = med_abnd, size = n_genes), alpha = 1) +
    scale_edge_color_manual(breaks =  c("positive", "negative"), values = c("positive" = "red", "negative" = "blue"), 
                            labels = c('> 0.5', '< -0.5')) +
    scale_color_gradient2(low = "blue", mid = 0, high = "red") +
    scale_edge_width(range = c(0.1, 1)) +
    scale_size_continuous(range = c(3, 10)) + 
    geom_label_repel(aes(label = name, x = x, y = y), 
                     size = 3, 
                     color = "black", 
                     box.padding = 0.35, 
                     point.padding = 0.5,
                     min.segment.length = 0.1,
                     segment.color = "black") +
    labs(title = nmf) + 
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
    guides(edge_width = "none", 
           edge_color = guide_legend(order = 1), 
           color = guide_colorbar(order = 2),
           size = guide_legend(order = 3))
  
  
  pdf(paste0(path_save, 'ALL_cell_line_top_modules_',nmf, '.pdf'), width = 6, height = 6)
  plot(graph_sub)
  dev.off()
  
})  


# Sample specific
nmf.net.abnd <- function(sample_id, top_scores_cor, unique_basins_profiles) {
  
  # lapply(names(nmf_cluster), function(sample_id) {
  
  # print(sample_id)
  nmf <- paste0('NMF', as.character(nmf_clusters[sample_id]))
  
  top_scores_cor_sub <- top_scores_cor[top_scores_cor$cluster == nmf, ]
  nmf_mods <- union(top_scores_cor_sub$Var1, top_scores_cor_sub$Var2)
  
  # Build edge list with aesthetics included
  edges <- top_scores_cor_sub %>%
    select(Var1, Var2, value, spec) %>%
    mutate(
      edge_color = ifelse(value > 0.5, "positive", ifelse(value < -0.5, "negative", 'else')),
      edge_width =  as.numeric(spec))
  
  top_scores_cor_sub$abnd1 <- unique_basins_profiles[as.character(top_scores_cor_sub$Var1), sample_id]
  top_scores_cor_sub$abnd2 <- unique_basins_profiles[as.character(top_scores_cor_sub$Var2), sample_id]
  
  # Build node list
  nodes1 <- top_scores_cor_sub %>%
    select(name = Var1, med_abnd = abnd1, n_genes = n_genes1) %>%
    mutate(across(c(med_abnd, n_genes), as.numeric))
  nodes2 <- top_scores_cor_sub %>%
    select(name = Var2, med_abnd = abnd2, n_genes = n_genes2) %>%
    mutate(across(c(med_abnd, n_genes), as.numeric))
  nodes <- bind_rows(nodes1, nodes2) %>%
    group_by(name) %>%
    summarise(
      med_abnd = mean(med_abnd, na.rm = TRUE),
      n_genes = mean(n_genes, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Build igraph
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  E(g)$weight <- ifelse(E(g)$edge_color == "positive", 2, 0.5)
  g_filtered <- delete_edges(g, which(E(g)$edge_color == "else"))
  
  set.seed(123)  # ensures reproducibility
  layout_coords <- layout_with_fr(g_filtered, weights = E(g_filtered)$weight)
  
  # Plot
  graph_sub <- ggraph(g_filtered, layout =layout_coords) +
    geom_edge_link(aes(width = edge_width, color = edge_color)) +
    geom_node_point(aes(color = med_abnd, size = n_genes), alpha = 1) +
    scale_edge_color_manual(values = c("positive" = "red", "negative" = "blue"), 
                            labels = c('> 0.5', '< -0.5')) +
    scale_color_gradient2(name = 'Abnd', low = "blue", mid = 0, high = "red") +
    scale_edge_width(range = c(0.1, 1)) +
    scale_size_continuous(range = c(3, 10)) + 
    geom_label_repel(aes(label = name, x = x, y = y), 
                     size = 3, 
                     color = "black", 
                     box.padding = 0.35, 
                     point.padding = 0.5,
                     min.segment.length = 0.1,
                     segment.color = "black") +
    labs(title = paste0(sample_id, ' - ', nmf)) + 
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
    guides(edge_width = "none", 
           edge_color = guide_legend(order = 1), 
           color = guide_colorbar(order = 2),
           size = guide_legend(order = 3))
  
  
  
}
## Sample specific
lapply(names(nmf_clusters), function(sample_id) {
  
  p <- nmf.net.abnd (sample_id = sample_id, top_scores_cor = top_scores_cor, unique_basins_profiles = unique_basins_profiles)
  
  pdf(paste0(path_save, 'ALL_cell_line_net_top_modules_', sample_id, '.pdf'), width = 6, height = 6)
  plot(p)
  dev.off()
})

################################################################################################################################
## Complex Heatmap membership
colnames(h_membership) <- paste0('NMF', 1:K)
h_membership$class <- apply(h_membership, 1, function(i) colnames(h_membership)[which.max(i)])
h_membership$max_val <- apply(h_membership[, paste0('NMF', 1:K)], 1, function(i) max(i))

sample_order <- rownames(h_membership)[order(h_membership$class, -h_membership$max_val)]

ht_memb <- Heatmap(t(h_membership[sample_order, paste0('NMF', 1:K)]), 
                   col = circlize::colorRamp2(breaks = c(0,0.5, 1), colors = c('blue', 'white', 'red')),
                   column_split = h_membership[sample_order, 'class'],
                   column_title = NULL,
                   show_column_names = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   row_names_side = 'left',
                   show_row_dend = TRUE,
                   show_row_names = TRUE,
                   column_gap = unit(3, "mm"),
                   border = TRUE, 
                   heatmap_legend_param = list(title = 'membership'),
                   column_title_side = "top",
                   column_names_side = "bottom",
                   row_dend_width = unit(2, "cm"),
                   row_names_gp = gpar(fontsize = 12, col = 'black'),
                   height = unit(4, "cm") 
                
)

write.table(h_membership[sample_order, ], paste0(path_file, 'ALL_cell_line_per_sample_membership_non_outliers.txt'), sep = '\t')

ht_prof <- Heatmap(unique_basins_profiles[ top_score_mods_order$value, sample_order], 
                   # col = circlize::colorRamp2(breaks = c(0,1), colors = c('white', 'red')),
                   row_split = top_score_mods_order$Var2,
                   column_split = h_membership[sample_order, 'class'],
                   column_title = NULL,
                   row_title = NULL,
                   show_column_names = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   clustering_distance_rows = 'pearson',
                   clustering_method_rows = 'ward.D2',
                   show_row_dend = FALSE,
                   show_row_names = TRUE,
                   row_gap = unit(0, "mm"),
                   column_gap = unit(3, "mm"),
                   border = TRUE, 
                   heatmap_legend_param = list(title = 'log2 ratio'),
                   column_title_side = "top",
                   row_title_rot = 0,
                   column_names_side = "bottom",
                   row_names_gp = gpar(fontsize = 12, col = 'black'),
                   height = unit(22, "cm"))

ht_all <-  ht_memb %v% ht_prof

pdf(file = paste0(path_save, 'ALL_cell_line_NMFs_heatmap_non_outliers.pdf'), width = 15, height = 12)
draw(ht_all, padding = unit(c(2, 20, 2, 2), "mm"))
dev.off()

###