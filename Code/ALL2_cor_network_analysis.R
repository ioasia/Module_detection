### ALL cell line network reconstruction

## Libraries
source('~/Downloads/ALL2_cell_line_proj/code/exploratory/classifier_functions.R')
library(Biobase)
library(data.table)
library(Hmisc)
library(openxlsx)


path_file <- '~/Downloads/ALL2_cell_line_proj/data/'
path_save <- '~/Downloads/ALL2_cell_line_proj/code/figures/'

# Data
proteomics_quant <- read.delim(paste0(path_file, 'ALL_cell_line_proteomics_median_centered.txt'), sep = '\t')
idxtokeep <- apply(proteomics_quant, 1, function(i) sum(!is.na(i)) >= 0.7 * ncol(proteomics_quant))
proteomics_quant_sub <- proteomics_quant[idxtokeep, ]
boxplot(proteomics_quant_sub)

proteomics_meta <- read.delim(paste0(path_file, 'ALL_cell_line_meta_data_init.txt'), sep = '\t')
##############################################################################################################################################################################################
## CORUM paits
corum_pairs <- read.delim(paste0(path_file, 'CORUM_complex_pairs.txt'), sep = '\t')
##############################################################################################################################################################################################
## Igraph 

## Correlation-based
coreg_res <- rcorr(t(proteomics_quant_sub), type= "pearson")
cor_est <- reshape2::melt(coreg_res$r)
cor_n <- reshape2::melt(coreg_res$n)
cor_p <- reshape2::melt(coreg_res$P)


coreg <- cbind.data.frame(cor_est, 'n' = cor_n$value, 'p' = cor_p$value)
coreg <- coreg[as.character(coreg$Var1) > as.character(coreg$Var2), ]
coreg <- coreg[coreg$n >= 30, ]

colnames(coreg)[1:3] <- c('Protein1', 'Protein2', 'sim')

f_dist <- feature.dist(features = coreg$sim, feature_name = 'Pearson cor.')
ggsave(paste0(path_save, 'ALL_cell_lines_proteomics_feature_dist.pdf'), plot = f_dist, width = 6, height = 5)


res_protein <- annotate.mat(cormat = coreg, corum_pairs = corum_pairs, y_title = 'Pearson cor.')
cormat_annot <- res_protein$cormat_annot

ggsave(paste0(path_save, 'ALL_cell_line_corum_dist.pdf'), plot = res_protein$plot, width = 7, height = 5)
####################################################################################################################################
# Threshold metric 
protein_class_res <- cv.fold(cor_mat = cormat_annot, max_iterations = 100, sim_col = 'sim', class_col = 'class')

p_fdr <- plot.fdr(class_res = protein_class_res, titl = 'Proteomics')
ggsave(paste0(path_save, 'ALL_cell_line_proteomics_fdr_iter_method_comp.pdf'), plot = p_fdr, width = 7.5, height = 4)

p_auc <- plot.auc(class_res = protein_class_res, titl = 'Proteomics')
ggsave(paste0(path_save, 'ALL_cell_line_proteomics_auc_iter_method_comp.pdf'), plot = p_auc, width = 6, height = 4)

p_pres_rec <- plot.pres.recall(class_res = protein_class_res, titl = 'Proteomics')
ggsave(paste0(path_save, 'ALL_cell_line_proteomics_prec_recall_iter_method_comp.pdf'), plot = p_pres_rec, width = 5, height = 7)


protein_cormat_pred <- generate.pred(class_res = protein_class_res, cor_mat = cormat_annot)


# Probability distributions
p_dist <- prob.dist(features = protein_cormat_pred$sim, preds = protein_cormat_pred$pred, feature_name = 'Pearson cor.')
ggsave(paste0(path_save, 'ALL_cell_line_proteomics_prob_dist.pdf'), plot = p_dist, width = 5, height = 8)


# Size of connected network by threshold used
p_net <- net.size(prob_mat = protein_cormat_pred, min_thres = 0.7, max_thres = 0.9, step_thres = 10, pred_col = 'pred')
ggsave(paste0(path_save, 'ALL_cell_line_proteomics_graph_quality.pdf'), plot = p_net, width = 5, height = 8)

# Subset to significant
write.table(protein_cormat_pred, paste0(path_file, 'ALL_cell_line_proteomics_igraph_predictions.txt'), sep = '\t')
####################################################################################################################################
# Set threshold 
colnames(protein_cormat_pred)[colnames(protein_cormat_pred) %in% c('sim', 'pred')] <- paste0('protein_',  c('sim', 'pred'))

gc()

res_mat <- protein_cormat_pred
thres_sim <- 0.8

# Change name of columns
res_mat_sig <- as.data.frame(res_mat[res_mat$protein_pred > thres_sim, ])


igraph_snn_gephi <- graph_from_data_frame(res_mat_sig[, c("Protein1", "Protein2")], directed = FALSE)
igraph_snn_gephi <- connected.component(i = igraph_snn_gephi)

p_degree <- degree.hist(g = igraph_snn_gephi, thres = thres_sim)
ggsave(paste0(path_save, 'ALL_cell_line_network_degree_dist.pdf'), plot = p_degree$p, width = 6, height = 5)

p_top <- degree.scatter(g = igraph_snn_gephi, thres = thres_sim)
ggsave(paste0(path_save, 'ALL_cell_line_network_degree_top.pdf'), plot  = p_top$p, width = 6, height = 5)


igraph_snn_gephi_dt <- as.data.frame(as_edgelist(igraph_snn_gephi))
igraph_snn_gephi_dt <- as.data.frame(t(apply(igraph_snn_gephi_dt, 1, sort, decreasing = TRUE)))
res_mat_sig <- res_mat_sig[paste0(res_mat_sig$Protein1, res_mat_sig$Protein2) %in% 
                             paste0(igraph_snn_gephi_dt[, 1], igraph_snn_gephi_dt[, 2]), ]

write.table(res_mat_sig, paste0(path_file, 'ALL_cell_line_merged_igraph_predictions_sig.txt'), sep = '\t', row.names = FALSE)
#####################################################################################################################################
## Tissue specificity with lung, brain, breast and colon 
sig_pairs <- paste0(res_mat_sig$Protein1, '_', res_mat_sig$Protein2)
coexpr_dat <- list.files(path = '/Volumes/LaCie/Nat_biotech_coppi/(fig2a)_probabilities_cohorts_combined/', pattern = 'output.txt')
tumor_types <- sapply(coexpr_dat, function(i) gsub(".*(blood|brain|lung|breast|colon|kidney|liver|ovary|pancreas|stomach|throat).*", "\\1", i))

res_coexp <- as.data.frame(do.call(rbind, lapply(1:length(coexpr_dat), function(j) {
  
  print(coexpr_dat[j])
  
  coexpr <- fread(paste0('/Volumes/LaCie/Nat_biotech_coppi/(fig2a)_probabilities_cohorts_combined/', coexpr_dat[j]), sep = '\t')
  idx_disim <- which(coexpr$prot1 < coexpr$prot2)
  coexpr[idx_disim, c("prot1", "prot2") := .(prot2, prot1)]
  
  coexpr$pair <- paste0(coexpr$prot1, '_', coexpr$prot2)
  
  common_pairs <- intersect(coexpr$pair, sig_pairs)
  coexpr <- coexpr[coexpr$pair %in% common_pairs, ]
  coexpr_sig <- coexpr[coexpr[[3]] > 0.5]
  
  perc <- nrow(coexpr_sig)/nrow(coexpr)
  
  data.frame('dat' = coexpr_dat[j],
             'tissue' = tumor_types[j],
             'perc' = perc)
  
})))

res_coexp$tissue <- factor(res_coexp$tissue, unique(res_coexp$tissue[order(res_coexp$perc, decreasing = FALSE)]))

p_coexpr <- ggplot(res_coexp, aes(x = tissue, y = perc)) + 
  geom_point(stat = 'identity') + 
  geom_line(group = 1) + 
  geom_text(mapping = aes(label = round(perc, 3)), hjust = -0.1) + 
  scale_y_continuous(expand = c(0,0,0.1,0.1)) + 
  labs(x = '', y = 'Probability to be likely (score > 0.5)', title = 'Tissue specificity') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip()

ggsave(paste0(path_save, 'ALL_cell_line_interaction_coexpr_prob.pdf'), plot = p_coexpr, width = 6, height = 4)
####################################################################################################################################
## Pairwise associations
col_names <- c("prot1", "prot2")
colnames(corum_pairs)[colnames(corum_pairs) %in% c('V1', 'V2')] <-  col_names

string_pairs <- read.delim(paste0(path_file, 'stringdb_400.txt'), sep = '\t')
colnames(string_pairs)[colnames(string_pairs) %in% c("gene1",  "gene2")] <- col_names

reactome_pairs <- read.delim(paste0(path_file, 'reactome_gene_symbols.txt'), sep = '\t')
colnames(reactome_pairs)[colnames(reactome_pairs) %in% c("gene_symbol1",  "gene_symbol2")] <- col_names

signor_pairs <-  read.delim(paste0(path_file,'signorDB.txt'), sep = '\t')
colnames(signor_pairs)[colnames(signor_pairs) %in% c("ENTITYA", "ENTITYB")] <- col_names

bioplex_pairs <- read.delim(paste0(path_file, 'Bioplex_combined.txt'), sep = '\t')
colnames(bioplex_pairs) <-  col_names

huri_pairs <- read.delim(paste0(path_file, 'HuRI.txt'), sep = '\t')
colnames(huri_pairs)[colnames(huri_pairs) %in% c("gene_symbol1", "gene_symbol2")] <- col_names

funmap_pairs <- read.delim(paste0(path_file, 'FunMap.tsv'), sep = '\t', header = FALSE)
colnames(funmap_pairs) <- col_names

pairs_list <- list('CORUM' = corum_pairs, 
                   'STRING' = string_pairs, 
                   'Reactome' = reactome_pairs, 
                   'Signor' = signor_pairs, 
                   'BioPlex' = bioplex_pairs,
                   'HuRI' = huri_pairs,
                   'FunMap' = funmap_pairs
)


all_genes <- union(res_mat_sig$Protein1, res_mat_sig$Protein2)

pair_dat_overlap <- sapply(pairs_list, function(k) {
  
  # k <- pairs_list[[2]]
  k_genes <- unique(c(k$prot1, k$prot2))
  # Shared genes with ALL
  shared_genes <- intersect(all_genes, k_genes)
  
  k <- k[((k$prot1 %in% shared_genes) | (k$prot2 %in% shared_genes)), ]
  k_pairs_shared <- paste0(k$prot1, '_', k$prot2)
  
  res_mat_sig_shared <- res_mat_sig[(res_mat_sig$Protein1 %in% shared_genes) | (res_mat_sig$Protein2 %in% shared_genes), ]
  sig_pairs_shared <- paste0(res_mat_sig_shared$Protein1, '_', res_mat_sig_shared$Protein2)
  
  pairs_all <- which(k_pairs_shared %in% sig_pairs_shared)
  length(pairs_all)/length(k_pairs_shared)
  
})

pair_dat_overlap <- data.frame(origin = names(pair_dat_overlap), 
                               overlap =  pair_dat_overlap)

pair_dat_overlap$origin <- factor(pair_dat_overlap$origin, levels = rev(c(names(pairs_list))))

p_overlap_dat <- ggplot(pair_dat_overlap, aes(x = origin, y = overlap * 100, fill = origin)) + 
  geom_bar(stat = 'identity') +
  geom_text(mapping = aes(label = round(overlap, 2) * 100), hjust = -0.1)  +
  scale_fill_manual(values = c('CORUM' = '#9c5d9fff', "STRING" = '#9c5d9fff',  "Reactome" = '#9c5d9fff', "Signor" = '#9c5d9fff',
                               "BioPlex"= '#51a4d9ff', "HuRI" = '#51a4d9ff', 'FunMap' = '#51a4d9ff', 'Quantified' = 'grey')) + 
  theme_classic() + 
  # facet_wrap(.~thres, ncol = 1, labeller = labeller(thres = c("0.5" = "thres: 0.5",
  #                                                             "0.8" = "thres: 0.8"))) + 
  labs(x = '', y = 'Overlap (%)', title = 'ALL cell line network overlap') + 
  scale_y_continuous(expand = c(0,0, 0.1, 0.1)) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip() + 
  guides(fill = 'none')

ggsave(paste0(path_save, 'ALL_cell_line_interaction_pair_db_overlap.pdf'), plot = p_overlap_dat, width = 6, height = 4)
####################################################################################################################################
# Print number of nodes and edges
final_graph_top <- igraph_snn_gephi
final_graph_top <- connected.component(final_graph_top)
final_graph <- final_graph_top
vcount(final_graph)
ecount(final_graph)
summary(final_graph)

nodes_df <- data.frame(ID = c(1:vcount(final_graph)), NAME = V(final_graph)$name)

nodes_att <- nodes_df[, 'ID', drop = FALSE]
nodes_att$Label <- V(final_graph_top)$name

edges_df <- as.data.frame(get.edges(final_graph, c(1:ecount(final_graph))))

# And without edge weights
write.gexf(nodes = nodes_df, 
           edges = edges_df,
           nodesAtt = nodes_att,
           defaultedgetype = "undirected", output = paste0(path_save, "ALL_cell_line_correlation_net.gexf"))

###