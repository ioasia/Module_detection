### Protein co-expression network reconstruction
#
# Builds a protein-protein interaction network from a quantitative abundance
# matrix by computing pairwise Pearson correlations, training a logistic
# regression classifier against CORUM complex pairs, and applying a
# probability threshold to retain significant edges.
#
# Usage (run from repo root):
#   Rscript Code/network_reconstruction.R \
#     --proteomics=Data/my_proteomics.txt \
#     --data-dir=Data \
#     --figures-dir=Figures \
#     --prefix=my_dataset \
#     --pred-threshold=0.8
#
# Required argument:
#   --proteomics       TSV file, proteins x samples, first column = protein names
#
# Optional arguments:
#   --corum            path to CORUM complex pairs file
#                      (default: {data-dir}/CORUM_complex_pairs.txt)
#   --data-dir         directory for reference files and outputs (default: Data)
#   --figures-dir      directory for figures (default: Figures)
#   --prefix           prefix for all output filenames (default: dataset)
#   --pred-threshold   probability threshold for significant edges (default: 0.8)
#
# Outputs (in data-dir):
#   {prefix}_igraph_predictions.txt       all protein pairs with predicted probabilities
#   {prefix}_igraph_predictions_sig.txt   filtered significant network

source("Code/classifier_functions.R")
library(Biobase)
library(data.table)
library(Hmisc)
library(ggrepel)

# ---- argument parsing -------------------------------------------------------
args <- (function() {
  a <- commandArgs(trailingOnly = TRUE)
  d <- list(
    "data-dir"       = "Data",
    "figures-dir"    = "Figures",
    "prefix"         = "dataset",
    "pred-threshold" = "0.8",
    "corum"          = "",
    "early-stop"     = FALSE
  )
  for (x in a) {
    m <- regmatches(x, regexec("^--([^=]+)=(.+)$", x))[[1]]
    if (length(m) == 3) {
      d[[m[2]]] <- m[3]
    } else if (grepl("^--", x)) {
      d[[sub("^--", "", x)]] <- TRUE
    }
  }
  d
})()

if (is.null(args$proteomics)) stop("--proteomics is required")

data_dir    <- args[["data-dir"]]
figures_dir <- args[["figures-dir"]]
prefix      <- args[["prefix"]]
thres_sim   <- as.numeric(args[["pred-threshold"]])
corum_file  <- if (nchar(args[["corum"]]) > 0) args[["corum"]] else
                 file.path(data_dir, "CORUM_complex_pairs.txt")

# ---- data -------------------------------------------------------------------
proteomics_quant <- read.delim(args$proteomics, sep = '\t', row.names = 1)
idxtokeep        <- apply(proteomics_quant, 1, function(i) sum(!is.na(i)) >= 0.7 * ncol(proteomics_quant))
proteomics_quant_sub <- proteomics_quant[idxtokeep, ]

corum_pairs <- read.delim(corum_file, sep = '\t')

# ---- pairwise correlations --------------------------------------------------
coreg_res <- Hmisc::rcorr(t(proteomics_quant_sub), type = "pearson")
cor_est   <- reshape2::melt(coreg_res$r)
cor_n     <- reshape2::melt(coreg_res$n)

coreg <- cbind.data.frame(cor_est, 'n' = cor_n$value, 'p' = reshape2::melt(coreg_res$P)$value)
coreg <- coreg[as.character(coreg$Var1) > as.character(coreg$Var2), ]
coreg <- coreg[coreg$n >= 30, ]
colnames(coreg)[1:3] <- c('Protein1', 'Protein2', 'sim')

f_dist <- feature.dist(features = coreg$sim, feature_name = 'Pearson cor.')
ggsave(file.path(figures_dir, paste0(prefix, '_feature_dist.pdf')), plot = f_dist, width = 6, height = 5)

# ---- annotate with CORUM and train classifier -------------------------------
res_protein  <- annotate.mat(cormat = coreg, corum_pairs = corum_pairs, y_title = 'Pearson cor.')
cormat_annot <- res_protein$cormat_annot
ggsave(file.path(figures_dir, paste0(prefix, '_corum_dist.pdf')), plot = res_protein$plot, width = 7, height = 5)

protein_class_res <- cv.fold(cor_mat = cormat_annot, max_iterations = 100, sim_col = 'sim', class_col = 'class')

ggsave(file.path(figures_dir, paste0(prefix, '_fdr.pdf')),        plot = plot.fdr(protein_class_res, prefix),        width = 7.5, height = 4)
ggsave(file.path(figures_dir, paste0(prefix, '_auc.pdf')),        plot = plot.auc(protein_class_res, prefix),        width = 6,   height = 4)
ggsave(file.path(figures_dir, paste0(prefix, '_prec_recall.pdf')), plot = plot.pres.recall(protein_class_res, prefix), width = 5,   height = 7)

protein_cormat_pred <- generate.pred(class_res = protein_class_res, cor_mat = cormat_annot)

ggsave(file.path(figures_dir, paste0(prefix, '_prob_dist.pdf')),
       plot = prob.dist(protein_cormat_pred$sim, protein_cormat_pred$pred, 'Pearson cor.'),
       width = 5, height = 8)

p_net <- net.size(prob_mat = protein_cormat_pred, pred_col = 'pred',
                  min_thres = 0.7, max_thres = 0.9, step_thres = 10)
ggsave(file.path(figures_dir, paste0(prefix, '_graph_quality.pdf')), plot = p_net, width = 5, height = 8)

write.table(protein_cormat_pred,
            file.path(data_dir, paste0(prefix, '_igraph_predictions.txt')), sep = '\t')

if (isTRUE(args[["early-stop"]])) {
  message("--early-stop: FDR diagnostics written to ", figures_dir, "/\n",
          "Inspect ", prefix, "_graph_quality.pdf to choose --pred-threshold,\n",
          "then re-run without --until=fdr.")
  quit(save = "no", status = 0)
}

# ---- apply threshold and extract largest connected component ----------------
colnames(protein_cormat_pred)[colnames(protein_cormat_pred) %in% c('sim', 'pred')] <-
  paste0('protein_', c('sim', 'pred'))

res_mat     <- protein_cormat_pred
res_mat_sig <- as.data.frame(res_mat[res_mat$protein_pred > thres_sim, ])

igraph_net <- graph_from_data_frame(res_mat_sig[, c("Protein1", "Protein2")], directed = FALSE)
igraph_net <- connected.component(i = igraph_net)

ggsave(file.path(figures_dir, paste0(prefix, '_degree_dist.pdf')),
       plot = degree.hist(g = igraph_net, thres = thres_sim)$p, width = 6, height = 5)
ggsave(file.path(figures_dir, paste0(prefix, '_degree_scatter.pdf')),
       plot = degree.scatter(g = igraph_net, thres = thres_sim)$p, width = 6, height = 5)

igraph_dt  <- as.data.frame(as_edgelist(igraph_net))
igraph_dt  <- as.data.frame(t(apply(igraph_dt, 1, sort, decreasing = TRUE)))
res_mat_sig <- res_mat_sig[paste0(res_mat_sig$Protein1, res_mat_sig$Protein2) %in%
                             paste0(igraph_dt[, 1], igraph_dt[, 2]), ]

write.table(res_mat_sig,
            file.path(data_dir, paste0(prefix, '_igraph_predictions_sig.txt')),
            sep = '\t', row.names = FALSE)

# ---- database overlap (optional; uses files in data-dir if present) ---------
db_specs <- list(
  CORUM    = list(file = corum_file,                                        sep = '\t', header = TRUE),
  STRING   = list(file = file.path(data_dir, "stringdb_400.txt"),           sep = '\t', header = TRUE),
  Reactome = list(file = file.path(data_dir, "reactome_gene_symbols.txt"),  sep = '\t', header = TRUE),
  Signor   = list(file = file.path(data_dir, "signorDB.txt"),               sep = '\t', header = TRUE),
  BioPlex  = list(file = file.path(data_dir, "Bioplex_combined.txt"),       sep = '\t', header = TRUE),
  HuRI     = list(file = file.path(data_dir, "HuRI.txt"),                   sep = '\t', header = TRUE),
  FunMap   = list(file = file.path(data_dir, "FunMap.tsv"),                 sep = '\t', header = FALSE)
)

pairs_list <- Filter(Negate(is.null), lapply(names(db_specs), function(nm) {
  spec <- db_specs[[nm]]
  if (!file.exists(spec$file)) return(NULL)
  db <- read.delim(spec$file, sep = spec$sep, header = spec$header)
  if (!spec$header) colnames(db)[1:2] <- c("prot1", "prot2")
  db
}))
names(pairs_list) <- names(db_specs)[sapply(db_specs, function(s) file.exists(s$file))]

if (length(pairs_list) > 0) {
  all_genes <- union(res_mat_sig$Protein1, res_mat_sig$Protein2)

  pair_dat_overlap <- sapply(pairs_list, function(k) {
    p1 <- k[[1]]; p2 <- k[[2]]
    shared   <- intersect(all_genes, unique(c(p1, p2)))
    k_pairs  <- paste0(p1[p1 %in% shared | p2 %in% shared], '_', p2[p1 %in% shared | p2 %in% shared])
    sig_sub  <- res_mat_sig[(res_mat_sig$Protein1 %in% shared) | (res_mat_sig$Protein2 %in% shared), ]
    sig_pairs <- paste0(sig_sub$Protein1, '_', sig_sub$Protein2)
    length(which(k_pairs %in% sig_pairs)) / length(k_pairs)
  })

  pair_dt <- data.frame(origin = names(pair_dat_overlap), overlap = pair_dat_overlap)
  pair_dt$origin <- factor(pair_dt$origin, levels = rev(names(pairs_list)))

  p_overlap <- ggplot(pair_dt, aes(x = origin, y = overlap * 100, fill = origin)) +
    geom_bar(stat = 'identity') +
    geom_text(mapping = aes(label = round(overlap * 100, 1)), hjust = -0.1) +
    labs(x = '', y = 'Overlap (%)', title = paste0(prefix, ' network - database overlap')) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0.1)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_flip() +
    guides(fill = 'none')

  ggsave(file.path(figures_dir, paste0(prefix, '_db_overlap.pdf')), plot = p_overlap, width = 6, height = 4)
}

# ---- export to Gephi --------------------------------------------------------
final_graph <- connected.component(igraph_net)
nodes_df    <- data.frame(ID = seq_len(vcount(final_graph)), NAME = V(final_graph)$name)
nodes_att   <- nodes_df[, 'ID', drop = FALSE]
nodes_att$Label <- V(final_graph)$name
edges_df    <- as.data.frame(get.edges(final_graph, seq_len(ecount(final_graph))))

write.gexf(nodes = nodes_df, edges = edges_df, nodesAtt = nodes_att,
           defaultedgetype = "undirected",
           output = file.path(figures_dir, paste0(prefix, "_correlation_net.gexf")))

message("network_reconstruction.R done. Network written to ",
        file.path(data_dir, paste0(prefix, '_igraph_predictions_sig.txt')))
