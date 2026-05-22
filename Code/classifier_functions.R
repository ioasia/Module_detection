### Classifier and network utility functions
### Used by network_reconstruction.R

library(tidyverse)
library(rms)
library(patchwork)
library(caret)
library(ROCR)
library(data.table)
library(umap)
library(igraph)
library(rgexf)


annotate.mat <- function(cormat, corum_pairs, y_title) {

  common_genes <- intersect(unique(c(corum_pairs[,1], corum_pairs[,2])),
                            unique(c(cormat$Protein1, cormat$Protein2)))

  corum_pairs_sub <- corum_pairs[corum_pairs[,1] %in% common_genes & corum_pairs[,2] %in% common_genes, ]
  cormat$db <- ifelse(paste0(cormat$Protein1, '_', cormat$Protein2) %in%
                        paste0(corum_pairs_sub[,1], '_', corum_pairs_sub[,2]), 'TRUE', 'FALSE')

  cormat$db <- as.factor(cormat$db)
  freq_db <- as.data.frame(table(cormat$db))
  colnames(freq_db)[1] <- 'db'

  cormat$loc1 <- corum_pairs[match(cormat$Protein1, corum_pairs$V1), 4]
  cormat$loc2 <- corum_pairs[match(cormat$Protein2, corum_pairs$V2), 5]

  cormat$same_loc <- cormat$loc1 == cormat$loc2
  cormat$same_loc <- cormat$same_loc & cormat$db == 'FALSE'

  cormat$cat <- ifelse(cormat$same_loc, cormat$loc1, ifelse(cormat$db == 'TRUE', 'Same Complex', 'Random'))
  cormat$cat <- factor(cormat$cat, levels = c('Random', 'Mitochondria', 'Secretory', 'Cytosol', 'Nuclear', 'Same Complex'))
  cormat$cat2 <- ifelse(cormat$cat %in% c('Mitochondria', 'Secretory', 'Cytosol', 'Nuclear'), 'Same location', as.character(cormat$cat))

  freq_cat   <- cormat[!is.na(cormat$cat), ] %>% group_by(cat, cat2) %>% summarise(n = n())
  median_cat <- cormat[!is.na(cormat$cat), ] %>% group_by(cat, cat2) %>% summarise(m = median(sim, na.rm = TRUE))
  colnames(freq_cat)[1] <- 'cat'

  p_cor_distr <- ggplot(cormat[!is.na(cormat$cat), ], aes(x = cat, y = sim)) +
    geom_violin(alpha = 0.4) +
    geom_text(data = freq_cat, aes(label = paste0(n), y = -1), vjust = 0.2) +
    geom_text(data = median_cat, aes(label = paste0(round(m, 2)), y = 1.1), vjust = 0.2) +
    geom_boxplot(width = 0.2, outlier.color = NA) +
    labs(x = '', y = y_title) +
    labs(fill = "Class") +
    facet_grid(.~cat2, space = 'free', scales = 'free') +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, angle = 45), strip.background = element_blank(),
          strip.text = element_blank(), panel.grid.minor = element_blank()) +
    guides(fill = 'none')

  complex_size <- table(corum_pairs_sub$complex)
  thres_complex_size <- 100
  tokeep <- names(complex_size)[complex_size < thres_complex_size]

  corum_pairs_sub <- corum_pairs_sub[corum_pairs_sub$complex %in% tokeep, ]
  common_genes <- intersect(unique(c(corum_pairs_sub[,1], corum_pairs_sub[,2])),
                            unique(c(cormat$Protein1, cormat$Protein2)))

  cormat$db[!(cormat$Protein1 %in% common_genes | cormat$Protein2 %in% common_genes)] <- FALSE
  cormat$class <- cormat$db

  return(list('cormat_annot' = cormat, 'plot' = p_cor_distr))
}


calc_metrics <- function(conf_matrix) {
  conf_ma <- conf_matrix

  TP <- conf_ma$table[1,1]
  TN <- conf_ma$table[2,2]
  FP <- conf_ma$table[1,2]
  FN <- conf_ma$table[2,1]

  FDR       <- (FP / (TP + FP))
  FPR       <- (FP / (TN + FP))
  precision <- (TP / (TP + FP))
  recall    <- (TP / (TP + FN))
  accuracy  <- (TP + TN) / (TP + TN + FP + FN)

  return(list('FDR' = FDR, 'FPR' = FPR, 'precision' = precision,
              'recall' = recall, 'accuracy' = accuracy))
}


thres_optim <- function(res_pred_dt, testFold) {
  p_step <- (max(res_pred_dt$Pred) - min(res_pred_dt$Pred)) / 1000
  p_seq  <- seq(min(res_pred_dt$Pred), max(res_pred_dt$Pred), by = p_step)
  res    <- data.frame('s' = 0, 'fdr' = 0, 'precision' = 0, 'recall' = 0, 'threshold' = 0, 'n' = 0)

  for (s in 1:length(p_seq)) {
    threshold    <- p_seq[s]
    actual_pred  <- res_pred_dt$Pred > threshold
    res_pred_dt_new <- data.frame('Pred' = factor(actual_pred, c('TRUE', 'FALSE')),
                                  'Init' = factor(testFold$class, c('TRUE', 'FALSE')))
    conf_ma      <- confusionMatrix(res_pred_dt_new[["Pred"]], res_pred_dt_new[["Init"]], positive = "TRUE")
    current_metric <- calc_metrics(conf_matrix = conf_ma)
    res[s, c('s', 'fdr', 'precision', 'recall', 'threshold', 'n')] <-
      c(s, current_metric$FDR, current_metric$precision, current_metric$recall, threshold, conf_ma$table[1,1])
  }
  return(res)
}


cv.fold <- function(cor_mat, max_iterations, sim_col, class_col) {
  res_list_iter <- lapply(1:max_iterations, function(k) {
    print(k)
    set.seed(k)
    idx_true  <- sample(which(cor_mat[[class_col]] == 'TRUE'), 1000)
    idx_false <- sample(which(cor_mat[[class_col]] == 'FALSE'), 1000)
    cor_mat_train <- cor_mat[c(idx_true, idx_false), ]

    set.seed(k)
    folds      <- createDataPartition(cor_mat_train[[class_col]], times = 1, p = 0.7)
    trainFold  <- cor_mat_train[folds[[1]], ]
    testFold   <- cor_mat_train[setdiff(1:nrow(cor_mat_train), folds[[1]]), ]

    set.seed(k)
    formula_text <- paste0(class_col, ' ~ ', sim_col)
    lr_model  <- glm(eval(formula_text), family = binomial(link = 'logit'), data = trainFold)
    res_pred  <- predict(lr_model, newdata = testFold, type = "response")
    res_pred_dt <- data.frame('Pred' = res_pred, 'Init' = testFold[[class_col]])

    p_thres <- thres_optim(res_pred_dt = res_pred_dt, testFold = testFold)
    p_thres$k <- k

    pred <- prediction(res_pred_dt$Pred, res_pred_dt$Init)
    perf <- performance(pred, "tpr", "fpr")
    auc  <- performance(pred, "auc")@y.values[[1]]

    list('auc' = auc, 'thres_p' = p_thres, 'model' = lr_model)
  })

  all_auc   <- data.frame('iter' = 1:length(res_list_iter),
                          'auc'  = unlist(lapply(res_list_iter, function(j) j$auc)))
  median_auc <- round(median(all_auc$auc), 3)

  all_thres <- do.call(rbind, lapply(res_list_iter, function(j) j$thres_p))

  final_thres_idx <- which.min(abs(median_auc - all_auc$auc))
  final_model <- res_list_iter[[final_thres_idx]]$model

  return(list('auc' = all_auc$auc, 'thres' = all_thres, 'final_model' = final_model))
}


feature.dist <- function(features, feature_name) {
  n_size <- length(features)
  idx <- if (n_size > 10^6) { set.seed(1234); sample(1:n_size, 10^5) } else 1:n_size
  dt  <- data.frame('f' = features[idx])
  ggplot(dt, aes(x = f)) +
    scale_y_continuous(expand = c(0, 0.1)) +
    geom_histogram(col = 'black', fill = 'white', bins = 100) +
    labs(x = feature_name, y = 'Count') +
    theme_classic()
}


prob.dist <- function(features, preds, feature_name) {
  n_size <- length(features)
  idx    <- if (n_size > 10^6) sample(1:n_size, 10^5) else 1:n_size

  pred_dt <- data.frame('p' = preds[idx])
  p_pred  <- ggplot(pred_dt, aes(x = p)) +
    scale_y_continuous(expand = c(0, 0.1)) +
    geom_histogram(col = 'black', fill = 'white', bins = 100) +
    labs(x = 'probabilities', y = 'Count') +
    theme_classic()

  p_feature <- feature.dist(features = features, feature_name = feature_name) + coord_flip()

  dt_both <- data.frame('p' = preds[idx], 'f' = features[idx])
  p_both  <- ggplot(dt_both, aes(x = p, y = f)) +
    geom_point() +
    scale_y_continuous(expand = c(0, 0.1)) +
    labs(x = 'probabilities', y = feature_name) +
    theme_classic()

  p_pred + plot_spacer() + p_both + p_feature +
    plot_layout(nrow = 2, widths = c(0.7, 0.3), heights = c(0.3, 0.7), axis_titles = 'collect')
}


plot.fdr <- function(class_res, titl) {
  dt_fdr    <- class_res$thres
  fdr_thres <- c(0.01, 0.05, 0.1, 0.2, 0.25, 0.3)

  p_cont <- ggplot(dt_fdr, aes(x = threshold, y = fdr, size = n, group = k)) +
    geom_point(alpha = 0.1) +
    scale_y_continuous(limits = c(0, 1.01)) +
    labs(x = 'probability thres.', y = 'FDR', title = titl) +
    geom_hline(yintercept = fdr_thres, lty = 2) +
    geom_smooth() +
    theme_classic() +
    guides(size = 'none')

  iter <- unique(dt_fdr$k)
  dt_fdr_disc <- as.data.frame(do.call(rbind, lapply(iter, function(k) {
    dt_fdr_sub <- dt_fdr[dt_fdr$k == k, ]
    do.call(cbind, lapply(fdr_thres, function(l) {
      dt_fdr_sub[which.min(abs(dt_fdr_sub$fdr - l)), 'threshold']
    }))
  })))
  colnames(dt_fdr_disc) <- paste0('fdr_', fdr_thres)

  p_m <- reshape2::melt(apply(dt_fdr_disc, 2, median, na.rm = TRUE))
  p_m$variable <- rownames(p_m)
  dt_fdr_disc$iter <- iter
  dt_fdr_disc <- reshape2::melt(dt_fdr_disc, id.var = 'iter')

  p_disc <- ggplot(dt_fdr_disc, aes(x = variable, y = value)) +
    geom_boxplot() +
    geom_text(data = p_m, mapping = aes(x = variable, y = 1.1, label = round(value, 4)), size = 2.5) +
    scale_x_discrete(label = function(i) gsub('.*_', '', i)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0.1)) +
    labs(x = 'FDR', y = 'probability thres.') +
    theme_classic() +
    coord_flip()

  p_cont + p_disc + plot_annotation(title = paste0(titl)) +
    plot_layout(axis_titles = 'collect') & theme(plot.title = element_text(hjust = 0.5))
}


plot.auc <- function(class_res, titl) {
  dt_auc <- data.frame('auc' = class_res$auc, 'group' = 'Proteomics')
  ggplot(dt_auc, aes(x = group, y = auc)) +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(width = 0.2) +
    scale_y_continuous(limits = c(0, 1.01)) +
    annotate(x = Inf, y = round(median(dt_auc$auc), 3), geom = 'text',
             label = paste0('m = ', round(median(dt_auc$auc), 3)), vjust = 1) +
    annotate(x = Inf, y = -Inf, geom = 'text',
             label = paste0('n_iter = ', nrow(dt_auc)), vjust = 1, hjust = -1) +
    labs(x = '', y = 'AUC', title = titl) +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    coord_flip()
}


plot.pres.recall <- function(class_res, titl) {
  dt <- class_res$thres

  p_pres <- ggplot(dt, aes(x = threshold, y = precision, size = n, group = k)) +
    geom_point(alpha = 0.1) +
    scale_y_continuous(limits = c(0, 1.01)) +
    labs(x = 'probability thres.', y = 'Precision', title = titl) +
    geom_smooth() +
    theme_classic() +
    guides(size = 'none')

  p_rec <- ggplot(dt, aes(x = threshold, y = recall, size = n, group = k)) +
    geom_point(alpha = 0.1) +
    scale_y_continuous(limits = c(0, 1.01)) +
    labs(x = 'probability thres.', y = 'Recall') +
    geom_smooth() +
    theme_classic() +
    guides(size = 'none')

  p_pres / p_rec + plot_annotation(title = paste0(titl)) +
    plot_layout(axis_titles = 'collect') & theme(plot.title = element_text(hjust = 0.5))
}


generate.pred <- function(class_res, cor_mat) {
  rownames(cor_mat) <- paste(cor_mat$Protein1, cor_mat$Protein2, sep = '_')
  cor_mat_pred <- cor_mat
  cor_mat_pred[, c("Protein1", "Protein2", "db", "loc1", "loc2", "same_loc", "cat", "cat2")] <- NULL

  res_pred <- predict(class_res$final_model, newdata = cor_mat_pred, type = "response")
  cor_mat_pred$pred <- as.numeric(res_pred)

  cor_mat_pred[, c('Protein1', 'Protein2')] <- do.call(rbind, strsplit(rownames(cor_mat_pred), '_'))

  cbind(cor_mat_pred[, c('Protein1', 'Protein2', 'sim', 'pred')],
        cor_mat[, c("class", "db", "loc1", "loc2", "same_loc", "cat", "cat2")])
}


connected.component <- function(i) {
  comps    <- igraph::components(i)
  max_comp <- which.max(comps$csize)
  tokeep   <- names(comps$membership)[comps$membership == max_comp]
  subgraph(i, vids = tokeep)
}


degree.scatter <- function(g, thres) {
  deg          <- igraph::degree(g)
  deg_freq     <- table(deg)
  log_deg_freq <- as.numeric(log10(deg_freq))
  log_deg      <- log10(as.numeric(names(deg_freq)))

  fit       <- lm(log_deg_freq ~ log_deg)
  r_squared <- summary(fit)$r.squared

  dt_log <- data.frame(x = log_deg, y = log_deg_freq, type = 'linear')

  p <- ggplot() +
    geom_point(dt_log, mapping = aes(x = x, y = y), col = 'black') +
    geom_smooth(dt_log, mapping = aes(x = x, y = y, col = type), method = 'lm', se = FALSE) +
    scale_color_manual(name = 'Fit', values = c('linear' = 'steelblue')) +
    annotate(geom = 'text', x = Inf, y = Inf,
             label = paste0('R^2_linear = ', round(r_squared, 2)), hjust = 1, vjust = 4, size = 5, col = 'steelblue') +
    labs(x = 'Node degree - log10', y = 'Number of nodes - log10',
         title = paste0('(nodes: ', length(V(g)), ', edges: ', length(E(g)), ')'),
         subtitle = paste0('thres: ', thres)) +
    theme_classic() +
    theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    guides(col = 'none')

  list(p = p, r_squared = r_squared)
}


degree.hist <- function(g, thres) {
  degree_vec <- igraph::degree(g)
  degree_dt  <- data.frame(gene = names(degree_vec), 'degree' = as.numeric(degree_vec))
  m_degree   <- round(mean(degree_vec), 0)

  p_degree <- ggplot(degree_dt, aes(x = degree)) +
    geom_histogram(colour = 'black') +
    geom_vline(xintercept = m_degree, color = 'red', lty = 2) +
    annotate(geom = 'text', x = m_degree, y = Inf,
             label = paste0('avg.degree = ', m_degree), color = 'red', hjust = -0.1, vjust = 1) +
    theme_classic() +
    labs(x = 'Degree', y = '# Nodes', subtitle = paste0('thres: ', thres)) +
    scale_x_continuous(expand = c(0, 0.1)) +
    scale_y_continuous(expand = c(0, 0.1)) +
    theme(plot.subtitle = element_text(hjust = 0.5))

  list('p' = p_degree, 'avg_d' = m_degree)
}


net.size <- function(prob_mat, pred_col, min_thres, max_thres, step_thres) {
  prob_step <- (max_thres - min_thres) / (step_thres - 1)
  p_seq     <- seq(min_thres, max_thres, by = prob_step)

  res <- do.call(rbind, lapply(p_seq, function(i) {
    print(i)
    prob_mat_sub       <- prob_mat[prob_mat[[pred_col]] > i, ]
    prob_mat_sub_graph <- graph_from_data_frame(prob_mat_sub[, c('Protein1', 'Protein2')], directed = FALSE)
    prob_mat_sub_graph <- connected.component(prob_mat_sub_graph)

    if (length(V(prob_mat_sub_graph)) == 0) {
      return(data.frame('n_nodes' = 0, 'n_edges' = 0, 'avg_d' = NA, 'r' = NA,
                        'cluster_coef' = NA, 'modularity' = NA, 'avg_path' = NA, 'thres' = i))
    }

    degree_res   <- degree.scatter(g = prob_mat_sub_graph, thres = i)
    degree_hist  <- degree.hist(g = prob_mat_sub_graph, thres = i)
    cluster_coef <- transitivity(graph = prob_mat_sub_graph, type = 'average')
    cluster_leiden <- cluster_leiden(graph = prob_mat_sub_graph,
                                     weights = E(prob_mat_sub_graph)$pred,
                                     objective_function = 'modularity')
    modularity_val <- modularity(x = prob_mat_sub_graph, membership = cluster_leiden$membership)
    avg_path <- mean_distance(graph = prob_mat_sub_graph, weights = NULL, directed = FALSE,
                              unconnected = TRUE, details = FALSE)

    data.frame('n_nodes' = length(V(prob_mat_sub_graph)), 'n_edges' = length(E(prob_mat_sub_graph)),
               'avg_d' = degree_hist$avg_d, 'r' = degree_res$r_squared, 'cluster_coef' = cluster_coef,
               'modularity' = modularity_val, 'avg_path' = avg_path, 'thres' = i)
  }))

  p_nodes   <- ggplot(res, aes(x = thres, y = n_nodes))    + geom_point() + geom_line() + theme_classic() + labs(x = 'probability thres.', y = '# nodes')
  p_edges   <- ggplot(res, aes(x = thres, y = n_edges))    + geom_point() + geom_line() + theme_classic() + labs(x = 'probability thres.', y = '# edges')
  p_avg_d   <- ggplot(res, aes(x = thres, y = avg_d))      + geom_point() + geom_line() + theme_classic() + labs(x = 'probability thres.', y = 'Average degree')
  p_r       <- ggplot(res, aes(x = thres, y = r))          + geom_point() + geom_line() + theme_classic() + labs(x = 'probability thres.', y = 'R^2')
  p_clust   <- ggplot(res, aes(x = thres, y = cluster_coef)) + geom_point() + geom_line() + theme_classic() + labs(x = 'probability thres.', y = 'Average clustering coef.')
  p_mod     <- ggplot(res, aes(x = thres, y = modularity)) + geom_point() + geom_line() + theme_classic() + labs(x = 'probability thres.', y = 'Modularity')
  p_path    <- ggplot(res, aes(x = thres, y = avg_path))   + geom_point() + geom_line() + theme_classic() + labs(x = 'probability thres.', y = 'Average shortest path')

  (p_nodes / p_edges / p_avg_d / p_r / p_clust / p_mod / p_path) + plot_layout(axes = 'collect')
}
