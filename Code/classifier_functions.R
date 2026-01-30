### Classifier functions

### Libraries
library(tidyverse)
library(rms)
library(patchwork)
library(caret)
library(ROCR)
library(data.table)
library(umap) 
library(igraph)
library(rgexf)

## Annotate wth CORUM/SCBC2 distributions
annotate.mat <- function(cormat, corum_pairs, y_title) {
  
  # Compare the rows of protein expression files with a vector containing the true pairs of each database. Unique makes you get unique elements from the combination
  common_genes <- intersect(unique(c(corum_pairs[,1], corum_pairs[,2])),
                            unique(c(cormat$Protein1, cormat$Protein2)))
  
  #Subset by row the database and protein expression files based on common genes 
  corum_pairs_sub <- corum_pairs[corum_pairs[,1] %in% common_genes & corum_pairs[,2] %in% common_genes, ]
  cormat$db <-   ifelse(paste0(cormat$Protein1, '_', cormat$Protein2) %in% 
                          paste0(corum_pairs_sub[,1], '_', corum_pairs_sub[,2]), 'TRUE', 'FALSE')
  
  # Make database entries Factors and use Table function 
  cormat$db <- as.factor(cormat$db)
  freq_db <- as.data.frame(table(cormat$db))
  colnames(freq_db)[1] <- 'db' 
  
  cormat$loc1 <- corum_pairs[match(cormat$Protein1, corum_pairs$V1), 4]
  cormat$loc2 <- corum_pairs[match(cormat$Protein2, corum_pairs$V2), 5]
  # cor_mat$comp1 <- corum_pairs$complex[match(cor_mat$Var1, corum_pairs$V1)]
  # cor_mat$comp2 <- corum_pairs$complex[match(cor_mat$Var2, corum_pairs$V2)]
  
  cormat$same_loc <- cormat$loc1 == cormat$loc2
  
  # cor_mat$same_comp <- cor_mat$comp1 == cor_mat$comp2
  
  cormat$same_loc <- cormat$same_loc & cormat$db == 'FALSE'
  
  # cor_mat$same_loc_same_comp <- cor_mat$same_loc & cor_mat$same_comp
  
  cormat$cat <- ifelse(cormat$same_loc, cormat$loc1, ifelse(cormat$db == 'TRUE', 'Same Complex', 'Random'))
  cormat$cat <- factor(cormat$cat, levels = c( 'Random','Mitochondria', 'Secretory', 'Cytosol', 'Nuclear','Same Complex'))
  cormat$cat2 <- ifelse(cormat$cat %in% c('Mitochondria', 'Secretory', 'Cytosol', 'Nuclear'), 'Same location', as.character(cormat$cat))
  
  freq_cat <- cormat[!is.na(cormat$cat), ] %>% group_by(cat, cat2) %>%
    summarise(n = n())
  
  colnames(freq_cat)[1] <- 'cat' 
  
  median_cat <- cormat[!is.na(cormat$cat), ] %>% group_by(cat, cat2) %>%
    summarise(m = median(sim, na.rm = TRUE))
  
  p_cor_distr <- ggplot(cormat[!is.na(cormat$cat), ], aes(x = cat, y = sim)) +
    geom_violin(alpha=0.4) +
    geom_text(data = freq_cat, aes(label = paste0(n), y = -1), vjust = 0.2) +
    geom_text(data = median_cat, aes(label = paste0(round(m, 2)), y = 1.1), vjust = 0.2) +
    geom_boxplot(width = 0.2, outlier.color = NA) + 
    # scale_x_discrete(labels = c('Random', levels(cor_mat$cat)[-c(1,6)], 'Same_complex')) + 
    labs(x = '', y = y_title) +
    labs(fill = "Class") +
    facet_grid(.~cat2, space = 'free', scales = 'free') + 
    theme_bw() + 
    theme(axis.text.x = element_text(hjust = 1, angle = 45), strip.background = element_blank(), 
          strip.text = element_blank(), panel.grid.minor = element_blank()) + 
    guides(fill = 'none')
  
  
  # Remove large complexes
  complex_size <- table(corum_pairs_sub$complex)
  # thres_complex_size <- quantile(complex_size, 0.75)
  thres_complex_size <- 100
  tokeep <- names(complex_size)[complex_size < thres_complex_size]
  
  corum_pairs_sub <- corum_pairs_sub[corum_pairs_sub$complex %in% tokeep, ]
  common_genes <- intersect(unique(c(corum_pairs_sub[,1], corum_pairs_sub[,2])), 
                            unique(c(cormat$Protein1, cormat$Protein2)))
  
  cormat$db[!(cormat$Protein1 %in% common_genes | cormat$Protein2 %in% common_genes)] <- FALSE
  # protein_cormat_sub <- cormat[cormat$Protein1 %in% common_genes & cormat$Protein2 %in% common_genes, ]
  
  ## Subset to same-complex and same-location
  # protein_cormat_sub <- protein_cormat_sub[protein_cormat_sub$cat2 %in% c('Same Complex', 'Same location'), ]
  # table(protein_cormat_sub$cat2)
  ## Create training data
  # protein_cormat_sub$class <- protein_cormat_sub$cat2
  # protein_cormat_sub$class <- ifelse(protein_cormat_sub$class == 'Same Complex', 'TRUE', 'FALSE')
  # protein_cormat_sub$class <- factor(protein_cormat_sub$class, levels = c('FALSE', 'TRUE'))
  # rownames(protein_cormat_sub) <- paste(protein_cormat_sub$Protein1, protein_cormat_sub$Protein2, sep = '_')
  # protein_cormat_sub[, c("Protein1", "Protein2", "db", "loc1", "loc2","same_loc", "cat", "cat2")] <- NULL
  # 
  # table(protein_cormat_sub$class)
  # protein_cormat_sub$class <-  protein_cormat_sub$db
  
  cormat$class <-  cormat$db 
  
  
  return(list('cormat_annot' = cormat, 'plot' = p_cor_distr))
  
}


calc_metrics <- function(conf_matrix){
  conf_ma <- conf_matrix
  
  TP <- conf_ma$table[1,1]
  TN <- conf_ma$table[2,2]
  FP <- conf_ma$table[1,2]
  FN <- conf_ma$table[2,1]
  
  FDR <- (FP/(TP + FP))
  FPR <- (FP/(TN + FP))
  precision <- (TP/(TP + FP))
  recall <- (TP/(TP + FN))
  accuracy = (TP + TN)/(TP + TN + FP + FN)
  return(list('FDR' = FDR, 'FPR' = FPR, 'precision' = precision, 'recall' = recall, 
              'accuracy' = accuracy))
}


## threshold value 
# thres_optim <- function(res_pred_dt, testFold, metric_thres){
#   
#   min = 0
#   max = 1.0
#   m = 0
#   res <- data.frame('m'=0, 'fdr'= 0, 'threshold' = 0, 'n' = 0)
#   
#   while (max - min > 0.0001){
#     m = m + 1
#     # first threshold
#     threshold = (min + max)/2
#     
#     # Convert probabilities to class predictions based on threshold
#     actual_pred <- res_pred_dt$Pred > threshold
#     
#     # Get prediction and groun truth df ready
#     res_pred_dt_new <- data.frame('Pred' = factor(actual_pred,c('TRUE', 'FALSE')),
#                                   'Init' = factor(testFold$class, c('TRUE', 'FALSE')))
#     
#     # Create confusion matrix
#     conf_ma <- confusionMatrix(res_pred_dt_new[["Pred"]], res_pred_dt_new[["Init"]], positive = "TRUE")
#     
#     # Calculate metric
#     current_metric <- calc_metrics(conf_matrix = conf_ma)
#     current_metric <- current_metric$FDR
# 
#     
#     res[m, c('m', 'fdr', 'threshold', 'n')] <- c(m,  current_metric, threshold, conf_ma$table[1,1])
#     
#     # Update bounds based on FDR
#     if(is.nan(current_metric)) {
#       
#       max <- threshold
#       
#     } else if (current_metric > metric_thres) {
#       min <- threshold
#     } else {
#       max <- threshold
#     }
#   }
#   return(res)
#   # return(threshold)
# }

thres_optim <- function(res_pred_dt, testFold){
  # min = 0
  # max = 1.0
  # m = 0
  p_step <-  (max(res_pred_dt$Pred) - min((res_pred_dt$Pred)))/1000
  p_seq <- seq(min((res_pred_dt$Pred)), max(res_pred_dt$Pred), by = p_step)
  res <- data.frame('s'=0, 'fdr'= 0, 'precision'= 0, 'recall' = 0,  'threshold' = 0, 'n' = 0)
  for (s in 1:length(p_seq)) {
    # Convert probabilities to class predictions based on threshold
    threshold <- p_seq[s]
    actual_pred <- res_pred_dt$Pred > threshold
    # Get prediction and groun truth df ready
    res_pred_dt_new <- data.frame('Pred' = factor(actual_pred,c('TRUE', 'FALSE')),
                                  'Init' = factor(testFold$class, c('TRUE', 'FALSE')))
    # Create confusion matrix
    conf_ma <- confusionMatrix(res_pred_dt_new[["Pred"]], res_pred_dt_new[["Init"]], positive = "TRUE")
    # Calculate metric
    current_metric <- calc_metrics(conf_matrix = conf_ma)
    fdr_res <- current_metric$FDR
    precision_res <- current_metric$precision
    recall_res <- current_metric$recall
    # current_metric <- calc_metric_FPR(conf_matrix = conf_ma)
    
    res[s, c('s', 'fdr', 'precision', 'recall', 'threshold', 'n')] <- c(s,  fdr_res, precision_res, recall_res, threshold, conf_ma$table[1,1])
    # Update bounds based on FDR
  }
  return(res)
  # return(threshold)
}


### Logistic regression
cv.fold <- function(cor_mat, max_iterations, sim_col, class_col) {
  # Create equally sized lists
  #why do this 100 times. 
  
  res_list_iter <- lapply(1:max_iterations, function(k) {
    
    print(k)
    set.seed(k)
    idx_true <- sample(which(cor_mat[[class_col]]=='TRUE'), 1000)
    idx_false <- sample(which(cor_mat[[class_col]]=='FALSE'), 1000)
    
    cor_mat_train <- cor_mat[c(idx_true, idx_false), ]
    
    # split to train and validation set
    set.seed(k)
    folds <- createDataPartition(cor_mat_train[[class_col]], times = 1, p = 0.7)
    
    trainFold <- cor_mat_train[folds[[1]], ]
    testFold <- cor_mat_train[setdiff(1:nrow(cor_mat_train), folds[[1]]), ]
    
    # Check parameters
    set.seed(k)
    formula_text <- paste0(class_col, ' ~ ', sim_col)
    
    lr_model <- glm(eval(formula_text), family=binomial(link='logit'), data=trainFold) # logistic regression
    
    # SVM model
    # svm_model <- svm(class ~ ., data = trainFold, probability = TRUE)
    
    # # Prediction of the 600 validation test data relationship between cor and probability of complex
    res_pred <- predict(lr_model, newdata =  testFold,  type="response")
    # res_pred <- predict(svm_model, newdata = testFold,  probability  = TRUE)
    # res_pred <- attr(res_pred, "probabilities")
  
    
    res_pred_dt <- data.frame(
      'Pred' =res_pred,
      # 'Pred' = res_pred[,1],
      'Init' = testFold[[class_col]])
    
    
    # Calibrate 
    # rms::val.prob(res_pred_dt$Pred, as.numeric(res_pred_dt$Init)- 1)
    
    
    p_thres <- thres_optim(res_pred_dt = res_pred_dt, testFold = testFold)
    p_thres$k = k
    
    # Performance metrics - 1st classifier
    pred <- prediction(res_pred_dt$Pred, res_pred_dt$Init)
    perf <- performance(pred, "tpr", "fpr")
    #auc_val <- unlist(slot(performance(pred, "auc"), "y.values"))
    auc <- performance(pred,"auc")@y.values[[1]]
    
    
    roc_data <- data.frame('FPR' = unlist(perf@x.values),
                           'TPR' = unlist(perf@y.values),
                           'Prob' = unlist(perf@alpha.values))
    
    # pdf(paste0('../figures/ROC_curve_iter_', k, '.pdf'), width = 6, height = 5)
    # plot(perf, colorize=FALSE, main = 'ROC curve, 1000/1000')
    # dev.off()
    
    list('auc' = auc,
         'thres_p' = p_thres, 
         'model' = lr_model
         # 'model' = svm_model
    ) 
    
  })
  
  # if(length(res_list_iter) == 0) {
  #   
  #   return(NULL)
  #   
  # } else {
  
  # toremove <- unlist(lapply(res_list_iter, function(j) j$thres_p)) == "cannot reach specified metric threshold"
  # res_list_iter <- res_list_iter[!toremove]
  
  
  all_auc <- data.frame('iter' = 1:length(res_list_iter), 
                        'auc' = unlist(lapply(res_list_iter, function(j) j$auc)))
  
  # Median cor. across iterations
  median_auc <- round(median(all_auc$auc), 3)
  
  
  # # Plot
  # p_auc <- ggplot(all_auc, aes(x = '', y = auc)) +
  #   geom_boxplot() +
  #   geom_jitter(width = 0.1) +
  #   annotate(x = 1, y = median_auc, geom = 'text', label = paste0('m = ', median_auc), vjust = -3) +
  #   scale_y_continuous(limits = c(0, 1)) +
  #   labs(x = '# Iterations', y = 'AUC') +
  #   theme_classic()
  # 
  # pdf(paste0('../figures/pp_auc.pdf'), width = 4.5, height = 6)
  # plot(p_auc)
  # dev.off()
  
  # # Extract thresholds 
  # all_thres <- data.frame('iter' = 1:length(res_list_iter), 
  #                         'thres' = unlist(lapply(res_list_iter, function(j) j$thres_p)))
  # 
  # all_thres$thres <- as.numeric(all_thres$thres)
  # median_thres <- round(median(all_thres$thres), 3)
  
  
  # p_thres <- ggplot(all_thres, aes(x = '', y = thres)) +
  #   geom_boxplot() +
  #   geom_jitter(width = 0.1) +
  #   scale_y_continuous(limits = c(0, 1)) + 
  #   annotate(x = 1, y = median_thres, geom = 'text', label = paste0('m = ', median_thres), vjust = -3) + 
  #   labs(x = '# Iterations', y = 'Probability') +
  #   theme_classic()
  
  # pdf(paste0('../figures/pp_thres_cor.pdf'), width = 4.5, height = 6)
  # plot(p_thres)
  # dev.off()
  
  
  all_thres <- do.call(rbind, lapply(res_list_iter, function(j) j$thres_p))
  
  # all_thres <- all_thres %>% group_by(k) %>%
  #   arrange(fdr, .by_group = TRUE) %>%
  #   slice(1)
  # 
  # # Exclude not fdr attained
  # all_thres_filt <- all_thres[all_thres$fdr < metric_thres, ]
  # all_thres_filt$thres <- as.numeric(all_thres_filt$threshold)
  # median_thres <- round(median(all_thres_filt$thres), 3)
  # 
  # 
  ## Final threshold
  final_thres <- median_auc

  # final_thres_idx <-  all_thres_filt$k[which.min(abs(final_thres - all_thres_filt$threshold))]
  final_thres_idx <-  which.min(abs(final_thres - all_auc$auc))
  final_model <- res_list_iter[[final_thres_idx]]$model
  
  return(list(
              # 'auc' = all_auc[all_auc$iter %in% all_thres_filt$k, ],
              'auc' = all_auc$auc,
              'thres' = all_thres,
              # 'final_thres' = final_thres, 
              'final_model' = final_model))
}


feature.dist <- function(features, feature_name) {
  
  n_size <- length(features)
  if(n_size > 10^6) {
    set.seed(1234)
    idx <- sample(1:n_size, 10^5)
  } else {
    idx <- 1:n_size
  }
  
  dt <- data.frame('f' = features[idx])
  
  p_feature <- ggplot(dt, aes(x = f)) + 
    scale_y_continuous(expand = c(0,0.1)) + 
    geom_histogram(col = 'black', fill = 'white', bins = 100) + 
    labs(x = feature_name, y = 'Count') + 
    theme_classic() 
  
  p_feature
}


## Probability distribution
prob.dist <- function(features, preds, feature_name) {

  n_size <- length(features)
  if(n_size > 10^6) {
    idx <- sample(1:n_size, 10^5)
  } else {
    idx <- 1:n_size
    }
  # preds <- genept_cormat_pred$pred
  # features <- genept_cormat_pred$sim
  pred_dt <- data.frame('p' = preds[idx])
  p_pred <- ggplot(pred_dt, aes(x = p)) + 
    scale_y_continuous(expand = c(0,0.1)) + 
    geom_histogram(col = 'black', fill = 'white', bins = 100) + 
    labs(x = 'probabilities', y = 'Count') + 
    theme_classic() 
  
  p_feature <- feature.dist(features = features, feature_name = feature_name ) + coord_flip()
  
  dt_both <- data.frame('p' = preds[idx], 
                        'f' = features[idx])
  
  p_both <- ggplot(dt_both, aes(x = p, y = f)) + 
    geom_point() + 
    scale_y_continuous(expand = c(0,0.1)) + 
    # geom_histogram(col = 'black', fill = 'white', bins = 100) + 
    labs(x = 'probabilities', y = feature_name) + 
    theme_classic() 
  
  
  p_all <- p_pred + plot_spacer() + p_both + p_feature + plot_layout(nrow = 2, widths = c(0.7, 0.3), heights = c(0.3, 0.7), axis_titles = 'collect')
  
}


                  


## FDR plot
plot.fdr <- function(class_res, titl) {
  
  dt_fdr <- class_res$thres
  
  fdr_thres <- c(0.01, 0.05, 0.1, 0.2, 0.25, 0.3)
  
  
  p_cont <- ggplot(dt_fdr ,aes(x = threshold, y = fdr, size = n, group = k)) + 
    geom_point(alpha = 0.1) + 
    scale_y_continuous(limits = c(0, 1.01))+ 
    # annotate(x = Inf, y = Inf, geom = 'text', label = paste0('n_iter = ', nrow(dt_fdr)),  vjust = 1, hjust = 1) +
    labs(x = 'probability thres.', y = 'FDR', title = titl) + 
    geom_hline(yintercept = fdr_thres, lty = 2) + 
    geom_smooth() + 
    theme_classic() + 
    # theme(
    #   plot.title = element_text(hjust = 0.5))  + 
    guides(size = 'none')
  
  iter <- unique(dt_fdr$k)
  
  # fdr_thres <- range(dt_fdr$fdr, na.rm = TRUE)
  dt_fdr_disc <- as.data.frame(do.call(rbind, lapply(iter, function(k) {
    
    # Largest value to achieve a FDR
    dt_fdr_sub <- dt_fdr[dt_fdr$k == k, ]
    
    res <- do.call(cbind, lapply(fdr_thres, function(l) {
      
      idx <- which.min(abs(dt_fdr_sub$fdr - l))
      dt_fdr_sub[idx, 'threshold']
      
    }))
    
  })))
  colnames(dt_fdr_disc) <- paste0('fdr_', fdr_thres )
  
  p_m <- reshape2::melt(apply(dt_fdr_disc, 2, median, na.rm = TRUE))
  p_m$variable <- rownames(p_m)
  dt_fdr_disc$iter <- iter
  dt_fdr_disc <- reshape2::melt(dt_fdr_disc, id.var = 'iter')
  

  p_disc <- ggplot(dt_fdr_disc ,aes(x = variable, y = value)) + 
    geom_boxplot() + 
    geom_text(data = p_m, mapping = aes(x = variable, y = 1.1, label = round(value, 4)),  size = 2.5) + 
    scale_x_discrete(label = function(i) {gsub('.*_', '', i)}) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0,0.1)) +
    labs(x = 'FDR', y = 'probability thres.') + 
    theme_classic() + 
    coord_flip()
  
  
  p_res <- p_cont + p_disc + plot_annotation(title = paste0(titl)) + plot_layout(axis_titles = 'collect') & theme(plot.title = element_text(hjust = 0.5))
  
  
  
  # scale_y_continuous(expand = c(0, 0, 0.1, 0.1)) + 
}


## AUC
plot.auc <- function(class_res, titl) {
  
  dt_auc <- data.frame('auc' = class_res$auc)
  dt_auc$group <- 'Proteomics_retro'
  p_auc <- ggplot(dt_auc ,aes(x = group, y = auc)) + 
    geom_boxplot(outlier.colour = NA) + 
    geom_jitter(width = 0.2) + 
    scale_y_continuous(limits = c(0, 1.01))+ 
    annotate(x = Inf, y = round(median(dt_auc$auc), 3), geom = 'text', label = paste0('m = ', c(round(median(dt_auc$auc), 3))), vjust = 1) +
    annotate(x = Inf, y = -Inf, geom = 'text', label = paste0('n_iter = ', nrow(dt_auc)),  vjust = 1, hjust = -1) +
    labs(x = '', y = 'AUC', title =  titl) + 
    theme_classic() + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          plot.title = element_text(hjust = 0.5)) + 
    # scale_y_continuous(expand = c(0, 0, 0.1, 0.1)) + 
    coord_flip()
}


# Precision - recall
plot.pres.recall <- function(class_res, titl) {
  
  dt <- class_res$thres

  p_pres <- ggplot(dt ,aes(x = threshold, y = precision, size = n, group = k)) + 
    geom_point(alpha = 0.1) + 
    scale_y_continuous(limits = c(0, 1.01))+ 
    # annotate(x = Inf, y = Inf, geom = 'text', label = paste0('n_iter = ', nrow(dt_fdr)),  vjust = 1, hjust = 1) +
    labs(x = 'probability thres.', y = 'Precision', title = titl) + 
    geom_smooth() + 
    theme_classic() + 
    # theme(
    #   plot.title = element_text(hjust = 0.5))  + 
    guides(size = 'none')
  
  
  p_rec <- ggplot(dt ,aes(x = threshold, y = recall, size = n, group = k)) + 
    geom_point(alpha = 0.1) + 
    scale_y_continuous(limits = c(0, 1.01))+ 
    # annotate(x = Inf, y = Inf, geom = 'text', label = paste0('n_iter = ', nrow(dt_fdr)),  vjust = 1, hjust = 1) +
    labs(x = 'probability thres.', y = 'Recall') + 
    geom_smooth() + 
    theme_classic() + 
    # theme(
    #   plot.title = element_text(hjust = 0.5))  + 
    guides(size = 'none')
  
  
  p_res <- p_pres/p_rec + plot_annotation(title = paste0(titl)) + plot_layout(axis_titles = 'collect') & theme(plot.title = element_text(hjust = 0.5))
  
  
  
}

# Calculate probabilities for all data
generate.pred <- function(class_res, cor_mat) {
  
  rownames(cor_mat) <- paste(cor_mat$Protein1, cor_mat$Protein2, sep = '_')
  cor_mat_pred <- cor_mat
  cor_mat_pred[, c("Protein1", "Protein2", "db", "loc1", "loc2", "same_loc", "cat", "cat2") ] <- NULL
  final_model <- class_res$final_model
  # res_pred <- predict(final_model, newdata =  cor_mat,  probability = TRUE)
  res_pred <- predict(final_model, newdata =  cor_mat_pred, type="response")
  
  # res_pred <- attr(res_pred, "probabilities")
  # cor_mat$pred <-  res_pred[,1]
  cor_mat_pred$pred <-  as.numeric(res_pred)
  
  # final_thres <-  class_res$final_thres
  # cor_mat_pred$pred_thes <- cor_mat_pred$pred > final_thres
  
  # # Calculate threshold for negative correlations
  # thres_cor <- min(cor_mat_all$value[cor_mat_all$pred_thes])
  # cor_mat_all$cor_thes <- cor_mat_all$pred_thes
  # cor_mat_all$cor_thes[abs(cor_mat_all$value) > thres_cor] <- TRUE
  # 
  # cor_mat_all_sig <- cor_mat_all[cor_mat_all$cor_thes, ]
  # cor_mat_all_pos <- cor_mat_all_sig[cor_mat_all_sig$pred_thes, ]
  # 
  # rm(cor_mat_all)
  # gc()
  
  # Reformat
  cor_mat_pred[,  c('Protein1', 'Protein2')] <- do.call(rbind, strsplit(rownames(cor_mat_pred), '_'))
  
  res_mat <- cbind(cor_mat_pred[,   c('Protein1', 'Protein2', 'sim',  "pred")],
                   cor_mat[, c("class", "db", "loc1", "loc2", "same_loc", "cat", "cat2")])
  return(res_mat)
}


connected.component <- function(i) {
  
  comps <- igraph::components(i)
  max_comp <- which.max(comps$csize)
  tokeep <- names(comps$membership)[comps$membership == max_comp]
  res <- subgraph(i, vids = tokeep)
  
}



degree.scatter <- function(g, thres) {
  
  deg <- igraph::degree(g)
  deg_freq <- table(deg)
  
  log_deg_freq <- as.numeric(log10(deg_freq))
  log_deg <- log10(as.numeric(names(deg_freq)))
  
  fit <- lm(log_deg_freq ~ log_deg)
  r_squared <- summary(fit)
  r_squared <- r_squared$r.squared
  
  ## Exponential fit
  # Fit the exponential model
  # fit <- nls(log_deg_freq ~ a -exp(b * log_deg), start = list(a = 0.1, b = -0.1))
  # y_pred <- predict(fit)
  # 
  # # Calculate R-squared
  # y_observed <- log_deg_freq
  # y_mean <- mean(y_observed)
  # 
  # # Residual sum of squares
  # RSS <- sum((y_observed - y_pred)^2)
  # 
  # # Total sum of squares
  # TSS <- sum((y_observed - y_mean)^2)
  # 
  # # R-squared
  # r_squared_exp <- 1 - (RSS / TSS)
  # 
  # # Extract the fitted parameters
  # params <- coef(fit)
  # 
  # # Create a data frame with fitted values for plotting
  # dt_fitted <- data.frame(
  #   x = log_deg,
  #   y = params["a"] - exp(params["b"] * log_deg),
  #   type = 'exp'
  # )
  # 
  
  dt_log <- data.frame(x = log_deg, 
                       y = log_deg_freq, 
                       type = 'linear')
  
  dt_all <- rbind(
    # dt_fitted, 
    dt_log)
  
  p <- ggplot() + 
    geom_point(dt_all[dt_all$type == 'linear', ], mapping = aes(x = x, y = y), col = 'black', show.legend = FALSE) + 
    geom_smooth(dt_all[dt_all$type == 'linear', ], mapping = aes(x = x, y = y, col = type), method = 'lm', se = FALSE) + 
    # geom_point(dt_all[dt_all$type == 'exp', ], mapping = aes(x = x, y = y, col = type)) + 
    # geom_smooth(dt_all[dt_all$type == 'exp', ], mapping = aes(x = x, y = y, col = type), method = 'loess') + 
    scale_color_manual(name = 'Fit', values = c('exp' = 'indianred', 'linear' = 'steelblue'))+ 
    # annotate(geom = 'text', x = Inf, y = Inf, label = paste0('R^2_exp = ', round(r_squared_exp, 2)), hjust = 1, vjust = 2, size = 5, col = 'indianred')+ 
    annotate(geom = 'text', x = Inf, y = Inf, label = paste0('R^2_linear = ', round(r_squared, 2)), hjust = 1, vjust = 4, size = 5, col =  'steelblue')+ 
    
    labs(x = 'Node degree - log10', y = 'Number of nodes - log10', title = paste0('(nodes: ', length(V(g)), ', edges: ', length(E(g)), ')'), 
         subtitle = paste0('thres: ', thres) )+ 
    theme_classic() + 
    theme(text = element_text(size = 10), 
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5)) + 
    # guides(col = guide_legend(override.aes = list(size = 3, shape = 15, col = c('indianred', 'steelblue'), fill = 'white')))
    guides(col = 'none')
  
  
  list(p = p,
       r_squared = r_squared)
}



## Degree distribution
degree.hist <- function(g, thres) {
  
  degree_vec <- igraph::degree(g)
  degree_dt <- data.frame(gene = names(degree_vec), 'degree' = as.numeric(degree_vec))
  
  m_degree <- round(mean(degree_vec),0)
  
  p_degree <- ggplot(degree_dt, aes(x = degree)) + 
    geom_histogram(colour = 'black') + 
    geom_vline(xintercept = m_degree, color = 'red', lty = 2) + 
    annotate(geom = 'text', x = m_degree, y = Inf, label = paste0('avg.degree = ', m_degree), color = 'red', hjust = -0.1, vjust = 1)  +
    theme_classic() + 
    labs(x = 'Degree', y = '# Nodes',  subtitle = paste0('thres: ', thres)) + 
    scale_x_continuous(expand = c(0, 0.1)) + 
    scale_y_continuous(expand = c(0, 0.1)) + 
    theme( plot.subtitle = element_text(hjust = 0.5))
  
  
  list('p' = p_degree, 
       'avg_d' = m_degree)
  

}



net.size <- function(prob_mat,pred_col, min_thres, max_thres, step_thres) {
  
  prob_step <- (max_thres - min_thres)/(step_thres -1)
  p_seq <- seq(min_thres, max_thres, by = prob_step)
  
  res <- do.call(rbind, lapply(p_seq, function(i) {
    
    # i = p_seq[1]
    print(i)
    prob_mat_sub <- prob_mat[prob_mat[[pred_col]] > i, ]
    prob_mat_sub_graph <- graph_from_data_frame(prob_mat_sub[, c('Protein1', 'Protein2')], directed = FALSE)
    prob_mat_sub_graph <- connected.component(prob_mat_sub_graph)
    
    if(length(V(prob_mat_sub_graph)) == 0) {
      
      r_squared <- NA
      avg_d <- NA
      cluster_coef <- NA
      modularity_val <- NA
      avg_path <- NA
      
    } else{
      degree_res <- degree.scatter(g = prob_mat_sub_graph, thres = i)
      r_squared <- degree_res$r_squared
      degree_hist <- degree.hist(g = prob_mat_sub_graph, thres = i)
      avg_d <- degree_hist$avg_d
      
      cluster_coef <- transitivity(graph = prob_mat_sub_graph, type = 'average')
      cluster_leiden <- cluster_leiden(graph = prob_mat_sub_graph, weights = E(prob_mat_sub_graph)$pred, objective_function = 'modularity')
      # cluster_louvain <- cluster_louvain(graph = prob_mat_sub_graph, weights = NULL)
      modularity_val <- modularity(x = prob_mat_sub_graph, membership = cluster_leiden$membership)
      avg_path <- mean_distance(graph = prob_mat_sub_graph, weights = NULL,  directed = FALSE, unconnected = TRUE,details = FALSE)
      
    }
    
    
    
    # prob_mat_sub_graph <- connected.component(prob_mat_sub_graph)
    
    data.frame('n_nodes' = length(V(prob_mat_sub_graph)), 
               'n_edges' =length(E(prob_mat_sub_graph)), 
               'avg_d' = avg_d, 
               'r' = r_squared,
               'cluster_coef' = cluster_coef,
               'modularity' = modularity_val,
               'avg_path' = avg_path,
               'thres' = i)
    
  }))
  
  
  p_nodes <- ggplot(res, aes(x = thres, y = n_nodes)) + 
    geom_point() + 
    geom_line() + 
    theme_classic() + 
    labs(x = 'probability thres.', y = '# nodes')
  
  p_edges <- ggplot(res, aes(x = thres, y = n_edges)) + 
    geom_point() + 
    geom_line() + 
    theme_classic() + 
    labs(x = 'probability thres.', y = '# edges')
  
  p_avg_d <- ggplot(res, aes(x = thres, y = avg_d)) + 
    geom_point() + 
    geom_line() + 
    theme_classic() + 
    labs(x = 'probability thres.', y = 'Average degree') 
  
  p_r_squared <- ggplot(res, aes(x = thres, y = r)) + 
    geom_point() + 
    geom_line() + 
    theme_classic() + 
    labs(x = 'probability thres.', y = 'R^2') 
  
  p_cluster_coef <- ggplot(res, aes(x = thres, y = cluster_coef)) + 
    geom_point() + 
    geom_line() + 
    theme_classic() + 
    labs(x = 'probability thres.', y = 'Average clustering coef.') 
  
  p_module <- ggplot(res, aes(x = thres, y = modularity)) + 
    geom_point() + 
    geom_line() + 
    theme_classic() + 
    labs(x = 'probability thres.', y = 'Modularity') 
  
  p_path <- ggplot(res, aes(x = thres, y = avg_path)) + 
    geom_point() + 
    geom_line() + 
    theme_classic() + 
    labs(x = 'probability thres.', y = 'Average shortest path') 
  
  
  
  p_all <- p_nodes/p_edges/p_avg_d/p_r_squared/p_cluster_coef/p_module/p_path
  
  p_all + plot_layout(axes = 'collect')
  
  
  
}


# reduce.size <-  function(n_size) {
#   
#   if(n_size > 10^6) {
#     set.seed(1234)
#     idx <- sample(1:n_size, 10^5)
#   } else {
#     idx <- 1:n_size
#   }
# }

bootstap.enrich <- function(cormat_pred, links_db) {
  
  min_thres <- min(cormat_pred$pred)
  max_thres <-max(cormat_pred$pred)
  prob_step <- (max_thres - min_thres)/(100 -1)
  p_seq <- seq(min_thres, max_thres, by = prob_step)
  
  links_db$db_id <- paste(links_db$gene1, links_db$gene2, sep = '_')
  
  # common_links <- intersect(links_db$db_id, cormat_pred$V1)
  # links_db <- links_db[links_db$db_id %in% common_links, ]
  # links_db <- links_db[!duplicated(links_db$db_id), ]
  # 
  # cormat_pred_sub <- cormat_pred[cormat_pred$V1 %in% links_db$db_id]
  
  
  
  res <- do.call(rbind, lapply(p_seq, function(i) {
    
    
    print(i)
    
    all_idx <- nrow(cormat_pred)
    
    dt <- do.call(rbind, lapply(1:100, function(j) {
      
      set.seed(j)
      idx <- sample(1:all_idx, 100000)
      
      cormat_pred_sub <- cormat_pred[idx]
      cormat_pred_in <- cormat_pred_sub[pred >= i]
      cormat_pred_out <- cormat_pred_sub[pred < i]
      
      in_db <- length(intersect(cormat_pred_in$V1, links_db$db_id))
      in_nondb <-  length(setdiff(cormat_pred_in$V1, db_id))
      out_db <-  length(intersect(cormat_pred_out$V1, db_id))
      out_nondb <-  length(setdiff(cormat_pred_out$V1, db_id))
      
      mat <- matrix(data = c(in_db, in_nondb, out_db, out_nondb), nrow = 2, ncol = 2, byrow = TRUE)
      res <- fisher.test(mat)
      
      data.frame('iter' = j,
                 'p' = i, 
                 'odds' = res$estimate,
                 'p-value' = res$p.value)
      
    }))
    
  }))
  
  return(res)
}


###