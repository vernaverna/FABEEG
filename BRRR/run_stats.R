###########################################
# Script for run statistics (10-fold CV)  #
###########################################
setwd("/projects/FABEEG/BRRR/")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("cvms")
library("RColorBrewer")
library("corrplot")
library("viridis")
library("rstatix")

# get the CV results and perform k-fold cross-validated paired t-test
# 


compute_10fold_t_test <- function(population, psd_seq){
  
  fname <- paste0('results/10foldCV/NEW_K30_', population, '_', psd_seq, '.RData')
  load(fname)
  
  CV_models <- lapply(CV_results, `[`, 3)
  
  CV_scores <- lapply(CV_results, `[[`, 1) #unlisting 
  CV_scores <- lapply(CV_scores, `[[`, 2)
  
  CV_null_scores <- lapply(CV_results, `[[`, 4)
  CV_null_scores <- lapply(CV_null_scores, `[[`, 2)
  
  
  accs <- unlist(lapply(CV_scores, `[[`, 1))
  print("Average CV accuracy:")
  print(mean(accs))
  
  null_accs <- unlist(lapply(CV_null_scores, `[[`, 1))
  print("Average CV 0-accuracy:")
  print(mean(null_accs))
  
  
  #define difference in accuracies
  p <- unlist(CV_scores) - unlist(CV_null_scores)
  p_m <- mean(p)
  dev_score <- sqrt(sum((p - p_m)**2) / 9)
  
  t <- p_m*sqrt(10) / dev_score
  p_val <- pt(q=t, df=9, lower.tail=FALSE)
  
  print("T-test results:")
  print(paste("T:", t))
  print(paste("p-value:", p_val))
  #crit_value <- qt(p=0.025, df=9)
  
  return(list(t, p_val))
}

population <- "all"
psd_seq <- "N1AN1BN2C"

stat_results <- compute_10fold_t_test(population, psd_seq) 


stat_results <- data.frame(
  row_id = c(1:12),
  population = c(rep("all", 6), rep("o7", 6))
)

pos <- c(rep("o7", 6))
psd_seqs <- c("N1AN1BN2C", "N1AN1BN2AN2C", "N1AN2BN2C", "N1AN2BN2AN1B",
              "N2AN2BN2C", "N2AN2BN2CN2D")
#psd_seqs <- c("N1AN1BN2C", "N1AN1BN2AN2C", "N1AN2BN2C", "N1AN2BN1B", "N1AN2BN2AN1B",
#                             "N2AN2BN2C", "N2AN2BN1B", "N2AN2BN2CN2D")
for(i in 1:6){
  population <- pos[i]
  psd_seq <- psd_seqs[i]
  print(psd_seq)
  res <- compute_10fold_t_test(population, psd_seq) 
}









