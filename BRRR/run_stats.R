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

# idea: get the CV results, and
# 1. apply permutation test to determine the chance-level success rate
# 2. apply permutation test to check if BRRR vs CORR models have significantly 
#    different success rates across the test sets



get_cv_data <- function(population, psd_seq){
  
  fname <- paste0('results/10foldCV/NEW_K30_', population, '_', psd_seq, '.RData')
  load(fname)
  
  CV_models <- lapply(CV_results, `[`, 3)
  #get these from all folds, take the union?
  test_data <- CV_mod$fold_1$data$phenotypes 
  
  CV_scores <- lapply(CV_results, `[[`, 1) #unlisting stuff; looks ugly
  CV_scores <- lapply(CV_scores, `[[`, 2)
  
  CV_ptves <- lapply(CV_results, `[[`, 2) #get ptves

  
  CV_null_scores <- lapply(CV_results, `[[`, 4)
  CV_null_scores <- lapply(CV_null_scores, `[[`, 2)
  
  
  accs <- unlist(lapply(CV_scores, `[[`, 1))
  print("Average CV accuracy:")
  print(mean(accs))
  
  ptvs = unlist(lapply(CV_ptves, `[[`, 1))
  print("Average train-PTVE:")
  print( mean(ptvs) )
  
  null_accs <- unlist(lapply(CV_null_scores, `[[`, 1))
  print("Average CV 0-accuracy:")
  print(mean(null_accs))
  
}