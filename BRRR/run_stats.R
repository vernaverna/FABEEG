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
  
  fname <- paste0('results/10foldCV/unseen_data/NEW_K30_', population, '_', psd_seq, '.RData')
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

population <- "o7"
psd_seq <- "N1AN1BN2C"

stat_results <- compute_10fold_t_test(population, psd_seq) 


stat_results <- data.frame(
  row_id = c(1:12),
  population = c(rep("all", 6), rep("o7", 6))
)

pos <- c(rep("o7", 9)) #or o7
#psd_seqs <- c("N1AN1BN2C", "N1AN1BN2AN2C", "N1AN2BN2C", "N1AN2BN2AN1B",
#              "N2AN2BN2C", "N2AN2BN2CN2D")
psd_seqs <- c("N1AN1BN2C", "N1AN1BN2AN2C", "N1AN2BN2C", "N1AN2BN1B", "N1AN2BN2AN1B",
              "N2AN2BN2C", "N2AN2BN1B", "N2AN2BN2CN2D", "N2AN2BN2CN1B")
for(i in 1:9){
  population <- pos[i]
  psd_seq <- psd_seqs[i]
  print(psd_seq)
  res <- compute_10fold_t_test(population, psd_seq)
}





#############################################################
library("reshape2")
library("viridis")

# Plot run results as bar charts

fname <- '/dataToR/unseen_subj_res_table.csv'
res0 <- read.csv(paste0(getwd(),fname))
#
names(res0)[3] <- 'BRRR'
names(res0)[5] <- 'Corr'
ptve_vec <- res0$PTVE
res0$PTVE <- NULL

res <- melt(res0, id.vars=c("Group", "Input"),
            variable.name = 'Model', value.name ='SR')
#res$Test <- as.factor(res$Test)
res$Input <- as.factor(res$Input)
res$PTVE <- c(ptve_vec, rep(NA, 12))
#barplot(height=res$SR, names=res$Input)

res <- res[res$Test=='N1',]

p <- ggplot() +
      geom_bar(data=res, aes(x=Input, y=SR, fill=Input), position='dodge', stat='identity') +
      geom_line(data=res, aes(x=Input,y=PTVE), group=interaction(res$Group, res$Model), color='red') +
      facet_grid(Group~Model) +
      #facet_wrap(~Model) +
      scale_fill_viridis(discrete=T) +
      theme_minimal() + ylab("Success rate") + ylim(0.0,1.0) + xlab("") +
      theme(legend.text = element_text(size = 13),
            legend.title = element_text(size = 18),
            axis.text.x = element_blank(), axis.text.y = element_text(size=13),
            axis.title.x = element_text(size=18),
            axis.title.y = element_text(size=18),
            panel.border = element_blank(),
            panel.spacing = unit(55, 'points'),
            panel.grid.major.x = element_blank(),
            strip.text = element_text(face="bold", size = 15),
            panel.grid.minor.x = element_blank(), 
            panel.grid.major.y = element_line(colour = "grey80"),
            panel.grid.minor.y = element_line(colour = "grey80"),
        axis.line = element_line(colour = "black"))
p
ggsave(file="figures/unseen_subj_model_comparison_ptve.pdf", plot=p, width=8, height=8)




