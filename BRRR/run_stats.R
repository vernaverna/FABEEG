###########################################
# Script for run statistics (10-fold CV)  #
###########################################
setwd("/projects/FABEEG/BRRR/")
library("ggplot2")
library("ggpattern")
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

pos <- c(rep("all", 10)) #or o7
#psd_seqs <- c("N1AN1BN2C", "N1AN1BN2AN2C", "N1AN2BN2C", "N1AN2BN2AN1B",
#              "N2AN2BN2C", "N2AN2BN2CN2D")
psd_seqs <- c("N1AN1BN2C", "N1AN1BN2AN2C", "N1AN2BN2C", "N1AN2BN1B", "N1AN2BN2AN1B",
              "N2AN2BN2C", "N2AN2BN1B", "N2AN2BN2CN2D", "N2AN2BN2CN1B", "N1AN2BN2AN2D")
for(i in 1:10){
  population <- pos[i]
  psd_seq <- psd_seqs[i]
  print(psd_seq)
  res <- compute_10fold_t_test(population, psd_seq)
}





#############################################################
library("reshape2")
library("viridis")

modality <- 'data' #'data'

# Plot run results as bar charts
fname <- '/dataToR/unseen_subj_res_table.csv'
res0 <- read.csv(paste0(getwd(),fname))
#
names(res0)[3] <- 'BRRR'
names(res0)[5] <- 'Corr'

# choose which sleep stage model to use
res0 <- res0[res0$Type=='Across',]

ptve_vec <- res0$PTVE
res0$PTVE <- NULL

res <- melt(res0, id.vars=c("Group", "Input", "Type"),
            variable.name = 'Model', value.name ='SR')
#res$Test <- as.factor(res$Test)
res$Input <- as.factor(res$Input)
res$PTVE <- c(ptve_vec, rep(NA, 6))
#barplot(height=res$SR, names=res$Input)

#res <- res[res$Test=='N1',]
#res <- res[res$Input != "3N2", ]
#res <- res[res$Input != "2N1+N2", ]
res2 <- res %>% filter(!is.na(PTVE))

if(modality=='subj'){
  p <- ggplot() +
    geom_col_pattern(data=res, aes(pattern=ifelse(Model == "Corr", "stripe", "none"), x=Input, y=SR, fill=Input),
                     colour                   = 'black', 
                     pattern_density          = 0.35, 
                     pattern_spacing = 0.035,
                     pattern_angle = 45,
                     position='dodge') +
    scale_pattern_manual(values = c("none" = NA, "stripe" = "stripe")) +  # Define patterns
    scale_pattern_fill_manual(values = c("black")) +
    geom_line(data=res2, aes(x=Input,y=PTVE), 
              group=interaction(res2$Group, res2$Model), color='orangered1', linewidth=1) +
    #facet_grid(Group~Model) +
    facet_wrap(~Group) +
    #scale_fill_viridis(option='cividis', discrete=T) +
    scale_fill_viridis(discrete=T) +
    guides(pattern=guide_legend(title='Model')) +
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
  ggsave(file="figures/across_unseen_subj_model_comparison_ptve.pdf", plot=p, width=8, height=5)
  
  
} else if(modality=='data'){
  p2 <- ggplot() +
    geom_col_pattern(data=res, aes(pattern=ifelse(Model == "Corr", "stripe", "none"), x=Input, y=SR, fill=Input),
                     colour                   = 'black', 
                     pattern_density          = 0.35, 
                     pattern_spacing = 0.035,
                     pattern_angle = 45,
                     position='dodge') +
    scale_pattern_manual(values = c("none" = NA, "stripe" = "stripe")) +  # Define patterns
    scale_pattern_fill_manual(values = c("black")) +
    geom_line(data=res2, aes(x=Input,y=PTVE), 
              group=interaction(res2$Model, res2$Group), color='orangered1', linewidth=1) +
    facet_grid(Test~Group) +
    #facet_wrap(~Group) +
    #scale_fill_viridis(option='cividis', discrete=T) +
    scale_fill_viridis(discrete=T) +
    guides(pattern=guide_legend(title='Model')) +
    theme_minimal() + ylab("Success rate") + ylim(0.0,1.0) + xlab("") +
    theme(legend.text = element_text(size = 13),
          legend.title = element_text(size = 18),
          axis.text.x = element_text(size = 13), axis.text.y = element_text(size=13),
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
  p2
  ggsave(file="figures/within_unseen_data_model_comparison_ptve.pdf", plot=p2, width=8, height=8)
}


# Try with some splitting
res2 <- res[res$Group=='A',]


p2 <- ggplot() +
  geom_bar(data=res2, aes(x=Input, y=SR, fill=Input, alpha=0.6), position='dodge', stat='identity') +
  geom_line(data=res2, aes(x=Input,y=PTVE), group=res2$Model, color='orangered1', linewidth=2) +
  facet_wrap(~Model) +
  #scale_fill_viridis(option='cividis', discrete=T) +
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
p2


ggsave(file="figures/GrA_unseen_N1data_model_comparison_ptve.pdf", plot=p2, width=8, height=4)




