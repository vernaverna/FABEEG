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
  
  CV_scores <- lapply(CV_results, `[[`, 1) #un-listing stuff
  CV_ptves <- lapply(CV_results, `[[`, 2)
  
  CV_models <- lapply(CV_results, `[`, 3)
  
  CV_scores <- lapply(CV_results, `[[`, 1) #un-listing 
  CV_scores <- lapply(CV_scores, `[[`, 2)
  
  CV_null_scores <- lapply(CV_results, `[[`, 4)
  CV_null_scores <- lapply(CV_null_scores, `[[`, 2)
  
  
  accs <- unlist(lapply(CV_scores, `[[`, 1))
  print("Average CV accuracy:")
  print(mean(accs))
  
  null_accs <- unlist(lapply(CV_null_scores, `[[`, 1))
  print("Average CV 0-accuracy:")
  print(mean(null_accs))
  
  ptvs = unlist(lapply(CV_ptves, `[[`, 1))
  print("Average train-PTVE:")
  print( mean(ptvs) )
  
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
  print("-------------------------------------")
  
  return(list(t, p_val))
}



do_OOS_validation <- function(population, psd_seq){
  
  test_seq <- tail(psd_seq,2)
  train_seq <- setdiff(psd_seq, test_seq)
  psd_seq_list <- paste(psd_seq, collapse='')
  
  fname <- paste0('results/full/', population, '_', psd_seq_list, '_BRRR_K30.RData')
  load(fname)
  
  Y = res$input_data[[1]]
  Y2 = res$input_data[[5]] #validation set data
  X =  res$input_data[[3]]
  W = res$scaling
  
  oos_proj <- Y2%*%W
  is_proj <- Y%*%W
  
  #distance matrices 
  D <- matrix(NA, ncol(X), ncol(X), dimnames=list(colnames(X), c(paste0('other', colnames(X)) )) )  
  D0 <- matrix(NA, ncol(X), ncol(X), dimnames=list(colnames(X), c(paste0('other', colnames(X)) )) )  
  
  # compute distances in full data matrix (correlation + BRRR, minimum linkage) 
  # assuming paired observations
   for(testidx in 1:ncol(X)){ 
     for(m in 1:ncol(X)){     
       group_members <- rownames(oos_proj)[m]   
       idxs = which(row.names(oos_proj) %in% group_members) #should be 2!
       other_data <- oos_proj[idxs[2], 1:30]
       other_data0 <- Y2[idxs[2],]
       
       dist <- sum(abs(oos_proj[testidx,1:30] - other_data)) #L1-dist
       dist0 <- cor(Y2[testidx,1:247], other_data0, method='pearson') 
   
       D[testidx,m] <-  dist 
       D0[testidx, m] <- dist0
      }
   } 
   
   PROJ <- D*0 #initialize 'projection' matrix
   PROJ0 <- D*0
   
   for(r in 1:nrow(D)){ #assign group based on minimal distance in lat. space
     index <- which.min(D[r,])
     index_cor <- which.max(D0[r,])
     PROJ[r,index] <- 1+PROJ[r,index]
     PROJ0[r,index_cor] <- 1+PROJ0[r,index_cor]
   }
   colnames(PROJ) <- c(paste0("class",rownames(D)))
   colnames(PROJ0) <- c(paste0("class",rownames(D0)))
   
   print(paste("Train seq:", paste(train_seq, collapse = ' '), "Test seq:", paste(test_seq, collapse=' ')))
   accuracy = sum(diag(PROJ))/nrow(D)
   print(paste("Model accuracy:", accuracy)) 
   
   print(paste("Model PTVE:", res$model$ptve)) 
   
   null_accuracy = sum(diag(PROJ0))/nrow(D0)
   print(paste("Null model accuracy:", null_accuracy)) 

   print("------------------------------------")
}



#############################################################


population <- "o7"
psd_seq <- "N1AN1BN2C"

stat_results <- compute_10fold_t_test(population, psd_seq) 


stat_results <- data.frame(
  row_id = c(1:12),
  population = c(rep("all", 6), rep("o7", 6))
)

#############################################################

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

#create set of conditions to loop over as a list
conds3 = list(c("N1A","N1B","N2A","N2B"),         # 1, 2    | 3, 4  aX near
              c("N2A","N2B","N1A","N1B"),         # 3, 4    | 1, 2  aX near, permuted 
              c("N1A","N1B","N2A","N2D"),         # 1, 2    | 3, 6  aX far
              c("N2A","N2D","N1A","N1B"),         # 3, 6    | 1, 2  aX far, permuted
              c("N2C","N2D","N1A","N1B"),         # 5, 6    | 1, 2  aX far, permuted
              c("N2B","N2C","N2D","N1A"),         # 4, 5    | 6, 2  mixed farish 
              c("N1A","N2B","N2A","N1B"),         # 1, 4    | 3, 2  mixed near
              c("N1A","N2D","N2A","N2B"),         # 1, 6    | 3, 4  mixed in between
              #c("N2A","N2B","N2A","N2D"),         # 3, 4    | 1, 6  mixed in between, permuted
              c("N1A","N2A","N2C","N2D"),         # 1, 3    | 5, 6  mixed far
              c("N2C","N2D","N1A","N2A"),         # 5, 6    | 1, 3  mixed far, permuted             
              c("N2A","N2B","N2C","N2D"),         # 3, 4    | 5, 6  N2 near
              c("N2A","N2D","N2B","N2C"))         # 3, 6    | 4, 5  N2 between

for(psd_seq in conds3){
  do_OOS_validation(psd_seq = psd_seq, population = 'o7')
}





#############################################################
library("reshape2")
library("viridis")

modality <- 'oos' #'subj' 'data', 'oos'

# Plot run results as bar charts

fname <- sprintf('/dataToR/unseen_%s_res_table.csv', modality)
if(modality=='oos'){
  fname <-  sprintf('/dataToR/%s_res_table.csv', modality)
}
res0 <- read.csv(paste0(getwd(),fname))
#
names(res0)[4] <- 'BRRR'
names(res0)[6] <- 'Corr'

# choose which sleep stage model to use
res0 <- res0[res0$Type=='Mixed',]

ptve_vec <- res0$PTVE
res0$PTVE <- NULL

res <- melt(res0, id.vars=c("Group", "Input", "Type", "Test"), #Test
            variable.name = 'Model', value.name ='SR')
#res$Test <- as.factor(res$Test)
res$Input <- as.factor(res$Input)
res$PTVE <- c(ptve_vec, rep(NA, 10))
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
  
  
} else if(modality=='data' | modality=='oos'){
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
  if(modality=='data'){
    ggsave(file="figures/within_unseen_data_model_comparison_ptve.pdf", plot=p2, width=8, height=8)
    
  } else {
    ggsave(file="figures/mixed_oos_model_comparison_ptve.pdf", plot=p2, width=8, height=6) 
  }
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




