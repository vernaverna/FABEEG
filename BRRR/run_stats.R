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
library("rlang")

make_age_groups <- function(age_df){
  
  n_bins <- 10
  # Compute quantile breaks once
  q_breaks <- quantile(age_df$Age, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  q_breaks_rounded <- round(q_breaks, 1)
  
  # Apply inside mutate
  age_df <- age_df %>%
    mutate(
      Age_group = cut(
        Age,
        breaks = q_breaks,
        labels = paste0(q_breaks_rounded[-length(q_breaks_rounded)], " – ", q_breaks_rounded[-1]),
        include.lowest = TRUE,
        right = TRUE
      )
    )
  
  return(age_df)
}

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


do_validation <- function(population, psd_seq, OOS=TRUE){
  
  test_seq <- tail(psd_seq,2)
  train_seq <- setdiff(psd_seq, test_seq)
  psd_seq_list <- paste(psd_seq, collapse='')
  
  if(population=='o7|all'){
    pop='all'
  } else {
    pop=population
  }
  
  fname <- paste0('results/full/', pop, '_', psd_seq_list, '_BRRR_K30.RData')
  load(fname)
  
  Y = res$input_data[[1]]
  Y2 = res$input_data[[5]] #validation set data
  X =  res$input_data[[3]]
  W = res$scaling
  
  subj <- row.names(X)
  ages = read.csv('data/new_age_df.csv')
  ages <- ages[,-1]
  ages <- make_age_groups(ages)
  ages[ages==" "] <- NA #replace empty strings with NA-values
  ages = ages[ages$File %in% unique(subj)==TRUE, ] #to get rid of some extra subjects that should not be there
  ages['Sex'] <- as.factor(ages$Sex)
  names(ages)[names(ages) == 'File'] <- 'subj'
  
  if(population=='o7|all'){
    ages <- subset(ages, Age > 7) 
    keep_subj <- ages$subj
    X  <- X[rownames(X) %in% keep_subj, colnames(X) %in% keep_subj, drop=FALSE]
    Y  <- Y[rownames(X) %in% keep_subj, , drop=FALSE]
    Y2 <- Y2[rownames(X) %in% keep_subj, , drop=FALSE]
  }
  
  age_df <- na.omit(ages)
  oos_proj <- Y2%*%W
  is_proj <- Y%*%W
  
  if(OOS){
    proj <- oos_proj
    ref_Y <- Y2
  } else {
    proj <- is_proj
    ref_Y <- Y
  }
  
  #distance matrices 
  D <- matrix(NA, ncol(X), ncol(X), dimnames=list(colnames(X), c(paste0('other', colnames(X)) )) )  
  D0 <- matrix(NA, ncol(X), ncol(X), dimnames=list(colnames(X), c(paste0('other', colnames(X)) )) )  
  
  # compute distances in full data matrix (correlation + BRRR, minimum linkage) 
   for(testidx in 1:ncol(X)){ 
     for(m in 1:ncol(X)){     
       group_members <- rownames(proj)[m]   
       idxs = which(row.names(proj) %in% group_members) #should be 2!
       other_data <- proj[idxs[2], 1:30]
       other_data0 <- ref_Y[idxs[2],]
       
       dist <- sum(abs(proj[testidx,1:30] - other_data)) #L1-dist
       dist0 <- cor(ref_Y[testidx,1:247], other_data0, method='pearson') 
   
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
   if(OOS){
     print(paste("Train seq:", paste(train_seq, collapse = ' '), "Test seq:", paste(test_seq, collapse=' ')))
     
   } else {
     print(paste("Sleep model:", paste(train_seq, collapse=' ')) )

   }
   accuracy = sum(diag(PROJ))/nrow(D)
   print(paste("Model accuracy:", accuracy)) 
   
   print(paste("Model PTVE:", res$model$ptve)) 
   
   null_accuracy = sum(diag(PROJ0))/nrow(D0)
   print(paste("Null model accuracy:", null_accuracy)) 

   print("------------------------------------")
   
   res_dataframe <- data.frame(BRRR=diag(PROJ))
   res_dataframe$CORR <- diag(PROJ0)
   res_dataframe$subj <- rownames(PROJ)
   
   result_df <- merge(res_dataframe, ages, by="subj")
   result_df$Cap <- as.factor(result_df$Cap)
   
   return(result_df)
}

plot_by_grouping <- function(res_df, title, grouping='Age_group', method='BRRR'){
  
  success_rates <- res_df %>%
    group_by(!!sym(grouping)) %>%
    summarise(success_rate = mean(!!sym(method), na.rm = TRUE),n = n())
  
  if(method=='BRRR'){
    opt <- "D"
  } else {
    opt <- "G"
  }
  
  p <- ggplot(data=success_rates, aes(x=!!sym(grouping), y=success_rate, fill=!!sym(grouping))) + geom_col() + 
    scale_fill_viridis(discrete=T, option=opt) + theme_minimal() + 
    labs(
      x = "",
      y = "Success Rate",
      fill = grouping,
      title = title
    ) +
    ylim(0.0,1.0) +
    theme(legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          axis.text.x = element_blank(), axis.text.y = element_text(size=20),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=20),
          plot.title = element_text(size=20),
          panel.border = element_blank(),
          panel.spacing = unit(55, 'points'),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(face="bold", size = 24),
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_line(colour = "grey80"),
          panel.grid.minor.y = element_line(colour = "grey80"),
          axis.line = element_line(colour = "black"))
  
  return(list(success_rates, p))
  
}


plot_comparison <- function(df, title){
  # Reshape the data
  plot_df <- df %>%
    pivot_longer(
      cols = -c(Age_group, n),
      names_to = c("Method"),
      names_pattern = "([^:]+)",  # everything before ':' = Method, after = Data
      values_to = "Success rate"
    )
  g <- ggplot() +
    geom_col_pattern(data=plot_df, 
                     aes(pattern=ifelse(Method == "CORR", "stripe", "none"), x=Age_group, y=`Success rate`, fill=Age_group),
                     colour                   = 'black', 
                     pattern_density          = 0.25, 
                     pattern_spacing = 0.035,
                     pattern_angle = 45,
                     position='dodge',
                     show.legend = c(pattern = FALSE)) +
    scale_pattern_manual(values = c("none" = NA, "stripe" = "stripe")) +  # Define patterns
    scale_pattern_fill_manual(values = c("CORR", "black")) +
    scale_fill_viridis(discrete=T) +
    guides(pattern=guide_legend(title='Method')) +
    theme_minimal() + ylab("Success rate") + ylim(0.0,1.0) + xlab("") +
    theme(legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          axis.text.x = element_blank(), axis.text.y = element_text(size=18),
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
  g
  ggsave(file= paste0("figures/full", grouping, "_", title, "comparison.pdf"), plot=g, width=8, height=6)
  
  
}
#############################################################

pos <- c(rep("all", 10)) #or o7
psd_seqs <- c("N1AN1BN2C", "N1AN1BN2AN2C", "N1AN2BN2C", "N1AN2BN1B", "N1AN2BN2AN1B",
              "N2AN2BN2C", "N2AN2BN1B", "N2AN2BN2CN2D", "N2AN2BN2CN1B", "N1AN2BN2AN2D")
for(i in 1:10){
  population <- pos[i]
  psd_seq <- psd_seqs[i]
  print(psd_seq)
  res <- compute_10fold_t_test(population, psd_seq, OOS=FALSE)
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
              c("N1A","N2A","N2C","N2D"),         # 1, 3    | 5, 6  mixed far
              c("N2C","N2D","N1A","N2A"),         # 5, 6    | 1, 3  mixed far, permuted             
              c("N2A","N2B","N2C","N2D"),         # 3, 4    | 5, 6  N2 near
              c("N2A","N2D","N2B","N2C"))         # 3, 6    | 4, 5  N2 between

grouping='Age_group'
plotting=FALSE
#grouping='Sex'

for(psd_seq in conds3){
  res_df <- do_validation(psd_seq = psd_seq, population = 'o7|all', OOS = TRUE)
  res_df2 <- do_validation(psd_seq = psd_seq, population = 'o7|all', OOS = FALSE)
  lat_space_stages <- paste(psd_seq[1:2], collapse='+') #relevant only for BRRR
  test_stages <- paste(psd_seq[3:4], collapse='+') #where the fingerprint is calculated
  title1 <- paste("BRRR:", lat_space_stages, '->', test_stages)
  title1b <- paste("BRRR:", lat_space_stages)
  title2 <- paste("CORR:", test_stages)
  title2b <- paste("CORR:", lat_space_stages)
  res_df <- na.omit(res_df)
  res_df2 <- na.omit(res_df2)
  
  if(plotting){
    #OOS plots (BRRR)
    res_brrr <- plot_by_grouping(res_df, title1, grouping=grouping, method='BRRR')
    p=res_brrr[[2]]
    
    #WS PLOTS (BRRR and CORR)
    res2_brrr <- plot_by_grouping(res_df2, title1b, grouping=grouping, method='BRRR')
    p1=res2_brrr[[2]]
    res_corr <- plot_by_grouping(res_df, title2, grouping=grouping, method='CORR')
    q=res_corr[[2]]
    res2_corr <- plot_by_grouping(res_df2, title2b, grouping=grouping, method='CORR')
    q1=res2_corr[[2]]
    total_oos <- merge(res_brrr[1], res_corr[1], by=c(grouping, "n"))
    names(total_oos)[names(total_oos) == 'success_rate.x'] <- title1
    names(total_oos)[names(total_oos) == 'success_rate.y'] <- title2
    
    total_ws <- merge(res2_brrr[1], res2_corr[1], by=c(grouping, "n"))
    names(total_ws)[names(total_ws) == 'success_rate.x'] <- title1b
    names(total_ws)[names(total_ws) == 'success_rate.y'] <- title2b
    plot_comparison(total_oos, title1)
    plot_comparison(total_ws, title1b)
    
    total <- merge(total_oos, total_ws, by = c(grouping, "n"))
    write.csv(total, file = paste0("results/full", grouping, "_", title1, ".csv"), row.names = FALSE)
    
    ggsave(file=paste0("figures/full_OOS_BRRR_",paste(grouping, lat_space_stages, '->', test_stages),'.pdf'), plot=p, width=8, height=6)
    ggsave(file=paste0("figures/full_CORR_",paste(grouping, test_stages),'.pdf'), plot=q, width=8, height=6)
    ggsave(file=paste0("figures/full_WS_BRRR_",paste(grouping, lat_space_stages),'.pdf'), plot=p1, width=8, height=6)
    ggsave(file=paste0("figures/full_CORR_",paste(grouping, lat_space_stages),'.pdf'), plot=q1, width=8, height=6)
  }
 
}

#############################################################
library("reshape2")
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


################################################################################
#   Convergence checks for coefficient matrices and ptve traces                # 
################################################################################
library("rstan")

# load both models
#load("results/full/i1000_all_2N2_BRRR_K30.RData")
load("results/full/o7_2N2N1_BRRR_K30.RData")
res1 <- res
load("results/full/i5000_all_2N2_BRRR_K30.RData")
res2 <- res
remove(res)

# Traces... there are 500
traces_2 <- res2$traces
traces_1 <- res1$traces

tr_psi1 <- traces_1$Psi
tr_gamma1 <- traces_1$Gamma
ptves <- unlist(traces_1$tpve)
iter <- traces_1$iter

# because of the rotation invariance, the convergence needs to be studied w.r.t. to
# the standard regression coefficient matrix Psi*Gamma

# for all iterations, compute theta, esitmate 
# for the set of coefficient matrices..
# ... sample 200 or so parameters (Gillberg et al)
# ... and store their values in a matrix
# which should be plotted. 

# Compute the standard regression coefficient matrices
coefMats <- vector(mode="list", length=length(iter))
for(i in 1:length(iter)){
  coef_mat <- tr_psi1[[i]]%*%tr_gamma1[[i]]
  coefMats[[i]] <- coef_mat
}

n_subj <- nrow(coef_mat)
set.seed(121)
rows <- sample(n_subj, size=200) # choose 200 random element indices from the matrix
cols <- sample(247, size=200)
# pick iteration series for each of the params
param_list <- vector(mode="list", length = 200)
for(j in 1:200){
  param_series <- lapply(coefMats, "[", rows[j], cols[j])
  param_list[[j]] <- unlist(param_series)
}

param_matrix <- do.call(rbind, param_list)
# collect some stats from indices of param matrix
R_hats <- c()
ESS_bs <- c()
ESS_ts <- c()
for(r in 1:200){
  samples <- param_matrix[r,]
  R_hats <- c(R_hats, Rhat(samples))
  ESS_bs <- c(ESS_bs, ess_bulk(samples))
  ESS_ts <- c(ESS_ts, ess_tail(samples))
}

convergence_stats <- as.data.frame(cbind(R_hats, ESS_bs, ESS_ts))
mu<-apply(convergence_stats, 2, mean)
sd<-apply(convergence_stats, 2, sd)
rbind(mu, sd) #These are ok!

Rhat = Rhat(unlist(ptves))
ESS_b = ess_bulk(unlist(ptves))
ESS_t = ess_tail(unlist(ptves))
pdf("figures/i1000_ptve_traces.pdf", width=9, height=6)
plot(iter, ptves, xlab='iteration', 'l', col='darkorchid4', ylab="PTVE score")#, ylim = c(0.84312, 0.84355))
text(870, 0.8435, paste0("Rhat=", round(Rhat, 4), ", ESS_b=", round(ESS_b, 4),
                         ", ESS_t=", round(ESS_t, 4)))
dev.off()
#legend(950, 0.8432, legend='PTVE', col="yellow3", pch=16)

