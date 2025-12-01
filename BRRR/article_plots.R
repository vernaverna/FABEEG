setwd("/projects/FABEEG/BRRR/")
library("vegan")
library("penalizedLDA")
library("ggplot2")
library("sjPlot")
library("ggpubr")
library("RColorBrewer")
library("corrplot")
library("viridis")

library("dplyr")
library("tidyr")
library("reshape2")
library("cvms")

library("rstatix")
library("lattice")
library("rlang")
library("lme4")
library("lmtest")
library("sandwich")

# functions ---------------------------------------------------------------
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

# # # # # # # # # # # # # # # #
#   load visualization data   #
# # # # # # # # # # # # # # # #

get_viz_data <- function(fname, OOS=FALSE){
  
  data_list <- vector(mode="list")
  
  load(fname)
  if(!OOS){
    X =  res$data$genotypes
    Y = res$data$phenotypes
  } else {
    X =  res$input_data[[3]]
    Y = res$input_data[[1]]
  }
  subj <- row.names(X)
  
  ages = read.csv('data/new_age_df.csv')
  ages <- ages[,-1]
  ages <- make_age_groups(ages)
  ages[ages==" "] <- NA #replace empty strings with NA-values
  ages = ages[ages$File %in% unique(subj)==TRUE, ] #to get rid of some extra subjects that should not be there
  ages['Sex'] <- as.factor(ages$Sex)
  
  rms_data = read.csv('data/emg_RMS_values.csv')
  names(rms_data)[names(rms_data) == 'X'] <- 'File'
  rms_data[,2:ncol(rms_data)] <- log10(rms_data[,2:ncol(rms_data)])
  ages = merge(x = ages, y = rms_data, by = "File", all.x = TRUE)
  
  age_df <- na.omit(ages$Age)
  
  inv_G <- res$scaling #inv(average(Gamma))
  lat_space <- Y%*%inv_G #obs x latent components
  lat_map = as.data.frame(rbind(lat_space))  
  
  data_list$X <- X
  data_list$Y <- Y
  data_list$subj <- subj
  data_list$ages <- ages
  data_list$age_df <- age_df
  data_list$lat_space <- lat_space
  data_list$lat_map <- lat_map
  data_list$W <- inv_G
  
  if(OOS){
    Y2 = res$input_data[[5]] #validation set data
    lat_space2 <- Y2%*%inv_G #obs x latent components
    lat_map2 = as.data.frame(rbind(lat_space2))  
    data_list$Y2 <- Y2
    data_list$lat_space2 <- lat_space2
    data_list$lat_map2 <- lat_map2
  }

  return(data_list)
}

#the validation function 
#' Function for validating the model
#' 
#' @param OOS Logical. If TRUE, does out-of-sample validation (unseen data points for all individuals)
#' @param dis which distance measure to use. One of 'L1' or 'L2', 'corr' or 'cos'
#' @param pK vector, which components to use in distance calculation
#' @param X the subject identifiers
#' @param lat_map the latent mapping used
#' @param lat_map2 is using unseen data
# 

validation <- function(OOS=T, dis='L1', pK=c(1:K), X=X, lat_map, lat_map2=NULL){
  
  #distance matrix# 
  D <- matrix(NA, ncol(X), ncol(X), dimnames=list(colnames(X), c(paste0('other', colnames(X)) )) )  

  # Helper function to calculate distance betewen 2 vectors
  dist_func <- function(x, y){
    if(dis=='L1'){
      return(sum(abs(x-y)))
    } else if(dis=='L2'){
      return(sqrt( sum(x-y)**2) )
    } else if(dis=='cos'){
      return( sum(x*y) / (sqrt(sum(x**2))*sqrt(sum(y**2))) )
    } else if(dis=='corr'){
      return(cor(x,y,method="pearson"))
    } else {
      browser(text="wrong distance measure!")
    }
  }
  
  if(OOS){
    proj <- lat_map2
  } else {
    proj <- lat_map
  }
  
  # compute distances in full data matrix (correlation + BRRR, minimum linkage) 
  for(testidx in 1:ncol(X)){ 
    for(m in 1:ncol(X)){     
      group_members <- rownames(proj)[m]   
      idxs = which(row.names(proj) %in% group_members) #should be 2!
      other_data <- proj[idxs[2], pK]

      dist <- dist_func(proj[testidx,pK], other_data) #L1-dist
      D[testidx,m] <-  dist 
    }
  } 
  
  PROJ <- D*0 #initialize 'projection' matrix

  for(r in 1:nrow(D)){ #assign group based on minimal distance/max correlation in lat. space
    if(dis=='L1'){
      index <- which.min(D[r,])
    } else {
      index <-which.max(D[r,])
    }
    PROJ[r,index] <- 1+PROJ[r,index]
  }
  colnames(PROJ) <- c(paste0("class",rownames(D)))
  accuracy = sum(diag(PROJ))/nrow(D)
  print(paste("Model accuracy:", accuracy)) 
  
  print("------------------------------------")
  return(list(D, accuracy))
  
}


get_distances <- function(viz_data, OOS = FALSE, K=30){
  X <- viz_data$X
  Y <- viz_data$Y
  subj <- viz_data$subj
  ages <- viz_data$ages
  n <- length(subj)/2
  W = viz_data$W
  
  lat_space <- viz_data$lat_space
  lat_map <- viz_data$lat_map

  if(!OOS){
    lat_map2=NULL
    Y2=Y
  } else {
    Y2 <- viz_data$Y2
    lat_map2 <- viz_data$lat_space2
  }
  
  D_list <- validation(OOS=OOS, dis='L1', pK=c(1:K), 
                       lat_map=lat_space, lat_map2=lat_map2,X=X)
  D <- D_list[[1]]
  # compute 0-model (correlations in full data matrix) 
  D0_list <- validation(OOS=OOS, dis='corr',pK=c(1:247), 
                        lat_map=Y, lat_map2=Y2, X=X)
  D0 <- D0_list[[1]]
  
  return(list(D, D0))
}


compute_differentiability <- function(distmat_vec, type='z_score'){
  #compute and add the mean distance to other subjects
  D = distmat_vec[[1]]
  C = distmat_vec[[2]] #correlation matrix
  
  if(type=='z_score'){
    n <- nrow(D)
    z_scores <- numeric(n)
    z_scores_c <- numeric(n)
    dist_to_others <- numeric(n)
    corr_to_others <- numeric(n)
    
    for (i in 1:n) {
      self_d <- D[i, i]                     # self distance
      self_c <- C[i, i]                     # self correlation
      
      other_d <- D[i, -i]                   # distances to others
      other_c <- C[i, -i]                   # correlations to others
      
      mu <- mean(other_d, na.rm = TRUE)     # mean of others
      sigma <- sd(other_d, na.rm = TRUE)    # sd of others
      mu_c <- mean(other_c, na.rm = TRUE)     # mean of others
      sigma_c <- sd(other_c, na.rm = TRUE)    # sd of others
            
      z_scores[i] <- -((self_d - mu) / sigma)  # z-scored self-distance (mirrored) 
      z_scores_c[i] <- (self_c - mu_c) / sigma_c  # z-scored self-correltaion 
      
      dist_to_others[i] <- mu 
      corr_to_others[i] <- mu_c
      
    }
    return(list(dist_to_others, z_scores, corr_to_others, z_scores_c))
  }
 
  else{
    dist_to_others <- c()
    corr_to_others <- c()
    for(i in 1:nrow(D)){
      other_d <- D[i, -i]                   # distances to others
      other_c <- C[i, -i]                   # correlations to others
      mu <- mean(other_d, na.rm = TRUE)     # mean of others
      mu_c <- mean(other_c, na.rm = TRUE)     # mean of others
      dist_to_others[i] <- mu 
      corr_to_others[i] <- mu_c
    }
    diff <- diag(D) /  dist_to_others
    diff_c <- diag(C) / corr_to_others
  }
  
  return(list(dist_to_others, diff, corr_to_others, diff_c))
}


# plots 1 ---------------------------------------------------------------
use_RMS <- TRUE
viz_data1 <- get_viz_data(fname = "results/full/all_2N1_BRRR_30.RData", OOS=F)
viz_data2 <- get_viz_data(fname= "results/full/all_2N2_BRRR_30.RData",  OOS=F)
viz_data3 <- get_viz_data(fname= "results/full/all_N1N2_BRRR_30.RData",  OOS=F)

ages <- viz_data1$ages
N1_dists <- get_distances(viz_data1, OOS=F)
N2_dists <- get_distances(viz_data2, OOS=F)
mixed_dists1 <- get_distances(viz_data3, OOS=F)

# Creating a dataframe for plotting....
reorder_val_idx <- match(rownames(N1_dists[[1]]), ages$File)
ages <- ages[reorder_val_idx,] 

D_n1 <- N1_dists[[1]]
D0_n1 <- N1_dists[[2]]
D_n2 <- N2_dists[[1]]
D0_n2 <- N2_dists[[2]]
D_n1n2 <- mixed_dists1[[1]][reorder_val_idx,reorder_val_idx]
D0_n1n2 <- mixed_dists1[[2]][reorder_val_idx,reorder_val_idx]


dist_df = as.data.frame(diag(D_n1), row.names = rownames(D_n1))
names(dist_df)[1] <- 'D_n1'
dist_df['D_n1'] = diag(D_n1)
dist_df['D_n2'] = diag(D_n2)
dist_df['D_n1n2'] = diag(D_n1n2)
dist_df['C_n1'] = diag(D0_n1)
dist_df['C_n2'] = diag(D0_n2) # C for correlation
dist_df['C_n1n2'] = diag(D0_n1n2)

rm()
dist_df['age'] = ages$Age #get age data
dist_df['sex'] = ages$Sex
dist_df['cap'] = as.factor(ages$Cap)
dist_df['age group'] = ages$Age_group
if(use_RMS){
  dist_df['RMS N1a'] = ages$RMS.N1a
  dist_df['RMS N1b'] = ages$RMS.N1b
  dist_df['RMS N2a'] = ages$RMS.N2a
  dist_df['RMS N2b'] = ages$RMS.N2b
  dist_df['RMS N2c'] = ages$RMS.N2c
  dist_df['RMS N2d'] = ages$RMS.N2d
}

dist_df['I(BRRR) N1'] <- compute_differentiability(N1_dists)[[2]]
dist_df['I(BRRR) N2'] <- compute_differentiability(N2_dists)[[2]]
dist_df['I(BRRR) N1N2'] <- compute_differentiability(mixed_dists1)[[2]][reorder_val_idx]
dist_df['I(corr) N1'] <- compute_differentiability(N1_dists)[[4]]
dist_df['I(corr) N2'] <- compute_differentiability(N2_dists)[[4]]
dist_df['I(corr) N1N2'] <- compute_differentiability(mixed_dists1)[[4]][reorder_val_idx]

# ========================================================

# FIGURE 6
# make a stacked histogram with subject info
age_df <- ages %>% mutate(age_group = cut(Age, breaks=14))
age_df <- na.omit(age_df)
prop_df <- age_df %>% group_by(age_group, Sex) %>% summarise(n = n()) %>% mutate(prop=n/sum(n) )
prop_df['prop'] <- round(prop_df$prop, 2)
prop_df$prop[prop_df$Sex=='M'] <- NA
prop_df <- prop_df %>% mutate(total=sum(n) )
age_df <- merge(x=age_df,y=prop_df, by=c('age_group', 'Sex'), all=TRUE)

pg <- ggplot(age_df, aes(x=age_group, fill=Sex)) + 
  geom_bar(position='stack') + theme_bw() + ylab("Frequency") + xlab('Age') +
  geom_text(aes(y=total, label=prop), vjust=-0.4, color='gray26', size=5) + 
  scale_fill_manual(values=c("#DCE319FF", "#33638DFF")) + # scale_fill_hue(l=70, c=90) + 
  scale_x_discrete(labels=c("0","","","4","","","8","","","12","","","16","")) +
  theme(legend.key = element_rect(linewidth = 16), 
        legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold", size = 18),
        axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file='demographics.pdf', plot=pg, width=10, height=7) 


# =================================================================================

# FIGURE 5
# tSNE projection of the subspace colored by age and sex.

library("Rtsne")
D <- normalize_input(lat_space)
set.seed(191) 
tsne <- Rtsne(D, perplexity=90, theta=0.0, check_duplicates = FALSE, max_iter=2000) 
reorder_val_idx <- match(rownames(D), ages$File)
ages <- ages[reorder_val_idx,] 
nsubj <- length(unique(subj))
reps = 1
# adding together N2 mappings  plus some covariates
lat_map = as.data.frame(rbind(lat_space))
lat_map['X1'] = tsne$Y[,1]
lat_map['X2'] = tsne$Y[,2]
lat_map['spectra'] = c(rep('N2A', nsubj), rep('N2B', nsubj))
lat_map['log_age'] = rep(log10(ages$Age), reps) #get age data
lat_map['age'] = rep(ages$Age, reps) #get age data
lat_map['group'] = rep(round(ages$Age, 0), reps) #get age data
lat_map['sex'] = rep(ages$Sex, reps)
lat_map['cap'] = rep(ages$Cap, reps)
lat_map['subject'] = subj[1:nsubj]


lat_map <- na.omit(lat_map)
p <- ggplot(data=lat_map, aes(x=X1, y=X2, shape=spectra, colour=age)) + 
            geom_point(alpha=0.5, size=3) +
            ggtitle("Subjects in latent mapping ") + 
            scale_color_viridis() +
            xlab("t-SNE #1") + ylab("t-SNE #2")  + #xlim(-42,48) + ylim(-20,25) +
            theme(legend.position="none") + theme_bw() + 
            theme(panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black"))
p
ggsave("figures/NEW_latmap_all_2N2_Age_K30.pdf", width=5.4, height=4.2)


q <- ggplot(data=lat_map, aes(x=X1, y=X2, shape=spectra, colour=cap)) + 
  geom_point(alpha=0.65, size=3) +
  ggtitle("Subjects in latent mapping ") + 
  scale_color_manual(values=c("#DCE319FF", "#33638DFF")) +
  xlab("t-SNE #1") + ylab("t-SNE #2") + #xlim(-42,48) + ylim(-20,25) +
  theme(legend.position="none") + theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
q
ggsave("figures/NEW_latmap_all_2N2_Sex_K30.pdf", width=5.4, height=4.2)


# =================================================================================
# FIGURE 4C ---- Differentiabiliy scores in lat.space // FULL DATA ...
#               ... according to age, factored by sex.

# remove NA -rows before muting the vars
dist_df <- na.omit(dist_df)
if(use_RMS){
  dist_df2 <- dist_df[,c('age', 'cap', 'sex', 'RMS N1a', 'RMS N1b', 'RMS N2a', 
                         'RMS N2b', 'RMS N2c', 'RMS N2d',
                         'I(BRRR) N1', 'I(BRRR) N2', 'I(BRRR) N1N2',
                         'I(corr) N1', 'I(corr) N2', 'I(corr) N1N2')]
  
} else{
  dist_df2 <- dist_df[,c('age', 'cap', 'sex', 
                         'I(BRRR) N1', 'I(BRRR) N2', 'I(BRRR) N1N2',
                         'I(corr) N1', 'I(corr) N2', 'I(corr) N1N2')]
}

# Reshape the data
long_dist_df <- dist_df2 %>%
  pivot_longer(
    cols = starts_with("I("),    # select all the "I(...)" columns
    names_to = c("Method", "Data"),
    names_pattern = "I\\(([^)]+)\\)\\s*(.*)",
    values_to = "value",
    values_drop_na = TRUE
  )

long_dist_df2 <- long_dist_df %>% rename('Zscore'='value')

#make a paired density / histogram plot with mean Z score annotations
means_df <- long_dist_df2 %>%
  group_by(Data, Method) %>%
  summarise(mean_Z = mean(Zscore, na.rm = TRUE), .groups = "drop")
# Compute the density height for each group
densities <- long_dist_df2 %>%
  group_by(Data, Method) %>%
  summarise(
    ymax = max(density(Zscore)$y),
    .groups = "drop"
  )
# Merge with means
means_df <- left_join(means_df, densities, by = c("Data", "Method")) %>%
  mutate(
    ymin = ifelse(Method == "corr", 0, -ymax.x),  # start of line
    ymax_plot = ifelse(Method == "corr", ymax.x, 0) # end of line
  )
means_p <- long_dist_df2 %>%
  group_by(Data) %>%
  filter(n_distinct(Method) == 2) %>%   # keep only groups with both methods
  t_test(Zscore ~ Method) %>%
  add_significance()


# Final chart
p <- ggplot(long_dist_df2, aes(x=Zscore, fill=Method, color=Method) ) +
  geom_density(data = subset(long_dist_df2, Method == "corr"),
                aes(y = after_stat(density)), alpha = 0.7) +
  geom_density(data = subset(long_dist_df2, Method == "BRRR"),
                aes(y = -after_stat(density)), alpha = 0.7) +
  geom_segment(data = means_df, aes(x = mean_Z, xend = mean_Z,
                                    y = ymin, yend = ymax_plot, linetype=Method),
              linewidth = 0.7, color='orangered') +
  geom_label(data = means_df, aes(x = mean_Z-3, y = ifelse(Method == "corr", 0.3, -0.3),
        label = sprintf("Mean: %.2f", mean_Z), fill = Method), size = 4, alpha=0.3, color = 'black',
        vjust = ifelse(means_df$Method == "corr", 0, 1), show.legend = FALSE) +
  geom_text(data = means_p, aes(x = -2, y = 1.0, label=sprintf("t = %.2f, %s", statistic, p.signif)), 
            inherit.aes = FALSE, size = 6) +
  facet_wrap(~Data) + 
  scale_fill_viridis_d(option='viridis', begin = 0.25, end = 0.85) +
  scale_color_viridis_d(begin = 0.25, end = 0.85) +
  theme_bw() +
  ylab("Frequency") + xlab("Z-scores") + 
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        panel.grid.major = element_blank(),
        strip.text = element_text(face="bold", size = 15),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

p
ggsave(file="figures/Z_score_densities_BRRR_vs_CORR_OOS.pdf", plot=p, width=10, height=10)


if(use_RMS){
  long_dist_df2 <- long_dist_df2 %>%
    pivot_longer(
      cols=starts_with("RMS"),
      names_to = "segment",             # new column for N1a, N1b, etc.
      names_pattern = "RMS\\s*(.*)",    # regex: remove "RMS " prefix
      values_to = "RMS"                 # new column for numeric values
    )
}

long_dist_df <- long_dist_df2 %>% filter(Method == "corr")

anova <- welch_anova_test(`Zscore` ~ Data, data = long_dist_df) # analysis of variance
summary(anova)
pwc <- games_howell_test(`Zscore` ~ Data, data = long_dist_df) # + post_hoc
print(pwc)
pwc <- pwc %>% add_xy_position(x='Data')

q <- ggplot(long_dist_df, aes(x=Data, y=`Zscore`))+
      geom_violin(aes(fill=Data), trim = TRUE, alpha=0.50) + 
      geom_boxplot(width=0.35, fill=NA) +
      scale_fill_viridis(option='cividis', discrete = T) +
      #scale_y_continuous(limits = c(0, 3)) +
      stat_pvalue_manual(pwc, hide.ns = TRUE) + 
      labs(subtitle = get_test_label(anova, detailed = TRUE),
            caption = get_pwc_label(pwc)) + ylim(-2.5,7.5)+
      theme_minimal() + ylab("Dist(self)") + xlab("") +
      theme(legend.text = element_text(size = 13),
            legend.title = element_text(size = 18),
            axis.text.x = element_blank(), axis.text.y = element_text(size=13),
            axis.title.x = element_text(size=18),
            axis.title.y = element_text(size=18),
            panel.border = element_blank(),
            panel.spacing = unit(85, 'points'),
            panel.grid.major = element_blank(),
            strip.text = element_text(face="bold", size = 15),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"))
q

ggsave(file="figures/Z_distance_latspace_fingerprint_all_corr_anova.pdf", plot=q, width=9, height=9)

################################################################################
#                        Regression plots (FIGURE 6)                           #
################################################################################

fit_RMS_multi_regression <- function(long_dist_df, method='BRRR', data='N1', cats=FALSE, OOS=FALSE){
  
  if(OOS){
    if(data=='N1ab'){
      s1="N1a"
      s2="N1b"
    } else if(data=='N2ad'){
      s1="N2a"
      s2="N2d"
    } else if(data=='N2ab_mix'){
      s1="N2a" 
      s2="N2b"
    } else if(data=='N2bc'){
      s1="N2b" 
      s2="N2c"
    }
  } else {
    s1=paste0(data,'a')
    s2=paste0(data,'b')
    if(data=='N1N2'){
      s1="N1a" 
      s2="N2b"
    }
  }
  
  data_df <- long_dist_df %>% filter(Method == method, Data==data, segment==s1 | segment==s2)
  data_df_w <- tidyr::pivot_wider(data_df, names_from = 'segment', values_from = 'RMS')
  data_df_w$noise_diff <- abs( data_df_w[,7][[1]] - data_df_w[,8][[1]] )

  theme_set(theme_sjplot(base_size = 18) )
  
  if(!cats){
    res_fname <- sprintf("results/regressions/%s_noisediff_age_results_all_%s.txt", method, data)
    model <- lm(noise_diff~age, data = data_df_w)
    p <- plot_model(model, type = "pred", terms = c("age"), vcov.fun = "HC3", show.data=TRUE, jitter=0.05) + 
      scale_color_manual(values=c("#339999"))
    ggsave(file=sprintf("figures/regressions/%s_noisediff_age_resid_all_%s.pdf", method, data), 
           plot=p, width=9, height=7)
    #m1 <- lm(Zscore ~ I(N1a-N1b), data=data_df_w)
  } else {
    res_fname <- sprintf("results/regressions/%s_Zscore_multi_cat_results_all_%s.txt", method, data)
    model <- lm(Zscore ~ noise_diff*sex + age + cap, data = data_df_w)
    p <- plot_model(model, vcov.fun = "HC3", show.values = TRUE, value.offset = .3) + 
       scale_color_manual(values=c("#339999", "#660066"))
    p2 <- plot_model(model, type = "pred", terms = c("noise_diff", "sex"), vcov.fun = "HC3", show.data=TRUE, jitter=0.05) + 
      scale_color_manual(values=c("#DCE319FF", "#33638DFF"))+ scale_fill_manual(values=c("#DCE319FF", "#33638DFF"))
    p3 <- plot_model(model, type = "pred", terms = c("age", "sex"), vcov.fun = "HC3", show.data=TRUE, jitter=0.05) + 
      scale_color_manual(values=c("#339999", "#660066"))+ scale_fill_manual(values=c("#339999", "#660066"))
    
    ggsave(file=sprintf("figures/regressions/%s_Zscore_multi_cat_results_all_%s.pdf", method, data), 
           plot=p, width=9, height=9)
    ggsave(file=sprintf("figures/regressions/%s_Zscore_multi_cat_RMS_resid_all_%s.pdf", method, data), 
           plot=p2, width=9, height=7)
    ggsave(file=sprintf("figures/regressions/%s_Zscore_multi_cat_AGE_resid_all_%s.pdf", method, data), 
           plot=p3, width=9, height=7)
  }
  
  sink(res_fname)    # start redirecting console output to file
  cat("=== Model Summary (excluding coefficients) ===\n\n")
  summ <- summary(model)
  cat("Call:\n")
  print(summ$call)
  cat("\nResidual standard error:", summ$sigma, "\n")
  cat("Multiple R-squared:", summ$r.squared, "\n")
  cat("Adjusted R-squared:", summ$adj.r.squared, "\n")
  cat("F-statistic:", summ$fstatistic[1], "on", summ$fstatistic[2], "and", summ$fstatistic[3], "DF,  p-value:",
      pf(summ$fstatistic[1], summ$fstatistic[2], summ$fstatistic[3], lower.tail=FALSE), "\n\n")
  cat("=== Robust Coefficients (HC3) ===\n\n")
  print(coeftest(model, vcov = vcovHC(model, type = "HC3")))
  sink()
  
}


fit_and_plot_regression <- function(long_dist_df, res_var='Zscore', method='BRRR', data='N1', seg='a', exp_var='age', gr_var='sex'){
  
  if(!is.na(seg)){
    if(data!='N1N2'){
      rms_seg = paste0(data,seg)
    } else {
      rms_seg='N1a' #N1a, N2b
    }
    data_df <- long_dist_df %>% filter(Method == method, Data==data, segment==rms_seg)
  } else {
    data_df <- long_dist_df %>% filter(Method == method, Data==data)
  }
  if(res_var=='Zscore' && !is.na(gr_var)){
    stat_t <-  data_df %>% t_test(as.formula(paste(res_var, " ~ ", gr_var))) %>% add_significance()
    stat_t <- stat_t %>% add_xy_position(x = gr_var)
    
    p <- ggplot(data_df, aes(x=!!sym(gr_var), y=!!sym(res_var)))+ geom_violin(aes(fill=!!sym(gr_var)), alpha=0.65) +
      scale_fill_manual(values=c("#DCE319FF", "#33638DFF")) + 
      #stat_summary(fun.data = "mean_sdl", geom="pointrange", size=2, color="black") +
      geom_boxplot(width=0.25) +
      #scale_y_continuous(limits = c(0, 2.5)) + 
      stat_pvalue_manual(stat_t, tip.length = 0) +
      labs(subtitle = get_test_label(stat_t, detailed = TRUE)) +
      theme_bw() + ylab(res_var) + xlab("") +
      theme(legend.key = element_rect(linewidth = 16), 
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 18),
            axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),
            axis.title.x = element_text(size=18),
            axis.title.y = element_text(size=18),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"))
    p
    figname = sprintf("figures/%s_Z_Differentiability_all_%s_%s.pdf", method, data, gr_var)
    ggsave(file=figname, plot=p, width=7, height=7)
    
  }
  if(is.na(gr_var)){
    res_fname <- sprintf("results/regressions/%s_%s_%s_results_all_%s.txt", method, res_var, exp_var, data)
    # --- Fit model and compute robust SEs ---
    form <- reformulate(exp_var, response = res_var)
    model <- lm(form, data = data_df)
    
    # Compute robust variance-covariance matrix 
    robust_vcov <- vcovHC(model, type = "HC3")
    coefs_robust <- coeftest(model, vcov. = robust_vcov)
    sink(res_fname)
    print(coefs_robust)
    sink()    
    
    # --- Build prediction data with robust SEs ---
    pred_df <- data_df %>% arrange(!!sym(exp_var)) %>%
      mutate(fit = predict(model, newdata = data_df, se.fit = TRUE)$fit)
    
    # Compute robust SE of fitted values
    X <- model.matrix(model)
    fit_vals <- X %*% coef(model)
    fit_var <- diag(X %*% robust_vcov %*% t(X))
    fit_se <- sqrt(fit_var)
 
    pred_df <- data.frame(x = data_df[[exp_var]], fit = fit_vals, se_robust = fit_se)
    pred_df <- setNames(pred_df, c(exp_var, "fit", "se_robust"))
    pred_df <- pred_df %>% mutate(lower = fit - 1.96 * se_robust,
                                  upper = fit + 1.96 * se_robust)
    # --- Regression equation (using robust p-value) ---
    r2 <- summary(model)$r.squared
    a <- coefs_robust[1, 1]
    b <- coefs_robust[2, 1]
    pval_b <- coefs_robust[2, 4]
    
    eqs <- paste0("italic(y) == ", round(a,2), " + ", round(b,2), " %.% italic(x)",
      "*',' ~italic(r)^2~'='~", round(r2,3), "*','~italic(p)~'='~", signif(pval_b,3))
    # ---  text position ---
    x_pos <- mean(data_df[[exp_var]], na.rm = TRUE) * 1.5
    y_pos <- max(data_df[[res_var]], na.rm = TRUE) * 1.05
    
    # --- Plot with robust SE ribbon ---
    p2 <- ggplot() + geom_point(data = data_df, aes(x = !!sym(exp_var), y = !!sym(res_var)),
                                alpha = 0.15) +
      geom_ribbon(data = pred_df, aes(x = !!sym(exp_var), ymin = lower, ymax = upper),
                  fill = "#660066", alpha = 0.25, color = NA) +
      geom_line(data = pred_df, aes(x = !!sym(exp_var), y = fit),
                color = "#660066", size = 1.2) +
      annotate("text", x = x_pos, y = y_pos, label = eqs, parse = TRUE, size = 5) +
      theme_minimal() + ylab(res_var) + xlab(exp_var) +
      theme(
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      )
    
    p2
    
    figname = sprintf("figures/regressions/%s_%s_%s_regplot_all_%s.pdf", method, res_var, exp_var, data)
    ggsave(file=figname, plot=p2, width=11, height=7)
    
  } else {
    res_fname <- sprintf("results/regressions/%s_%s_%s_by%s_results_all_%s.txt", method, res_var, exp_var, gr_var, data)
    robust_df <- data_df %>% group_by(!!sym(gr_var)) %>%
      group_modify(~{
        model <- lm(as.formula(paste(res_var, "~", exp_var)), data = .x)
        robust_vcov <- vcovHC(model, type = "HC3")
        # Prediction grid for this group
        x_seq <- seq(min(.x[[exp_var]], na.rm = TRUE),
                     max(.x[[exp_var]], na.rm = TRUE),
                     length.out = 100)
        newdata <- data.frame(x_seq)
        names(newdata) <- exp_var
        # Predictions with robust SEs
        pred <- predict(model, newdata = newdata, se.fit = TRUE, vcov. = robust_vcov)
        out <- tibble(
          !!exp_var := x_seq,
          fit = pred$fit,
          se_robust = pred$se.fit,
          lower = pred$fit - 1.96 * pred$se.fit,
          upper = pred$fit + 1.96 * pred$se.fit,
          !!gr_var := unique(.x[[gr_var]])[1]
        )
        out
      }) %>% ungroup()
    
    sink(res_fname)
    eqs_df <- data_df %>% group_by(!!sym(gr_var)) %>%
      group_modify(~{
        group_name <- unique(.x[[gr_var]])
        cat("====================================\n")
        cat("Group:", group_name, "\n")
        cat("====================================\n\n")
        model <- lm(as.formula(paste(res_var, "~", exp_var)), data = .x)
        robust_vcov <- vcovHC(model, type = "HC3")
        coefs <- coeftest(model, vcov. = robust_vcov)
        # --- Print results to file ---
        cat("OLS Summary:\n")
        print(summary(model))
        cat("\nRobust Coefficients:\n")
        print(coefs)
        cat("\n\n")
        
        # slope and intercept
        a <- coefs[1,1]
        b <- coefs[2,1]
        se_b <- coefs[2,2]
        pval_b <- coefs[2,4]
        r2 <- summary(model)$r.squared
        
        # Build equation string for ggplot parse
        eq_str <- paste0("italic(y) == ", round(a,2), " + ", round(b,2), " %.% italic(x)",
          "*',' ~italic(r)^2~'='~", round(r2,3), "*','~italic(p)~'='~", signif(pval_b,3) )
        # Determine x and y position for label dynamically
        x_pos <- mean(.x[[exp_var]], na.rm = TRUE) * 1.5
        y_pos <- max(.x[[res_var]], na.rm = TRUE) * 1.05  # slightly above max
        
        tibble(!!gr_var := unique(.x[[gr_var]]),
               x = x_pos,
               y = y_pos,
               label = eq_str)
      }) %>% ungroup()
    sink()
    
    p2 <- ggplot() + geom_point(data = data_df, aes(x = !!sym(exp_var), y = !!sym(res_var)),
                                alpha = 0.15) +
      geom_line(data = robust_df, aes(x = !!sym(exp_var), y = fit, color = !!sym(gr_var)),
                size = 1.2) +
      geom_ribbon(data = robust_df, aes(x = !!sym(exp_var), ymin = lower, ymax = upper, fill = !!sym(gr_var)),
                  alpha = 0.25, color = NA, inherit.aes = FALSE) +
      geom_text(data = eqs_df, aes(x = x, y = y, label = label),
                parse = TRUE, inherit.aes = FALSE, size = 6) +
      facet_grid(reformulate(gr_var)) +
      scale_fill_manual(values = c("#DCE319FF", "#33638DFF")) +
      scale_color_manual(values = c("#DCE319FF", "#33638DFF")) +
      theme_minimal() + xlab(exp_var) + ylab(res_var) +
      theme(
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      )
    
    p2
    
    figname = sprintf("figures/regressions/%s_%s_%s_by%s_regplot_all_%s.pdf", method, res_var, exp_var, gr_var, data)
    ggsave(file=figname, plot=p2, width=11, height=5)
  }
}

fit_and_plot_regression(long_dist_df2, res_var='RMS', method='corr', data='N1', seg='a', exp_var='age', gr_var='sex')
if(use_RMS){
  fit_RMS_multi_regression(long_dist_df2, method='corr', data='N2', cats=T, OOS=F)
}

################################################################################
#                Plot distance heatmaps for BRRRR model (fig  3)               #
################################################################################

## convert to tibble, add row identifier, and shape "long"
dat2 <- melt(C[30:50, 30:50])

p3 <- ggplot(dat2, aes(Var1, Var2, fill=value)) +
        geom_tile() +
        labs(fill='Corr') +
        #scale_fill_viridis(direction = -1) + xlab('') + ylab('') +
        scale_fill_viridis() + xlab('') + ylab('') +
        theme_bw() +
        theme(legend.key = element_rect(linewidth = 16), 
              legend.text = element_text(size = 13),
              legend.title = element_text(size = 18),
              axis.text.x = element_text(size=0), axis.text.y = element_text(size=0),
              axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18),
              axis.ticks=element_blank(), 
              panel.border = element_blank())
p3
ggsave(file="figures/small4_distmat_2N1_corr.pdf", plot=p3, width=11, height=10)


################################################################################
#           Mantell test for latent map matrices                               #
################################################################################
library("vegan")
viz_data <- get_viz_data(fname= "results/full/o7_N2AN2DN1AN1B_BRRR_K30.RData", OOS=T)
viz_data2 <- get_viz_data(fname= "results/full/o7_N2AN2BN2CN2D_BRRR_K30.RData", OOS=T)

K=30
D_list <- get_distances(viz_data, OOS=F)
D_list2 <- get_distances(viz_data2, OOS=F)
Dn1 <- D_list[[1]]
Cn1 <- D_list[[2]]
Dn2 <- D_list2[[1]]
Cn2 <- D_list2[[2]]

reorder_val_idx <- match(rownames(Dn1), rownames(Dn2))

#stats <- mantel.test(D1, D2, graph=TRUE)
print("BRRR matrices")
mantel_test <- mantel(xdis=Dn1, ydis=Dn2[reorder_val_idx, reorder_val_idx])
print(mantel_test)
print("----------------------------------")
print("Corr matrices")
mantel_test <- mantel(xdis=Cn1, ydis=Cn2[reorder_val_idx, reorder_val_idx])
print(mantel_test)
print("----------------------------------")
print("BRRR vs Corr matrices")
mantel_test <- mantel(xdis=-Dn2, ydis=Cn2)
print(mantel_test)
