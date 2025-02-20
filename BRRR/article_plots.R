# collection of plots for teh article.

setwd("/projects/FABEEG/BRRR/")
library("vegan")
library("penalizedLDA")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("cvms")
library("RColorBrewer")
library("corrplot")
library("viridis")
library("rstatix")

# For viridis color scheme, see e.g. 
#   https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

# # # # # # # # # # # # # # # #
#   load visualization data   #
# # # # # # # # # # # # # # # #

get_viz_data <- function(fname){
  
  data_list <- vector(mode="list")
  
  load(fname)
  Y <- res$data$phenotypes
  X <- res$data$genotypes
  subj <- row.names(X)
  ages = read.csv('data/new_age_df.csv')
  ages <- ages[,-1]
  ages[ages==" "] <- NA #replace empty strings with NA-values
  ages = ages[ages$File %in% unique(subj)==TRUE, ] #to get rid of some extra subjects that should not be there
  ages['Sex'] <- as.factor(ages$Sex)
  
  age_df <- na.omit(ages)
  
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
  return(data_list)
}

#the validation function 
#' Function for validating the model
#' 
#' @param within_sample Logical. If FALSE, does out-of-sample validation (unseen data points for all individuals)
#' @param dis which distance measure to use. One of 'L1' or 'L2', 'corr' or 'cos'
#' @param pK vector, which components to use in distance calculation
#' @param Xt the subject identifiers
#' @param lat_map the latent mapping used
#' @param lat_map2 is using unseen data
#' @param linkage either 'average' or 'minimum', for deciding how to calculate distances
# 

validation <- function(within_sample=FALSE, dis='L1', pK=c(1:K), Xt=X, lat_map, lat_map2=NULL, 
                       linkage='minimum'){
  
  #distance matrix# 
  D <- matrix(NA, ncol(Xt), ncol(Xt), dimnames=list(colnames(Xt), c(paste0('other', colnames(Xt)) )) )  
  
  # permute lat_map to have same row_names as Xt - safety measure
  # TODO: check this??? do i mess it up. 
  if(!is.null(lat_map2)){
    reorder_val_idx <- match(colnames(Xt), rownames(lat_map2))
    lat_map2 <- lat_map2[reorder_val_idx,]    
  }
  
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
  
  #calculates the distances between individuals in lat.space   
  for(testidx in 1:ncol(Xt)){ 
    for(m in 1:ncol(Xt)){     
      group_members <- rownames(lat_map)[m]   
      idxs = which(row.names(lat_map) %in% group_members) #should be at least 2!
      
      if(within_sample){
        
        if(length(idxs) == 3){
          if(linkage=='average'){
            other_data <- colMeans(lat_map[tail(idxs,2), pK])
          } else if(linkage=='minimum'){
            other_data <- lat_map[tail(idxs,2), pK]
            dists <- apply(other_data, 1, function(x) dist_func(lat_map[testidx,pK], x))
          }
        } else {
          if(linkage=='average'){
            other_data <- lat_map[tail(idxs,1), pK]
          } else if(linkage=='minimum'){
            other_data <- lat_map[tail(idxs,1), pK]
            dists <- dist_func(lat_map[testidx,pK],other_data) 
          }
        }
        # Finally, determine the distances in the latent coordinates
        if(linkage=='average'){ 
          D[testidx,m] <- dist_func(lat_map[testidx,pK],other_data)
        } else if (linkage=='minimum'){
          D[testidx,m] <- if(dis!='corr') min(dists) else max(dists)
        }
        
        
      } else { #OOS validation (unseen data)
        
        if(length(idxs)>1){
          if(linkage=='average'){
            other_data <- colMeans(lat_map[idxs, pK])
          } else if(linkage=='minimum'){
            other_data <- lat_map[idxs, pK]
            dists <- apply(other_data, 1, function(x) dist_func(lat_map2[testidx,pK], x))
          }
          
        } else {
          if(linkage=='average'){
            other_data <- lat_map[idxs, pK]
          } else if(linkage=='minimum'){
            other_data <- lat_map[idxs, pK]
            dists <- dist_func(lat_map2[testidx,pK],other_data)
          }
        }
        if(linkage=='average'){
          D[testidx,m] <- dist_func(lat_map2[testidx,pK],other_data)
        } else if(linkage=='minimum'){
          D[testidx,m] <- if(dis!='corr') min(dists) else max(dists)
        }
      }
    } 
  }
  
  PROJ <- D*0 #initialize 'projection' matrix
  
  for(r in 1:nrow(D)){ #assign group based on minimal distance in lat. space
    if(dis!='corr'){
      index <- which.min(D[r,])
    } else {
      index <- which.max(D[r,]) #in case of correlation, we are of course looking for maximum correlation between subs
    }
    
    PROJ[r,index] <- 1+PROJ[r,index]
  }
  colnames(PROJ) <- c(paste0("class",rownames(D)))
  
  accuracy = sum(diag(PROJ))/nrow(D)
  print(paste("Model accuracy:", accuracy)) #alright this is how I like it ;)
  
  rankings <- apply(D, 2, rank)
  avg_rank = sum(diag(rankings))/nrow(rankings)
  
  return(list(D, accuracy, avg_rank))
  
}

########

double_subs <- read.csv("longitudinal_subset.csv")

viz_data1 <- get_viz_data(fname = "results/full/all_2N1_BRRR_30.RData")
viz_data2 <- get_viz_data(fname= "results/full/all_2N2_BRRR_30.RData")
viz_data3 <- get_viz_data(fname= "results/full/all_N1N2_BRRR_30.RData")

ages <- viz_data1$ages

get_distances <- function(viz_data, viz_data2, within_sample = T, K=30, half='first'){
  X <- viz_data$X
  Y <- viz_data$Y
  subj <- viz_data$subj
  ages <- viz_data$ages
  n <- length(subj)/2
  W = viz_data$W
  
  lat_space <- viz_data$lat_space
  lat_map <- viz_data$lat_map

  if(within_sample){
    lat_map2=NULL
    Y2=Y
  } else {
    if(half=='second'){
      Y2 <- viz_data2$Y[(n+1):(2*n),] #788:,? 
    } else {
      Y2 <- viz_data2$Y[1:n,] #788:,? 
    }
    lat_map2 <- Y2%*%W
  }
  
  D_list <- validation(within_sample = within_sample, dis='L1', pK=c(1:K), 
                       lat_map=lat_space, lat_map2=lat_map2,Xt=X)
  D <- D_list[[1]]
  # compute 0-model (correlations in full data matrix) 
  D0_list <- validation(within_sample = within_sample, dis='corr',pK=c(1:247), 
                        lat_map=Y, lat_map2=Y2, Xt=X)
  D0 <- D0_list[[1]]
  
  return(list(D, D0))
}


compute_differentiability <- function(distmat_vec){
  #compute and add the mean distance to other subjects
  D = distmat_vec[[1]]
  D0= distmat_vec[[2]]
  
  dist_to_others <- c()
  dist_to_others_corr <- c()
  for(i in 1:nrow(D)){
    to_include <- setdiff((1:nrow(D)), i) #take all distances except to self
    to_include2 <- setdiff((1:nrow(D0)), i)
    avg_dist <- mean(D[i, to_include])
    avg_dist_corr <- mean(D0[i, to_include2])
    dist_to_others <- c(dist_to_others, avg_dist)
    dist_to_others_corr <- c(dist_to_others_corr, avg_dist_corr)
  }
  diff <- diag(D) /  dist_to_others
  diff_c <- -diag(D0) / dist_to_others_corr
  
  return(list(dist_to_others, diff, dist_to_others_corr, diff_c))
}

N1_dists <- get_distances(viz_data1, viz_data2, within_sample = T)
N2_dists <- get_distances(viz_data2, viz_data2=viz_data1, within_sample = T)
mixed_dists1 <- get_distances(viz_data3, viz_data2, within_sample = T)
mixed_dists2 <- get_distances(viz_data3, viz_data2=viz_data1, within_sample = F, half='second')

# Creating a dataframe for plotting....
reorder_val_idx <- match(rownames(N1_dists[[1]]), ages$File)
ages <- ages[reorder_val_idx,] 

D_n1 <- N1_dists[[1]]
D0_n1 <- N1_dists[[2]]
D_n2 <- N2_dists[[1]]
D0_n2 <- N2_dists[[2]]
D_n1n2 <- mixed_dists1[[1]]
D0_n1n2 <- mixed_dists1[[2]]
D_n1n2_ <- mixed_dists2[[1]]
D0_n1n2_ <- mixed_dists2[[2]]

dist_df = as.data.frame(diag(D_n1), row.names = rownames(D_n1))
names(dist_df)[1] <- 'D_n1'
dist_df['D_n2'] = diag(D0_n2)
dist_df['D_n1n2'] = diag(D0_n1n2)
dist_df['D0_n1'] = diag(D0_n1)
dist_df['D0_n2'] = diag(D0_n2)
dist_df['D0_n1n2'] = diag(D0_n1n2)
dist_df['D_n1n2_'] = diag(D_n1n2_)
dist_df['D_n1n2_'] = diag(D0_n1n2_)
dist_df['age'] = ages$Age #get age data
#dist_df['group'] = rep(round(ages$Age, 0), reps) #get age data
dist_df['sex'] = ages$Sex
dist_df['cap'] = ages$Cap


dist_df['WS-dist N1'] <- compute_differentiability(N1_dists)[[2]]
dist_df['WS-dist N2'] <- compute_differentiability(N2_dists)[[2]]
dist_df['WS-dist N1N2'] <- compute_differentiability(mixed_dists1)[[2]]
#dist_df['WS-dist N1N2_N1'] <- compute_differentiability(mixed_dists1)[[2]]

dist_df['BS-dist N1'] <- compute_differentiability(N1_dists)[[1]]
dist_df['BS-dist N2'] <- compute_differentiability(N2_dists)[[1]]
dist_df['BS-dist N1N2'] <- compute_differentiability(mixed_dists1)[[1]]
#dist_df['BS-dist N1N2_N1'] <- compute_differentiability(mixed_dists2)[[1]]

### for general use?

X <- viz_data1$X
Y <- viz_data1$Y
subj <- viz_data1$subj
ages <- viz_data1$ages
age_df <- viz_data1$age_df
lat_space <- viz_data1$lat_space
lat_map <- viz_data1$lat_map
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
# NONNIH.

# =================================================================================

# FIGURE 5
# tSNE projection of the subspace colored by age and sex.

## UPDATE: NOPE - not happening.

library("Rtsne")
D <- normalize_input(lat_space)
set.seed(191) #20230505
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
  #scale_color_viridis(discrete=TRUE) +
  scale_color_manual(values=c("#DCE319FF", "#33638DFF")) +
  xlab("t-SNE #1") + ylab("t-SNE #2") + #xlim(-42,48) + ylim(-20,25) +
  theme(legend.position="none") + theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
q
ggsave("figures/NEW_latmap_all_2N2_Sex_K30.pdf", width=5.4, height=4.2)


for(comp in c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')){
  print(paste0('Correlation between age and ', comp, ':'))
  rho=cor(lat_map$age, lat_map[comp], method='pearson', use = "complete.obs")
  print(rho)
  
  print(paste0('Correlation between log-age and ', comp, ':'))
  rho2=cor(lat_map$log_age, lat_map[comp], method='pearson', use = "complete.obs")
  print(rho2)
}

# =================================================================================
# FIGURE N ---- Intra-subject distances in lat.space // FULL DATA ...
#               ... according to age, factored by sex.
#  possibilities: lm-plot, violinplot
library("reshape2")
library("tidyr")
library("dplyr")
#library("ggpattern")

# remove NA -rows before muting the vars
dist_df <- na.omit(dist_df)
#dist_df['log_age'] = log10(dist_df$age) #get age data
dist_df2 <- dist_df[,c('age', 'sex', 'WS-dist N1', 'WS-dist N2', 'WS-dist N1N2',# 'WS-dist N1N2_N1',
                       'BS-dist N1', 'BS-dist N2', 'BS-dist N1N2')]
#dist_df2 <- dist_df[,c('age', 'sex', 'diff N1', 'diff N2', 'diff N1N2_N1', 'diff N1N2_N2')]
#names(dist_df2)[3:6] <- c('N1', 'N2', 'N1N2_1', 'N1N2_2') 
#long_dist_df <- melt(dist_df2, id.vars=c("age", "sex"))

# Reshape the data
long_dist_df <- dist_df2 %>%
  pivot_longer(
    cols = starts_with("WS-dist") | starts_with("BS-dist"),   
    names_to = c("Type", "Data"),       # Create new columns `type` and `group`
    names_pattern = "([WB])_?(.*)",     # Regex to capture "W" or "B" and "aa", "bb", "cc"
    values_drop_na = TRUE       # Optional: drop any NA values if they exist
  )

long_dist_df$Data = sub("^S-dist ", "", long_dist_df$Data)
long_dist_df <- long_dist_df %>% rename('Distance'='value')

# add info on test segments, if applicable
test <- c()
for(i in 1:nrow(long_dist_df)){
  input <- long_dist_df$Segments[i]
  
  if(input %in% c('N2', 'N1N2_1')){
    test <- c('N1', test)
  } else {
    test <- c('N2', test)
  }
}
long_dist_df$Test <- as.factor(test)

# AHHHH I think I'm gonna still use only W or B
long_dist_df <- long_dist_df %>% filter(Type == "W")

# analysis of variance
anova <- welch_anova_test(Distance ~ Data, data = long_dist_df)
summary(anova)
# post_hoc
pwc <- games_howell_test(Distance ~ Data, data = long_dist_df)
print(pwc)
pwc <- pwc %>% add_xy_position(x='Data')

q <- ggplot(long_dist_df, aes(x=Data, y=Distance))+
      geom_violin(aes(fill=Data), trim = TRUE, alpha=0.50) + 
      geom_boxplot(width=0.35, fill=NA) +
      scale_fill_viridis(option='cividis', discrete = T) +
      scale_y_continuous(limits = c(0, 3)) +
      stat_pvalue_manual(pwc, hide.ns = TRUE) + 
      labs(subtitle = get_test_label(anova, detailed = TRUE),
            caption = get_pwc_label(pwc)) + #ylim(0.0,2.0)+
      theme_minimal() + ylab("Within-Subject Distance") + xlab("") +
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


ggsave(file="figures/WS_distance_latspace_fingerprint_all_BRRR_anova.pdf", plot=q, width=9, height=9)

###############################################################################


###############################################################################
# longitudinal subset
# check distance to later meas. sessions
longitudinal_dists <- c()

for(s in 1:nrow(double_subs)){
  sub1 = double_subs[s,1]
  sub2 = double_subs[s,2]
  
  if(sub1 %in% rownames(D) & sub2 %in% rownames(D)){
    longitudinal_dist <- D[sub1, paste0('other',sub2)]
    print(paste(sub1, sub2))
    print("distance to later instance")
    print(longitudinal_dist)
    longitudinal_dists <- c(longitudinal_dist, longitudinal_dists)
    print(paste("minimum distance:", sub1))
    print(min(D[sub1,]))
    print(paste("minimum distance:", sub2))
    print(min(D[sub2,]))
  } else {
    print("Both subjects not found in the data")
  }

}
  

################################################################################

dist_df['differentiability'] = dist_df['WS-dist N2']

# compare if significant difference between F vs M using t-test
stat_t <- dist_df %>% t_test(differentiability~sex) %>% add_significance()
stat_t <- stat_t %>% add_xy_position(x = "sex")

#y=`diag(D)`
p <- ggplot(dist_df, aes(x=sex, y=differentiability))+ geom_violin(aes(fill=sex), alpha=0.65) +
     scale_fill_manual(values=c("#DCE319FF", "#33638DFF")) + 
     #stat_summary(fun.data = "mean_sdl", geom="pointrange", size=2, color="black") +
     geom_boxplot(width=0.25) +
     scale_y_continuous(limits = c(0, 2.5)) + 
     stat_pvalue_manual(stat_t, tip.length = 0) +
     labs(subtitle = get_test_label(stat_t, detailed = TRUE)) +
     theme_bw() + ylab("Self-distance") + xlab("") +
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
ggsave(file="figures/Differentiability_all_N2_sex_K30.pdf", plot=p, width=7, height=7)


# print also summary stats
dist_df <- dist_df %>% rename('WS-distance'='diag(D)')
dist_df['WS-distance^2'] = dist_df['WS-distance']**2
dist_df['differentiability^2'] = dist_df['differentiability']**2
stat_df <- dist_df[c('differentiability', 'age', 'sex')] %>% 
             group_by(sex) %>%
              summarise(mean(differentiability), sd(differentiability), mad(differentiability))
print(stat_df)

stat_df2 <- dist_df[c('differentiability', 'age', 'sex')] %>% 
  group_by(sex) %>%
  summarise(mean(differentiability), sd(differentiability), mad(differentiability))
print(stat_df2)


# then creating lm plot for checking the age-dep of WS-distances or differentiability measures
lm_eqn <- function(df){
  #m <- lm(age ~ `WS-distance`, df); #standard regression model
  #m <- lm(age ~ differentiability + `differentiability^2`, df)
  m <- lm(age ~ differentiability, df)
  #m <- lm(log_age ~ `WS-distance`, df); #exp. model
  #m <- lm(age~ `WS-distance`+ `WS-distance^2`,df)
  eq <- substitute(italic(y) == a + b %.% italic(x) + c %.% italic(x)^2*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~pvalue,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        c = format(unname(coef(m)[3]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        pvalue = format(summary(m)$coefficients[2,4], digits = 2)))
  as.character(as.expression(eq));
}

eqs <- dist_df %>% group_by(sex) %>% group_map(~lm_eqn(.x))
egs <- unlist(eqs)
eqs_df <- data.frame(matrix(eqs, nrow=length(eqs), byrow=TRUE), 
                     stringsAsFactors=FALSE)
colnames(eqs_df) <- "V1"
eqs_df['sex'] <- as.factor(c('F', 'M'))

p2 <- ggplot(dist_df, aes(x=age, y=differentiability)) +
      geom_point(alpha=0.15) +
      geom_smooth(data=dist_df, aes(color=sex, fill=sex), method=lm, 
                  #formula=y~x+I(x^2), se=TRUE) +
                  formula=y~x, se=TRUE) +
      facet_grid(~sex) +
      scale_fill_manual(values=c("#DCE319FF", "#33638DFF")) +
      scale_color_manual(values=c("#DCE319FF", "#33638DFF")) +
      geom_text(data=eqs_df, aes(x=8, y=1.2, label = V1), parse = TRUE, inherit.aes=FALSE) +
      theme_minimal() + ylab("Differentiability") + ylim(0.05,1.5) + xlab("Age") +
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"))
p2


ggsave(file="figures/Differentiability_age_sex_regplot_all_N1N2_K30.pdf", plot=p2, width=11, height=5)

################################################################################
#                Plot distance heatmaps for BRRRR model (fig  3)               #
################################################################################
library("reshape2")

# Note: will be needing a subset, o7 model for these visualizations

#p3 <- heatmap(D, Colv=NA, Rowv=NA)
## convert to tibble, add row identifier, and shape "long"
dat2 <- melt(D0[30:50, 30:50])
#D0 <- 1-D0
#dat2 <- melt(D0[50:70, 50:70])

p3 <- ggplot(dat2, aes(Var1, Var2, fill=value)) +
        geom_tile() +
        labs(fill='Corr') +
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
ggsave(file="figures/small4_distmat_o7_2N2_corr.pdf", plot=p3, width=11, height=10)



################################################################################
#           Mantell test for latent map matrices                               #
################################################################################
library("vegan")
mat1 <- viz_data1$lat_space #N1 sleep
mat2 <- viz_data2$lat_space #N2 sleep

K=30
D_list <- validation(within_sample = T, dis='L1', pK=c(1:K), lat_map=mat1, Xt=viz_data1$X)
D1 <- D_list[[1]]
# compute 0-model (correlations in full data matrix) 
D0_list <- validation(within_sample = T, dis='L1',  pK=c(1:K), lat_map=mat2, Xt=viz_data2$X)
D2 <- D0_list[[1]]  

# match subjects for uneven data sets
common_subjects = intersect(rownames(D1), rownames(D2)) 
to_keep = which(rownames(D1)%in%common_subjects)
D1 <- D1[to_keep, to_keep]
reorder_idx <- match(rownames(D1), rownames(D2))
D2 <- D2[reorder_idx, reorder_idx]


#stats <- mantel.test(D1, D2, graph=TRUE)
mantel_test <- mantel(xdis=D1, ydis=D2)
print(mantel_test)




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

# choose 200 random element indices from the matrix
set.seed(121)
rows <- sample(n_subj, size=200)
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











