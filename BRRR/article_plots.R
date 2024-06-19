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
  
  return(data_list)
}


viz_data1 <- get_viz_data(fname = "results/full/all_2N1_BRRR_K30.RData")
viz_data2 <- get_viz_data(fname= "results/full/all_2N1_BRRR_K30.RData")



X <- viz_data2$X
Y <- viz_data2$Y
subj <- viz_data2$subj
ages <- viz_data2$ages
age_df <- viz_data2$age_df
lat_space <- viz_data2$lat_space
lat_map <- viz_data2$lat_map

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
tsne <- Rtsne(D, perplexity=90, theta=0.0, check_duplicates = FALSE, max_iter=2000) 

nsubj <- length(unique(subj))
reps = 2
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

lat_map <- na.omit(lat_map)

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
            dists <- apply(other_data, 1, function(x) dist_func(lat_map[testidx,pK], x))
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


K=30
D_list <- validation(within_sample = T, dis='L1', pK=c(1:K), lat_map=lat_space, Xt=X)
D <- D_list[[1]]
# compute 0-model (correlations in full data matrix) 
D0_list <- validation(within_sample = T, dis='corr',  pK=c(1:247), lat_map=Y, Xt=X)
D0 <- D0_list[[1]]  
# Creating a dataframe for plotting....
reorder_val_idx <- match(rownames(D), ages$File)
ages <- ages[reorder_val_idx,] 

dist_df = as.data.frame(diag(D), row.names = rownames(D))
dist_df['age'] = ages$Age #get age data
#dist_df['group'] = rep(round(ages$Age, 0), reps) #get age data
dist_df['sex'] = ages$Sex
dist_df['cap'] = ages$Cap

# remove NA -rows before muting the vars
dist_df <- na.omit(dist_df)
dist_df['log_age'] = log10(dist_df$age) #get age data

# compare if significant difference between F vs M using t-test
stat_t <- dist_df %>% t_test(`diag(D)`~sex) %>% add_significance()
stat_t <- stat_t %>% add_xy_position(x = "sex")


p <- ggplot(dist_df, aes(x=sex, y=`diag(D)`)) + geom_violin(aes(fill=sex), alpha=0.65) +
     scale_fill_manual(values=c("#DCE319FF", "#33638DFF")) + 
     #stat_summary(fun.data = "mean_sdl", geom="pointrange", size=2, color="black") +
     geom_boxplot(width=0.1, fill="white") +
     stat_pvalue_manual(stat_t, tip.length = 0) +
     labs(subtitle = get_test_label(stat_t, detailed = TRUE)) +
     theme_bw() + ylab("W-S distance") + xlab("") +
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
ggsave(file="figures/WS_distances_all_2N2_sex_K30.pdf", plot=p, width=7, height=7)


# print also summary stats
dist_df <- dist_df %>% rename('WS-distance'='diag(D)')
dist_df['WS-distance^2'] = dist_df['WS-distance']**2
stat_df <- dist_df[c('WS-distance', 'age', 'sex')] %>% 
             group_by(sex) %>%
              summarise(mean(`WS-distance`), sd(`WS-distance`), mad(`WS-distance`))
print(stat_df)

# then creating lm plot for checking the age-dep of WS-distances
lm_eqn <- function(df){
  #m <- lm(age ~ `WS-distance`, df); #standard regression model
  #m <- lm(log_age ~ `WS-distance`, df); #exp. model
  m <- lm(age~ `WS-distance`+ `WS-distance^2`,df)
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

p2 <- ggplot(dist_df, aes(x=age, y=`WS-distance`)) +
      geom_point(alpha=0.15) +
      geom_smooth(data=dist_df, aes(color=sex, fill=sex), method=lm, 
                  formula=y~x+I(x^2), se=TRUE) +
      facet_grid(~sex) +
      scale_fill_manual(values=c("#DCE319FF", "#33638DFF")) +
      scale_color_manual(values=c("#DCE319FF", "#33638DFF")) +
      geom_text(data=eqs_df, aes(x=6, y=7, label = V1), parse = TRUE, inherit.aes=FALSE) +
      theme_bw() + ylab("W-S distance") + ylim(0.5,7.5) +
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"))
p2


ggsave(file="figures/WS_age_sex_quat_regplot_all_2N1_K30.pdf", plot=p2, width=11, height=5)

################################################################################
#                Plot distance heatmaps for BRRRR model                        #
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
load("results/full/i1000_all_2N2_BRRR_K30.RData")
res1 <- res
load("results/full/i5000_all_2N2_BRRR_K30.RData")
res2 <- res
remove(res)

# is it legal to compare averages of over mcmc samples?
#coefmat_1 <- averagePsi(res1)%*%averageGamma(res1)
#coefmat_2 <- averagePsi(res2)%*%averageGamma(res2)

# Traces... there is either 100 (for i=1000) or 500 (for i=5000)
traces_1 <- res1$traces
traces_2 <- res2$traces

tr_psi1 <- traces_1$Psi
tr_gamma1 <- traces_1$Gamma
ptves <- traces_1$tpve
iter <- traces_1$iter

# because of the rotation invariance, the convergence needs to be studied w.r.t. to
# the standard regression coefficient matrix Psi*Gamma

# pseudocode

# for all iterations, compute theta esitmate (Lapply???)
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

# choose 200 random element indices from the matrix
set.seed(121)
rows <- sample(788, size=200)
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
plot(iter, ptves, xlab='iteration', 'l', col='darkorchid4', ylab="PTVE score", ylim = c(0.84312, 0.84355))
text(870, 0.8435, paste0("Rhat=", round(Rhat, 4), ", ESS_b=", round(ESS_b, 4),
                         ", ESS_t=", round(ESS_t, 4)))
dev.off()
#legend(950, 0.8432, legend='PTVE', col="yellow3", pch=16)











