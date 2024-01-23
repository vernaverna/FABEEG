# collection of plots for teh article.

setwd("/projects/FABEEG/BRRR/")
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

load("results/full/all_2N1_BRRR_12.RData")
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
tsne <- Rtsne(D, perplexity=200, theta=0.0, check_duplicates = FALSE, max_iter=2000) 

nsubj <- length(unique(subj))
reps = 2
# adding together N2 mappings  plus some covariates
lat_map = as.data.frame(rbind(lat_space))
lat_map['X1'] = tsne$Y[,1]
lat_map['X2'] = tsne$Y[,2]
lat_map['spectra'] = c(rep('N1B', nsubj), rep('N1A', nsubj))
lat_map['log_age'] = rep(log10(ages$Age), reps) #get age data
lat_map['age'] = rep(ages$Age, reps) #get age data
lat_map['group'] = rep(round(ages$Age, 0), reps) #get age data
lat_map['sex'] = rep(ages$Sex, reps)
lat_map['cap'] = rep(ages$Cap, reps)
lat_map['subject'] = subj[1:nsubj]



p <- ggplot(data=lat_map, aes(x=V2, y=V3, shape=spectra, colour=age)) + 
            geom_point(alpha=0.3, size=2) +
            ggtitle("Subjects in latent mapping ") + 
            scale_color_viridis() +
            xlab("t-SNE #1") + ylab("t-SNE #2") + #xlim(-42,48) + ylim(-20,25) +
            theme(legend.position="none") + theme_bw() + 
            theme(panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black"))
p
ggsave("figures/NEW_latmap_all_2N2_Age.pdf", width=7.4, height=5.2)

q <- ggplot(data=lat_map, aes(x=V2, y=V3, shape=spectra, colour=sex)) + 
  geom_point(alpha=0.5, size=3) +
  ggtitle("Subjects in latent mapping ") + 
  scale_color_viridis(discrete=TRUE) +
  xlab("t-SNE #1") + ylab("t-SNE #2") + #xlim(-42,48) + ylim(-20,25) +
  theme(legend.position="none") + theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
q
ggsave("figures/NEW_latmap_all_2N2_Sex.pdf", width=7.4, height=5.2)


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
                       linkage='average'){
  
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


K=12
D_list <- validation(within_sample = T, dis='L1', pK=c(1:K), lat_map=lat_space, Xt=X, linkage = 'average')
D <- D_list[[1]]
# compute 0-model (correlations in full data matrix) 
D0_list <- validation(within_sample = T, dis='corr',  pK=c(1:247), lat_map=Y, Xt=X)
  
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
ggsave(file="figures/WS_distances_all_2N2_sex.pdf", plot=p, width=7, height=7)



# then creating lm plot for checking the age-dep of WS-distances
lm_eqn <- function(df){
  m <- lm(age ~ `diag(D)`, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~pvalue,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
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

p2 <- ggplot(dist_df, aes(x=age, y=`diag(D)`)) +
      geom_point(alpha=0.15) +
      geom_smooth(data=dist_df, aes(color=sex, fill=sex), method=lm, se=TRUE) +
      facet_grid(~sex) +
      scale_fill_manual(values=c("#DCE319FF", "#33638DFF")) +
      scale_color_manual(values=c("#DCE319FF", "#33638DFF")) +
      geom_text(data=eqs_df, aes(x=10, y=5, label = V1), parse = TRUE, inherit.aes=FALSE) +
      theme_bw() + ylab("W-S distance") +
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"))
p2


ggsave(file="figures/WS_age_sex_regplot_all_2N1.pdf", plot=p2, width=11, height=5)


