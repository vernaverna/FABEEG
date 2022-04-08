#This script has post hoc analysis methods etc.
# TODO: clean this up
setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("dplyr")
library("cvms")
library("RColorBrewer")
library("corrplot")


## Test component simlarities ##

load("results/full/over7_2N2_BRRR_K12.RData") #try with the best ptve models tho!
comp1 <- res$scaling
res2 <- load("results/full/over7_2N1_BRRR_K12.RData")
comp2 <- res$scaling
rm(res)

cossim <- function(x,y){ #calculate cossimilarity between vectors
  return( sum(x*y) / (sqrt(sum(x**2))*sqrt(sum(y**2))) )
}

for(i in 1:ncol(comp1)){
  print(paste0("cosine similarity of K",i))
  print( cossim(comp1[,i], comp2[,i]) )
}

M <- matrix(NA, nrow=ncol(comp1), ncol=ncol(comp1))
for(i in 1:ncol(comp1)){
  for(j in 1:ncol(comp2)){
    M[i,j] <- cossim(comp1[,i], comp2[,j])
  }
}
rownames(M) <- c(paste0('A_K', c(1:12))) #rows=res1
colnames(M) <- c(paste0('B_K', c(1:12)))

#maybe better to use abs values with gradient palette
heatmap(abs(M),  Colv = NA, Rowv = NA, col=brewer.pal(8,"BuPu"), 
        xlab="model B",ylab="model A",main="Heatmap of cosine distances") 

P=cor(Ys, method = 'spearman')
rownames(P) <- plotLabels
colnames(P) <- plotLabels
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(P, method = "color", col=col(200), type = "lower",
         title="Correlations between global PSDs of one subject",
         addCoef.col = "black", mar=c(0,0,1,0))

  ## WORKING WITH FULL DATA MODEL ##  
load("results/full/over1_ind_2N2_BRRR_12.RData")
Y <- res$data$phenotypes
X <- res$data$genotypes
subj <- row.names(X)
ages = read.csv('data/new_age_df.csv')
ages <- ages[,-1]
ages = ages[ages$File %in% unique(subj)==TRUE, ] #to get rid of some extra subjects that should not be there



inv_G <- res$scaling #inv(average(Gamma))
inv_PG <- res$scaling2 #inv(average(Gamma)%*%average(Psi))

p_G <- res$model$brr$context$Gamma #posterior Gamma
p_P <- res$model$brr$context$Psi #posterior Psi

Y_hat <- X%*%p_P%*%p_G
lat_space <- Y%*%inv_G #obs x latent components
lat_proj <- Y%*%inv_PG #has same dimensions as X



source("visNetwork.R")
net <- list(omitMag=omitMag,Y=res$scaling,keepFeat=keepFeat,
            penLDA=4,omitN=omitN,l1d=l1d)
net$lambda <- lambda

net$freq <- freq
net$fname <- fname

save(net,file=fname)
visNetwork(net,onlyPdf=TRUE)
visNetwork(net,onlyPdf=TRUE,levelplot=TRUE)



## 2. COMPONENT PLOTTING ##

nsubj <- length(unique(subj))
# adding together N2 mappings  plus some covariates
lat_map = as.data.frame(rbind(lat_space))
lat_map['spectra'] = c(rep('N2A', nsubj), rep('N2B', nsubj))
lat_map['age'] = rep(ages[which(ages$File%in%subj),]$Age, 2) #get age data
lat_map['group'] = rep(round(ages[which(ages$File%in%subj),]$Age, 0), 2) #get age data
lat_map['sex'] = rep(ages[which(ages$File%in%subj),]$Sex, 2)
lat_map['cap'] = rep(ages[which(ages$File%in%subj),]$Cap, 2)
lat_map['subject'] = subj[1:nsubj]

gender <- ifelse(lat_map$sex == "F",1,0)

#Create a latent space mapping, highlighting only certain individuals

#get colorblind-friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


lat_map_mini=lat_map[c(1:8, 592:599), ] #make highlight points

ggplot(data=as.data.frame(lat_map), 
       aes(lat_map[,1], lat_map[,2], shape=spectra)) + 
  geom_point(alpha=0.2) + geom_point(data=lat_map_mini, 
                                     aes(lat_map_mini[,1], lat_map_mini[,2], shape=spectra, colour=subject),
                                     size=3) +
  ggtitle("Subjects in latent mapping ") + scale_colour_manual(values=cbbPalette)+
  xlab("Component #1") + ylab("Component #2") + xlim(-45,55) + ylim(-25,30) +
  theme(legend.position="none") + theme_bw()


# trying 3D plots

library("scatterplot3d")
library("viridis")

shapes <- c(16, 17)
shapes <- shapes[factor(lat_map$condition)]
colors <- viridis_pal(option = "D")(length(unique(lat_map$group)))
colors <- colors[factor(lat_map$group)]

scatterplot3d(lat_map[,3:5], pch=shapes, color=colors, angle=20)

#trying t-SNE

library("Rtsne")
age_col = round(ages[which(ages$File%in%subj),]$Age, 0)
sex_col = ages[which(ages$File%in%subj),]$Sex

colors = rainbow(length(unique(age_col)))
names(colors) = unique(age_col)

tsne <- Rtsne(D, dims=2, is_distance = T,
              perplexity=4, verbose=TRUE, max_iter = 500, check_duplicates = F)
#plot(tsne$Y, t='n', main="tsne")
#text(tsne$Y, labels=x, col=colors[x])

plot(tsne$Y,col=colors, asp=1)



## 3. GETTING AGE OUT ##

#trying distance matrix with GGplot
# TODO: something is wrong
data_df <- data.frame(tsne$Y)
rownames(data_df) <- rownames(D)
data_df$age <- ages[which(ages$File%in%subj),]$Age
data_df$sex <- ages[which(ages$File%in%subj),]$Sex
data_df$group <- round(ages[which(ages$File%in%subj),]$Age, 0)

ggplot(data_df, aes(x=X1, y=X2, color=age, shape=sex)) +
  geom_point(size=4) + theme_minimal()


#how distance correlates with age?
# -> make a matrix with shape identical to D, fill with age diffs and 
#    find correlation between these matrices (or each row vector I guess?)


G <- D*0
D_ages <- ages[ages$File %in% rownames(D),]
#match the ordering that is mixed due to random sampling
reorder_indexes <- match(rownames(D), D_ages$File)
D_ages <- D_ages[reorder_indexes,]

for(r in rownames(D)){
  age_diffs = abs(D_ages[D_ages$File==r,]$Age - D_ages$Age) #age difference
  G[r,] <- age_diffs
}
C=cor(t(D), t(G)) # calculates correlation between rows of G and D
heatmap(C)

cors=c()
for(r in 1:nrow(G)){cors=c(cors, cor(G[r,], D[r,]))}
D_ages$cors <- cors


#try hierarchical clustering with dissimilarity matrix
library("dendextend")
heatmap(S)

colnames(S) <- D_ages$Age
rownames(S) <- D_ages$Age

hc <- hclust(as.dist(S), method="average")
plot(hc)
dend <- as.dendrogram(hc)
dend %>% set("branches_k_color", k = 5) %>% plot(main = "Nice defaults")



# Or how about k-means?
install.packages("ClusterR")
install.packages("cluster")




# Plotting dependence of train-PTVE, accuracy and lat. space dimension

#K=c(6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)
#Accuracy=c(0.5026919,0.6054551,0.6409343,0.6762854,0.6796895,0.7216066,0.7083464,0.7184304,
#           0.7436832,0.7470019,0.7739353,0.750349,0.770631,0.7841903,0.7993876,0.780715,
#           0.7808432,0.8077339,0.780715,0.8061672,0.8010255,0.8110668,0.8026207)
#PTVE=c(0.783122,0.7989688,0.8105054,0.8198825,0.8272988,0.8337677,0.8387554,0.8431643,
#       0.8468616,0.8499344,0.8528079,0.8554512,0.8575301,0.8595379,0.8612431,0.8628522,
#       0.8643616,0.8657661,0.8670517,0.868357,0.8694523,0.870628,0.871647)

P<-c(80,180,280,380,480,580,680,780)
Accuracy<-c(0.3375,0.5592593,0.5396825,0.552924,0.537326,0.5421456,0.5395425,0.4837607) #for the different Ns; K=12
PTVE <-c(0.7839801,0.7560055,0.7771796,0.7729471,0.772513275,0.7694992,0.768224,0.7740776)


library(reshape2)
K <- c(1:12)
data=as.data.frame(cbind(K, N1, N2))
melt(data, id.vars = "K", measure.vars = c("N1", "N2"))
colnames(data) <- c("K", "model", "correlation")

ggplot(data=data, aes(x=K, y=correlation, color=model) ) + geom_point(size=4) + 
  scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
  #geom_point(data=data, aes(x=K, y=N2), color='darkorange', size=4) +
  geom_text(
    label=data$correlation, 
    nudge_x = 0.1, nudge_y = 0.05, 
    check_overlap = T
  ) + 
 # xlim(0, 12) +
  theme_minimal()

plot(Accuracy~P, type='l', col='darkorange', ylim=c(0.30, 0.8),
     ylab="", lwd=3, bty='n')
lines(PTVE~P, col='darkcyan', lwd=3, bty='n')
grid(nx = NA,
     ny = NULL,
     lty = 2, col = "gray", lwd = 0.8)

legend('bottomright', legend = c('Accuracy', 'PTVE'),  
       col=c('darkorange', 'darkcyan'), pch=19, bty = "n", 
       pt.cex = 1.8, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , inset = c(0.07, 0.07) )

