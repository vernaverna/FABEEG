#This script has post hoc analysis methods etc.
# TODO: clean this up
setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("cvms")

## 1. CONFUSION MATRIX ##

#convert model matrices to factors and plot confusion matrix
P_factor <- as.factor(colnames(PROJ))[PROJ %*% 1:ncol(PROJ)]
X_factor <- as.factor(colnames(X))[X %*% 1:ncol(X)]


#results are somewhat catastrophic
cmat = confusion_matrix(X_factor, P_factor)
print("Accurcy:")
print(cmat$`Overall Accuracy`)
print("Specificity:")
print(cmat$Sensitivity)


png("figures/K6full_confmat_77.png")
plot_confusion_matrix(cmat$`Confusion Matrix`[[1]])
dev.off()  




## 2. COMPONENT PLOTTING ##

nsubj <- length(unique(subj))
# adding together N2 mappings  plus some covariates
lat_map = as.data.frame(rbind(lat_map_n2[1:nsubj,], lat_map))
lat_map['condition'] = c(rep('test', nsubj), rep('train', 2*nsubj))
lat_map['age'] = rep(ages[which(ages$File%in%subj),]$Age, 3) #get age data
lat_map['group'] = rep(round(ages[which(ages$File%in%subj),]$Age, 0), 3) #get age data
lat_map['sex'] = rep(ages[which(ages$File%in%subj),]$Sex, 3)

subj_number <- factor(c(x,seq(1,nsubj)))

ggplot(data=as.data.frame(lat_map), aes(lat_map[,6], lat_map[,7], shape=condition, col=factor(sex))) + 
  geom_point(aes(size=age)) + ggtitle("Subjects in latent mapping ") + 
  xlab("Component #1") + ylab("Component #2")


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


