setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("cvms")

#every value of the variance parameter of Ω can
#be immediately interpreted as the percentage of variance explained by the noise model as compared to the covariates

# Choosing the channels and frequencies
chs <- 1:19
ex="N2"
omitFreq = c(0,60)

# Reading data into workspace
loadfile <- paste0("data/",ex,"spectrum.RData")

#extract the age groups 
ages = read.csv('data/new_age_df.csv')
ages <- ages[,-1]
ages$Age.group <- round(ages$Age, 0)

# Group +15 year-olds together in group 15
old_idx = which(ages$Age.group >= 15)
ages$Age.group[old_idx] <- 15

use_all=TRUE #perform LOO-CV or use all data in training?

# Check if the file exists
if(file.exists(loadfile)) {
  load(loadfile)
  ages = subset(ages, ages$File %in% names(Y)) #remove subjects not in the dataset
  
  # Omits frequencies
  frequencies <- which(freq[,2] >= omitFreq[1] & freq[,2] <= omitFreq[2])
  freq <- freq[frequencies, ]
  
  # Searching the names/identifiers of the subjects
  subjects <- unlist(individuals)
  S <- length(subjects)
  A <- vector("list",length=S)
  names(A) <- subjects
  
  corrupted = c()
  
  # Saving the data into matrix A
  for(s in names(Y)) {
    
    if(s %in% ages$File){ #temporal solution to get over the filename hassle
      tmp <- t(Y[[s]]) #transposed
      #tmp <- log10(tmp) #TODO: is this necessary???
      if(any(is.na(tmp))) browser()
      A[[s]] <- tmp[frequencies,chs]
      
      
    } else {
      corrupted = c(corrupted, s)
    } 
  }
  
  #corrupted = c("FLE13539", "FLE141089") #this is bubblegum solution to remove singular groups
  
  if(length(corrupted)>0){
    obs <- subjects[subjects %in% corrupted == FALSE]
    A = A[names(A) %in% corrupted == FALSE] 
    ages = ages[ages$File %in% corrupted == FALSE,]
  } else { 
    obs <- subjects
  }
  
} else {
  print(paste0("File '",loadfile,"' does not exist, returning!"))
  return(0)
}

#Sampling even age groups
group_counts = table(ages$Age.group)
min_groupsize = group_counts[which.min(group_counts)] #find minimal age group


# The number of classes and storing the subjects; skip too small groups
group_counts = group_counts[c(which(group_counts >= min_groupsize))]

groups = names(group_counts)
M = length(groups)


for(m in groups){
  n_subj = as.numeric(group_counts[m]) #how many subjects are there
  age_group = ages[ages$Age.group==as.numeric(m),]
  
  if(n_subj-min_groupsize > 0){ #to skip the minimum group(s)
    extras = age_group$File[1:(n_subj-min_groupsize)] #excess subjects to remove
    obs <- obs[obs %in% extras == FALSE]
    A = A[names(A) %in% extras == FALSE] 
    ages = ages[ages$File %in% extras == FALSE,]
  }
}

S <- length(obs)

print("Data dimension per subject:")
print(dim(A[[1]]))
keepFeat <- NA

# Making MEG matrix X and response vec y and filling them
Y <- matrix(NA,length(A),prod(dim(A[[1]])),dimnames=list(names(A),c()))
x = matrix(NA, length(A), dimnames=list(names(A)))
colnames(Y) <- c(outer(paste0("s",1:nrow(A[[1]]),"."),1:ncol(A[[1]]),paste0))


for(i in obs) { 
  tmp <- A[[i]]
  if(nrow(tmp)==ncol(tmp)) tmp[lower.tri(tmp)] <- 0
  Y[i,] <- c(tmp)
}

keepFeat <- which(apply(Y,2,var,na.rm=T)>0)
Y <- Y[,keepFeat]

print("LDA matrix dimension:")
print(dim(Y))

# Scaling the Y for BRRR 
tmp <- scale(Y[!is.na(Y[,1]),],center=T,scale=T)
Y <- scale(Y, attr(tmp,"scaled:center"), attr(tmp,"scaled:scale"))


# Saving classes to x
x <- rep(NA,S); names(x) <- obs
for(s in 1:S){
  x[obs[s]] <- ages[ages$File==obs[s],"Age.group"]
}

# When choosing BRRR x is modified into X. If there are just 2 classes, N x 1 matrix containing 
# identifier values +-1 is created

source("brrr.R")
X <- matrix(0,S,M,dimnames=list(obs,paste0("class",0:(M-1)))) 
for(i in 1:length(x)) if(!is.na(x[i])) X[i,x[i]+1] <- 1 #classes go from 0 to 15, shift by one
if(M==2) {
  print("Condensing two classes into one +-1 covariate.")
  X <- X[,1,drop=FALSE] 
  X[X==0] <- -1
}


# Create LOO-CV for the small data set
# TODO: change pred to Yhat?
looCV <- function(X, Y, model){
  pred <- Y*NA
  #distance matrix to save the L1 distances in latent space Y%*%inv(Gamma)
  D <- matrix(NA, nrow(X), ncol(X), dimnames=list(names(x), c()))
  
  for(testidx in 1:nrow(X)){
    print("-------------")
    print(paste0("CV ",testidx,"/",nrow(X)))
    
    if(model=="brrr"){
      res <- brrr(X=X[-testidx,], Y=Y[-testidx,], K=3, n.iter=500, thin=5, fam=x[-testidx]) #res==mcmc.output
      pred[testidx,] <- X[testidx,,drop=F]%*%res$model$brr$context$Psi%*%res$model$brr$context$Gamma
      
      lat_space=Y%*%ginv(averageGamma(res)) #only makes sense to return this when trained with all data?
      
      for(testidx in 1:nrow(lat_space)){ #calculates the distances between individual and group mean in lat.space
        for(m in 1:M){
          group_members <- ages[ages$Age.group==m,]$File
          idxs = which(row.names(lat_space) %in% group_members)
          group_mean <- colMeans(lat_space[idxs,]) #mean over all individuals, vector of length K
          D[testidx,m] <- sum(abs(lat_space[testidx,]-group_mean)) 
        }
      }
      
      #TODO: distances for penLDA  
    } else if(model=="penlda"){
      class <- x[-testidx] #Class required to run from 1 to max
      class <- match(class,unique(class))
      xte <- Y[obs,] #should data and xte be passed vice versa?
      res2 <- PenalizedLDA(Y[-testidx,],class,xte=xte,
                          lambda=0,K=6,standardized=TRUE) #changed standardized to TRUE
      res2$scaling <- res2$discrim
      
      pred[testidx,] <- res2$xteproj%*%t(res2$discrim)
      lat_space = res2$xteproj #scaling is different!
    }
    
  } 
  
  return(list(pred, D, lat_space))
  
}


#' Function for visualizing a matrix
#' 
#' @param x The matrix to visualize
#' @param ... Further arguments passed to levelplot
visMatrix <- function(x, ...) {
  M <- max(abs(x))
  cols <- colorRampPalette(c("red","white","blue"))
  print(levelplot(x, col.regions=cols, at=seq(-M, M, length=50), aspect=1, ...))
}






######## BRRR #########

if(use_all==TRUE){
  pred <- X*NA
  res <- brrr(X=X,Y=Y, K=10,n.iter=500,thin=5,init="LDA", fam = x) #fit the model
  res$scaling <- ginv(averageGamma(res))
  W <- res$scaling
  lat_space=Y%*%W
  
  pred <- X%*%res$model$brr$context$Psi%*%res$model$brr$context$Gamma #X%*%Psi%*%Gamma
  D <- matrix(NA, nrow(X), ncol(X), dimnames=list(names(x), c())) #distance matrix
  

  for(testidx in 1:nrow(lat_space)){ #calculates the distances between individual and group mean in lat.space
    for(m in groups){
      group_members <- ages[ages$Age.group==m,]$File
      idxs = which(row.names(lat_space) %in% group_members)
      group_mean <- colMeans(lat_space[idxs,]) #mean over all individuals, vector of length K
      
      ix <- as.integer(m) + 1
      D[testidx,ix] <- sum(abs(lat_space[testidx,]-group_mean)) #L1 distance
    }
  }
  
  
  
  PROJ <- D*0
  
  for(r in 1:nrow(D)){ #assign age groups based on lat. space distances
    index <- which.min(D[r,])
    PROJ[r,index] <- 1
  }
  colnames(PROJ) <- c(paste0("class",0:(M-1)))
  
  #convert model matrices to factors and plot confusion matrix
  P_factor <- as.factor(colnames(PROJ))[PROJ %*% 1:ncol(PROJ)]
  X_factor <- as.factor(colnames(X))[X %*% 1:ncol(X)]

  
  #results are somewhat catastrophic
  png("figures/K6full_confmat_77.png")
  cmat = confusion_matrix(X_factor, P_factor)
  plot_confusion_matrix(cmat$`Confusion Matrix`[[1]])
  dev.off()  
  
  png("figures/K6full_Yinv(G)_77.png")
  heatmap(lat_space) # lat_space; these two should be the 
  dev.off()
  png("figures/K6full_XPsi_77.png")
  heatmap(X%*%res$model$brr$context$Psi)
  dev.off()
  #however, they are absolutely NOT
  
  #check the correlation:
  library("ggcorrplot")
  corrmat <- cor(as.vector(lat_space), as.vector(X%*%res$model$brr$context$Psi))
  print("Correlation matrix:")
  print(corrmat)
  ggcorrplot(corrmat)
  
  
  print("Overall correlation:")
  print(cor(c(lat_space), c(X%*%res$model$brr$context$Psi)))
  heatmap(corrmat)
  
  
  
} else if(use_all=FALSE){
  results <- looCV(X=X, Y=Y, model="brrr") #perform LOO-CV
  distmat <- results[[2]]
  PROJ <- distmat*0
  
  for(r in 1:nrow(distmat)){ #assign age groups based on lat. space distances
    index <- which.min(distmat[r,])
    PROJ[r,index] <- 1
  }
  colnames(PROJ) <- c(paste0("class",1:M))
  
  #convert model matrices to factors and plot confusion matrix
  P_factor <- as.factor(colnames(PROJ))[PROJ %*% 1:ncol(PROJ)]
  X_factor <- as.factor(colnames(X))[X %*% 1:ncol(X)]
  
  
  #results are somewhat catastrophic
  png("figures/LOO_confmat_316.png")
  cmat = confusion_matrix(X_factor, P_factor)
  plot_confusion_matrix(cmat$`Confusion Matrix`[[1]])
  dev.off()
}





heatmap(Y)
heatmap(results[[1]])
heatmap(Y%*%W) #ONLY WHEN TRAINED WITH FULL DATA
heatmap(Yhat%*%W)


######### LDA ###########

#using PenalizedLDA.cv instead of own LOO-function to tune model params
class=match(x, unique(x)) #this is for full data
xte=Y[obs,]
#cv_results <- PenalizedLDA.cv(Y, class, nfold=6) 
res2 <- PenalizedLDA(Y, class, xte=xte, lambda=0, K=14)
#png("figures/K6full_penLDA_confmat_77.png")
cmat = confusion_matrix(res2$y, res2$ypred[,14]) #results with 6 discriminant vectors used
plot_confusion_matrix(cmat$`Confusion Matrix`[[1]])
dev.off()


# testing testing 