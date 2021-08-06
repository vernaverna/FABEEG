setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")

#every value of the variance parameter of Ω can
#be immediately interpreted as the percentage of variance explained by the noise model as compared to the covariates

# Choosing the channels and frequencies
chs <- 1:19
ex="N2"
omitFreq = c(0,60)

# Reading data into workspace
loadfile <- paste0("data/",ex,"spectrum.RData")

#extract the age groups
ages = read.csv('data/age_df.csv')
ages <- ages[,-1]


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
  
  # Saving the data into matrix A
  for(s in names(Y)) {
    
    if(s %in% ages$File){ #temporal solution to get over the filename hassle
      tmp <- t(Y[[s]]) #transposed
      tmp <- log10(tmp) #TODO: is this necessary???
      if(any(is.na(tmp))) browser()
      A[[s]] <- tmp[frequencies,chs]
      
      corrupted=NULL
      
    } else {
      corrupted = s
    } 
  }
  
  corrupted = c("FLE13539", "FLE141089") #this is bubblegum solution to remove singular groups
  
  if(!is.null(corrupted)){
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



# The number of classes and storing the subjects
M <- length(unique(ages$Age.group)) 
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
# TODO: write warnings for case with singular groups!
source("brrr.R")
X <- matrix(0,S,M,dimnames=list(obs,paste0("class",3:(M+2)))) #NOTE! shifted because 2 first age groups removed
for(i in 1:length(x)) if(!is.na(x[i])) X[i,x[i]-2] <- 1
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
      res <- brrr(X=X[-testidx,], Y=Y[-testidx,], K=3, n.iter=100, thin=5, fam=x[-testidx]) #res==mcmc.output
      pred[testidx,] <- X[testidx,,drop=F]%*%res$model$brr$context$Psi%*%res$model$brr$context$Gamma
      
      lat_space=Y%*%ginv(averageGamma(res))
      
      for(m in 1:M){
        D[testidx,m] <- mean(abs(lat_space[testidx,]-res$model$brr$context$Psi[m,]))
      }
    
      
      #TODO: distances for penLDA  
    } else if(model=="penlda"){
      class <- x[-testidx] #Class required to run from 1 to max
      class <- match(class,unique(class))
      xte <- Y[obs,] #should data and xte be passed vice versa?
      res2 <- PenalizedLDA(Y[-testidx,],class,xte=xte,
                          lambda=0,K=3,standardized=TRUE) #changed standardized to TRUE
      res2$scaling <- res2$discrim
      
      pred[testidx,] <- res2$xteproj%*%t(res2$discrim)
      lat_space = res2$xteproj #scaling is different!
    }
    
  } 
  
  return(list(pred, D))
  
}




######## BRRRR #########

pred <- X*NA
res <- brrr(X=X,Y=Y, K=3,n.iter=500,thin=5,init="LDA", fam = x) #fit the model
res$scaling <- ginv(averageGamma(res))
W <- res$scaling

pred <- X%*%res$model$brr$context$Psi%*%res$model$brr$context$Gamma #X%*%Psi%*%Gamma


results <- looCV(X=X, Y=Y, model="brrr") #perform LOO-CV
distmat <- results[[2]]
PROJ <- distmat*0

for(r in 1:nrow(distmat)){ #assign age groups based on lat. space distances
  index <- which.min(distmat[r,])
  PROJ[r,index] <- 1
}
colnames(PROJ) <- c(paste0("class",3:(M+2)))

#TODO: construct a confusion matrix
#results are somewhat catastrophic

heatmap(Y)
heatmap(Yhat)
heatmap(Y%*%W) #ONLY WHEN TRAINED WITH FULL DATA
heatmap(Yhat%*%W)


######### LDA ###########



Ydot <- looCV(X=X, Y=Y, model="penlda")





