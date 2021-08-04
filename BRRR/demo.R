setwd("/projects/FABEEG/BRRR")

# ## Create artificial data
# N <- 500  #500 samples
# classes <- sample(c(-1,0,1), N, replace=TRUE) #Three classes
# W <- c(.2, .6, -.3) #Mapping from classes to observations
# x <- outer(classes, W) #The observations; outer product 
# x <- x[,] + rnorm(prod(dim(x))) #Noise
# 
# 
# # Learn BRRR model
# source("brrr.R")
# res <- brrr(matrix(classes, ncol=1), x, K=1,n.iter=100,init="penLDA")
# 
# # The true and inferred mapping (note: scaling is different due to Psi in BRRR)
# mapping <- data.frame(actual=W, final_posterior_sample=c(res$model$brr$context$Gamma))
# print(mapping)
# 
# #Posterior samples for Gamma. n.iter should be long enough so that these are converged
# # -- initialization will have effect as well
# GammaPosterior <- simplify2array(res$traces$Gamma)


# !!!!!!!!
#every value of the variance parameter of Ω can
#be immediately interpreted as the percentage of variance explained by the noise model as compared to the covariates

# Choosing the channels and frequencies
chs <- 1:19
ex="N2"

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
  
  if(!is.null(corrupted)){
    obs <- subjects[-which(subjects==corrupted)]
    A = A[names(A) %in% corrupted == FALSE]  
  } else { 
    obs <- subjects
  }
  
} else {
  print(paste0("File '",loadfile,"' does not exist, returning!"))
  return(0)
}



# The number of classes and storing the subjects
M <- length(unique(ages$Age.group)) 

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
if(model == "brrr") {
  tmp <- scale(Y[!is.na(Y[,1]),],center=T,scale=T)
  Y <- scale(Y, attr(tmp,"scaled:center"), attr(tmp,"scaled:scale"))
}

# Saving classes to x
x <- rep(NA,S); names(x) <- subjects

for(s in 1:S){
  x[individuals[[s]]] <- ages[ages$File==individuals[[s]],"Age.group"]
}

# When choosing BRRR x is modified into X. If there are just 2 classes, N x 1 matrix containing 
# identifier values +-1 is created
source("brrr.R")
X <- matrix(0,S,M,dimnames=list(subjects,paste0("class",1:M)))
for(i in 1:length(x)) if(!is.na(x[i])) X[i,x[i]] <- 1
if(M==2) {
  print("Condensing two classes into one +-1 covariate.")
  X <- X[,1,drop=FALSE] 
  X[X==0] <- -1
}


pred <- X*NA
res <- brrr(X=X,Y=Y, K=discTop,n.iter=n.iter,thin=5,init="LDA", fam = x) #fit the model
res$scaling <- ginv(averageGamma(res))
W <- res$scaling

pred <- X%*%res$model$brr$context$Psi%*%res$model$brr$context$Gamma #X%*%Psi%*%Gamma

#examine matrices
heatmap(Y)
heatmap(pred)
heatmap(Y%*%W)























for(i in obs) { 
  tmp <- A[[i]]
  if(nrow(tmp)==ncol(tmp)) tmp[lower.tri(tmp)] <- 0
  X[i,] <- c(tmp)
  y[i,] <- ages[ages$File==i,]$Age #get the age of the obs
}

keepFeat <- which(apply(X,2,var,na.rm=T)>0)
X <- X[,keepFeat]

print("LDA matrix dimension:")
print(dim(X))

# Fit BRRR model
res <- brrr(X=, Y=, K=6, n.iter=500, init="LDA")

final_posterior_gamma = res$model$brr$context$Gamma
final_posterior_psi = res$model$brr$context$Psi

#actual fit; weird?
psi = t(final_posterior_gamma)
gamma = t(final_posterior_psi)
omega = res$model$brr$context$Omega

Y_pred = (X%*%psi + omega)%*%gamma
print(Y_pred)
