setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("cvms")



prepare_data <- function(ex){
  
  # Choosing the channels and frequencies
  chs <- 1:19
  omitFreq = c(0,60)
  
  # Reading data into workspace
  loadfile <- paste0("data/",ex,"spectrum.RData")
  
  #extract the age groups
  ages = read.csv('data/age_df.csv')
  ages <- ages[,-1]
  
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
        tmp <- log10(tmp) #TODO: is this necessary???
        if(any(is.na(tmp))) browser()
        A[[s]] <- tmp[frequencies,chs]
        
        
      } else {
        corrupted = c(corrupted, s)
      } 
    }

    
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
  min_groupsize = nrow(ages[ages$Age.group==1,])
  
  # The number of classes and storing the subjects
  M <- length(unique(ages$Age.group)) 
  
  for(m in 1:M){
    age_group = ages[ages$Age.group==m,] #choose age group
    n_subj = nrow(age_group) #how many subjects are there
    
    if(n_subj-min_groupsize > 0){ #to skip the minimum group
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
  
  # Scaling & centering the Y for BRRR 
  tmp <- scale(Y[!is.na(Y[,1]),],center=T,scale=T)
  Y <- scale(Y, attr(tmp,"scaled:center"), attr(tmp,"scaled:scale"))
  
  
  # Saving classes to x
  x <- rep(NA,S); names(x) <- obs
  for(s in 1:S){
    x[obs[s]] <- ages[ages$File==obs[s],"Age.group"]
  }
  
  # Making design matrix out of classes
  X <- matrix(0,S,M,dimnames=list(obs,paste0("class",1:M))) 
  for(i in 1:length(x)) if(!is.na(x[i])) X[i,x[i]] <- 1
  if(M==2) {
    print("Condensing two classes into one +-1 covariate.")
    X <- X[,1,drop=FALSE] 
    X[X==0] <- -1
  }
  
  return(list(Y, x, X, M))
}


n1_data <- prepare_data(ex="N1")
n2_data <- prepare_data(ex="N2")

Y1 = n1_data[[1]]
Y2 = n2_data[[1]]

x = n2_data[[2]]
X = n2_data[[3]]
M = n2_data[[4]]
subj = dimnames(Y1)[[1]]


### TRAINING ###

# The model is trained using N2 data, but the performance is evaluated with N1 data

source("brrr.R")
pred <- X*NA
res <- brrr(X=X,Y=Y1, K=6,n.iter=500,thin=5,init="LDA", fam = x) #fit the model
res$scaling <- ginv(averageGamma(res))
W <- res$scaling
lat_space=Y2%*%W #validate the results using N1 data

pred <- X%*%res$model$brr$context$Psi%*%res$model$brr$context$Gamma #X%*%Psi%*%Gamma
D <- matrix(NA, nrow(X), ncol(X), dimnames=list(names(x), c())) #distance matrix

for(testidx in 1:nrow(lat_space)){ #fills the distance matrix
  for(m in 1:M){ #D[i,j] = avg. L1 dist between subject i and group j 
    D[testidx,m] <- mean(abs(lat_space[testidx,]-res$model$brr$context$Psi[m,])) 
  }
}


PROJ <- D*0 #initialize 'projection' matrix

for(r in 1:nrow(D)){ #assign age groups based on lat. space distances
  index <- which.min(D[r,])
  PROJ[r,index] <- 1
}
colnames(PROJ) <- c(paste0("class",1:M))

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




























