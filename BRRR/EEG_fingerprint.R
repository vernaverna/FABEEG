# operate the process here...
setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("cvms")


## TODO: use less individuals! try with old and young kids as well

# Function for preparing the data for analysis
# ex = N1 or N2, depending which data is read in
prepare_data <- function(ex){
  
  # Choosing the channels and frequencies
  chs <- 1:19
  omitFreq = c(0,60)
  
  # Reading data into workspace
  loadfile <- paste0("data/",ex,"spectrum.RData")
  
  #extract the age groups
  ages = read.csv('data/age_df.csv')
  ages <- ages[,-1]
  
  use_all=FALSE #should we use all subjects in training?
  min_age=0 #exclude some of the younger children? 
  #TODO: change to range?
  
  # Check if the file exists
  if(file.exists(loadfile)) {
    load(loadfile)
    
    if(use_all==F){
      set.seed(121)
      individuals = sample(individuals, 40)
      Y = Y[names(Y) %in% individuals]
      
    }
    
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
        age_data = ages[ages$File==s,]
        if(age_data$Age >= min_age){
          tmp <- t(Y[[s]]) #transposed
          tmp <- log10(tmp) #TODO: is this necessary???
          if(any(is.na(tmp))) browser()
          A[[s]] <- tmp[frequencies,chs] 
        } else {
          corrupted = c(corrupted, s)
        }
      } else {
        corrupted = c(corrupted, s)
      } 
    }
    
    
    if(length(corrupted)>0){
      obs <- subjects[subjects %in% corrupted == FALSE]
      A = A[names(A) %in% corrupted == FALSE] 
    } else { 
      obs <- subjects
    }
    
  } else {
    print(paste0("File '",loadfile,"' does not exist, returning!"))
    return(0)
  }
  
  
  
  # The number of classes and storing the subjects
  M <- length(obs) #nclasses
  S <- length(obs) #nsubj
  
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
    x[obs[s]] <- s
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


n2_data <- prepare_data(ex="N2")
n2b_data <- prepare_data(ex="N2B")

validation_data <- prepare_data(ex="N1")
Y3 <- validation_data[[1]]

Y1 = n2_data[[1]]
Y2 = n2b_data[[1]]

Y = rbind(Y1, Y2)

x = c(n2_data[[2]], n2b_data[[2]])
X = rbind(n2_data[[3]], n2b_data[[3]])
M = n2_data[[4]]
subj = dimnames(Y)[[1]]


### TRAINING ###

# The model is trained using two sets of N2 data, and the performance is evaluated using third set

source("brrr.R")
pred <- X*NA
res <- brrr(X=X,Y=Y, K=6,n.iter=500,thin=5,init="LDA", fam = x) #fit the model
res$scaling <- ginv(averageGamma(res))
W <- res$scaling

save(res, file = "results/full/over7_indN2_BRRR_K6.RData")
lat_map <- Y%*%W
lat_map_n2 <- Y3%*%W #mapping to latent space with N2_C data!# 

D <- matrix(NA, nrow(lat_map_n2), ncol(X), dimnames=list(unique(names(x)), c())) #distance matrix# # # 

for(testidx in 1:nrow(lat_map_n2)){ #calculates the distances between individual and group mean in lat.space   
  for(m in 1:M){     
    group_members <- rownames(lat_map)[m]   
    idxs = which(row.names(lat_map) %in% group_members)#     
    group_mean <- colMeans(lat_map[idxs,]) #mean over all individuals, vector of length K
    D[testidx,m] <- sum(abs(lat_map_n2[testidx,]-group_mean)) #L1 distance#   
  } 
}



#### VALIDATION ####

PROJ <- D*0 #initialize 'projection' matrix

for(r in 1:nrow(D)){ #assign age groups based on lat. space distances
  index <- which.min(D[r,])
  PROJ[r,index] <- 1
}
colnames(PROJ) <- c(paste0("class",1:M))

accuracy = sum(diag(PROJ))/M
print(paste("Model accuracy:", accuracy))

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

