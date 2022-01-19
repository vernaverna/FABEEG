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
  ages = read.csv('data/new_age_df.csv')
  ages <- ages[,-1]
  ages$Age.group <- round(ages$Age, 0)
  
  # Group +15 year-olds together in group 15
  old_idx = which(ages$Age.group >= 15)
  ages$Age.group[old_idx] <- 15
  to_exclude=c(0,1,2,3) #c(0,1,2,3) #which age groups to exclude
  use_all=F #use all individuals?
  stabilize=T #get even age groups?
  
  # Check if the file exists
  if(file.exists(loadfile)) {
    load(loadfile)
    
    if(use_all==F){
      set.seed(121)
      individuals = sample(individuals, 500)
      Y = Y[names(Y) %in% individuals]
      
    }
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
        group = ages[ages$File==s,]$Age.group
        
        if(!(group %in% to_exclude)){
          tmp <- t(Y[[s]]) #transposed
          tmp <- log10(tmp) 
          if(any(is.na(tmp))) browser() #TODO: add to corrupted or something? in python??
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
  
  if(stabilize == T){
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
  X <- matrix(0,S,M,dimnames=list(obs,paste0("class",4:max(x) )) ) #TODO: I do not want to change these manually!
  for(i in 1:length(x)) if(!is.na(x[i])) X[i,x[i]-3] <- 1 #classes go from 4 to 15, shift by one
  if(M==2) {
    print("Condensing two classes into one +-1 covariate.")
    X <- X[,1,drop=FALSE] 
    X[X==0] <- -1
  }
  
  return(list(Y, x, X, M, groups, ages))
}



### READING IN THE DATA ###
n2_data <- prepare_data(ex="N2A")
n2b_data <- prepare_data(ex="N2B") #validation set

Y1 = n2_data[[1]]
Y2 = n2b_data[[1]]

#Y = rbind(Y1, Y2)

x = n2_data[[2]]
X = n2_data[[3]]
M = n2_data[[4]]
subj = dimnames(Y1)[[1]]
groups = n2_data[[5]]
ages = n2_data[[6]]


### TRAINING ###

L1=FALSE

#cosine similarity
cossim <- function(x, y){
  sum(x*y) / (sqrt(sum(x**2))*sqrt(sum(y**2)) )
}

# The model is trained using N2 data, the performance is evaluated with another data
omegas = c(1e+05,1e+04,1e+03,1e+02,1e+01,1,1e-01,1e-02,1e-03,1e-04,1e-05)
source("brrr.R")
for(o in omegas){
  res <- brrr(X=X,Y=Y1, K=5, n.iter=1000,thin=5,init="LDA", fam = x, omg=o) #fit the model  
}



res$scaling <- ginv(averageGamma(res))
W <- res$scaling
#save(res, file="results/full/age_over4_N2_BRRR_K6.RData")

lat_map <- Y2%*%W #mapping to latent space with N2 sleep
lat_comp <- X%*%res$model$brr$context$Psi + res$model$brr$context$Omega #latent space à la BRRR

D <- matrix(NA, nrow(X), ncol(X), dimnames=list(names(x), c())) #distance matrix


for(testidx in 1:nrow(lat_map)){ #calculates the distances between individual and group mean in lat.space
  for(m in groups){
    group_members <- ages[ages$Age.group==m,]$File
    idxs = which(row.names(lat_map) %in% group_members)
    group_mean <- colMeans(lat_map[idxs,]) #mean over all individuals, vector of length K
    
    ix <- as.integer(m) -3 #+1
    #D[testidx,ix] <- cossim(lat_map[testidx,],group_mean) #cosine similarity
    D[testidx,ix] <- sum(abs(lat_map[testidx,]-group_mean)) #L1 distance
  }
}



PROJ <- D*0

for(r in 1:nrow(D)){ #assign age groups based on lat. space distances
  if(L1){
    index <- which.min(D[r,])
  } else {
    index <- which.max(D[r,])
  }
  PROJ[r,index] <- 1
}
colnames(PROJ) <- c(paste0("class",4:max(x)))

#convert model matrices to factors and plot confusion matrix
P_factor <- as.factor(colnames(PROJ))[PROJ %*% 1:ncol(PROJ)]
X_factor <- as.factor(colnames(X))[X %*% 1:ncol(X)]



#results are somewhat catastrophic
cmat = confusion_matrix(X_factor, P_factor)
print("Accurcy:")
print(cmat$`Overall Accuracy`)

#create more informative classification metric
# IDEA: hard-code error rates and then MAE per class?
Tab=cmat$Table[[1]]

# ja nyt vähän pseudokoodia
penalty_vec = rep(0, ncol(Tab))
names(penalty_vec) <- colnames(Tab)
for(i in 1:ncol(Tab)){
  target = as.integer(sub("class","",colnames(Tab)[i]))
  predlabs = unlist(lapply(names(Tab[,i]), function(x) {as.integer(sub("class","", x))} ))
  
  print(mean(predlabs))
  predvec = Tab[,i]
  penalties = -(0.1)*abs(predlabs-target)
  classpred = t(rbind(predlabs, predvec, penalties))
  
  which(classpred$predlabs==classpred$target)
  
  total_penalty = predvec*penalties
  penalty_vec[i] = sum(10^(penalties)*predvec)/13
}
penalty_vec


png("figures/K15full_confmat_208.png")
plot_confusion_matrix(cmat$`Confusion Matrix`[[1]])
dev.off()  


#' Function for visualizing a matrix
#' 
#' @param x The matrix to visualize
#' @param ... Further arguments passed to levelplot
visMatrix <- function(x, ...) {
  M <- max(abs(x))
  cols <- colorRampPalette(c("red","white","blue"))
  print(levelplot(x, col.regions=cols, at=seq(-M, M, length=50), aspect=1, ...))
}


visMatrix(res$model$brr$context$Psi, xlab="Feature #", ylab="Component #", main="Inferred matrix Psi")

visMatrix(res$model$brr$context$Gamma, xlab="Component #", ylab="Output #", main="Inferred matrix Gamma")
visMatrix(t(W), xlab="Component #", ylab="Output #", main="Averaged matrix inv-Gamma")

visMatrix(Y2%*%W, xlab="Component #", ylab="Output #", main="Projection to latent space")



#### SANITY CHECKS ####
visMatrix(lat_comp) #these are a-ok
visMatrix(lat_map)

model_diff <- abs(lat_map - lat_comp) #scaling is different, why?
print(paste0("this should be about 0: ", sum(model_diff))) #L1 distance again; scaling causes some error??

# adding together ec & eo mappings
lat_map = as.data.frame(rbind(lat_map_ec[1:nsubj,], lat_map_eo))
lat_map['condition'] = c(rep('ec', nsubj), rep('eo', nsubj))


ggplot(data=as.data.frame(lat_map), aes(lat_map[,1], lat_map[,2], shape=condition, col=factor(x))) + 
  geom_point(size=3) + ggtitle("Subjects in latent mapping ") + 
  xlab("Component #1") + ylab("Component #2")



######### LDA ###########

#using PenalizedLDA.cv instead of own LOO-function to tune model params
class=match(x, unique(x)) #this is for full data
xte=Y1[subj,]
#cv_results <- PenalizedLDA.cv(Y, class, nfold=6) 
res2 <- PenalizedLDA(Y1, class, xte=xte, lambda=0, K=6)
#png("figures/K6full_penLDA_confmat_77.png")
cmat = confusion_matrix(res2$y, res2$ypred[,6]) #results with 6 discriminant vectors used
plot_confusion_matrix(cmat$`Confusion Matrix`[[1]])
print("Accurcy:")
print(cmat$`Overall Accuracy`)

dev.off()



#TODO: try t-SNE/ Sammons' mapping for the same data!

library("Rtsne")
colors = rainbow(length(unique(x)))
names(colors) = unique(x)

tsne <- Rtsne(lat_map[,1:10], dims=2, perplexity=5, verbose=TRUE, max_iter = 500, check_duplicates = F)
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=x, col=colors[x])

plot(tsne$Y,col=lat_map$group, asp=1)











