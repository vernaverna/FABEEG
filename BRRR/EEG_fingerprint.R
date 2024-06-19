# operate the process here...
setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("cvms")



#' Function for preparing the EEG spectra for analysis
#' 
#' @param spectra list, tells which spectra are used for the model: e.g. c('N1C', 'N2B')
#' @param validation_set which spectra to use in validation, e.g. 'N2C'
#' @param n_inds is using only a subset of data, how many individuals to use?
#' @param group_by_spectra boolean, classify spectra rather than individuals?
#' @param data_type either 'spectra' or 'PSD'
#' 
# ex = N1 or N2, depending which data is read in
prepare_data <- function(spectra, validation_set, n_inds=180, 
                         group_by_spectra = F, data_type='spectrum',seed=191){
  
  
  # Choosing the channels and frequencies
  chs <- 1:19
  omitFreq = c(0,60) #so takes in all
  
  #extract the age groups
  ages = read.csv('data/new_age_df.csv')
  ages <- ages[,-1]
  ages[ages==" "] <- NA #replace empty strings with NA-values
  
  use_all=T #should we use all subjects in training?
  age_gap=c(0,19) #ages to include
  #Cap='-'
  
  data_Y = vector(mode='list',length=length(spectra)) #containers for targets Y and covariates X
  data_X = vector(mode='list',length=length(spectra))
  
  names(data_Y) <- spectra
    
  #Go through each data file
  for(i in 1:length(spectra)){
    # Reading data into workspace
    spec = spectra[i]
    loadfile <- paste0("data/new_",spec, data_type,".RData")
    
    # Check if the file exists
    if(file.exists(loadfile)){
      load(loadfile)
      if(use_all==F){
        set.seed(seed)
        #load("/projects/FABEEG/BRRR/ref_subjects.RData")
        #individuals=obs
        individuals = sample(individuals, n_inds)
        Y = Y[names(Y) %in% names(individuals)]
        
      }
      set.seed(seed) #for reproducibility.
      if(use_all==T){
        individuals = sample(individuals, length(individuals)) #shuffle the data, just in case. 
        Y = Y[names(Y) %in% names(individuals)]
      }

      
      # Omits frequencies
      if(data_type=='spectrum'){
        frequencies <- which(freq[,2] >= omitFreq[1] & freq[,2] <= omitFreq[2])
        freq <- freq[frequencies, ]        
      } else {
        freq <- freq$V1
        frequencies=1:length(freq)
      }

      
      
      # Searching the names/identifiers of the subjects
      #subjects <- unlist(individuals)
      subjects <- names(individuals)
      S <- length(subjects)
      A <- vector("list",length=S)
      names(A) <- subjects
      
      corrupted = c() #container for files that are dropped
      
      # Saving the data into matrix A
      for(s in names(Y)) {
        
        if(s %in% ages$File){ #temporal solution to get over the filename hassle
          age_data = ages[ages$File==s,]
          if(!is.na(age_data$Age)){
            if(age_data$Age >= age_gap[1] && age_data$Age <= age_gap[2]){
              #if(age_data$Cap == Cap){
              tmp <- t(Y[[s]]) #transposed
              if(data_type!= 'spectrum'){
                tmp <- Y[[s]]
              }
              tmp <- log10(tmp) #TODO: is this necessary?
              if(any(is.na(tmp))){ #browser()
                corrupted = c(corrupted, s)
              } 
              A[[s]] <- tmp[frequencies,chs] 
            } else {
              corrupted = c(corrupted, s)
            }
            
          } else { #include the non-aged anyway?
            tmp <- t(Y[[s]]) #transposed
            if(data_type!='spectrum'){
              tmp <- Y[[s]]
            }
            tmp <- log10(tmp) #TODO: is this necessary?
            A[[s]] <- tmp[frequencies,chs] 
            
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
    M <- length(obs) #nclasses #alright this is how I like it ;)
    S <- length(obs) #nsubj
    
    print("Data dimension per subject:")
    print(dim(A[[1]]))
    keepFeat <- NA
    
    # Making MEG matrix X and response vec y and filling them
    Y <- matrix(NA,length(A),prod(dim(A[[1]])),dimnames=list(names(A),c()))
    x = matrix(NA, length(A), dimnames=list(names(A)))
    colnames(Y) <- c(outer(paste0("s",1:nrow(A[[1]]),"."),1:ncol(A[[1]]),paste0))
    
    
    for(j in obs) { 
      tmp <- A[[j]]
      if(nrow(tmp)==ncol(tmp)) tmp[lower.tri(tmp)] <- 0
      Y[j,] <- c(tmp)
    }
    
    keepFeat <- which(apply(Y,2,var,na.rm=T)>0)
    Y <- Y[,keepFeat]
    
    print("LDA matrix dimension:")
    print(dim(Y))
    
    data_Y[[i]] <- Y
    
  }

  # Combine data into 1 big matrix

  val=which(names(data_Y)==validation_set)
  validation_data = data_Y[[val]] #extract validation data
  data_Y[[val]] <- NULL
  
  data_Y <- do.call(rbind, data_Y) #make into matrices
  # TODO: make sure that data_Y has the same elements!!!!
  
  #choose those subjects who are both in the validation set and train set
  common_subjects = intersect(rownames(data_Y), rownames(validation_data)) 
  to_keep = which(rownames(data_Y)%in%common_subjects)
  data_Y <- data_Y[to_keep,]
  
  #Center and scale Y & validation data
  tmp1 <- scale(validation_data[!is.na(validation_data[,1]),],center=T,scale=T)
  Y_val <- scale(validation_data, attr(tmp1,"scaled:center"), attr(tmp1,"scaled:scale"))
  if(use_all){
    reorder_idx <- match(rownames(Y_val), rownames(data_Y))
    if(length(spectra)==3){
      ordering <- c(reorder_idx, (length(common_subjects)+reorder_idx))
    } else if(length(spectra)==4){
      ordering <- c(reorder_idx, (length(common_subjects)+reorder_idx), (2*length(common_subjects)+reorder_idx))
    }
    data_Y <- data_Y[ordering,]
    
  } else {
    reorder_idx <- match(rownames(data_Y), rownames(Y_val))
    Y_val <- Y_val[reorder_idx,]
  }
  
  
  
  
  
  tmp <- scale(data_Y[!is.na(data_Y[,1]),],center=T,scale=T) #scale after re-ordering
  Y1 <- scale(data_Y, attr(tmp,"scaled:center"), attr(tmp,"scaled:scale"))
  
  #Construct identifier matrix (covariates)
  subj = dimnames(Y1)[[1]]
  
  if(group_by_spectra){
    classes <- setdiff(spectra, validation_set)
    M <- 2 #either N1 or N2 spectra
    S = length(subj) #number of observations
    
    x <- rep(NA,S); names(x) <- subj
    # TODO: automate for M
    group_idx <- seq(1, length(x), length(x)/M)
    x[group_idx[1]:group_idx[2]] <- 1
    x[(group_idx[2]+1):length(x)] <- 2
    
    X <- matrix(0,S,M,dimnames=list(subj,classes) )
    for (i in 1:length(x)) if (!is.na(x[i])) X[i, x[i]] <- 1
    
    print("Condensing two classes into one +-1 covariate.")
    X <- X[, 1, drop = FALSE]
    X[X == 0] <- -1
    
    
  } else {
    M = length(unique(subj)) #number of groups (=subjects) 
    S = length(subj) #number of observations
    
    x <- rep(NA,M); names(x) <- unique(subj)
    for(s in 1:M){
      x[subj[s]] <- s
    }
    x=rep(x, length(spectra)-1) 
    
    X <- matrix(0,S,M,dimnames=list(subj,unique(subj)))
    for(i in 1:length(x)) if(!is.na(x[i])) X[i,x[i]] <- 1
  }

  
  #Clean up the subject info dataframe
  ages = ages[ages$File %in% unique(subj)==TRUE, ] #to get rid of some extra subjects that should not be there
  
  
  #Trying with added Z
  Z1 = ages[ages$File %in% unique(subj),]
  reorder_idx <- match(rownames(Y1), Z1$File)
  Z1 <- Z1[reorder_idx,]

  Z = matrix(c(Z1$Age, Z1$Age, Z1$Age))
  
  data <- list(Y1, x, X, ages, Y_val, Z=Z)
  names(data) <- c("Y", "x", "X", "ages", "Y2", "Z")
  return(data)
}


#### VALIDATION ####

#' Function for validating the model
#' 
#' @param within_sample Logical. If FALSE, does out-of-sample validation (unseen data points for all individuals)
#' @param dis which distance measure to use. One of 'L1' or 'L2' or 'cos'
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
            dists <- apply(other_data, 1, function(x) dist_func(lat_map2[testidx,pK], x))
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



#### CROSS-VALIDATION ####

#' Function for performing cross-validation: calls 
#' 
#' @param data list containing the data  
#' @param n_folds how many folds to use
#' @param K BRRR rank 
#' @param iter how many iterations given to BRRR
#' @param dis which distance measure to use in validation
#' @param validation_scheme 'subject' if dropping subjects, 'unseen_data' if validating on new spectra


do_CV <- function(data, n_folds=5, K=20, iter=500, dis='L1', validation_scheme='subject') {
  
  #library(foreach)    # install.packages('foreach') options(repos = c(CRAN = "http://cran.rstudio.com"))
  #library(doParallel) # install.packages('doParallel')
  #registerDoParallel(makeCluster(4)) # Use 4 cores for parallel CV
  Ds <- vector(mode = "list", length = n_folds) #list containing results from each CV fold
  names(Ds) <- c(paste0('fold_', seq(1:n_folds)) )
  
  
  subjects <- unique(names(data$x))
  cvId <- c()
  S <- ncol(data$X) #number of UNIQUE subjects - leaving subjects out, not data points 
  
  
  while (length(cvId) < S) cvId <- c(cvId, sample(1:n_folds))
  # if the number of subjects/observations is not a multiple of folds, truncate the labels
  cvId <- cvId[1:S]
  #cvId <- c(cvId, cvId)
  
  for(fold in 1:n_folds){
    
    Y <- data$Y #extract these here to avoid racing problems
    X <- data$X
    x <- data$x
    Y2 <- data$Y2
    
    testsubj <- subjects[cvId==fold] #extract testsubj
    train_Y <- Y[which(!rownames(Y)%in%testsubj),]
    train_X <- X[which(!rownames(X)%in%testsubj) , which(!colnames(X)%in%testsubj) ]
    train_fam <- x[which(!names(x)%in%testsubj)]
    
    #fit the model for training data
    source("brrr.R")
    res <- brrr(X=train_X,Y=train_Y,K=K,Z=NA,n.iter=iter,thin=5,init="LDA",fam=train_fam)  
    
    #res$scaling2 <- ginv(averagePsi(res)%*%averageGamma(res)) # i have seen this as well
    res$scaling <- MASS::ginv(averageGamma(res))
    W <- res$scaling
    
    #lat_map <- train_Y%*%W
    
    
    #call validation function
    if(validation_scheme=='subject'){
      test_Y <- Y[which(rownames(Y)%in%testsubj),]
      test_X <- X[which(rownames(X)%in%testsubj) , which(colnames(X)%in%testsubj) ]
      lat_map <- test_Y%*%W #mapping to latent space with unseen individuals
      
      #is not actually within-sample but oh well--- works now
      D <- validation(within_sample = T, dis=dis, pK=c(1:K), lat_map=lat_map, 
                      Xt=test_X)
      # compute 0-model (correlations in full data matrix) 
      D0 <- validation(within_sample = T, dis='corr',  pK=c(1:247), lat_map=test_Y,   
                       Xt=test_X)
      
    } else if(validation_scheme=='unseen_data'){ #validate on the test spectra
      test_Y <- Y2[which(rownames(Y2)%in%testsubj),] #used to be ! before rownames
      test2_Y <- Y[which(rownames(Y)%in%testsubj),]
      #test_X <- train_X[1:nrow(test_Y),]
      test_X <- X[which(rownames(X)%in%testsubj) , which(colnames(X)%in%testsubj) ]
      lat_map2 <- test_Y%*%W
      lat_map <- test2_Y%*%W #used to be train_Y
      
      #is not actually within-sample but oh well--- works now
      D <- validation(within_sample = F, dis=dis, pK=c(1:K), lat_map=lat_map, lat_map2=lat_map2,  
                      Xt=test_X)
      
      # compute 0-model (correlations in full data matrix) 
      D0 <- validation(within_sample = F, dis='corr', pK=c(1:247), lat_map=test_Y, lat_map2=test2_Y,  
                      Xt=test_X)
      #calculate test MSE & PTVE
      # PRED <- test_X%*%averagePsi(res)%*%averageGamma(res)
      # test_MSE <- mean((test_Y-PRED)^2)
      # test_PTVE <- 1 - sum(apply((test_Y-PRED),2,var)) / sum(apply(test_Y, 2, var))
      # res$testMSE <- test_MSE
      # res$testPTVE <- test_PTVE
      # 
      # print(paste0("Test MSE: \n", test_MSE))
      # print(paste0("Test PTVE: \n", test_PTVE))
      # 
    }
    
    Ds[[fold]] <- list(D, res$model$ptve, res, D0) #append results
    
  }
   
  return(Ds) 
}
    
    

  


### TRAINING / RUNS ###


Ns <- seq(480, 780, by=30)
accuracies <- c()
ptve <-c()
null_accus <- c()
rank <- c()

#create set of conditions to loop over as a list
conds = list(c("N1A","N1B","N2C"),
             c("N1A","N1B","N2A","N2C"),
             c("N1A","N2B","N2C"),
             #c("N1A","N2B","N1B"),
             c("N1A","N2B","N2A","N1B"),
             c("N2A","N2B","N2C"),
             #c("N2A","N2B","N1B"),
             c("N2A","N2B","N2C","N2D"))

for(n in 1:length(conds)){
  
  spectra_list = unlist(conds[n])
  # read in the data
  n2_data <- prepare_data(spectra = spectra_list, validation_set = tail(spectra_list,1), 
                          n_inds=792, data_type='spectrum')
  Y = n2_data[[1]]
  X = n2_data[[3]]
  x = n2_data[[2]]
  ages = n2_data[[4]]
  # 
  Y2 = n2_data[[5]] #validation set data
  Z = n2_data[[6]]
  # 
  
  # The model is trained using two sets of N2 data, and the within-sample performance is evaluated using
  # MSE, PTVE and accuracy (L1 distances in the projection)
  source("brrr.R")
  
  CV_results = do_CV(data=n2_data, n_folds=10, K=20, iter=1000, validation_scheme='subject')
  
  CV_scores <- lapply(CV_results, `[[`, 1) #unlisting stuff; looks ugly
  #CV_ranks <- lapply(CV_scores, `[[`, 3)
  CV_scores <- lapply(CV_scores, `[[`, 2)
  
  CV_ptves <- lapply(CV_results, `[[`, 2) #get ptves
  #CV_2 <- lapply(CV_ptves, `[[`, 1)
  #CV_ptve <- lapply(CV_2, `[[`, 3)
  
  CV_null_scores <- lapply(CV_results, `[[`, 4)
  CV_null_scores <- lapply(CV_null_scores, `[[`, 2)
  
  
  accs <- unlist(lapply(CV_scores, `[[`, 1))
  print("Average CV accuracy:")
  print(mean(accs))
  
  ptvs = unlist(lapply(CV_ptves, `[[`, 1))
  print("Average train-PTVE:")
  print( mean(ptvs) )
  
  null_accs <- unlist(lapply(CV_null_scores, `[[`, 1))
  print("Average CV 0-accuracy:")
  print(mean(null_accs))
  
  # ranks = unlist(lapply(CV_ranks, `[[`, 1))
  # print("Average self-identification ranks:")
  # print( mean(ranks) )
  
  # rank = c(rank,mean(ranks))
  save(CV_results, file=paste0('results/', 10, 'foldCV/NEW_K30_o7_',paste(spectra_list, collapse=''), '.RData'))
  #write.csv(x=c(n, mean(accs),mean(ptvs),mean(ranks)), file=paste0("result_N1_all.csv"))
  
  accuracies = c(accuracies,mean(accs))
  ptve = c(ptve,mean(ptvs))
  null_accus = c(null_accus, mean(null_accs))
}



# loop thru results and do unseen_data validation
accs = c()
for(f in 1:length(CV_results)){
  res <- CV_results[[f]][[3]]
  res$scaling <- ginv(averageGamma(res))
  Xt <- res$data$genotypes
  lat_map <- res$data$phenotypes
  lat_map2 <- Y2[which(rownames(Y2)%in%colnames(Xt)),]
  D <- validation(within_sample = T, dis='L1', pK=c(1:K), 
                  lat_map=lat_map, lat_map2=lat_map2, Xt=Xt)
  
  accs = c(D[[2]], accs)
}
  
CV_results = do_CV(n_folds=10, K=12, iter=1000, validation_scheme='unseen_data')  
save(CV_results, file=paste0('results/', 10, 'foldCV/K12all_N1.RData'))

CV_scores <- lapply(CV_results, `[[`, 1)

n=lapply(CV_scores, sapply, mean)
accs <- unlist(lapply(n, `[[`, 2))
print("Average CV accuracy:")
print(mean(accs))



## Training with all data
source("brrr.R")
K=30
res <- brrr(X=X,Y=Y,K=K,Z=NA,n.iter=1000,thin=5,init="LDA",fam=x, omg=0.01) 
res$scaling2 <- ginv(averagePsi(res)%*%averageGamma(res)) # i have seen this as well
res$scaling <- ginv(averageGamma(res))

# ptve = res$factor_variance/sum(res$factor_variance)
# Ks <- c(1:K)
# plot(Ks,ptve, 'l', col='firebrick', ylab="ptve %", bty="n")

save(res, file = paste0("results/full/o7_N1N2_BRRR_",K, ".RData") )
W <- res$scaling
lat_map <- Y%*%W
lat_map2 <- Y2%*%W #mapping to latent space with unseen N2_D data! (HOLDOUT METHOD)

D1 <- validation(within_sample = F, dis='L1', pK=c(1:K), lat_map=lat_map, lat_map2=lat_map2, Xt=X)

#baseline: calculate distances on all data
nsubj <- length(unique(row.names(Y)))

#lat_map <- Y[1:(2*nsubj),] 
#lat_map2 <- Y[(2*nsubj+1):(3*nsubj),]
lat_map <- Y
lat_map2 <- Y2
D <- validation(within_sample = F, dis='L1', pK=c(1:247), lat_map=lat_map, lat_map2=lat_map2, Xt=X)


### Loop over different number of components
Ks <- seq(from=2,to=248, by=4)
Ns <- seq(100, 790, by=20)
ptves <- c()
accs_1 <- c()
accs_2 <- c()
accs_3 <- c()
N_s <- c() 

for(n in Ns){
  n2_data <- prepare_data(spectra = c("N2A", "N2B", "N2C"), validation_set = "N2C", 
                          n_inds=n, data_type='spectrum')
  Y = n2_data[[1]]
  X = n2_data[[3]]
  x = n2_data[[2]]
  ages = n2_data[[4]]
  # 
  Y2 = n2_data[[5]] #validation set data
  Z = n2_data[[6]]
  
  N_s = c(N_s, length(x))
  
  source("brrr.R") 
  K=30
  lkg='minimum'
  res <- brrr(X=X,Y=Y,K=K,Z=NA,n.iter=1000, burnin=0.5, thin=1, init="LDA",fam=x) #fit the model
  res$scaling <- ginv(averageGamma(res))
  save(res, file = paste0("results/full/all_2N2_BRRR_K",K,".RData") )

  res$scaling2 <- ginv(averagePsi(res)%*%averageGamma(res)) # i have seen this as well
  res$scaling <- ginv(averageGamma(res))
  W <- res$scaling

  lat_map <- Y%*%W
  lat_map2 <- Y2%*%W #mapping to latent space with unseen N2_D data! (HOLDOUT METHOD)
  
  D1 <- validation(within_sample = T, dis='L1', pK=c(1:K), lat_map=lat_map, 
                   lat_map2=NULL, Xt=X, linkage=lkg)
  D2 <- validation(within_sample = F, dis='L1', pK=c(1:K), lat_map=lat_map, 
                   lat_map2=lat_map2, Xt=X, linkage=lkg)
  D3 <- validation(within_sample = F, dis='corr', pK=c(1:247), lat_map=Y, 
                   lat_map2=Y2, Xt=X, linkage=lkg)
  
  ptves = c(ptves, res$model$ptve)
  accs_1 = c(accs_1, D1[[2]])
  accs_2 = c(accs_2, D2[[2]])
  accs_3 = c(accs_3, D3[[2]])
  #D <- validation(within_sample = T, dis='L1', pK=k, lat_map=lat_map, lat_map2=NULL, Xt=X)
  #if((n%%100)==0){
  #  save(res, file = paste0("results/full/all_2N2_BRRR_K",K,".RData") )
  #} 
}

N_s <- N_s/2
data_list <- list(N_s, ptves, accs_1, accs_2, accs_3)

save(data_list, file = "N-dependencies.RData")

#dataK <- load("K-dependencies.RData")
Ks <- data_list[[1]]
ptves <- data_list[[2]] 
accs_1 <- data_list[[3]] 
accs_2 <- data_list[[4]] 
accs_3 <- data_list[[5]]

pdf("figures/N2_all_latspace_dim_K.pdf", width=7, height=6)
plot(Ks, accs_1, xlab='K', 'l', col='yellow3', ylab="Score", bty="n", lwd=2, ylim=c(0,0.99), axes=FALSE)
lines(Ks, accs_3, xlab='K', 'l', col='darkorchid4', lty=2, lwd=2, axes=FALSE)
lines(Ks, ptves, xlab='K', 'l', col='deepskyblue4', lwd=2, axes=FALSE)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9))
axis(1, at=c(0,100,200,300,400,500,600,700,750))
legend(20,0.20, legend=c('BRRR', 'Null model', 'PTVE'), 
       col=c("yellow3", "darkorchid4", "deepskyblue4"), pch=c(16,15)) #topright?
dev.off()
# grid(nx = NULL, ny = NULL,
#      lty = 2,      # Grid line type
#      col = "gray", # Grid line color
#      lwd = 1)      # Grid line width

#### VISUALIZATIONS ####

Y <- res$data$phenotypes
X <- res$data$genotypes

x <- rep(NA, times=nrow(X))
names(x) <- rownames(X)
for(i in 1:length(x)) x[i] <- which(X[i,]==1)
subj <- names(x)


nsubj <- length(unique(subj))
# adding together N2 mappings  plus some covariates
lat_map = as.data.frame(rbind(lat_map2[1:nsubj,], lat_map))
lat_map['condition'] = c(rep('test', nsubj), rep('train', 1*nsubj))
lat_map['age'] = rep(ages[which(ages$File%in%subj),]$Age, 3) #get age data
lat_map['group'] = rep(round(ages[which(ages$File%in%subj),]$Age, 0), 3) #get age data
lat_map['sex'] = rep(ages[which(ages$File%in%subj),]$Sex, 3)

subj_number <- factor(c(x,seq(1,nsubj)))

ggplot(data=as.data.frame(lat_map), aes(lat_map[,1], lat_map[,2], shape=condition, col=factor(sex))) + 
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
#TODO: test for larger data set


library("Rtsne")
D_ages <- ages[ages$File %in% rownames(D),]
#match the ordering that is mixed due to random sampling
reorder_indexes <- match(rownames(D), D_ages$File)
D_ages <- D_ages[reorder_indexes,]


tsne <- Rtsne(D, dims=2, is_distance = T,
              perplexity=3, verbose=TRUE, max_iter = 500, check_duplicates = F)

#plot(tsne$Y,col=colors, asp=1)
#scatterplot3d(tsne$Y, pch=4, color=D_ages$Age, angle=60)



#trying distance matrix with GGplot
# TODO: something is wrong
D_ages$Y1 <- tsne$Y[,1]
D_ages$Y2 <- tsne$Y[,2]


ggplot(D_ages, aes(x=Y1, y=Y2, color=Sex)) +
  geom_point(size=4) + theme_minimal()




# Sammon mapping

#https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/sammon.html



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

#create another distance matrix from test data
#TODO: calculate per component and do the clustering analysis
library("dendextend")

for(k in 1:1){
  S <- D*0 #S is dissimilarity matrix
  colnames(S) <- rownames(S)
  for(testidx in 1:nrow(D)){ #calculates the distances between individuals   
    for(m in 1:M){     
      group_members <- rownames(lat_map_n2)[m]   
      idxs = which(row.names(lat_map_n2) %in% group_members) #     
      #group_data <- lat_map_n2[idxs,][1:k]
      group_data <- lat_map_n2[c(idxs), 1:k]
      S[testidx,m] <- sum(abs(lat_map_n2[testidx,1:k]-group_data)) #L1 distance#  
    } 
  }
  
  colnames(S) <- D_ages$Cap
  rownames(S) <- D_ages$Cap
  
  hc <- hclust(as.dist(S), method="ward.D")
  plot(hc, main=paste0("Clustering on component #", 1, "-", k))
  
}


heatmap(S)

dend <- as.dendrogram(hc)
dend %>% set("branches_k_color", k = 5) %>% plot(main = "Nice defaults")




### quick plots for results

#over 1y olds, N2A and N2B within-model
MSES <- c(0.2429,0.2165,0.2007,0.1893,0.1802,0.1729,0.1665,0.1615,0.1573,0.1535,0.1504,0.1475,0.1450,0.1428)
ptves <- c(0.7569,0.7833,0.7991,0.8106,0.8196,0.8270,0.8334,0.8384,0.8426,0.8463,0.8494,0.8254,0.8549,0.8571)
accs <- c(0.3926,0.5736,0.7530,0.8426,0.8426,0.8934,0.9154,0.9171,0.9323,0.9238,0.9205,0.9323,0.9425,0.9391)
K <- c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30)


#would be good on all data too!!!
#### so use dat 

plot(Ns, accs, type = 'l', col='red')
lines(Ns,ptves)


########## PENLDA for sleep stage differences #################


######### LDA ###########

#using PenalizedLDA.cv instead of own LOO-function to tune model params

#sample randomly
rand <- sample(x, 300)
class=match(rand, unique(rand)) #this is for full data

xte=Y[names(rand),]
#cv_results <- PenalizedLDA.cv(Y, class, nfold=6) 
res2 <- PenalizedLDA(Y[names(rand),], class, xte=xte, type="ordered", lambda=0, lambda2=0,K=1)
#png("figures/K6full_penLDA_confmat_77.png")
cmat = confusion_matrix(res2$y, res2$ypred) #results with 6 discriminant vectors used
plot_confusion_matrix(cmat$`Confusion Matrix`[[1]])
print("Accurcy:")
print(cmat$`Balanced Accuracy`)

dev.off()



