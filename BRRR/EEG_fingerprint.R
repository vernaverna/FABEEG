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
prepare_data <- function(spectra, validation_set, n_inds=180, start_age=0,
                         group_by_spectra = F, data_type='spectrum', seed=191){
  
  
  # Choosing the channels and frequencies
  chs <- 1:19
  omitFreq = c(0,60) #so takes in all
  
  #extract the age groups
  ages = read.csv('data/new_age_df.csv')
  ages <- ages[,-1]
  ages[ages==" "] <- NA #replace empty strings with NA-values
  
  use_all=T #should we use all subjects in training?
  age_gap=c(start_age,19) #ages to include
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

  # Choose common rownames
  rowname_list <- lapply(data_Y, rownames)
  common_rows <- Reduce(intersect, rowname_list)
  data_Y <- lapply(data_Y, function(mat) mat[common_rows, , drop = FALSE])
  
  # Extract validation data
  val=which(names(data_Y)==validation_set)
  validation_data = data_Y[val] 
  data_Y[val] <- NULL
  
  data_Y <- do.call(rbind, data_Y) #make into single matrix
  if(length(val) > 1){
    validation_data <- do.call(rbind, validation_data)
  } else {
    validation_data <- validation_data[1] #unlist the matrix
  } 

  #Center and scale Y & validation data
  tmp1 <- scale(validation_data[!is.na(validation_data[,1]),],center=T,scale=T)
  Y_val <- scale(validation_data, attr(tmp1,"scaled:center"), attr(tmp1,"scaled:scale"))
  
  
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
    x=rep(x, length(spectra)-length(val)) 
    
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
#' @param validation_scheme 'subject' if dropping subjects, 'unseen_data' if validating on additional spectra


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
    
    
################################################################################

### TRAINING / RUNS ###


#Ns <- seq(480, 780, by=30)
accuracies <- c()
ptve <-c()
null_accus <- c()
rank <- c()

#create set of conditions to loop over as a list
conds = list(c("N1A","N1B","N2C"),               # 1, 2    | 5
             c("N1A","N1B","N2A","N2C"),         # 1, 2, 3 | 5
             c("N1A","N2B","N2C"),               # 1, 4    | 5
             c("N1A","N2B","N1B"),               # 1, 4    | 2
             c("N1A","N2B","N2A","N1B"),         # 1, 3, 4 | 2
             c("N2A","N2B","N2C"),               # 3, 4    | 5
             c("N2A","N2B","N1B"),               # 3, 4    | 2
             c("N2A","N2B","N2C","N2D"),         # 3, 4, 5 | 6
             c("N1A","N2B","N2A","N2D"),         # 1, 3, 4 | 6
             c("N2A","N2B","N2C","N1B"))         # 3, 4, 5 | 2

# try with these combinations as well
conds2 = list(c("N1A","N1B","N2A"),              # 1, 2    | 3
             c("N1A","N1B","N2A","N2B"),         # 1, 2, 3 | 4
             c("N1A","N2A","N2C"),               # 1, 3    | 5
             c("N1A","N2A","N1B"),               # 1, 3    | 2
             c("N1A","N2C","N2D","N1B"),         # 1, 5, 6 | 2
             c("N2A","N2B","N2D"),               # 3, 4    | 6
             c("N2A","N2D","N1A"),               # 3, 6    | 1
             c("N2A","N2B","N2D","N2C"),         # 3, 4, 6 | 5
             c("N1B","N2B","N2A","N2C"),         # 2, 3, 4 | 5
             c("N2B","N2C","N2D","N1A"))         # 4, 5, 6 | 2


for(n in 1:length(conds2)){
  
  spectra_list = unlist(conds2[n])
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
  
  CV_results = do_CV(data=n2_data, n_folds=10, K=30, iter=1000, validation_scheme='unseen_data')
  
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
  save(CV_results, file=paste0('results/', 10, 'foldCV/unseen_data/NEW_K30_o7_',paste(spectra_list, collapse=''), '.RData'))
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

################################################################################

## Training with all data - OOS 
accuracies <- c()
ptve <-c()
null_accus <- c()


#create set of conditions to loop over as a list
conds3 = list(c("N1A","N1B","N2A","N2B"),         # 1, 2    | 3, 4  aX near
             c("N2A","N2B","N1A","N1B"),         # 3, 4    | 1, 2  aX near, permuted 
             c("N1A","N1B","N2A","N2D"),         # 1, 2    | 3, 6  aX far
             c("N2A","N2D","N1A","N1B"),         # 3, 6    | 1, 2  aX far, permuted
             c("N2C","N2D","N1A","N1B"),         # 5, 6    | 1, 2  aX far, permuted
             c("N2B","N2C","N2D","N1A"),         # 4, 5    | 6, 2  mixed farish 
             c("N1A","N2B","N2A","N1B"),         # 1, 4    | 3, 2  mixed near
             c("N1A","N2D","N2A","N2B"),         # 1, 6    | 3, 4  mixed in between
             c("N2A","N2B","N2A","N2D"),         # 3, 4    | 1, 6  mixed in between, permuted
             c("N1A","N2A","N2C","N2D"),         # 1, 3    | 5, 6  mixed far
             c("N2C","N2D","N1A","N2A"),         # 5, 6    | 1, 3  mixed far, permuted             
             c("N2A","N2B","N2C","N2D"),         # 3, 4    | 5, 6  N2 near
             c("N2A","N2D","N2B","N2C"))         # 3, 6    | 4, 5  N2 between


source("brrr.R")
for(n in 1:length(conds3)){
  
  spectra_list = unlist(conds3[n])
  # read in the data
  n2_data <- prepare_data(spectra = spectra_list, validation_set = tail(spectra_list,2), 
                          n_inds=792, start_age=0, data_type='spectrum')
  
  Y = n2_data[[1]]
  X = n2_data[[3]]
  x = n2_data[[2]]
  ages = n2_data[[4]]
  # 
  Y2 = n2_data[[5]] #validation set data
  Z = n2_data[[6]]
  
  K=30
  res <- brrr(X=X,Y=Y,K=K,Z=NA,n.iter=1000,thin=5,init="LDA",fam=x, omg=0.01) 
  #res$scaling2 <- ginv(averagePsi(res)%*%averageGamma(res)) # i have seen this as well
  res$scaling <- ginv(averageGamma(res))
  res$lat_map <- Y%*%res$scaling
  res$lat_map2 <-Y2%*%res$scaling
  res$input_data <- n2_data #needs to be saved if validation is done at later point
  
  #TODO: needs a new validation function! do also the o7 tomorrow.


  save(res, file = paste0("results/full/all_",paste(spectra_list, collapse=''),"_BRRR_K",K, ".RData") )

}


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
D <- validation(within_sample = F, dis='corr', pK=c(1:247), lat_map=lat_map, lat_map2=lat_map2, Xt=X)

###############################################################################
#  FIGURE 8 (Appendix)

### Loop over different number of components
Ks <- seq(from=2,to=248, by=4)
Ns <- seq(100, 790, by=20)
ptves <- c()
accs_1 <- c()
accs_2 <- c()
accs_3 <- c()
N_s <- c() 

for(n in Ns){
  n2_data <- prepare_data(spectra = c("N1A","N2B","N2A","N2C"), validation_set = c("N2A","N2C"), 
                          n_inds=n, data_type='spectrum')
  Y = n2_data[[1]]
  X = n2_data[[3]]
  x = n2_data[[2]]
  ages = n2_data[[4]]
  # 
  Y2 = n2_data[[5]] #validation set data
  Z = n2_data[[6]]
  
  #N_s = c(N_s, length(x))
  
  source("brrr.R") 
  K=30
  lkg='minimum'
  res <- brrr(X=X,Y=Y,K=K,Z=NA,n.iter=1000, burnin=0.5, thin=1, init="LDA",fam=x) #fit the model
  res$scaling <- ginv(averageGamma(res))
  save(res, file = paste0("results/full/o7_2N2N1_BRRR_K",K,".RData") )

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

# Plotting - repeat for K and N separately!

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

