library(MASS)
library(lattice)
library(penalizedLDA)
library(abind)

# NOTE: in BRRR Y contains MEG and X classes

#' A function for running BRRR (or LDA) for a class set (siblings, control vs patient, test subjects etc)
#'Important parameters
#' 
#' @param network If true, the model is run using the full data. Else cross-validation is done, and accuracies returned
#' @param discTop The number of components (K) to use
#' @param model Which model to use? "brrr" or "penlda"
#' @param n.iter Number of brrr iterations
#' @param omitMag logical value indicating whether to omit magnetometers or not
#' @param omitFreg The lowest and highest frequency values (as a vector) to take into account
#' @param ex if there are several __spectrum files, this is the prefix
#' @param saving if TRUE, saves D into RData
#' @param CVfolds Number of cross-validation folds
#' 
acc <- function(network=F,discTop=15,model="brrr",lambda=0,l1d=200,omitN=2000, n.iter = 100,
                omitFreq = c(0,100), ex='N2', saving=T,CVfolds=3,seed=NA, network_vis=F) {
  
  if(network){if(model == 'penlda'){stop('Choose brrr as model if you want to plot network figures')}}
  
  # Sets the seed
  if(is.na(seed)) set.seed(1) else set.seed(seed)
  
  # Saving data into files
  fname <- paste0(paste(model, sep="_"),"_K",discTop)
  #fname2 <- paste0(paste(prop, model, sep="_"),"_DK",discTop)
  #fname3 <- paste0(paste(prop, model, sep="_"),"_xteK",discTop)
  if(network) {
    fname <- paste0("results/full/",fname,".RData")
    #fnameD <- paste0("results/full/",fname2,".RData")
    #fnameP <- paste0("results/full/",fname3,".RData")
  } else {
    fname <- paste0("results/",CVfolds,"foldCV/",fname,".RData")
    #fnameD <- paste0("results/",CVfolds,"foldCV/",fname2,".RData")
    #fnameP <- paste0("results/",CVfolds,"foldCV/",fname3,".RData")
  }
  print(fname)
  dir.create(dirname(fname),showWarnings=F,recursive=T)
  #dir.create(dirname(fnameD),showWarnings=F,recursive=T)
  #dir.create(dirname(fnameP),showWarnings=F,recursive=T)
  
  #extract the age groups
  ages = read.csv('data/new_age_df.csv')
  ages <- ages[,-1]
  
  # Helper function for preparing the data for analysis
  #' @param ex if there are several __spectrum files, this is the prefix
  #' @param type eithrt 'ind' for individual data or 'age' for age groups
  prepare_data <- function(ex, type="ind"){
    
    # Choosing the channels and frequencies
    chs <- 1:19
    omitFreq = c(0,60)
    
    # Reading data into workspace
    loadfile <- paste0("data/",ex,"spectrum.RData")
    
    use_all=F #should we use all subjects in training?
    age_gap=c(1,19) #exclude some of the younger children? 
    #TODO: change to range?
    
    # Check if the file exists
    if(file.exists(loadfile)) {
      load(loadfile)
      
      if(use_all==F){
        set.seed(11)
        individuals = sample(individuals, 160)
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
          if(age_data$Age >= age_gap[1] && age_data$Age <= age_gap[2]){
            tmp <- t(Y[[s]]) #transposed
            tmp <- log10(tmp) #TODO: is this necessary???
            if(any(is.na(tmp))){ #browser()
              corrupted = c(corrupted, s)
            } 
            A[[s]] <- tmp[frequencies,chs] 
          } else {
            corrupted = c(corrupted, s)
          }
        } else {
          corrupted = c(corrupted, s)
        } 
      }
      #TODO: even age groups!!!!!
      
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
  
  # select 3 data points 
  n2_data <- prepare_data(ex="N2B")
  n2b_data <- prepare_data(ex="N2A")
  n2c_data <- prepare_data(ex="N2D")
  
  # choose valudation data
  #validation_data <- prepare_data(ex="N2C")
  #Y4 <- validation_data[[1]]
  
  Y1 = n2_data[[1]]
  Y2 = n2b_data[[1]]
  Y3 = n2c_data[[1]]
  
  #some of the subjects do not have all spectras: drop them out
  if(length(rownames(Y3)) != length(rownames(Y2))){
    to_del= which(!rownames(Y2) %in% rownames(Y3))
    Y2 <- Y2[-to_del,]
    Y1 <- Y1[-to_del,]
    #Y4 <- Y4[-to_del,]
  }

  
  Y = rbind(Y1, Y2, Y3)
  
  x = c(n2c_data[[2]], n2c_data[[2]], n2c_data[[2]])
  X = rbind(n2c_data[[3]], n2c_data[[3]], n2c_data[[3]]) #row binding
  M = n2c_data[[4]]
  subj = dimnames(Y)[[1]]
  
  #Trying with added Z
  #TODO: put inside prepare-data!
  Z1 = ages[ages$File %in% rownames(Y1),]
  reorder_idx <- match(rownames(Y1), Z1$File)
  Z1 <- Z1[reorder_idx,]
  
  # Scaling and centering Z for BRRR
  Z = matrix(c(Z1$Age, Z1$Age, Z1$Age)) 
  Z <- scale(Z)
  
  # Scaling & centering the Y for BRRR 
  tmp <- scale(Y[!is.na(Y[,1]),],center=T,scale=T)
  Y <- scale(Y, attr(tmp,"scaled:center"), attr(tmp,"scaled:scale"))
  
  S <- nrow(X)
  # If network figures are created, only the model is taught with all data and no CV is done. 
  if(network) {
    CVfolds <- 1
    #D <- array(NA,dim=c(S,M,discTop),dimnames=list(subj,unique(subj),paste0("K",1:discTop))) #CHANGE!
  } else {
    D <- array(NA,dim=c(M,M,discTop),dimnames=list(subjects,paste0("class",1:M),paste0("K",1:discTop)))
  }
  
  # Creating CV folds from data 
  cvId <- c()
  while(length(cvId)<M) cvId <- c(cvId, sample(1:CVfolds)) # S-> M; cvID per subject, not data point?
  cvId <- cvId[1:M]
  
  obs <- subj #TODO: separate observations and subjects clearly!
  subjects <- rownames(Y1) #unique subjects
  classes <- list(subjects, colnames(X)) #subject-class mapping
  
  # Going through cross validation: different procedures for penLDA and brrr
  source("brrr.R")
  Dk <- array(NA,dim=c(M,M,discTop),dimnames=list(subjects,subjects,c()))
  for(fold in 1:CVfolds) {
    
    print("-------------")
    print(paste0("CV ",fold,"/",CVfolds))
    
    #Dk changed to record data points rather than subjects M -> S
    #Dk <- array(NA,dim=c(M,M,discTop),dimnames=list(subjects,subjects,c()))
    testsubj <- if(CVfolds==1) c() else subjects[cvId==fold] #subj->subjects
    
    if(model=="penlda") {
      class <- x[obs[!obs%in%testsubj]] #Class required to run from 1 to max
      class <- match(class,unique(class))
      xte <- Y[obs,]
      
      res <- PenalizedLDA(Y[obs[!obs%in%testsubj],],class,xte=xte,
                          lambda=lambda,K=discTop,standardized=TRUE) #changed standardized to TRUE
      res$scaling <- res$discrim
      
      
    } else if(model=="brrr") {
      #pred <- X*NA
      if(CVfolds == 1){
        res <- brrr(X=X,Y=Y, K=discTop,n.iter=n.iter,thin=5,init="LDA", fam = x)
      } else {
        res <- brrr(X=X[obs[!obs%in%testsubj],,drop=FALSE],Y=Y[obs[!obs%in%testsubj],],
                    K=discTop,n.iter=n.iter,thin=5,init="LDA", fam=x[obs[!obs%in%testsubj]])
      }
      res$scaling <- ginv(averageGamma(res))
      #PROJ <- res$model$brr$context$Gamma
    }
    
    W <- res$scaling
    if(model=="brrr") {
      xte <- Y[obs,]
      res$xteproj <- xte%*%W
      P <- res$xteproj
    }
    
    # TODO: CHECK THESE
    if(length(testsubj)>0) {
      for(k in 1:discTop) {
        for(s1 in testsubj) {
          for(s2 in subjects[!subjects%in%testsubj]) {
            Dk[s1,s2,k] <- mean(abs(res$xteproj[which(obs==s1),1:k]-res$xteproj[which(obs==s2),1:k]))
          }
          
          #for(m in 1:M) {
          #  D[s1,m,k] <- mean(Dk[s1,subjects[m],k], na.rm=T) #TODO: check this! overrides something
          #} 
        }
      } 
    }
      
    cat(".")
    
  }
  
  
    # Checking if network figures were applied and if so, creating them
  if(network_vis) {
    source("visNetwork.R")
    net <- list(omitMag=omitMag,Y=res$scaling,keepFeat=keepFeat,
                penLDA=4,omitN=omitN,l1d=l1d)
    net$lambda <- lambda
    
    net$freq <- freq
    net$fname <- fname
    
    save(net,file=fname)
    visNetwork(net,onlyPdf=TRUE)
    visNetwork(net,onlyPdf=TRUE,levelplot=TRUE)
    
    save(D, file=fnameD) #and the distances
    save(P, file=fnameP) #and the individual projection
    return(net)
    
  } else {
    if(saving) {
      save(W,file=fname) #saves the projection
      save(D, file=fnameD) #and the distances
      save(P, file=fnameP) #and the individual projection
    }
    return(list(D, W))
  }
}