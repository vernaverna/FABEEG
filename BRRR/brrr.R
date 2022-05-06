#
# copyright by the authors
#
#################################################################
#                                                               #
# This is a demo about latent-noise and the method presented in # 
# "Multiple Output Regression with Latent Noise", accepted for  #
# publication at the Journal of Machine Learning Research by    #
# Gillberg et al.                                               #
#                                                               #
#################################################################

library(coda)
library(MASS)
library(penalizedLDA)
source("var.R")

# the following function source_directory is by 
# Mehmet Gonen (mehmet.gonen@gmail.com)

#Helper functions to compute average matrices from MCMC traces

averageGamma <- function(res) {
  G <- res$traces$Gamma[[1]]*0
  ps <- length(res$traces$Gamma)
  for(i in 1:ps) {
    Gamma <- res$traces$Gamma[[i]]/sqrt(apply(res$traces$Gamma[[i]]^2,1,mean)) #Move the scale to Psi
    G <- G + Gamma/ps
  }
  return(G)
}

averagePsi <- function(res) {
  P <- res$traces$Psi[[1]] * 0
  ps <- length(res$traces$Psi)
  for (i in 1:ps) {
    Psi <- res$traces$Psi[[i]] / sqrt(apply(res$traces$Psi[[i]]^2, 1, mean)) # Move the scale to Psi
    P <- P + Psi / ps
  }
  return(P)
}


compute.factorwise.variance <- function(data, Psi, Gamma) {
  
  
  #how about eigenvalue decomposition ???
  n.patients <- nrow(data$genotypes)
  genotype.cov <-  data$crossprod.genotypes / (n.patients-1)
  aux.Gamma <- tcrossprod(Gamma)  # same as Gamma %*% t(Gamma)
  aux.Psi <- t(Psi) %*% genotype.cov %*% Psi
  
  total <- aux.Gamma * aux.Psi # covariance matrix of sorts?
  
  
  # Total explained variance is obtained by summing elements
  # in total. Note that off-diagonal elements are added twice,
  # as should.
  return (diag(total))
  
}

brrr <- function(X=NULL, Y=NULL, K=NULL, Z=NA, n.iter=500, burnin=0.5, thin=1, init="LDA", fam=NULL,
                 seed=1, snip=TRUE, pruned=TRUE, snpScale=1/10, omg=1e-2) {
  
  source_directory <- function(path) {
    files <- sort(dir(path, "\\.[rR]$", full.names = TRUE))
    lapply(files, source, chdir = TRUE)
  }
  ################################
  # load functions and libraries #
  ################################
  current.path <- getwd() 
  if (current.path == 'addMyPathHere') {
    stop('Add the path to the folder containing this file as current.path')
  }
  code.path <- paste0(current.path, '/codes')
  
  
  # load codes
  source_directory(paste0(code.path, '/simulation_study/'))
  source_directory(paste0(code.path, '/full_low_rank_brr/'))
  source_directory(paste0(code.path, '/infinite_brr/'))
  source_directory(paste0(code.path, '/sparse_fa/'))
  source_directory(paste0(code.path, '/common/'))
  
  X <- X[,colSums(abs(X))>0,drop=FALSE]
  data <- list(genotypes=X,phenotypes=Y,confounders=Z) #include confounders (=matrix Z)? 
  data$crossprod.genotypes <- crossprod(data$genotypes) #same as t(X) %*% X
  
  n.pheno <- ncol(data$phenotypes)
  n.train <- nrow(data$phenotypes)
  n.snps <- ncol(data$genotypes) #change to n.covariates?
  
  
  ############################
  # Learn latent-noise BRRR  #
  ############################
  
  
  # rank is sampled with NULL
  brr.rank <- K
  
  # latent signal-to-noise-ratio (latent SNR)
  # Omega.coef: 1/(latent SNR)
  # we are dealing with weak effects, this parameterization seems 
  # convenient
  # So here most of the data is explained by the covariates, not the noise?
  #Omega.coef <- 1e-03 #1e-6 #0.1 #7.5
  Omega.coef <- omg
  
  
  
  
  learnt.params <- list()

  # BRRR rank used for inititialization
  if (is.null(brr.rank)) brr.rank.init <- 3 else brr.rank.init <- brr.rank
  rank <- brr.rank.init
  
  if(!any(is.na(Z))){
    n.confounders=ncol(data$confounders)
  } else n.confounders=NA
  
  ######################
  #  initialize model  #
  ######################
  # n.pheno = K (target variables)
  # n.snps = P (covariates)
  # n.patients = N (observations)
  # brr.rank = rank of the regression/latent-noise part
  # n.confounders=ncol(data$confounders)
  init.model <- initialize.from.prior(n.pheno= n.pheno, n.snps=n.snps, n.patients=n.train, n.confounders=n.confounders,
                                      fa.rank=3, brr.rank=brr.rank.init,  a.sigma = 2.2, b.sigma = 0.5)
  Gamma <- init.model$brr$context$Gamma
  Psi <- init.model$brr$context$Psi
  if(init=="PCA") {
    pcaInit <- prcomp(Y)
    Gamma <- t(pcaInit$rotation[,1:rank])*apply(pcaInit$x[,1:rank],2,sd)
    Psi[,] <- rnorm(n.snps*brr.rank.init, sd=1)
    
  } else if(init=="LDA") {
    fam <- match(fam,sort(unique(fam)))
    res <- PenalizedLDA(Y,fam,lambda=0,K=rank,standardized=TRUE) ###standardized true ?
    Gamma <- t(res$discrim)
    Psi[,] <- 0
    if(exists("Xg") && ncol(Xg)>0) { #Genotype in covariates - initialize to fit the data
      #f <- lm(Y%*%ginv(Gamma) ~ Xg+0)
      #tmp <- f$coefficients[-1,]
      #Psi[1:ncol(Xg)+nrow(Psi)-ncol(Xg),] <- ginv(Xg)%*%Y%*%ginv(Gamma)
    } else {
      Xg <- matrix(NA,1,1); colnames(Xg) <- NA
    }
    # Rest of variance to be explained with the families
    if(prod(dim(Psi))>1) {
      for(k in 1:rank) {
        Yres <- Y - X %*% Psi %*% Gamma #residual
        for(f in 1:max(fam))
          Psi[f,k] <- mean(Yres[fam==f,]%*%Gamma[k,]/sum(Gamma[k,]^2))
      }
      #Generalized inverse probably off worse
      #Psi <- ginv(X)%*%Y%*%ginv(Gamma)
      Ksd <- apply(Psi,2,sd)
      Psi <- sweep(Psi, MARGIN=2, Ksd, "/")
      Gamma <- Gamma*Ksd
    } else {
      tmp <- optimize(function(z) mean((Y-z*X%*%Gamma)^2), interval=c(1e-4, 1e4))
      Psi[1,1] <- tmp$minimum
    }
    
  } else if(init=="random") {
    # random initialization with excess variance
    Psi[,] <- rnorm(n.snps*brr.rank.init, sd=1)
    Gamma[,] <- rnorm(brr.rank.init*n.pheno, sd=1)
  }
  Yhat <- X %*% Psi %*% Gamma
  init.model$ptve <- 1-mean((Y-Yhat)^2)/mean(Y^2)
  
  init.model$brr$context$Gamma <- Gamma
  init.model$brr$context$Psi <- Psi
  
  
  ## Set very small noise residuals
  init.model$fa$context$variances <- rep(0.01, n.pheno)
  
  # do not sparsify Gamma
  # -> set shrinkage parameter to non-shrinkage
  init.model$brr$context$Gamma.local.shrinkage[,] <- 1
  init.model$brr$context$star.deltas <- rep(1.1, length(init.model$brr$context$star.deltas)) #1
  
  # parameters related to sampling the rank: how long should it be updated
  init.model$brr$prior$alpha0 <- -2
  init.model$brr$prior$alpha1 <- (log(0.1) - init.model$brr$prior$alpha0) / (n.iter*burnin)
  
  
  # set BRRR shrinkage parameters. Note that in the notation of the
  # paper these are a1 and a2 (and vice versa)
  init.model$brr$context$a3a4 <- c(2, 2) #c(10, 4.1) #
  
  # Testing! More variance to Psi, to allow for different scales in X!
  init.model$brr$prior$psiPrec <- 1 #10000 #
  
  #########################
  # latent-noise variance #
  #########################
  
  # Set the variance of the latent noise to match the prior assumption
  # about the a priori latent signal-to-noise ratio
  
  # first compute prior variance of X\Psi, then set
  # prior variance of \Omega to match it when multiplied with
  # a factor of "Omega.coef"
  # for how to compute var(X %*% Psi), see e.g.
  # http://en.wikipedia.org/wiki/Variance#Basic_properties
  
  
  
  # each component has same variance
  # do following var calculation 1 comp at a time and sum
  prior.var.X.Psi <- 0
  prior.vars.Psi <- diag(ncol(data$genotypes))
  
  for (rank.tmp in 1:init.model$brr$context$brr.rank) {
    
    prior.var.X.Psi <- prior.var.X.Psi + sum(diag(data$genotypes %*% prior.vars.Psi %*% t(data$genotypes)))
  }
  
  
  # all elements of Omega assumed to have the same variance
  init.model$brr$context$latent.noise.var <- Omega.coef * prior.var.X.Psi / (n.train * brr.rank.init)
  init.model$brr$context$Omega <- matrix(rnorm(brr.rank.init*n.train,
                                               sd = sqrt(init.model$brr$context$latent.noise.var)),
                                         nrow=n.train, ncol=brr.rank.init)
  
  ptm <- proc.time()
  mcmc.output <- gibbs.full.low.rank.brr(model=init.model, data=data, n.iter=n.iter,
                                         thin=thin, fixed.brr.rank=brr.rank,
                                         fa.vars.to.record=c('variances','local.shrinkage',
                                                             'rank','a1a2','deltas'),
                                         brr.vars.to.fix=c('Gamma.local.shrinkage','a3a4')) #,'star.deltas')) #,
  #'Psi.local.shrinkage'))
  
  test.scores <- compute.prediction.error(data, mcmc.output, burnin = burnin)
  print(test.scores) #current performance metrics are MSE and PTVE
  
  # to use the independent-noise BRRR model
  #
  #  in the initialization, set:
  #    init.model$brr$context$Omega <- NULL
  #
  #  give the following arguments to gibbs.full.low.rank.brr:
  #    ind.struct.noise = TRUE
  #    latent.noise = FALSE
  
  
  
  print('model learnt')	
  
  
  # study MCMC chain: with large number of samples and data points, 
  # this will converge to the model used to generate the data
  if(burnin>0)
    mcmc.output <- remove.burnin(mcmc.output=mcmc.output, burnin=burnin*round(n.iter/thin))
  mcmc.output$init <- init.model
  mcmc.output$data <- data
  mcmc.output$runtime <- proc.time() - ptm
  
  #Use check_mcmc_result.R to study convergence 
  #context <- list()
  #context$Psi <- res$model$brr$context$Psi #what should be put into context?
  #context$Gamma <- res$model$brr$context$Gamma
  #name='coefMat'
  #plot.path = paste0(getwd(), "/figures/mcmc_")
  #plot.title = "coefMat=Psi*Gamma"
  #result = check.mcmc.result(context = context, mcmc.output = mcmc.output, name=name, plot.path = plot.path, plot.title = plot.title)
  
  #  From Gillberg. et. al. (2016):
  #
  #     Averaged effective sample sizes (ESS) and potential scale reduction
  #     factors (PSRF) were computed for 200 randomly selected parameters of the regression
  #     coefficient matrix (Gelman et al., 2004). 
  #
  
  
  #ADDED: PTVE per each component 
  #factor_variance <- compute.factorwise.variance(data=data, Psi=averagePsi(mcmc.output),
  #                                                Gamma=averageGamma(mcmc.output))
  #
  # Psi = t(ginv(averagePsi(mcmc.output)) or just averagePsi(mcmc.output)? Depends?
  #factor_variance <- compute.factorwise.variance(data=data, Psi=t(ginv(averagePsi(mcmc.output))),
  #                                               Gamma=t(ginv(averageGamma(mcmc.output))))
  #mcmc.output$factor_variance <- factor_variance #total variation explained =sum
  #print(factor_variance)

  
  return(mcmc.output)
  
}


