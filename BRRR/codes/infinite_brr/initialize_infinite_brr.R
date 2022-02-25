# copyright by the authors

#' A function for initializing the sparse infinite
#' Bayesian reduced rank regression.
#'
#' @param local.shrinkage.nu see initialize_fa_from_prior
#' @param a3.shape dito
#' @param a3.rate dito
#' @param a3.lower.bound dito
#' @param a4.shape dito
#' @param a4.rate dito
#' @param a4.lower.bound dito
#' @param brr.factor.relevance.cutoff  dito
#' @param alpha0 dito
#' @param alpha1 dito
#' @param a.sigma If a.sigma and b.sigma are given (non-NA), variances
# 	will be simulated using these values. Otherwise,
# 	the variances will be set to NA.
#' @param b.sigma dito
#' @param brr.rank the rank assumed
#' @param n.covariates number of covariates
#' @param n.targets number of targets
#' @param step.size ?
#'
#' @return context: a list with fields Psi, Gamma,
#'       Psi.local.shrinkage, Gamma.local.shrinkage,
#'       star.deltas, a3a4, brr.rank
#'
#'   prior: a list fields local.shrinkage.nu,
#'       a3.shape, a3.rate, a3.lower.bound,
#'       a4.shape, a4.rate, a4.lower.bound,
#'       brr.factor.relevance.cutoff, alpha0
#'       alpha1, variances
#'   (Note: variances is actually not a hyperparameter;
#'   however, it is required when updating the low-rank
#'   regression coefficent matrix, but it is not updated
#'   itself.)
#' @export
#'

initialize.infinite.brr <- function(local.shrinkage.nu=3, a3.shape=18, a3.rate=2, a3.lower.bound=2, a4.shape=18, a4.rate=2, a4.lower.bound=3, brr.factor.relevance.cutoff=0.01, alpha0=-1, alpha1=-5E-4, a.sigma=1, b.sigma=1, brr.rank=3, n.snps=50, n.pheno=10, step.size=10) {

	## Hyperparameters:    
	if (!is.na(a.sigma)) {
		precisions <- rgamma(n=n.pheno, shape=a.sigma, rate=b.sigma)
		variances <- 1/precisions
	} else {
		variances <- NA
	}
	
	prior <- list()
	prior$local.shrinkage.nu <- local.shrinkage.nu
	prior$a3.shape <- a3.shape
	prior$a3.rate <- a3.rate
	prior$a3.lower.bound <- a3.lower.bound
	prior$a4.shape <- a4.shape
	prior$a4.rate <- a4.rate
	prior$a4.lower.bound <- a4.lower.bound
	prior$brr.factor.relevance.cutoff <- brr.factor.relevance.cutoff
	prior$alpha0 <- alpha0
	prior$alpha1 <- alpha1
	prior$variances <- variances
	prior$step.size <- step.size
	
	

	
	## Variables to update:
	a3 <- -1
	while (a3 < a3.lower.bound) {
		a3 <- rgamma(n=1, shape=a3.shape, rate=a3.rate)
	}
	a4 <- -1
	while (a4 < a4.lower.bound) {
		a4 <- rgamma(n=1, shape=a4.shape, rate=a4.rate)
	}
	

	star.deltas <- rep(0, brr.rank)
	star.deltas[1] <- rgamma(n=1, shape=a3, rate=1)
	if (brr.rank>1) {
		star.deltas[2:brr.rank] <- rgamma(n=brr.rank-1, shape=a4, rate=1)
	}
	star.taus <- cumprod(star.deltas)
	 
	
	Psi <- matrix(rnorm(n=n.snps*brr.rank, mean=0, sd=1), nrow=n.snps, ncol=brr.rank)
	
	
	
	Gamma.local.shrinkage <- matrix(rgamma(n=brr.rank*n.pheno, shape=local.shrinkage.nu/2, local.shrinkage.nu/2), nrow=brr.rank, ncol=n.pheno)

	
	Gamma.precisions <- Gamma.local.shrinkage * (star.taus^2)
	
	Gamma.sd <- 1/sqrt(Gamma.precisions)
	Gamma <- matrix(rnorm(n=brr.rank*n.pheno, mean=0, sd=Gamma.sd), nrow=brr.rank, ncol=n.pheno)
		
	context <- list()
	context$Psi <- Psi
	context$Gamma <- Gamma
	context$Gamma.local.shrinkage <- Gamma.local.shrinkage
	context$star.deltas <- star.deltas
	context$a3a4 <- c(a3,a4)
	context$brr.rank <- brr.rank
	
	return(list(context=context, prior=prior))
}

