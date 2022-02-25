# copyright by the authors
#' Function for initializing the inputs for "sparse.fa.gibbs"
#'
#' Hyperparameter of the model. These will be collected
#'   in a single list.
#' @param a1.shape parameters for the regression priors
#' @param a1.rate parameters for the regression priors
#' @param a1.lower.bound parameters for the regression priors
#' @param a2.shape parameters for the regression priors
#' @param a2.rate parameters for the regression priors
#' @param a2.lower.bound parameters for the regression priors
#' @param a.sigma parameters of the residual noise parameter
#'                   distribution
#' @param b.sigma dito
#' @param local.shrinkage.nu element-wise shrinkage parameter for Lambda
#' @param factor.relevance.cutoff ?
#' @param alpha0 ?
#' @param alpha1 ?
#' @param rank the rank of the model
#' @param n.patients the number of observations
#' @param n.targets the number of target variables
#' @param SHRINKAGE define local shrinkage?
#' @param only.variances if the noise model has only diagonal terms set to true
#' @param prior.var.Eta the variance of the latent factors for the
#'                independent-noise model
#'
#' @return list of two:
#' 	context: contains initial values for all variables that
#'            will be updated by the "sparse.fa.gibbs" algorithm:
#'            "Lambda", "variances", "Eta", "local.shrinkage", '
#'            "deltas", "a1a2", "rank".
#'
#' 	prior: contains hyperparameters of the model:
#'          "local.shrinkage.nu", "a.sigma", "b.sigma",
#'          "a1.shape", "a1.rate", "a1.lower.bound",
#'          "a2.shape", "a2.rate", "a2.lower.bound",
#'          "factor.relevance.cutoff", "alpha0", "alpha1"
#' @export


initialize.fa.from.prior <- function(a1.shape=18, a1.rate=2, a1.lower.bound=2, a2.shape=18, a2.rate=2, a2.lower.bound=3, a.sigma=2.2, b.sigma=0.3, local.shrinkage.nu=3, factor.relevance.cutoff=0.001, alpha0=-1, alpha1=-0.005, rank=20, n.patients=200, n.pheno=10, SHRINKAGE = TRUE, only.variances = FALSE, prior.var.Eta = 1) {
	#
	# Function for initializing the inputs for "sparse.fa.gibbs". 
	#
	# Inputs:
	#	Hyperparameter of the model. These will be collected
	#   in a single list.
	#
	#
	# Outputs:
	#	context: contains initial values for all variables that
	#            will be updated by the "sparse.fa.gibbs" algorithm: 
	#            "Lambda", "variances", "Eta", "local.shrinkage", '
	#            "deltas", "a1a2", "rank".
	#
	#	prior: contains hyperparameters of the model:
	#          "local.shrinkage.nu", "a.sigma", "b.sigma",
	#          "a1.shape", "a1.rate", "a1.lower.bound",
	#          "a2.shape", "a2.rate", "a2.lower.bound",
	#          "factor.relevance.cutoff", "alpha0", "alpha1"
	
	if (rank<2) {
		stop('Rank must be >= 2')
	}
	
	
	prior <- list(local.shrinkage.nu=local.shrinkage.nu, a.sigma=a.sigma, b.sigma=b.sigma, a1.shape=a1.shape, a1.rate=a1.rate, a1.lower.bound=a1.lower.bound, a2.shape=a2.shape, a2.rate=a2.rate, a2.lower.bound=a2.lower.bound, alpha0=alpha0, alpha1=alpha1, factor.relevance.cutoff=factor.relevance.cutoff, prior.var.Eta=prior.var.Eta)

	context <- list()

	context$variances <- 1/rgamma(n=n.pheno, shape=a.sigma, rate=b.sigma)

	# when noise model only has diagonal terms, the remaining things are not needed
	if (!only.variances) {
		
		if (SHRINKAGE) {
			context$local.shrinkage <- matrix(rgamma(n=n.pheno*rank, shape=local.shrinkage.nu/2, rate=local.shrinkage.nu/2), nrow=n.pheno, ncol=rank)

		} else {
			context$local.shrinkage <- matrix(1, nrow=n.pheno, ncol=rank)			
		}
		
		
		
		context$rank <- rank
		
		context$Eta <- matrix(rnorm(n=n.patients*rank, mean=0, sd=sqrt(prior.var.Eta)), nrow=n.patients, ncol=rank)
		
		a1 <- -1
		while (a1 < a1.lower.bound) {
			a1 <- rgamma(n=1, shape=a1.shape, rate=a1.rate)
		}
		a2 <- -1
		while(a2 < a2.lower.bound) {
			a2 <- rgamma(n=1, shape=a2.shape, rate=a2.rate)
		}
		context$a1a2 <- c(a1,a2)
		
		context$deltas <- rep(0, rank) 
		context$deltas[1] <- rgamma(n=1, shape=a1, rate=1)
		context$deltas[-1] <- rgamma(n=rank-1, shape=a2, rate=1)
		
		taus <- cumprod(context$deltas)
		Lambda.sd <- t(t(1/context$local.shrinkage)*1/taus)
		context$Lambda <- matrix(rnorm(n=n.pheno*rank, mean=0, sd=Lambda.sd), nrow=n.pheno, ncol=rank)
		
		
	}
	
	return(list(context=context, prior=prior))
}
