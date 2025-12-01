# copyright by the authors

#' Initialization of the model
#'
#' @param n.targets the number of target variables (paper notation: K)
#' @param n.covariates the number of covariates (paper notation: P)
#' @param n.patients the number of observations (paper notation: N)
#' @param fa.rank rank of the independent-noise part
#' @param brr.rank rank of the regression/latent-noise part
#' @param alpha0 define the update schedule of the rank
#' @param alpha1 see alpha0
#' @param fa.relevance.cutoff threshold used in shutting down model components for fa part
#' @param brr.relevance.cutoff threshold used in shutting down model components for brr part
#' @param local.shrinkage.nu element-wise shrinkage parameter for Gamma and
#'                       Lambda parameters
#' @param n.confounders the number of observed confounders in the model, NA = no confounders
#' @param a3.shape parameters of the regression
#'                                      weight variance distribution
#'                    NOTE: a3 and a4 are a1  and a2 in paper notation and vice versa!
#' @param a3.rate dito
#' @param a3.lower.bound dito
#' @param a4.shape dito
#' @param a4.rate dito
#' @param a4.lower.bound dito
#' @param independent.noise include the terms for the independent structured
#'              noise (H,  Lambda)
#' @param a.sigma parameters of the residual noise parameter
#'                   distribution
#' @param b.sigma dito
#' @param prior.var.Eta the variance of the latent factors for the
#'                independent-noise model (H)
#'
#' @return the initialized model
#' @export

initialize.from.prior <- function(n.pheno=10, n.snps=15, n.patients=200, fa.rank=10, brr.rank=3, alpha0=-1, alpha1=-0.0005, fa.relevance.cutoff=0.01, brr.relevance.cutoff=0.01, local.shrinkage.nu=3, n.confounders=1, a3.shape=3, a3.rate=0.28, a3.lower.bound=2.1, a4.shape=4.1, a4.rate=0.31, a4.lower.bound=3.1, independent.noise = TRUE, a.sigma = 2.2, b.sigma=0.3, prior.var.Eta=1) {
  
	model <- list()

	
	if (independent.noise) {
		# Simulate the low-rank covariance:
			
		fa <- initialize.fa.from.prior(a1.shape=3, a1.rate=1, a1.lower.bound=2.1, a2.shape=3.1, a2.rate=1, a2.lower.bound=3.1, a.sigma=a.sigma, b.sigma=b.sigma, local.shrinkage.nu=local.shrinkage.nu, factor.relevance.cutoff=fa.relevance.cutoff, alpha0=alpha0, alpha1=alpha1, rank=fa.rank, n.patients=n.patients, n.pheno=n.pheno, SHRINKAGE = TRUE, prior.var.Eta=prior.var.Eta)	
		
		
	} else {
		
		# only simulate the variances
		fa <- initialize.fa.from.prior(a1.shape=18, a1.rate=2, a1.lower.bound=2, a2.shape=18, a2.rate=2, a2.lower.bound=3, a.sigma=2.2, b.sigma=0.3, local.shrinkage.nu=local.shrinkage.nu, factor.relevance.cutoff=fa.relevance.cutoff, alpha0=alpha0, alpha1=alpha1, rank=fa.rank, n.patients=n.patients, n.pheno=n.pheno, only.variances = TRUE)
		
		
	}
	

	
	model$fa$context <-fa$context
	model$fa$prior <- fa$prior	
	
	# Simulate the reduced-rank regression coefficients:
	inf.brr <- initialize.infinite.brr(local.shrinkage.nu=local.shrinkage.nu, a3.shape=a3.shape, a3.rate=a3.rate, a3.lower.bound=a3.lower.bound, a4.shape=a4.shape, a4.rate=a4.rate, a4.lower.bound=a4.lower.bound, brr.factor.relevance.cutoff=brr.relevance.cutoff, alpha0=alpha0, alpha1=alpha1, a.sigma=NA, b.sigma=NA, brr.rank=brr.rank, n.snps=n.snps, n.pheno=n.pheno)
	

	model$brr$context <- inf.brr$context
	model$brr$prior <- inf.brr$prior
	# NOTE: a.sigma is NA, therefore, variances are not simulated
	# (they are already in model$fa$context)


	# Simulate the regression coefficients for the confounders:
	# Just some random initialization.
	# Formally, A has improper prior.

	if (!is.na(n.confounders)) {
		model$A <- matrix(rnorm(n=n.confounders*n.pheno, mean=0, sd=1), nrow=n.confounders, ncol=n.pheno)  	
	}
	
	
	return(model)
}

