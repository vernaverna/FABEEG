# copyright by the authors

#' creates a new row for Gamma
#'
#' @param star.deltas hierarchical priors for global shrinkage parameter tau
#' @param local.shrinkage.nu hierarchical priors (rate parameters) for local shrinkage parameter psi?
#' @param a4 rate parameter for gamma prior distribution used to build deltas 
#' @param n.targets number of targets
#' @param n.covariates number of covariates
#'
#' @return list containing simulated draws from prior(?) distributions
#' @export 
#'

simulate.new.brr.factor <- function(star.deltas, local.shrinkage.nu, a4, n.pheno, n.snps) {
	
  # Column to be added is never the first.
	new.star.delta <- rgamma(n=1, shape=a4, rate=1)  
	new.star.tau <- prod(star.deltas) * new.star.delta
	
	# Psi
	
	
	# Psi is not shrunk
	new.Psi.col <- rnorm(n=n.snps, mean=0, sd=1)
		
	
	
	
	#n.pheno.clusters <- length(unique(output.clustering))
	
	new.Gamma.local.shrinkage.parameters <- rgamma(n=n.pheno, shape=local.shrinkage.nu/2, rate=local.shrinkage.nu/2)
	
	new.Gamma.local.shrinkage.row <- new.Gamma.local.shrinkage.parameters#[output.clustering]
	

		# Gamma is doubly shrunk
		new.Gamma.row <- rnorm(n=n.pheno, mean=0, sd=1/sqrt(new.Gamma.local.shrinkage.row * (new.star.tau)^2))		
		
	
	
	return(list(new.star.delta=new.star.delta, new.star.tau=new.star.tau, new.Psi.col=new.Psi.col, new.Gamma.local.shrinkage.row=new.Gamma.local.shrinkage.row, new.Gamma.row=new.Gamma.row))
	
}
