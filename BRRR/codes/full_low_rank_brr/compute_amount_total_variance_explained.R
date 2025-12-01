# copyright by the authors

#' Function for computing the amount of total variation
#' explained by the reduced rank part of the model, given
#' fixed parameters Psi and Gamma.
#'
#' @param covariates matrix of centered covariates
#' @param Psi   low rank representation of the
#' 		coefficient matrix such that Psi \%*\% Gamma is equivalent to the regression coefficients
#' @param Gamma see Psi
#'
#' @return total variance explained
#' @export

compute.amount.total.variance.explained <- function(genotypes, Psi, Gamma) {

	aux <- genotypes %*% Psi %*% Gamma
	amount.total.var.explained <- sum(apply(aux,2,var))

}


