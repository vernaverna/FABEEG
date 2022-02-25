# copyright by the authors
#' Function for fitting the Bayesian multivariate
#' linear regression model with conjugate priors:
#'
#' @param X matrix with covariates (n*q)
#' @param y vector of responses (n*1)
#' @param noise.var noise variance (scalar)
#' @param prior.var prior variances (q*1)
#' @param crossprod.X optional parameter for passing in a precomputed crossproduct of the covariates
#'                    if NULL the crossproduct will be computed
#'
#' It is possible that y is a matrix (n*p),
#' in which case noise.var must be a vector
#' of length p specifying the variances of
#' the different target variables. prior.var must
#' in this case be q*p
#'
#' NOTE: It is assumed that the prior mean
#' is equal to zero, and the prior covariance
#' is diagonal.
#'
#'
#' @return (list with 2 elements):
#' 	posterior.cov
#' 	posterior.mean
#' @export
#'

fit.bayes.lm.diag <- function(X, y, noise.var, prior.var, crossprod.X=NULL) {

	
	n.geno <- ncol(X)

	if (is.null(crossprod.X)) {
		crossprod.X <- crossprod(X)
	}

	if (is.matrix(y)) {

		posterior.mean <- list()
		posterior.cov <- list()
	
		p <- ncol(y)
		
		crossprod.X.y <- crossprod(X,y)  # Calculate outside the loop for speed.

		
		
		for (i in 1:p) {
			posterior.cov[[i]] <- chol2inv(chol( 1/prior.var[,i]*diag(n.geno) + 1/noise.var[i]*crossprod.X ))

			posterior.mean[[i]] <- posterior.cov[[i]] %*% ( 1/noise.var[i] * crossprod.X.y[,i] )
		}

	} else {
		posterior.cov <- chol2inv(chol( 1/prior.var*diag(n.geno) + 1/noise.var*crossprod.X ))
		posterior.mean <- posterior.cov %*% ( 1/noise.var * t(X) %*% y )
	}

	return(list(posterior.mean=posterior.mean, posterior.cov=posterior.cov))
}
