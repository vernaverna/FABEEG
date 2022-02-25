# copyright by the authors
#' random values from a multivariate normal distribution
#'
#' @param mu the mean of the distribution
#' @param Sigma Sigma
#'
#' @return the random values
#' @export
#'

mvr.norm.own <- function(mu, Sigma) {
	d <- nrow(Sigma)
	A <- try(chol(Sigma), silent=T)
	iter <- 0
	while(inherits(A, 'try-error')) {
		diag(Sigma) <- diag(Sigma) + 10^-(8-iter)
		A <- try(chol(Sigma), silent = TRUE)
		iter <- iter + 1		
	}
	if (iter > 0) {
		print('numeric stabilization: mvr.norm.own ')	
	}
	
	A <- t(A)
	z <- rnorm(n=d)
	val <- mu + A %*% z
	return(val)
}
