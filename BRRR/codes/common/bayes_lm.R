test.bayes.lm <- function(n.patients=500, p=200) {
	noise.var <- 0.1
	X <- matrix(rnorm(n=p*n.patients, mean=0, sd=1), nrow=n.patients, ncol=p)
	crossprod.X <- crossprod(X)

	v.true <- rnorm(n=p, mean=0, sd=2)

	patient.means <- X %*% v.true
	y <- rnorm(n=n.patients, mean=patient.means, sd = sqrt(noise.var))

	prior.cov <- diag(p)
	time.start1 <- proc.time()[3]
	#blm <- fit.bayes.lm(X, y, noise.var, prior.cov)
 	blm <- fit.bayes.lm.ref.prior(X, matrix(rep(y,2), ncol=2), rep(noise.var,2), crossprod.X=crossprod.X)
	time.finish1 <- proc.time()[3]
	print(paste('Time 1: ', time.finish1-time.start1, sep=''))
	
	time.start2 <- proc.time()[3]
	blm2 <- fit.bayes.lm.diag(X, matrix(rep(y,2), ncol=2), rep(noise.var,2), 1E10, crossprod.X=crossprod.X)
	time.finish2 <- proc.time()[3]
	print(paste('Time2: ', time.finish2-time.start2, sep=''))
	
	print(max(abs(blm$posterior.mean[[2]] - blm2$posterior.mean[[1]])))
	print(max(abs(blm$posterior.cov[[2]] - blm2$posterior.cov[[1]])))
	#print('Post mean: ')
	print(mean(abs(blm$posterior.mean[[1]])))
	#print('Post cov: ')
	print(mean(abs(blm$posterior.cov[[1]])))
	#print('Post sd: ')
	#print(sqrt(diag(blm$posterior.cov)))
}

fit.bayes.lm <- function(X, y, noise.var, prior.cov) {
	# Function for fitting the Bayesian multivariate
	# linear regression model with conjugate priors:
	# 
	# Inputs:
	#	X: matrix with covariates (n*p)
	#	y: vector of responses (n*1)
	#	noise.var: noise variance (scalar)
	#	prior.cov: prior covariance matrix (p*p)
	# 
	# NOTE: It is assumed that the prior mean
	# is equal to zero.
	#
	# Outputs (list with 2 elements): 
	#	posterior.cov
	#	posterior.mean
	
	
	n.patients <- length(y)
	
	posterior.cov <- chol2inv(chol( chol2inv(chol(prior.cov)) + 1/noise.var*crossprod(X) ))
	
	posterior.mean <- posterior.cov %*% ( 1/noise.var * crossprod(X,y) )
	
	return(list(posterior.mean=posterior.mean, posterior.cov=posterior.cov))
}


fit.bayes.lm.multivariate <- function(X, Y, noise.var, prior.prec, crossprod.X= NULL) {
	# Function for fitting the Bayesian multivariate
	# linear regression model with conjugate priors:
	# 
	# Inputs:
	#	X: matrix with covariates (n*p)
	#	y: matrix of responses (n*q)
	#	noise.var: 	target variable -specific noise variances (vector)
	#	prior.cov: prior covariance matrix (p*p)
	# 
	# NOTE: It is assumed that the prior mean
	# is equal to zero.
	#
	# Outputs (list with 2 elements): 
	#	posterior.cov
	#	posterior.mean
	
	latent.dim <- dim(X)[2]
	n.patients <- dim(Y)[1]
	XTy.full <- as.vector(apply(Y,2,function(yp){crossprod(X,yp)}))
	if (is.null(crossprod.X)) {
		XTX <- crossprod(X)
	} else {
		XTX <- crossprod.X
	}

	posterior.cov <- chol2inv(chol( prior.prec + diag(1/noise.var)%x%XTX ) )
	
 	posterior.mean <- posterior.cov %*% ( rep(1/noise.var, each=latent.dim) * XTy.full )
	
	return(list(posterior.mean=posterior.mean, posterior.cov=posterior.cov))
}

fit.bayes.lm.ref.prior <- function(X, y, noise.var, crossprod.X=NULL) {
	# Function for fitting the Bayesian multivariate
	# linear regression model with conjugate priors:
	# 
	# Inputs:
	#	X: matrix with covariates (n*q)
	#	y: vector of responses (n*1)
	#	noise.var: noise variance (scalar)
	# 
	# It is possible that y is a matrix (n*p),
	# in which case noise.var must be a vector
	# of length p specifying the variances of
	# the different phenotypes.
	#
	# NOTE: the noninformative reference prior 
	# corresponding to infinite variance for the
	# regression coefficients is assumed.
	#
	# Outputs (list with 2 elements): 
	#	posterior.cov
	#	posterior.mean
	
	if (is.null(crossprod.X)) {
		crossprod.X <- crossprod(X)
	}

	n.patients <- length(y)

	if (is.matrix(y)) {

		posterior.mean <- list()
		posterior.cov <- list()
		p <- ncol(y)
		crossprod.X.y <- crossprod(X,y)
		for (i in 1:p) {
			posterior.cov[[i]] <- chol2inv( chol( 1/noise.var[i] * crossprod.X ) )
			posterior.mean[[i]] <- posterior.cov[[i]] %*% ( 1/noise.var[i] * crossprod.X.y[,i])
		}
	} else {
		posterior.cov <- chol2inv( chol( 1/noise.var * crossprod.X ) )
		posterior.mean <- posterior.cov %*% ( 1/noise.var * crossprod(X, y) )
	}

	return(list(posterior.mean=posterior.mean, posterior.cov=posterior.cov))
}


fit.bayes.lm.diag <- function(X, y, noise.var, prior.var, crossprod.X=NULL) {
	# Function for fitting the Bayesian multivariate
	# linear regression model with conjugate priors:
	# 
	# Inputs:
	#	X: matrix with covariates (n*q)
	#	y: vector of responses (n*1)
	#	noise.var: noise variance (scalar)
	#	prior.var: prior variances (q*1)
	# 
	# It is possible that y is a matrix (n*p),
	# in which case noise.var must be a vector
	# of length p specifying the variances of
	# the different phenotypes. prior.var must
	# in this case be q*p
	#
	# NOTE: It is assumed that the prior mean
	# is equal to zero, and the prior covariance
	# is diagonal.
	#
	# Outputs (list with 2 elements): 
	#	posterior.cov
	#	posterior.mean
	
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

fit.bayes.lm.diag.unequal.noise <- function(X, y, noise.vars, prior.var) {
	# Function for fitting the Bayesian multivariate
	# linear regression model with conjugate priors:
	# 
	# Inputs:
	#	X: matrix with covariates (n*p)
	#	y: vector of responses (n*1)
	#	noise.vars: noise variances (n*1)
	#	prior.var: prior variances (p*1)
	# 
	# NOTE: It is assumed that the prior mean
	# is equal to zero, and the prior covariance
	# is diagonal.
	#
	# Outputs (list with 2 elements): 
	#	posterior.cov
	#	posterior.mean
	
	n.patients <- length(y)

	posterior.cov <- chol2inv(chol( diag(1/prior.var) + crossprod(X, 1/noise.vars*X) ) )
	posterior.mean <- posterior.cov %*% ( crossprod(1/noise.vars*X, y) )
	
	return(list(posterior.mean=posterior.mean, posterior.cov=posterior.cov))
}


fit.bayes.lm.unequal.noise.ref.prior <- function(X, y, noise.vars) {
	# Function for fitting the Bayesian multivariate
	# linear regression model with conjugate priors:
	# 
	# Inputs:
	#	X: matrix with covariates (n*p)
	#	y: vector of responses (n*1)
	#	noise.vars: noise variances (n*1)
	# 
	# NOTE: the noninformative reference prior 
	# corresponding to infinite variance for the
	# regression coefficients is assumed.
	#
	# Outputs (list with 2 elements): 
	#	posterior.cov
	#	posterior.mean
	
	n.patients <- length(y)	

	posterior.cov <- chol2inv( chol( crossprod(X, 1/noise.vars*X) ) )
	posterior.mean <- posterior.cov %*% ( crossprod(1/noise.vars*X, y) )
	
	return(list(posterior.mean=posterior.mean, posterior.cov=posterior.cov))
}