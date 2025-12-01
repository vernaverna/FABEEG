# copyright by the authors
#' this function computes the prediction error for the low-rank BRRR model
#'
#' @param test.data the data used to generate the model
#' @param mcmc.output the output of the MCMC modelling step
#' @param burnin the proportion of the MCMC chain to consider as burnin phase
#'
#' @return mse and ptve between input data and predicted data
#' @export
#'

compute.prediction.error <- function(test.data, mcmc.output, burnin = 0.5) {
	# this function computes the prediction error for the low-rank BRRR model


	# remove burn-in
	n.iter <- length(mcmc.output$traces$Gamma)	
	tmp.output <- remove.burnin(mcmc.output, burnin = round(burnin * n.iter))
	tmp.output <- remove.burnin(mcmc.output, burnin = 1) #why the reformulation?

	
	
	# temporary matrix initialization
	preds <- test.data$phenotypes
	preds[, ] <- 0
	
	
	# compute predictions
	for (i in 1:length(tmp.output$traces$Psi)) { #for each MCMC trace?
		
		# prediction
		preds	<- preds + 1/n.iter * test.data$genotypes %*% tmp.output$traces$Psi[[i]] %*% tmp.output$traces$Gamma[[i]]
		
	}
	
	
		
	# subtract prediction mean, compute variance
	tmp <- (test.data$phenotypes - preds)
	test.data.ptve <- 1 - sum(apply(tmp,2,var)) / sum(apply(test.data$phenotypes, 2, var))
		
		
		
	# traditional MSE
	MSE <- mean((tmp)^2) 
	
	return(list(MSE=MSE, test.data.ptve = test.data.ptve))
	
}
