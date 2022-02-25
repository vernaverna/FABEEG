# copyright by the authors
#' this function removes the burnin from a trace list
#'
#' @param mcmc.output the original (full) mcmc.output from
#'                gibbs.full.low.rank.brr
#' @param burnin the number of burnin samples

#'
#' @return the mcmc.output without the first n samples where n = burnin
#' @export
#'


remove.burnin <- function(mcmc.output, burnin, n.to.use = NULL) {

  
	if (is.null(n.to.use)) {
		for (name in names(mcmc.output$traces)) {
			mcmc.output$traces[[name]] <- mcmc.output$traces[[name]][-seq(1:burnin)]
		}
	} 


	return(mcmc.output)
}

