# copyright by the authors

#' a function to check an MCMC result, by computing a comparison and plotting it to file (pdf)
#'
#' @param context  contains the true variable values
#' @param mcmc.output contains the "traces" variable
#' @param name the name of the variable/ currently the only option is "coefMat"
#' @param plot.path the path to store the plot files at, if NULL, plots will not be saved, but output to an x11 window
#' @param plot.title the plots title
#'
#' @return a list containing mse                              the mean squared error between the correct and estimated mean
#'                           prop.in.interval                 the proportion of values used to generate the model that is withing 2 standard deviations of the predicted mean
#'                           correct, est.mean, est.std       the input data
#'                           correlation                      the correlation between the values use to generate the model and the models prediction
#' @export


check.mcmc.result <- function(context, mcmc.output, name, plot.path = NULL, plot.title = NULL) {

  
  if (name=='coefMat') {
    
    # Compare the coefficient matrices (obtainable from Psi and Gamma)
    
    res <- compute.coef.matrix.mcmc.estimate(mcmc.output)
    coef.matrix.mean <- res$coef.matrix.mean
    coef.matrix.std <- res$coef.matrix.std
    
    correct <- context$Psi %*% context$Gamma
    comp <- compare.matrices(correct=correct, est.mean=coef.matrix.mean, est.std=coef.matrix.std)
    
    
    plot.matrix.comparison(correct=correct, est.mean=coef.matrix.mean, plot.path = plot.path, plot.title)
    
    
  } 
  return(comp)
}
