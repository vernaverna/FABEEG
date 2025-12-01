# copyright by the authors
library("viridis")
#' this function is used for comparing the true model parameters
#' to an estimate by plotting both.
#'
#' @param correct the parameters (Psi \%*\% Gamma) used to generate the data
#' @param est.mean the estimated mean
#' @param plot.path path to save the resulting plots, if null, opens an x11 connection
#' @param plot.title title for the plot
#'
#' @export
#'

plot.matrix.comparison <- function(correct, est.mean=NULL, plot.path = NULL, plot.title) {
  
  n.col <- ncol(correct)
  n.row <- nrow(correct)
  
  #cols=rich.colors(30)
  cols = viridis(n=30)
  if (is.null(plot.path)) {
    x11()
  } else {
    pdf(file = paste(plot.path, 'true_params', '.pdf', sep = ''))
  }
  image(x=seq(1,n.row), y=seq(1,n.col), z=correct, col=cols, axes=FALSE, xlab='', ylab='',main=paste('correct', plot.title))
  image.legend( n.row/5*4, n.col/10, zlim = range(correct), bg='white', yju=0, col=cols, at.z=range(correct), legnd=range(round(correct, digits=1)))
  if (!is.null(plot.path)) {
    dev.off()
  }
  
  if (!is.null(est.mean)) {
    
    if (is.null(plot.path)) {
      x11()
    } else {
      pdf(file = paste(plot.path, 'estimated_params', '.pdf', sep = ''))
    }
    image(x=seq(1,n.row), y=seq(1,n.col), z=est.mean, col=cols, axes=FALSE, xlab='', ylab='',main=paste('estimated', plot.title))
    image.legend( n.row/5*4, n.col/10, zlim = range(est.mean), bg='white', yju=0, col=cols, at.z=range(est.mean), legnd=range(round(est.mean, digits=1)))
    
    if (!is.null(plot.path)) {
      dev.off()
    }
  }
}
