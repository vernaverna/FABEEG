# copyright by the authors
#' A function for simulating independent SNP genotypes
#'
#' @param n.patients number of patients (observations)
#' @param n.snps number of genotypes (covariates)
#' @param mafs minor allele frequencies, if not given will be drawn
#' from a uniform random distribution
#'
#' @return unnormalized genotypes
#' @export
#'


simulate.genotypes <- function(n.patients, n.snps, mafs = NULL) {

  
  if (is.null(mafs)) {
    mafs <- runif(n=n.snps, min=0.05, 0.5)
  }
  
  prob.minor.homozygote <- mafs^2
  prob.heterozygote <- 2* mafs * (1-mafs)
  prob.major.homozygote <- (1-mafs)^2
  
  table <- matrix(runif(n=n.patients*n.snps, min=0, max=1),nrow=n.patients, ncol=n.snps)
  
  genotypes <- matrix(2, nrow=n.patients, ncol=n.snps)
  heterozygote <- t(t(table)>prob.minor.homozygote & t(table)<prob.minor.homozygote+prob.heterozygote)
  genotypes[heterozygote]=1
  major.homozygote <- t(t(table)>prob.minor.homozygote+prob.heterozygote)
  genotypes[major.homozygote]=0
  
  #return(list(genotypes=genotypes, mafs=mafs))
  return(genotypes)
}
