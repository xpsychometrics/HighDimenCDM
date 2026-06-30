#'@title Batch variance calculation
#'
#' @description Computes variance of parameter estimates across batches in the stochastic EM algorithm.
#'
#' @param plist A list of parameter estimates collected across batches.
#' @param n Number of batches
#'
#' @return A vector of batch variance estimates.
#'
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#'


batch.var <- function(plist,n){
  phi.bar <- sapply(plist,rowMeans) # of par x n
  phi.hat <- rowMeans(Reduce(cbind,plist))
  rowMeans((phi.bar-phi.hat)^2)/(n-1)
}
