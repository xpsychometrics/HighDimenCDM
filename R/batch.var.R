#'@title Batch variance calculation
#'
#' @description Internal function used by autoGDINA.stEM.
#'
#' @param plist A list of parameter estimates collected across batches.
#' @param n Number of batches
#'
#' @return A vector of batch variance estimates.
#'
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#'
#' @examples
#'\dontrun{
#' dat <- realdata_ECPE$dat
#' Q <- realdata_ECPE$Q
#' fit <- GDINA(dat = dat, Q = Q, model = "GDINA")
#' fit
#' CA(fit)
#'
#'
#'
#' }


batch.var <- function(plist,n){
  phi.bar <- sapply(plist,rowMeans) # of par x n
  phi.hat <- rowMeans(Reduce(cbind,plist))
  rowMeans((phi.bar-phi.hat)^2)/(n-1)
}
