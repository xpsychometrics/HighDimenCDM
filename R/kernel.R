#'@title Gibbs Sampling Kernel for stEM Algorithm
#'
#' @description Performs Gibbs sampling updates of attribute profiles and item parameters within the stochastic EM algorithm.
#'
#' @param dat A required \eqn{N \times J} data matrix of N examinees to J items. Missing values are currently not allowed.
#' @param Q A required \eqn{J \times K} item-attribute association matrix, where J is the number of items and K is the number of attributes.
#' @param J The number of items.
#' @param reduced.profiles A list of reduced attribute profiles for each item.
#' @param alpha The initial latent attribute profiles.
#' @param B The number of Gibbs sampling iterations.
#' @param ip.only If \code{TRUE}, only item parameter estimates are returned.
#'
#'
#' @return a list with elements
#' \describe{
#' \item{ip}{estimated item parameters from Gibbs sampling iterations}
#' \item{alpha}{latent attribute profiles from the final Gibbs sampling iteration}
#' \item{a}{saved latent attribute profiles across Gibbs sampling iterations}
#' \item{item.parm}{final item parameter estimates}
#' }
#'
#' @author Wenchao Ma, The University of Minnesota, \email{wma@umn.edu}
#'
#'


kernel <- function(dat,Q,J,reduced.profiles,alpha,B,ip.only=TRUE){

  item.parm <- a <- ip <- list()

  for(b in 1:B){
    # item parameter estimation
    for(j in 1:J){
      col.loc <- which(Q[,j]!=0)
      if(length(col.loc)==1){
        if(all(alpha[col.loc,]==0)){
          item.parm[[j]] <- c(max(mean(dat[j,],na.rm = TRUE),.3),min(mean(dat[j,],na.rm = TRUE)+.4,.8))
        }else if(all(alpha[col.loc,]==1)){
          item.parm[[j]] <- c(max(mean(dat[j,],na.rm = TRUE)-.4,.1),max(mean(dat[j,],na.rm = TRUE),.6))
        }else{
          z <- aggregate(dat[j,],list(alpha[col.loc,]),mean,na.rm = TRUE)$x
          z[z<1e-4] <- 1e-4
          z[z> 1 - 1e-4] <- 1- 1e-4

          if(z[1]>z[2]){
            z[1] <- .4
            z[2] <- .6
          }
          item.parm[[j]] <- z

        }

      }else{
        lg <- c(gdina_match_matrix(t(reduced.profiles[[j]]),t(alpha[col.loc,])))
        if(length(unique(lg))!=ncol(reduced.profiles[[j]])){
          temp <- aggregate(dat[j,],list(lg),mean,na.rm = TRUE)
          item.parm[[j]] <- rep(0.5,ncol(reduced.profiles[[j]]))
          item.parm[[j]][temp$Group.1] <- temp$x
        }else{
          item.parm[[j]] <- aggregate(dat[j,],list(lg),mean,na.rm = TRUE)$x
        }
        for(k in 1:length(item.parm[[j]])){
          if(item.parm[[j]][k]<1e-4){
            item.parm[[j]][k] <- 1e-4
          }else if(item.parm[[j]][k]> 1- 1e-4){
            item.parm[[j]][k] <- 1- 1e-4
          }
        }
      }
    }

    ip[[b]] <- unlist(item.parm)
    # alpha sampling
    alpha <- seqGibbs(alpha,Q,gdina_l2m(item.parm),dat,reduced.profiles)
    if(!ip.only){
      a[[b]] <- alpha
    }
  }

  list(ip=do.call(cbind,ip),alpha=alpha,a=a,item.parm=item.parm)

}
