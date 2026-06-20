#'@title Stochastic EM Algorithm for High-Dimensional Cognitive Diagnosis Models
#'
#' @description This function estimates the model parameters of CDM using the stochastic EM algorithm with sequential Gibbs sampling.
#'
#' @param dat A required \eqn{N \times J} data matrix of N examinees to J items. Missing values are currently not allowed.
#' @param Q A required binary item and attribute association matrix.A required \eqn{J \times K} item-attribute association matrix, where J is the number of items and K is the number of attributes.
#' @param item.parm A initial item parameter estimates.
#' @param maxitr The maximum iterations allowed.
#' @param eps1 A positive convergence criterion used during the burn-in stage of the stochastic EM algorithm.
#' @param eps2 A positive convergence criterion used after the burn-in.
#' @param frac1 Proportion of the first part of the Markov chain used in the Geweke diagnostic.
#' @param frac2 Proportion of the last part of the Markov chain used in the Geweke diagnostic.
#'
#' @return a list with elements
#' \describe{
#' \item{catprob.parm}{estimated item category response probabilities}
#' \item{ip}{estimated final parameter obtained by averaging from each batches}
#' \item{alpha}{estimated attribute mastery profiles}
#' \item{alpha.list}{latent attribute profile estimates from each batches}
#' \item{total.number.of.batch}{total number of retained batches used for estimation}
#' \item{final.chain}{final length of the Markov chain}
#' \item{burn.in.size}{the number of burn-in iterations discarded}
#' \item{plist}{estimated item parameters from each batch}
#' }
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
#' }
#'@export

autoGDINA.stEM <- function(dat, Q, item.parm = NULL, maxitr = 100,
                           eps1=2,eps2=0.4,frac1=.1,frac2=.5){
  Rcpp::sourceCpp("seqGibbs.cpp")

  dat0 <- as.matrix(dat)
  Q0 <- as.matrix(Q)
  dat <- t(dat0)
  Q <- t(Q0)

  K <- nrow(Q)

  Kj <- colSums(Q > 0)  # vector with length of J
  Lj <- 2^Kj
  J <- nrow(dat)
  N <- ncol(dat)
  reduced.profiles <- lapply(Kj,function(x) t(GDINA::attributepattern(x)))

  which.k.required <- list()
  for(j in 1:J){
    which.k.required[[j]] <- which(Q[,j]==1)
  }
  L <- 2^K  # The number of latent classes

  if(is.null(item.parm)){
    if(any(is.na(dat))){
      alpha <- (matrix(runif(K*N),ncol = N) > .5)*1
    }else{
      alpha <- (Q%*%dat)/rowSums(Q)
      alpha <- 1 * (alpha>=matrix(runif(K*N),ncol = N))
    }
    item.parm <- list()
  }else{
    prior <- rep(1/L, L)
    logprior <- log(prior)
    parloc <- GDINA::LC2LG(Q=Q0)
    estep <- GDINA:::LikNR(as.matrix(GDINA:::l2m(item.parm)),
                           as.matrix(dat0), as.matrix(logprior), rep(1,N),
                           as.matrix(parloc), rep(1,N), FALSE)

    lc <- apply(exp(estep$logpost),1,function(x)sample(1:L,1,replace = TRUE,prob=x))
    alpha <- t(AlphaPattern[lc,])
  }


  npar <- sum(Lj)
  ###########################
  #
  # Warm-up phase - to be discarded
  #
  ###########################

  x <- kernel(dat,Q,J,reduced.profiles,alpha,20,ip.only=TRUE)

  item.parm <- x$item.parm
  alpha <- x$alpha # most recent alpha sample

  ###########################
  #
  # burn-in phase
  #
  ###########################


  alist <- plist <- list()

  total.number.of.batch <- 0
  M <- 10
  B <- 20
  burn.in.size <- 0
  # ip.output <- matrix(NA,nrow = npar,ncol = maxitr)
  # the initial MxB iterations
  cat("\nBurn-in phase:")
  for(m in 1:M){

    cat("\n  # of batch = ",m," with batch size of ", B)

    x <- kernel(dat,Q,J,reduced.profiles,alpha,B,ip.only=F)

    plist[[m]] <- x$ip
    alist[[m]] <- x$a
    alpha <- x$alpha # most recent alpha sample

    # if(m>=maxitr) break

  }

  #terminates burn-in?
  mcmc.par <- coda::mcmc(t(Reduce(cbind,plist)))
  z <- coda::geweke.diag(mcmc.par,frac1 = frac1,frac2 = frac2)$z
  z <- z[!is.na(z)]
  #m == M
  while(sum(z^2)/length(z)>=eps1&&m<60){
    #m<=60 means at most 50 batch (50x20=1000 iterations) will be burn-in; 10 batch or 200 iterations will be kept
    ## additional runs for burn-in?
    if(m>M){
      cat("\n  # of batch = ",m," sum z^2 / npar = ",sum(z^2)/npar, "eps 1 criterion = ",eps1)
    }

    x <- kernel(dat,Q,J,reduced.profiles,alpha,B,ip.only=F)

    plist <- plist[-1] # remove first burn-in batch
    alist <- alist[-1]
    plist[[M]] <- x$ip # add the new batch to the Mth
    alist[[M]] <- x$a
    alpha <- x$alpha # most recent alpha sample

    #terminates burn-in?
    mcmc.par <- coda::mcmc(t(Reduce(cbind,plist)))
    z <- coda::geweke.diag(mcmc.par,frac1 = frac1,frac2 = frac2)$z
    z <- z[!is.na(z)]

    burn.in.size <- burn.in.size + B # num of iterations
    m <- m+1 #number of batch
  }
  total.number.of.batch <- m
  cat("\n End of Burn-in phase summary:")
  if(sum(z^2)/length(z)<eps1){
    cat("\n  burn-in ends because sum z^2 / npar = ",sum(z^2)/length(z),"< eps 1 criterion = ",eps1)
  }else{
    cat("\n burn-in ends because total # of batch = ",total.number.of.batch,">= max allowed = 60")
  }
  cat("\n Total burn-in iterations = ",burn.in.size)


  ###########################
  #
  # determining chain length
  #
  ###########################
  n <- M #M=10 batch kept
  d.hat <- batch.var(plist,n=n)
  cat("\nAfter burn-in phase:")
  while(max(d.hat*N)>=eps2&&n<50){
    #n<=50 means adding at most 40 batch, resulting at most 50 batches or 1000 iterations
    cat("\n  # of valid batch after burn-in = ",n," max delta hat = ",max(d.hat)*N, "eps 2 criterion = ",eps2)
    n <- n+1

    x <- kernel(dat,Q,J,reduced.profiles,alpha,B,ip.only=F)

    plist[[n]] <- x$ip # add the new batch to the Mth
    alist[[n]] <- x$a
    alpha <- x$alpha # most recent alpha sample

    d.hat <- batch.var(plist,n=n)
  }
  cat("\n  # of valid batch after burn in = ",n," max delta hat = ",max(d.hat)*N, " criterion = ",eps2)

  cat("\nLenth of final MC chain = ",n*B)
  total.number.of.batch <- total.number.of.batch + n - M

  phat <- rowMeans(Reduce(cbind,plist))

  list.alpha <- do.call(c,alist)
  est.alpha2 <- t(Reduce("+",list.alpha))/length(list.alpha)

  return(list(catprob.parm=x$item.parm,ip=phat,alpha=1 * (est.alpha2 > 0.5),alpha.list=alist,total.number.of.batch=total.number.of.batch,
              final.chain=n*B,burn.in.size=burn.in.size,plist=plist))

}
