#######  Undocumented "tpx" utility functions #########

## ** Only referenced from topics.R


## check counts (can be an object from tm, slam, or a simple co-occurance matrix)
CheckCounts <- function(counts){
  if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
  if(is.null(dimnames(counts)[[1]])){ dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
  if(is.null(dimnames(counts)[[2]])){ dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }
  empty <- row_sums(counts) == 0
  if(sum(empty) != 0){
    counts <- counts[!empty,]
    cat(paste("Removed", sum(empty), "blank documents.\n")) }
  return(as.simple_triplet_matrix(counts))
}


## theta initialization
smash.tpxinit <- function(X, initheta, K1, alpha, verb, nbundles=1,
                    use_squarem=FALSE, init.adapt){
## initheta can be matrix, or c(nK, tmax, tol, verb)

  if(is.matrix(initheta)){
    if(ncol(initheta)!=K1){ stop("mis-match between initheta and K.") }
    if(prod(initheta>0) != 1){ stop("use probs > 0 for initheta.") }
    return(smash.normalizetpx(initheta, byrow=FALSE)) }

  if(is.matrix(alpha)){
    if(nrow(alpha)!=ncol(X) || ncol(alpha)!=K1){ stop("bad matrix alpha dimensions; check your K") }
    return(smash.normalizetpx(alpha, byrow=FALSE)) }

  if(is.null(initheta)){ ilength <- K1-1 }else{ ilength <- initheta[1] }
  if(ilength < 1){ ilength <- 1 }

  ## set number of initial steps
  if(length(initheta)>1){ tmax <- initheta[2] }else{ tmax <- 3 }
  ## set the tolerance
  if(length(initheta)>2){ tol <- initheta[3] }else{ tol <- 0.5 }
  ## print option
  if(length(initheta)>3){ verb <- initheta[4] }else{ verb <- 0 }


  if(verb){ cat("Building initial topics")
            if(verb > 1){ cat(" for K = ") }
            else{ cat("... ") } }

  nK <- length( Kseq <-  unique(ceiling(seq(2,K1,length=ilength))) )

  if(!init.adapt){
  initheta <- smash.tpxThetaStart(X, matrix(col_sums(X)/sum(X), ncol=1), matrix(rep(1,nrow(X))), K1)
#  return(initheta)
  } else{
  initheta <- smash.tpxThetaStart(X, matrix(col_sums(X)/sum(X), ncol=1), matrix(rep(1,nrow(X))), 2)

  if(verb > 0)
    { cat("\n")
      print(list(Kseq=Kseq, tmax=tmax, tol=tol)) }

  ## loop over topic numbers
  for(i in 1:nK){

    ## Solve for map omega in NEF space
    fit <- smash.tpxfit(X=X, theta=initheta, alpha=alpha, tol=tol, verb=verb,
                  admix=TRUE, method_admix=1, grp=NULL, tmax=tmax, wtol=-1, qn=-1,
                  nbundles = nbundles, use_squarem = FALSE, light=FALSE)
    if(verb>1){ cat(paste(Kseq[i],",", sep="")) }

    if(i<nK){ initheta <- smash.tpxThetaStart(X, fit$theta, fit$omega, Kseq[i+1]) }else{ initheta <- fit$theta }
  }
  if(verb){ cat("done.\n") }
#  return(initheta)
  }
  return(initheta)
}


## ** called from topics.R (predict) and smash.tpx.R
## Conditional solution for topic weights given theta
smash.tpxweights <- function(n, p, xvo, wrd, doc, start,
                             theta, verb=FALSE, nef=TRUE, wtol=10^{-5},
                             tmax=1000)
{
  K <- ncol(theta)
  start[start == 0] <- 0.1/K
  start <- start/rowSums(start)
  omega <- .C("Romega",
              n = as.integer(n),
              p = as.integer(p),
              K = as.integer(K),
              doc = as.integer(doc),
              wrd = as.integer(wrd),
              X = as.double(xvo),
              theta = as.double(theta),
              W = as.double(t(start)),
              nef = as.integer(nef),
              tol = as.double(wtol),
              tmax = as.integer(tmax),
              verb = as.integer(verb),
              PACKAGE="smashtpx")
  return(t(matrix(omega$W, nrow=ncol(theta), ncol=n))) }

## ** Called only in smash.tpx.R

smash.tpxsquarEM <- function(param_vec_in, X, m, K,
                       alpha, admix, method_admix, grp){
 omega_in <- inv.logit(matrix(param_vec_in[1:(nrow(X)*K)], nrow=nrow(X), ncol=K));
#  omega_in <- matrix(param_vec_in[1:(nrow(X)*K)], nrow=nrow(X), ncol=K);
 theta_in <- inv.logit(matrix(param_vec_in[-(1:(nrow(X)*K))], nrow=ncol(X), ncol=K))
#  theta_in <- matrix(param_vec_in[-(1:(nrow(X)*K))], nrow=ncol(X), ncol=K);
 out <- smash.tpxEM(X, m, theta_in, omega_in, alpha, admix, method_admix, grp);
 param_vec_out <- c(as.vector(logit(out$omega)),as.vector(logit(out$theta)))
# param_vec_out <- c(as.vector(out$omega),as.vector(out$theta))
 return(param_vec_out)
}


## Quasi Newton update for q>0



## log marginal likelihood
smash.tpxML <- function(X, theta, omega, alpha, L, dcut, admix=TRUE, grp=NULL){
  ## get the indices
  K <- ncol(theta)
  p <- nrow(theta)
  n <- nrow(omega)

  theta[theta==1] <- 1 - 1e-14;
  theta[theta==0] <- 1e-14;
  theta <- smash.normalizetpx(theta, byrow = FALSE)

  omega[omega==1] <- 1 - 1e-14;
  omega[omega==0] <- 1e-14;
  omega <- smash.normalizetpx(omega, byrow = TRUE)


  ## return BIC for simple finite mixture model
  if(!admix){
    qhat <- smash.tpxMixQ(X, omega, theta, grp, qhat=TRUE)$qhat
    ML <- sum(X$v*log(row_sums(qhat[X$i,]*theta[X$j,])))
    return( ML - 0.5*( K*p + (K-1)*n )*log(sum(X)) ) }

  ML <- L  + lfactorial(K) # lhd multiplied by label switching modes

  ## block-diagonal approx to determinant of the negative log hessian matrix
  q <- smash.tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)
  D <- smash.tpxHnegDet(X=X, q=q, theta=theta, omega=omega, alpha=alpha)

  D[D < dcut] <- dcut
  ML <- ML - 0.5*sum( D  )   # -1/2 |-H|

  ML <- ML + (K*p + sum(omega>0.01))*log(2*pi)/2  # (d/2)log(2pi)
  if(is.null(nrow(alpha))){ # theta prior normalizing constant
    ML <- ML + K*( lgamma(p*(alpha+1)) - p*lgamma(alpha+1) )  }
  else{ ML <- ML + sum(lgamma(col_sums(alpha+1)) - col_sums(lgamma(alpha+1))) } # matrix version
  ## omega(phi) prior normalizing constant number of parameters
  ML <- ML +  sum(D[-(1:p)]>dcut)*( lfactorial(K) - K*lgamma( 1+1/K ) ) #

  return(ML) }

## find residuals for X$v
smash.tpxResids <- function(X, theta, omega, grp=NULL, nonzero=TRUE)
{
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }

  m <- row_sums(X)
  K <- ncol(theta)
  n <- nrow(X)
  phat <- sum(col_sums(X)>0)
  d <- n*(K-1) + K*( phat-1 )

  if(nrow(omega) == nrow(X)){
    qhat <- smash.tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)
    xhat <- qhat*m[X$i]
  } else{
    q <- smash.tpxMixQ(X=X, omega=omega, theta=theta, grp=grp, qhat=TRUE)$qhat
    qhat <- row_sums(q[X$i,]*theta[X$j,])
    xhat <- qhat*m[X$i] }

  if(nonzero || nrow(omega) < nrow(X)){
    ## Calculations based on nonzero counts
    ## conditional adjusted residuals
    e <- X$v^2 - 2*(X$v*xhat - xhat^2)
    s <- qhat*m[X$i]*(1-qhat)^{1-m[X$i]}
    r <- sqrt(e/s)
    df <- length(r)*(1-d/(n*phat))
    R <- sum(r^2)
  }
  else{
    ## full table calculations
    e <- (X$v^2 - 2*X$v*m[X$i]*qhat)
    s <- m[X$i]*qhat*(1-qhat)
    fulltable <- .C("RcalcTau",
                 n = as.integer(nrow(omega)),
                 p = as.integer(nrow(theta)),
                 K = as.integer(ncol(theta)),
                 m = as.double(m),
                 omega = as.double(omega),
                 theta = as.double(theta),
                 tau = double(1), size=double(1),
                 PACKAGE="smashtpx" )
    tau <- fulltable$tau
    R <- sum(e/s) + tau
    df <-  fulltable$size - phat  - d
    r <- suppressWarnings(sqrt(e/s + tau))
    r[is.nan(r)] <- 0 ## should not happen, but can theoretically
  }

  ## collect and output
  sig2 <- R/df
  rho <- suppressWarnings(pchisq(R, df=df, lower.tail=FALSE))
  D <- list(dispersion=sig2, pvalue=rho, df=df)
  return( list(s=s, e=e, r=r, D=D) ) }


## fast initialization functions for theta (after increasing K) and omega (given theta)
smash.tpxThetaStart <- function(X, theta, omega, K)
  {
    R <- smash.tpxResids(X, theta=theta, omega=omega, nonzero=TRUE)
    X$v <- R$e*(R$r>3) + 1/ncol(X)
    Kpast <- ncol(theta)
    Kdiff <- K-Kpast
    if(Kpast != ncol(omega) || Kpast >= K){ stop("bad K in smash.tpxThetaStart") }
    initheta <- smash.normalizetpx(Kpast*theta+rowMeans(theta), byrow=FALSE)
    n <- nrow(X)
    ki <- matrix(1:(n-n%%Kdiff), ncol=Kdiff)
    for(i in 1:Kdiff){ initheta <- cbind(initheta, (col_sums(X[ki[,i],])+1/ncol(X))/(sum(X[ki[,i],])+1)) }
    return( initheta )
  }

smash.tpxOmegaStart <- function(X, theta)
  {
    if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
    omega <- try(tcrossprod_simple_triplet_matrix(X, solve(t(theta)%*%theta)%*%t(theta)), silent=TRUE )
    if(inherits(omega,"try-error")){ return( matrix( 1/ncol(theta), nrow=nrow(X), ncol=ncol(theta) ) ) }
    omega[omega <= 0] <- .5
    return( smash.normalizetpx(omega, byrow=TRUE) )
  }


## fast computation of sparse P(X) for X>0
smash.tpxQ <- function(theta, omega, doc, wrd){

  theta[theta==1] <- 1 - 1e-14;
  theta[theta==0] <- 1e-14;
  theta <- smash.normalizetpx(theta, byrow = FALSE)

  omega[omega==1] <- 1 - 1e-14;
  omega[omega==0] <- 1e-14;
  omega <- smash.normalizetpx(omega, byrow = TRUE)

  if(length(wrd)!=length(doc)){stop("index mis-match in tpxQ") }
  if(ncol(omega)!=ncol(theta)){stop("theta/omega mis-match in smash.tpxQ") }

  out <- .C("RcalcQ",
            n = as.integer(nrow(omega)),
            p = as.integer(nrow(theta)),
            K = as.integer(ncol(theta)),
            doc = as.integer(doc-1),
            wrd = as.integer(wrd-1),
            N = as.integer(length(wrd)),
            omega = as.double(omega),
            theta = as.double(theta),
            q = double(length(wrd)),
            PACKAGE="smashtpx" )

  return( out$q ) }

## model and component likelihoods for mixture model
smash.tpxMixQ <- function(X, omega, theta, grp=NULL, qhat=FALSE){
  if(is.null(grp)){ grp <- rep(1, nrow(X)) }

  theta[theta==1] <- 1 - 1e-14;
  theta[theta==0] <- 1e-14;
  theta <- smash.normalizetpx(theta, byrow = FALSE)

  omega[omega==1] <- 1 - 1e-14;
  omega[omega==0] <- 1e-14;
  omega <- smash.normalizetpx(omega, byrow = TRUE)

  K <- ncol(omega)
  n <- nrow(X)
  mixhat <- .C("RmixQ",
               n = as.integer(nrow(X)),
               p = as.integer(ncol(X)),
               K = as.integer(K),
               N = as.integer(length(X$v)),
               B = as.integer(nrow(omega)),
               cnt = as.double(X$v),
               doc = as.integer(X$i-1),
               wrd = as.integer(X$j-1),
               grp = as.integer(as.numeric(grp)-1),
               omega = as.double(omega),
               theta = as.double(theta),
               Q = double(K*n),
               PACKAGE="smashtpx")
  ## model and component likelihoods
  lQ <- matrix(mixhat$Q, ncol=K)
  lqlhd <- log(row_sums(exp(lQ)))
  lqlhd[is.infinite(lqlhd)] <- -600 # remove infs
  if(qhat){
    qhat <- exp(lQ-lqlhd)
    ## deal with numerical overload
    infq <- row_sums(qhat) < .999
    if(sum(infq)>0){
      qhat[infq,] <- 0
      qhat[n*(apply(matrix(lQ[infq,],ncol=K),1,which.max)-1) + (1:n)[infq]] <- 1 }
  }
  return(list(lQ=lQ, lqlhd=lqlhd, qhat=qhat)) }

## negative log hessian block diagonal matrix for theta & omega
smash.tpxHnegDet <- function(X, q, theta, omega, alpha){
  K <- ncol(theta)
  n <- nrow(omega)

  ## sparse Xij/Qij^2
  Xq <- X
  Xq$v <- Xq$v/q^2

  ## negative 2nd derivitive matrices for theta
  HT <- tcrossprod_simple_triplet_matrix(t(Xq), apply(omega, 1, function(v) v%o%v ) )
  HT[,K*(0:(K-1))+1:K] <- HT[,K*(0:(K-1))+1:K] + alpha/theta^2 # will break for alpha<=1
  DT <- apply(HT, 1, smash.tpxlogdet)

  ## ditto for omega
  HW <- matrix(.C("RnegHW",
                  n = as.integer(nrow(omega)),
                  p = as.integer(nrow(theta)),
                  K = as.integer(K-1),
                  omeg = as.double(omega[,-1]),
                  thet = as.double(theta[,-1]),
                  doc = as.integer(X$i-1),
                  wrd = as.integer(X$j-1),
                  cnt = as.double(X$v),
                  q = as.double(q),
                  N = as.integer(length(q)),
                  H = double(n*(K-1)^2),
                  PACKAGE="smashtpx")$H,
               nrow=(K-1)^2, ncol=n)
  DW <- apply(HW, 2, smash.tpxlogdet)
  return( c(DT,DW) )  }

## functions to move theta/omega to and from NEF.
smash.tpxToNEF <- function(theta, omega){
  n <- nrow(omega)
  p <- nrow(theta)
  K <- ncol(omega)
  return(.C("RtoNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=double((p-1)*K + n*(K-1)),
            theta=as.double(theta), tomega=as.double(t(omega)),
            PACKAGE="smashtpx")$Y)
}

## 'From' NEF representation back to probabilities
smash.tpxFromNEF <- function(Y, n, p, K){
  bck <- .C("RfromNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=as.double(Y), theta=double(K*p), tomega=double(K*n),
            PACKAGE="smashtpx")
  return(list(omega=t( matrix(bck$tomega, nrow=K) ), theta=matrix(bck$theta, ncol=K)))
}

## utility log determinant function for speed/stabilty
smash.tpxlogdet <- function(v){
    v <- matrix(v, ncol=sqrt(length(v)))
    if( sum(zeros <- colSums(v)==0)!=0 ){
      cat("warning: boundary values in laplace approx\n")
      v <- v[-zeros,-zeros] }

    return(determinant(v, logarithm=TRUE)$modulus)
}

