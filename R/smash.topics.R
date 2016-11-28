##### Estimation for Topic Models ######

## intended main function; provides defaults and selects K via marginal lhd
smash.topics <- function(counts,
                   K,
                   shape=NULL,
                   initopics=NULL,
                   tol=0.1,
                   bf=FALSE,
                   kill=2,
                   ord=TRUE,
                   verb=1,
                   admix=TRUE,
                   nbundles=1,
                   use_squarem=FALSE,
                   init.adapt=FALSE,
                   light=1,
                   sample_init=TRUE,
                   init.method = "taddy",
                   smash_gap = 100,
                   smash_method = "gaussian",
                   use_flash_ti=FALSE,
                   with_qn = FALSE,
                   K_flash = 5,
                   reflect = FALSE,
                   wtol=10^{-4},
                   qn=100,
                   tmax=10000,...)
  ## tpxselect defaults: tmax=10000, wtol=10^(-4), qn=100, grp=NULL,
  ## nonzero=FALSE, dcut=-10, top_genes=100, burn_in=5
{

  ###  add more features to make number of features a power of 2
  ###  if reflect is FALSE, we just add zeros as features
  ###  if reflect is TRUE, we reflect back the data for last few features


  ceil <- ceiling(log(dim(counts)[2])/log(2));
  if(log(dim(counts)[2])%%log(2)!=0) {
    cat(sprintf("number of features not a power of 2"));
    if(reflect){
      fcounts <- cbind(counts, counts[,dim(counts)[2]-(1:(2^{ceil}-dim(counts)[2]))]);
    }
    if(!reflect){
      fcounts <- cbind(counts, matrix(0, dim(counts)[1], 2^{ceil}-dim(counts)[2]));
    }}else{
      fcounts <- counts;
    }

  levels <- ceil+1;
  X <- CheckCounts(fcounts)

  library(smashr)
  library(ashr)


  p <- ncol(X)
  n <- nrow(X)

  if(verb>0)
    cat(sprintf("\nEstimating on a %d document collection.\n", nrow(X)))

  ## check the prior parameters for theta
  if(prod(shape>0) != 1){ stop("use shape > 0\n") }

  ## check the list of candidate K values
  if(prod(K>1)!=1){ stop(cat("use K values > 1\n")) }
  K <- sort(K)

  ## Null model log probability
  sx <- sum(X)
  qnull <- col_sums(X)/sx
  null <- sum( X$v*log(qnull[X$j]) ) - 0.5*(n+p)*(log(sx) - log(2*pi))

  ## initialize


  if(init.method=="taddy"){
      index_init <- 1:min(ceiling(nrow(X)*.05),100);
      if(sample_init==TRUE){
            samp_length <- length(index_init);
            index_init <- sample(1:nrow(X),samp_length);
       }

  ## initialize
  if(init.adapt==FALSE){

  initopics <- smash.tpxinit(X[index_init,], initopics, K[1],
                       shape, verb, nbundles=1, use_squarem=FALSE,
                       init.adapt)
    #initopics <- t(gtools::rdirichlet(4, rep(1+ 1/K*B, B)))
  }else{
 #   if(change_start_points){
 #      initopics <- smash.tpxinit(X[1:min(ceiling(nrow(X)*.05),100),], initopics, K[1]+3,
 #                          shape, verb, nbundles=1, use_squarem=FALSE, init.adapt)
 #      initopics <- initopics[,sort(sample(1:(K[1]+2), K[1], replace=FALSE))];
 #   }else{
      initopics <- smash.tpxinit(X[index_init,], initopics, K[1],
                           shape, verb, nbundles=1, use_squarem=FALSE,
                           init.adapt)
 #    }
  }}

  if(init.method=="kmeans"){
    kmeans.init=kmeans(fcounts, K, nstart=5, iter.max=100)
    phi0 = kmeans.init$centers;
    phi0 = t(apply(phi0, 1, function(x) return(x/sum(x))))
    theta <- t(phi0);
    omega = rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))
  }

   alpha <- 1/(K*p)
   if(is.matrix(alpha)){ if(nrow(alpha)!=p || ncol(alpha)!=K){ stop("bad matrix alpha dimensions") }}

   L <- smash.tpxlpost(y=y, theta=theta, omega=omega,
                      alpha=alpha, admix=admix, grp=grp);
   iter <- 0
   dif <- tol+1+qn
   update <- TRUE
   row_total <- tapply(X$v, X$i, sum);
   if(verb){ cat(paste("Fitting the",K,"topic model.\n")) }
   if(verb>0){
     cat("log posterior increase: " )
     digits <- max(1, -floor(log(tol, base=10))) }
   y <- as.matrix(X);
   Y <- NULL;

  while( update  && iter < tmax ){

    omega <- smash.normalizetpx(omega + 1e-15, byrow=TRUE);
    theta <- smash.normalizetpx(theta + 1e-15, byrow=FALSE);

    plot(theta[,1], type="l", col="red");
    lines(theta[,2], col="blue");
    lines(theta[,3], col="green")

    m <- row_sums(X);
    moveEM <- smash.tpxEM(y=y, m=m, theta=theta, omega=omega,
                        alpha=alpha, admix=admix, grp=grp)
    lambda <- moveEM$lambda;
    lambda.unsmoothed <- moveEM$lambda;

    if(with_qn){
       moveQN <- list(theta = moveEM$theta, omega = moveEM$omega);
       QNup <- smash.tpxQN(move=moveQN, Y=Y, y=y, alpha=alpha, verb=verb,
                            admix=admix, grp=grp, doqn=qn-dif)
       move <- QNup$move
       lambda <- t(move$theta)*moveEM$lscale;
    }

    #  L_new <- smash.tpxlpost(y=y, theta=move$theta, omega=move$omega, alpha=alpha, admix=admix, grp=grp)
    #  QNup <- list("move"=move, "L"=L_new, "Y"=NULL)
    ## quasinewton-newton acceleration
    # moveQN <- list(theta = moveEM$theta, omega = moveEM$omega);
    # QNup <- smash.tpxQN(move=moveQN, Y=Y, y=y, alpha=alpha, verb=verb,
    #                      admix=admix, grp=grp, doqn=qn-dif)
    # move <- QNup$move
    # lambda <- t(move$theta)*moveEM$lscale;

    if(iter %% smash_gap==0){

      if(smash_method=="poisson"){
          if(use_flash_ti){
            ti_tab <- smashr::TI_table_construct(lambda, cxx=TRUE, K_flash=K_flash)
            lambda=smooth.lambda(lambda, optional_ti_table = ti_tab)
          }else{
            lambda=smooth.lambda(lambda, optional_ti_table = NULL)
          }
          lambda[is.na(lambda)]=lambda.unsmoothed[is.na(lambda)]
          phi_smoothed=lambda/moveEM$lscale
          move <- list(theta=t(phi_smoothed), omega=moveEM$omega)
          QNup_L <-  smash.tpxlpost(y, theta = move$theta,
                                  omega = move$omega, alpha=alpha,
                                  admix=admix, grp=grp)
        }else if(smash_method=="gaussian"){
            z_leaf_est <- moveEM$theta
            z_leaf_smoothed <- do.call(cbind, lapply(1:dim(z_leaf_est)[2],
            function(k)
            {
              if(sum(z_leaf_est[,k])>0){
                out <- suppressMessages(smashr::smash.gaus(z_leaf_est[,k],
                            ashparam = list(control=list(maxiter=50))))
                out[ out < 0] = 0
                return(out)
                }else{
                  return(z_leaf_est[,k])
              }
            }))
          theta_smoothed <- smash.normalizetpx(z_leaf_smoothed+1e-06,
                                               byrow=FALSE)
          move <- list(theta=theta_smoothed, omega=moveEM$omega)
          QNup_L <-  smash.tpxlpost(y, theta = move$theta,
                                    omega = move$omega,
                                    alpha=alpha, admix=admix, grp=grp)
          }
      }

      dif <- abs(QNup_L-L)

      L <- QNup_L


      ## check convergence
      if(abs(dif) < tol){
        if(sum(abs(theta-move$theta)) < tol){ update = FALSE } }

      ## print
      if(verb>0 && (iter-1)%%ceiling(1/verb)==0 && iter>0){
          ##if(verb>0 && iter>0){
        cat( paste( round(dif,digits), #" (", sum(abs(theta-move$theta)),")",
                      ", ", sep="") ) }

     ## heartbeat for long jobs
      if(((iter+1)%%1000)==0){
          cat(sprintf("p %d iter %d diff %g\n",
                    nrow(move$theta), iter+1,round(dif))) }

      ## iterate
      iter <- iter+1
      theta <- move$theta
      omega <- move$omega

  }

  if(smash_method=="poisson"){
    lambda=smooth.lambda(moveEM$lambda)
    lambda[is.na(lambda)]=lambda.unsmoothed[is.na(lambda)]
    phi_smoothed=lambda/moveEM$lscale
    theta_smoothed <- t(phi_smoothed);
    move <- list(theta=theta_smoothed, omega=omega)
    L <-  smash.tpxlpost(y, theta = move$theta,
                              omega = move$omega, alpha=alpha,
                              admix=admix, grp=grp)
  }

  if(smash_method=="gaussian"){
    z_leaf_est <- move$theta
    z_leaf_smoothed <- do.call(cbind, lapply(1:dim(z_leaf_est)[2],
    function(k)
    {
      if(sum(z_leaf_est[,k])>0){
        out <- suppressMessages(smashr::smash.gaus(z_leaf_est[,k]))
        out[ out < 0] = 0
        return(out)
      }else{
        return(z_leaf_est[,k])
      }
    }))
    theta_smoothed <- smash.normalizetpx(z_leaf_smoothed+1e-10, byrow=FALSE)
    L <-  smash.tpxlpost(y, theta = theta_smoothed, omega = omega,
                         alpha=alpha, admix=admix, grp=grp)
  }

  ## final log posterior
  ## summary print
  if(verb>0){
    cat("done.")
    if(verb>1) { cat(paste(" (L = ", round(L,digits), ")", sep="")) }
    cat("\n")
  }

  tpx <- list(theta=move$theta,
              omega=move$omega,
              K=K,
              alpha=alpha,
              L=L,
              iter=iter)

  K <- tpx$K

  ## clean up and out
  if(ord){ worder <- order(col_sums(tpx$omega), decreasing=TRUE) } # order by decreasing usage
  else{ worder <- 1:K }
  ## Main parameters
  theta=matrix(tpx$theta[,worder], ncol=K, dimnames=list(phrase=dimnames(X)[[2]], topic=paste(1:K)) )
  omega=matrix(tpx$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(nrow(omega)==nrow(X)){ dimnames(omega)[[1]] <- dimnames(X)[[1]] }
  theta = theta[1:dim(counts)[2],];
  ## topic object
  out <- list(K=K, theta=theta, omega=omega, X=X)
  class(out) <- "topics"
  invisible(out)

}


## single EM update. two versions: admix and mix
smash.tpxEM <- function(y, m, theta, omega, alpha, admix, grp)
{

  n <- nrow(y)
  p <- ncol(y)
  K <- ncol(theta)

  phi <- t(theta);
  pi <- omega

  gamma=pi[rep(1:n,each=p),]*t(phi)[rep(1:p,n),]
  gamma=gamma/rowSums(gamma)
  gamma[is.na(gamma)]=1/K
  gammab=(as.vector(t(y))%o%rep(1,K))*gamma
  pi.num=t(apply(array(gammab,dim=c(p,n,K)),2,colSums))
  pi=pi.num/(rowSums(y)%o%rep(1,K))
  ybt=t(apply(array(gammab,dim=c(p,n,K)),1,colSums))
  theta=ybt/(rep(1,p)%o%colSums(gammab))

  lscale=((colSums(ybt)/colSums(pi))%o%rep(1,p))
  lambda=t(theta)*lscale

  omega <- smash.normalizetpx(pi, byrow=TRUE)
  theta <- smash.normalizetpx(theta, byrow=FALSE);

  return(list(theta=theta, omega=omega, lscale=lscale, lambda=lambda))
}

## Quasi Newton update for q>0
smash.tpxQN <- function(move, Y, y, alpha, verb, admix, grp, doqn)
{
  move$theta[move$theta==1] <- 1 - 1e-14;
  move$omega[move$omega==1] <- 1 - 1e-14;
  move$omega[move$omega==0] <- 1e-14;
  move$theta[move$theta==0] <- 1e-14;
  move$theta <- smash.normalizetpx(move$theta, byrow = FALSE)
  move$omega <- smash.normalizetpx(move$omega, byrow = TRUE)

  ## always check likelihood
  L <- smash.tpxlpost(y=y, theta=move$theta, omega=move$omega,
                      alpha=alpha, admix=admix, grp=grp)

  if(doqn < 0){ return(list(move=move, L=L, Y=Y)) }

  ## update Y accounting
  Y <- cbind(Y, smash.tpxToNEF(theta=move$theta, omega=move$omega))
  if(ncol(Y) < 3){ return(list(Y=Y, move=move, L=L)) }
  if(ncol(Y) > 3){ warning("mis-specification in quasi-newton update; please report this bug.") }

  ## Check quasinewton secant conditions and solve F(x) - x = 0.
  U <- as.matrix(Y[,2]-Y[,1])
  V <- as.matrix(Y[,3]-Y[,2])
  sUU <- sum(U^2)
  sVU <- sum(V*U)
  Ynew <- Y[,3] + V*(sVU/(sUU-sVU))
  qnup <- smash.tpxFromNEF(Ynew, n=nrow(move$omega),
                           p=nrow(move$theta), K=ncol(move$theta))
  ## check for a likelihood improvement
  Lqnup <- try(smash.tpxlpost(y=y, theta=qnup$theta, omega=qnup$omega,
                              alpha=alpha, admix=admix, grp=grp), silent=TRUE)

  if(inherits(Lqnup, "try-error")){
    if(verb>10){ cat("(QN: try error) ") }
    return(list(Y=Y[,-1], move=move, L=L)) }

  if(verb>10){ cat(paste("(QN diff ", round(Lqnup-L,3), ")\n", sep="")) }

  if(Lqnup < L){
    return(list(Y=Y[,-1], move=move, L=L)) }
  else{
    L <- Lqnup
    Y <- cbind(Y[,2],Ynew)
    return( list(Y=Y, move=qnup, L=L) )
  }
}

smash.tpxlpost_squarem <- function(param_vec_in,  y, m, K,
                                   alpha, admix=TRUE, method_admix, grp=NULL)
{
  omega_in <- inv.logit(matrix(param_vec_in[1:(nrow(X)*K)], nrow=nrow(X), ncol=K));
  #  omega_in <- matrix(param_vec_in[1:(nrow(X)*K)], nrow=nrow(X), ncol=K);
  theta_in <- inv.logit(matrix(param_vec_in[-(1:(nrow(X)*K))], nrow=ncol(X), ncol=K))
  #  theta_in <- matrix(param_vec_in[-(1:(nrow(X)*K))], nrow=ncol(X), ncol=K);
  return(smash.tpxlpost(y, theta_in, omega_in, alpha, admix, grp))
}


## unnormalized log posterior (objective function)
smash.tpxlpost <- function(y, theta, omega, alpha, admix=TRUE, grp=NULL)
{
  theta[theta==1] <- 1 - 1e-10;
  omega[omega==1] <- 1 - 1e-10;
  omega[omega==0] <- 1e-10;
  theta[theta==0] <- 1e-10;
  theta <- smash.normalizetpx(theta, byrow = FALSE)
  omega <- smash.normalizetpx(omega, byrow = TRUE)
  K <- ncol(theta)

  L.ini=log(omega%*%t(theta))
  yL=y*L.ini
  yL[is.na(yL)]=0
  L = sum(yL)
  return(L) }

normalize=function(x){
  #if(sum(abs(x))!=0){
  return(x/sum(x))
  #}else{
  #  return(rep(0,length(x)))
  #}
}

smooth.lambda = function(lambda, optional_ti_table){
  #return(t(apply(lambda,1,ashsmooth.pois,cxx = FALSE)))
  if(is.null(optional_ti_table)){
   return(t(apply(lambda,1,smashr::smash.poiss, cxx = TRUE, optional_ti_table=optional_ti_table)))
  }else{
    out <- do.call(rbind, lapply(1:dim(lambda)[1], function(l) {
      tmp <- smashr::smash.poiss(lambda[l,], cxx = TRUE, optional_ti_table=optional_ti_table[l,]);
      return(tmp)
    }))
    return(out)
  }
}

library(slam)
CheckCounts <- function(counts){
  if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
  if(is.null(dimnames(counts)[[1]])){ dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
  if(is.null(dimnames(counts)[[2]])){ dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }
  empty <- slam::row_sums(counts) == 0
  if(sum(empty) != 0){
    counts <- counts[!empty,]
    cat(paste("Removed", sum(empty), "blank documents.\n")) }
  return(as.simple_triplet_matrix(counts))
}


## S3 method predict function
predict.topics <- function(object, newcounts, loglhd=FALSE, ...)
  ## tpxweights optional arguments and defauls are verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000
{
  if(is.vector(newcounts)){ newcounts <- matrix(newcounts, nrow=1) }
  if(class(newcounts)[1] == "TermDocumentMatrix"){ newcounts <- t(newcounts) }
  X <- as.simple_triplet_matrix(newcounts)

  if(!(class(object)%in%c("topics","matrix"))){ stop("object class must be `topics' or 'matrix'.") }

  if(class(object)=="topics"){
    theta <- object$theta
    if(nrow(theta) != ncol(X)){ stop("Dimension mismatch: nrow(theta) != ncol(X)") }
    if(nrow(object$X) != nrow(object$omega)) # simple mixture
      { Q <- matrix(tpxMixQ(X, omega=object$omega, theta=theta, ...)$lQ, ncol=ncol(theta))
        return( (1:ncol(theta))[apply(Q,1,which.max)] ) }
  }
  else{ theta <- object }

  start <- tpxOmegaStart(X=X, theta=theta)

  ## re-order xvec in doc-blocks, and build indices
  doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))
  xvo <- X$v[order(X$i)]
  wrd <- X$j[order(X$i)]-1

  W <- tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc, start=start, theta=theta, ...)

  if(loglhd){
    L <- sum( X$v*log(tpxQ(theta=theta, omega=W, doc=X$i, wrd=X$j)) )
    return(list(W=W, L=L)) }
  else{ return(W) }
}

