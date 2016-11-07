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
                   method_admix=1,
                   sample_init=TRUE,
                   smooth_gap = 10,
                   smooth_method = c("gaussian", "poisson"),
                   reflect = FALSE,
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


  index_init <- 1:min(ceiling(nrow(X)*.05),100);
  if(sample_init==TRUE){
    samp_length <- length(index_init);
    index_init <- sample(1:nrow(X),samp_length);
  }

  ## initialize
  if(init.adapt==FALSE){

  initopics <- tpxinit(X[index_init,], initopics, K[1],
                       shape, verb, nbundles=1, use_squarem=FALSE, init.adapt)
    #initopics <- t(gtools::rdirichlet(4, rep(1+ 1/K*p, p)))
  }else{
 #   if(change_start_points){
 #      initopics <- tpxinit(X[1:min(ceiling(nrow(X)*.05),100),], initopics, K[1]+3,
 #                          shape, verb, nbundles=1, use_squarem=FALSE, init.adapt)
 #      initopics <- initopics[,sort(sample(1:(K[1]+2), K[1], replace=FALSE))];
 #   }else{
      initopics <- tpxinit(X[index_init,], initopics, K[1],
                           shape, verb, nbundles=1, use_squarem=FALSE,
                           init.adapt)
 #    }
  }

  initopics[initopics==1] <- 1 - 1e-14;
  initopics[initopics==0] <- 1e-14;
  initopics <- normalizetpx(initopics, byrow = FALSE)

 # initopics <- initopics[,sort(sample(1:(K[1]+2), K[1], replace=FALSE))];
 # initopics <- initopics[,1:K[1]];
  ## either search for marginal MAP K and return bayes factors, or just fit
  tpx <- smash.tpxSelect(X, K, bf, initopics,
                         alpha=shape, tol, kill, verb, nbundles,
                         use_squarem,
                         light, tmax, admix=TRUE,
                         method_admix=1,
                         sample_init=TRUE,
                         smooth_gap = smooth_gap,
                         smooth_method = smooth_method,
                         grp=NULL, wtol=10^{-4}, qn=100,
                         nonzero=FALSE, dcut=-10, top_genes=150, burn_in=5)

  K <- tpx$K

  ## clean up and out
  if(ord){ worder <- order(col_sums(tpx$omega), decreasing=TRUE) } # order by decreasing usage
  else{ worder <- 1:K }
  ## Main parameters
  theta=matrix(tpx$theta[,worder], ncol=K, dimnames=list(phrase=dimnames(X)[[2]], topic=paste(1:K)) )
  omega=matrix(tpx$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(nrow(omega)==nrow(X)){ dimnames(omega)[[1]] <- dimnames(X)[[1]] }

  ## topic object
  out <- list(K=K, theta=theta, omega=omega, BF=tpx$BF, D=tpx$D, X=X)
  class(out) <- "topics"
  invisible(out) }


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

