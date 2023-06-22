euspca = function(X,is.data.mat=TRUE,k,lambda,scale=FALSE,...) {
  
  ## data setup
  X = as.matrix(X)
  if (is.data.mat==FALSE && eigen(X)$values[nrow(X)] < 0)
    stop("The covariance (or corrleation) matrix X is not a positive semi-definite matrix.")
  
  ## If p<n, switch to covariance matrix 
  if(is.data.mat == TRUE){
    n = nrow(X); p = ncol(X)
    if(p<n) { X = t(X) %*% X / n
              is.data.mat = FALSE 
    }
  }  
  
  ## run euspca
  if(is.data.mat==FALSE){
    outlist = euspca_sig(Sig=X, k=k, lambda=lambda, scale=scale, ...)
  } 
  
  if(is.data.mat==TRUE){
    outlist = euspca_dat(dat=X, k=k, lambda=lambda, scale=scale, ...)
  }

  outlist
}
