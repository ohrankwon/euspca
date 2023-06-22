## When n<p

euspca_dat = function(dat, k, scale, lambda,
  eps1 = 0.05, eps2 = 0.0001, eps.sub = 0.005,
  max.iter.outer = 1000, max.iter.inner = 2000, track=5,
  parm.outer = list(sig = 0.2, tau = 1.1), 
  parm.inner = list(eta = 10, gam = 10^-4, 
  M = 5, beta.min = 10^-15, beta.max = 10^100) ){

  ## set parameters
  sig=parm.outer$sig; tau=parm.outer$tau
  eta=parm.inner$eta; gam=parm.inner$gam; M=parm.inner$M; beta.min=parm.inner$beta.min; beta.max=parm.inner$beta.max
  
  
  if(scale==TRUE){ X = scale(dat,center=TRUE,scale=scale)/sqrt(nrow(dat)-1) } else { # X^T X is an empirical correlation matrix, that is normalized version of empirical covariance matrix
    X = scale(dat,center=TRUE,scale=scale)/sqrt(nrow(dat)) } # X^T X is an empirical covariance matrix
  
  p = ncol(X)
  
  ## initial feasible point 
  if( !is.null(track) ){ print("Finding a feasible point.") }
  feasinit = feasP(X,num=k,eps=.Machine$double.eps,track=NULL)
  V_feas = t(feasinit$V)

  val_perturb = 10^-10 ## a tiny value to perturb $\Sigma_n$ for numerical stability
  r = 2 * sqrt(k/val_perturb)
  
  V = V_feas
  
  rho = 1 
  Lamb = matrix(0,k,k)
  Ups = max( aug_dat(V,lambda,Lamb,rho,X,r=r,val_perturb),
             aug_dat(V=V_feas,lambda,Lamb=matrix(0,k,k),rho=0,X,r=r,val_perturb)
  )
  
  ## outer iterations
  for(t_iter in 1:max.iter.outer){
    
    if(t_iter==1) {
      L = Ups 
    } else { 
      L = aug_dat(V,lambda,Lamb,rho,X,r=r,val_perturb)
      if(L > Ups) { 
        V = V_feas 
        L = aug_dat(V,lambda,Lamb,rho,X,r=r,val_perturb)
      } 
    }
    
    
    ## inner iterations
    for(u_iter in 1:max.iter.inner){ 
      
      if(u_iter==1) {
        beta = 1
        dL = d_aug1_dat(V,Lamb,rho,X,val_perturb)
      } else {
        beta = sum(dif_V*(dL$dV-dL_pre$dV)) / dif
        beta = min(beta,beta.max)
        beta = max(beta,beta.min)
      }
      while(1){
        ## solve sub-problem
        Zv = V - 1/beta * dL$dV
        Vup = V_sub(Zv,lambda,beta,r=r)
        Lup = aug_dat(Vup,lambda,Lamb,rho,X,r=r,val_perturb)
        dif_V = Vup-V
        dif = sum( dif_V^2 )
        if( Lup <= max(L[1:M],na.rm=T) - gam/2*dif ){ break }
        beta = eta * beta
      }
      
      prec1 = sign(Vup) * lambda # $ell_1$-norm is convex. It is okay to use subdifferential instead of limiting subdifferential.
      
      dL_pre = dL
      dL = d_aug1_dat(Vup,Lamb,rho,X,val_perturb)
      
      prec2 = prec1 + dL$dV 
      prec2_zero = sapply( abs(prec2[which(Vup==0)])-lambda, function(x) max(0,x) )
      prec2_notzero = prec2[which(Vup!=0)]
      if(length(which(Vup==0))>1) { prec2 = c(prec2_zero,prec2_notzero) } else { prec2 = prec2_notzero }
      
      if( sqrt(sum(prec2^2))/((1+abs(Lup))*sqrt(p*k)) <= max(exp(-t_iter-1),eps.sub) ) {      
        V = Vup
        break 
      }

      V = Vup
      L = c(Lup,L)
    }
    
    if(u_iter == max.iter.inner) { print(paste0("Inner algorithm did not converge at ",t_iter," outer iteration"))  }
    
    prec2 = prec1 + d_aug1_dat(V,Lamb,rho=0,X,val_perturb)$dV
    prec2_zero = sapply( abs(prec2[which(Vup==0)])-lambda, function(x) max(0,x) )
    prec2_notzero = prec2[which(Vup!=0)]
    if(length(which(Vup==0))>1) { prec2 = c(prec2_zero,prec2_notzero) } else { prec2 = prec2_notzero }
    
    fV = aug_dat(V,lambda,Lamb=matrix(0,k,k),rho=0,X,r=r,val_perturb)
    
    norm_const = max(abs( V%*% t(X) %*% X %*%t(V) + val_perturb* V%*% t(V) - diag(k)))
    if( sqrt(sum(prec2^2)) / ((1+abs(fV))*sqrt(p*k)) <= eps1 && norm_const <= eps2 ) { break }
    
    Lamb_pre = Lamb
    rho_pre = rho
    
    Lamb = Lamb + rho * ( V%*% t(X) %*% X %*%t(V) + val_perturb* V%*% t(V) - diag(k) )
    rho = max( tau*rho, sum(Lamb^2)^((1+sig)/2) )
    
    if( !is.null(track) && t_iter%%track==1 ) { print(paste0("Iteration: ", t_iter, ", Objective: ",fV,
                ", Constraint: ", norm_const
        )) }
    
  }
  
  if(t_iter == max.iter.outer) { print("Outer algorithm did not converge.") }

  ## orgV
  orgV = V
  
  ## normalize V
  V = V/sqrt(rowSums(V^2))
  
  ## adjvar total variance
  svdobj = svd(X)
  totalvariance = sum((svdobj$d)^2)
  u = X %*% t(V)
  R = qr.R(qr(u))
  pev = diag(R^2)/totalvariance
  
  outlist = list(loadings = V, V=orgV, lambda = lambda, k=k,
                 p.nz = 100*length(which(V!=0))/length(V),
                 pc.cor = diag(1/sqrt(diag(V %*% t(X)%*%X %*% t(V)))) %*% (V %*% t(X)%*%X %*% t(V)) %*% diag(1/sqrt(diag(V %*% t(X)%*%X %*% t(V)))) ,
                 p.ev = 100*sum(pev)
  )
  
  class(outlist) = "euspca"
  outlist
  
} #end function








## Functions 

#Phi: Part of Lagrangian function
aug1_dat <- function(V,Lamb,rho,X,val_perturb){
  temp = ( V %*% t(X) ) %*% ( X %*% t(V) )  +  val_perturb * ( V %*% t(V) ) ## V (\Sigma+\epsI) V^T slightly perturbed; 
  diag(temp) = diag(temp) - 1 ## V (\Sigma+\eps) V^T - I
  return( - sum( (( V %*% t(X) ) %*% X + val_perturb*V)^2 ) ## -tr( {V(\Sigma+\eps)}^T V(\Sigma+\epsI) )
          + sum( Lamb * temp ) 
          + rho / 2 * sum( temp^2 )  )
}


#L: Lagrangian function
aug_dat <- function(V,lambda,Lamb,rho,X,r,val_perturb){
  if( max(abs(V)) > r ){ chi = .Machine$double.xmax # if infinity
    } else { chi=0 }
  return( 
    aug1_dat(V,Lamb,rho,X,val_perturb) + lambda * sum( abs(V) ) + chi
  ) 
}

#deri Phi: derivative of Phi
# Note that Lamb is symmetric: t(Lamb) + Lamb = 2*Lamb
d_aug1_dat <- function(V,Lamb,rho,X,val_perturb){
  temp = -2*V%*%t(X)%*%X -2*val_perturb*V + 2*Lamb%*%V - 2*rho*V + 2*rho*V%*%t(X)%*%X%*%t(V)%*%V + 2*val_perturb*rho*V%*%t(V)%*%V  
  return(
    list( 
      dV = temp %*% t(X) %*% X + val_perturb*temp
    )
  )
}


#Gsoft for single element
sol <- function(z,lambda,beta,r){
  if( abs(z) > r+lambda/beta ){ v = r*sign(z) } else { v = sign(z) * max( c(abs(z) - lambda/beta,0) ) }
  return(v)
}


#Gsoft for a matrix
V_sub <- function(Zv,lambda,beta,r){
  return( apply(Zv, 1:2, sol, beta=beta, lambda=lambda, r=r) )
}




## Derive the "first" eigen-value and the corresponding eigen-vector.
power1 <- function(X,eps=.Machine$double.eps,track=NULL){
  
  u = rep(1,ncol(X))/sqrt(ncol(X))
  u = t(X) %*% (X %*% u)
  i=0
  
  while(1){
    
    x = u / sqrt(sum(u^2))
    u = t(X) %*% (X %*% x)
    lam = t(x) %*% u
    
    temp = sum( (u-c(lam)*x)^2 )
    if( temp < eps) {break}
    if( !is.null(track) && i%%track==0 ) { print(temp) }
    
    i=i+1
  }
  
  return(list(value=c(lam),vector=x)) 
}

## Derive the "second" eigen-value and the corresponding eigen-vector given the first eigen-value and -vector.
power2 <- function(X,lam1,v1,eps=.Machine$double.eps,track=NULL){
  
  u = rep(1,ncol(X))/sqrt(ncol(X))
  u = t(X) %*% (X %*% u) - lam1*sum(v1*u)*v1
  i = 0
  
  while(1){
    
    x = u / sqrt(sum(u^2))
    u = t(X) %*% (X %*% x) - lam1*sum(v1*x)*v1
    lam = t(x) %*% u
    
    temp = sum( (u-c(lam)*x)^2 )
    if( temp < eps) {break}
    if( !is.null(track) && i%%track==0 ) { print(temp) }
    
    i = i+1
  }
  
  return(list(value=c(lam),vector=x)) 
}

## Derive the k-th eigen-value and the corresponding eigen-vector given 1st to (k-1)th eigen-values and -vectors.
power3 <- function(X,Lam,V,eps=.Machine$double.eps,track=NULL){ 
  
  k = dim(V)[2] 
  u = rep(1,ncol(X))/sqrt(ncol(X))
  
  u = t(X) %*% (X %*% u) 
  for(i in 1:k){
    u = u - Lam[i]*sum(V[,i]*u)*V[,i] # O(k^2)
  }
  
  i = 0
  
  while(1){
    
    x = u / sqrt(sum(u^2))
    u = t(X) %*% (X %*% x)
    for(i in 1:k){
      u = u - Lam[i]*sum(V[,i]*x)*V[,i]
    }
    lam = t(x) %*% u
    
    temp = sum( (u-c(lam)*x)^2 )
    if( temp < eps) {break}
    if( !is.null(track) && i%%track==0 ) { print(temp) }
    
    i = i+1
  }
  
  return(list(value=c(lam),vector=x)) 
}


feasP = function(X,num,eps=.Machine$double.eps,track=NULL){
  pow1 = power1(X,eps=eps,track=track)
  if(num==1){ return(list(V=cbind(1/sqrt(pow1$value)*pow1$vector),Lam=pow1$value)) }
  
  if(num>1){
    pow2 = power2(X,lam1=pow1$value,v1=pow1$vector,eps=eps,track=track) 
    
    Lam = c(pow1$value, pow2$value)
    V = cbind(pow1$vector,pow2$vector)
    V_feas = cbind(1/sqrt(pow1$value)*pow1$vector,1/sqrt(pow2$value)*pow2$vector) }
  
  if(num==2){ return(list(V=V_feas,Lam=Lam,eigenV=V)) }
  
  if(num>2){
    for(i in 3:num){
      pow3 = power3(X,Lam,V,eps=eps,track=track)
      Lam = c(Lam,pow3$value)
      V = cbind(V,pow3$vector)
      V_feas = cbind(V_feas,1/sqrt(pow3$value)*pow3$vector) }
    return(list(V=V_feas,Lam=Lam,eigenV=V))
  }
  
}

