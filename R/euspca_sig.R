## When n>p or Sigma is used

euspca_sig = function(Sig, k, scale, lambda,
  eps1 = 1e-4, eps2 = 1e-5, eps.sub = 1e-4,
  max.iter.outer = 1000, max.iter.inner = 2000, track=5,
  parm.outer = list(sig = 0.2, tau = 1.1), 
  parm.inner = list(eta = 10, gam = 1e-4, 
                     M = 5, beta.min = 1e-15, beta.max = 1e+100) ){
  
  ## set parameters
  sig=parm.outer$sig; tau=parm.outer$tau
  eta=parm.inner$eta; gam=parm.inner$gam; M=parm.inner$M; beta.min=parm.inner$beta.min; beta.max=parm.inner$beta.max
  
  val_perturb = 10^-10 ## a tiny value to perturb $\Sigma_n$ for numerical stability
  
  diag(Sig) = diag(Sig) + val_perturb
  
  if(scale==TRUE){ Sig = diag(diag(Sig)^{-1/2}) %*% Sig %*% diag(diag(Sig)^{-1/2}) }# X^T X is an empirical correlation matrix, that is normalized version of empirical covariance matrix

  p = ncol(Sig)
  
  ## initial feasible point 
  if( !is.null(track) ){ print("Finding a feasible point.") }
  
  dec_Sig = eigen(Sig) # Note Sig == vec_Sig %*% diag(val_Sig) %*% t(vec_Sig)
  val_Sig = dec_Sig$values
  vec_Sig = dec_Sig$vectors 
  r = 2*sqrt(k/val_Sig[length(val_Sig)])
  
  Sig_hlf = vec_Sig %*% diag(sqrt(val_Sig)) %*% t(vec_Sig) ## Sigma^{1/2}
  Sig_twi = vec_Sig %*% diag(val_Sig^2) %*% t(vec_Sig) ## Sigma^2
  Sig_nhlf = vec_Sig %*% diag(1/sqrt(val_Sig)) %*% t(vec_Sig) ## Sigma^{-1/2}
  
  W_feas = t(matrix(vec_Sig[,1:k],p,k))  
  V_feas = W_feas %*% Sig_nhlf

  V = V_feas
  
  rho = 1 
  Lamb = matrix(0,k,k)
  Ups = - sum( diag( t(V_feas) %*% V_feas %*% Sig_twi ) ) + lambda * sum( (abs(V) ) )
  
  
  ## outer iterations
  for(t_iter in 1:max.iter.outer){
    
    if(t_iter==1) {
      L = Ups 
    } else { 
      L = aug_sig(V,lambda,Lamb,rho,Sig,Sig_twi,k,r)
      if(L > Ups) { 
        V = V_feas 
        L = aug_sig(V,lambda,Lamb,rho,Sig,Sig_twi,k,r)
      } 
    }
    
    
    ## inner iterations
    for(u_iter in 1:max.iter.inner){ 
      
      if(u_iter==1) {
        beta = 1
        dL = d_aug1_sig(V,Lamb,rho,Sig,Sig_twi,k)
      } else {
        beta = abs(sum(dif_V * (dL$dV-dL_pre$dV))) / dif
        beta = min(beta,beta.max)
        beta = max(beta,beta.min)
      }
      while(1){
        ## solve sub-problem
        Zv = V - 1/beta * dL$dV
        Vup = V_sub(Zv,lambda,beta,r)
        Lup = aug_sig(Vup,lambda,Lamb,rho,Sig,Sig_twi,k,r)
        dif_V = Vup-V
        dif = sum( dif_V^2 )
        if( Lup <= max(L[1:M],na.rm=T) - gam/2*dif ){ break }
        beta = eta * beta
      }
      
      prec1 = sign(Vup) * lambda # $ell_1$-norm is convex. It is okay to use subdifferential instead of limiting subdifferential.
      
      dL_pre = dL
      dL = d_aug1_sig(Vup,Lamb,rho,Sig,Sig_twi,k)
      
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
    
    prec2 = prec1 + d_aug1_sig(V,Lamb,rho=0,Sig,Sig_twi,k)$dV
    prec2_zero = sapply( abs(prec2[which(Vup==0)])-lambda, function(x) max(0,x) )
    prec2_notzero = prec2[which(Vup!=0)]
    if(length(which(Vup==0))>1) { prec2 = c(prec2_zero,prec2_notzero) } else { prec2 = prec2_notzero }
    
    fV = aug_sig(V,lambda,Lamb=matrix(0,k,k),rho=0,Sig,Sig_twi,k,r)
    
    norm_const = max( abs( V %*% Sig %*% t(V) - diag(k) ) )
    if( sqrt(sum(prec2^2)) / ((1+abs(fV))*sqrt(p*k)) <= eps1 && norm_const <= eps2 ) { break }
    
    Lamb_pre = Lamb
    rho_pre = rho
    
    Lamb = Lamb + rho * ( V%*%Sig%*%t(V)-diag(k) )
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
  diag(Sig) = diag(Sig) - val_perturb ## not perturbed Sigma
  Sig_hlf = vec_Sig %*% diag(sqrt(val_Sig - val_perturb)) %*% t(vec_Sig) ## not perturbed Sigma^{-1/2}
  svdobj = svd(Sig_hlf)
  totalvariance = sum((svdobj$d)^2)
  u = Sig_hlf %*% t(V)
  R = qr.R(qr(u))
  pev = diag(R^2)/totalvariance
  
  outlist = list(loadings = V, V=orgV, lambda = lambda, k=k,
                p.nz = 100*length(which(V!=0))/length(V),
                pc.cor = diag(1/sqrt(diag(V %*% Sig %*% t(V)))) %*% (V %*% Sig %*% t(V)) %*% diag(1/sqrt(diag(V %*% Sig %*% t(V)))) ,
                p.ev = 100*sum(pev) # Percentage of explained variance
  )
  
  class(outlist) = "euspca"
  outlist
  
} #end function








## Functions 

#Phi: Part of Lagrangian function
aug1_sig <- function(V,Lamb,rho,Sig,Sig_twi,k){
  temp1 = V %*% Sig_twi %*% t(V)
  temp2 = V %*% Sig %*% t(V) - diag(k)
  return( - sum(diag( temp1 ))
          + sum(diag( t(Lamb) %*% temp2 ))
          + rho / 2 * sum( temp2^2 )
  )
}


#L: Lagrangian function
aug_sig <- function(V,lambda,Lamb,rho,Sig,Sig_twi,k,r){
  if( max(abs(V)) > r ){ chi = 100^100 } else { chi=0 }
  return( 
    aug1_sig(V,Lamb,rho,Sig,Sig_twi,k) + lambda * sum( (abs(V) ) ) + chi
  )
}

#deri Phi: derivative of Phi 
d_aug1_sig <- function(V,Lamb,rho,Sig,Sig_twi,k){
  return(
    list( 
      dV = -2 * V %*% Sig_twi + ( t(Lamb) + Lamb - 2 * rho * diag(k) + 2 * rho * V %*% Sig %*% t(V) ) %*% V %*% Sig
    )
  )
}



#Gsoft for single element
sol <- function(z,lambda,beta,r){ #note apply(matrix, 1:2, sol_log, beta=, lambda=)
  if( abs(z) > r+lambda/beta ){ v = r*sign(z) } else { v = sign(z) * max( c(abs(z) - lambda/beta,0) ) }
  return(v)
}


#Gsoft for a matrix
V_sub <- function(Zv,lambda,beta,r){
  return( apply(Zv, 1:2, sol, beta=beta, lambda=lambda, r=r) )
}

