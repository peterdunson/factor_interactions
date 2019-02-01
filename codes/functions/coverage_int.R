coverage_int = function(beta,Omega,inf = 0.025,sup = 0.975){
  
  p = min(ncol(beta),nrow(beta))
  beta_c = apply(beta, 2, function(x) quantile(x,probs = c(inf,sup)))
  Omega_c = apply(Omega, c(2,3), function(x) quantile(x,probs = c(inf,sup)))
  
  beta_cov = rep(0,p)
  for (j in 1:p){
    if(sign(beta_c[1,j])==sign(beta_c[2,j])){
      beta_cov[j] = 1
    }
  }
  
  Omega_cov = matrix(0,p,p)
  for (j in 1:p){
    for(k in 1:p){
      if(sign(Omega_c[1,j,k])==sign(Omega_c[2,j,k])){
        Omega_cov[j,k] = 1
      }
    }
  }
  
  #if beta[j] = 1, it means that zero is not contained in the 
  # credible interval
  return(list(beta_cov,Omega_cov))
}
