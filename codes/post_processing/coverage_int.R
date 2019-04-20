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
  
  p = dim(Omega)[2]; q = dim(Omega)[3]
  Omega_cov = matrix(0,p,q)
  for (j in 1:p){
    for(k in 1:q){
      if(sign(Omega_c[1,j,k])==sign(Omega_c[2,j,k])){
        Omega_cov[j,k] = 1
      }
    }
  }
  
  #if beta[j] = 1, it means that zero is not contained in the 
  # credible interval
  return(list(beta_cov,Omega_cov))
}


coverage_y = function(y_pred,y,alpha = 0.05){
   
   inf = alpha/2
   sup = 1 - alpha/2
   n = length(y)
   y_c = apply(y_pred, 2, function(x) quantile(x,probs = c(inf,sup)))

   y_cov = rep(0,n)
   for (i in 1:n){
      if(y[i] < y_c[2,i] & y[i] > y_c[1,i]){
         y_cov[i] = 1
      }
   }
   
   #if beta[j] = 1, it means that zero is not contained in the 
   # credible interval
   return(y_cov)
}
