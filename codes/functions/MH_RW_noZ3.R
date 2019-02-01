MH_RW_noZ3 = function(phi,phi3,sigmasq_y,Lambda,ps,k,X,y,eta,acp,s,burn,Psi,delta_rw){
   #update of MH s.t. the posterior of eta_i is normal and the Covariance matrix is
   #the same for every eta_i
   
   n = nrow(X)
   p = ncol(X)
   
   eta3 = eta^3
   
   for (i in 1:n){
      
      eta_star = as.numeric(rmvnorm(1,eta[i,],diag(k)*delta_rw))
      eta_star3 = eta_star^3
      
      # MH update
      logr = log_posterior_eta3(X[i,],y[i],eta_star,eta_star3,phi,phi3,Lambda,ps,k,sigmasq_y,Psi) - 
         log_posterior_eta3(X[i,],y[i],eta[i,],eta3[i,],phi,phi3,Lambda,ps,k,sigmasq_y,Psi) 
      
      logu = log(runif(1))
      
      if (logr > logu){
         eta[i,] = eta_star
         if(s>burn){
            acp[i] = acp[i] + 1
         }
      }
   }
   return(list(eta,acp))
}
