MH_fast_new_mean2_noZ = function(phi,sigmasq_y,Lambda,ps,k,X,y,eta,acp,s,burn,Psi){
   #update of MH s.t. the posterior of eta_i is normal and the Covariance matrix is
   #the same for every eta_i
   
   n = nrow(X)
   p = ncol(X)
   
   Vn = solve(phi%*%t(phi)/sigmasq_y+t(Lambda)%*%diag(ps)%*%Lambda+diag(rep(1,k)))
   eta0 = MB_update_eta_noZ(phi,sigmasq_y,Lambda,ps,k,X,y,eta)
   
   for (i in 1:n){
      eta0_i = eta[i,]
      a0 = as.vector(Psi%*%eta0_i)
      eta0_sq = as.vector(t(eta0_i)%*%Psi%*%eta0_i)
      mun = Vn%*%(t(Lambda)%*%diag(ps)%*%X[i,] + 
                     (-2*eta0_sq + as.vector(2*y[i]/sigmasq_y-
                                                4*t(eta0_i)%*%phi/sigmasq_y))*a0 - 
                     2/sigmasq_y*eta0_sq*phi + 
                     phi*(y[i])/sigmasq_y)
      
      
      eta_star = as.numeric(rmvnorm(1,mun,Vn))
      
      
      # MH update
      logr = log_posterior_eta(X[i,],y[i],eta_star,phi,Lambda,ps,k,sigmasq_y,Psi) - 
         log_posterior_eta(X[i,],y[i],eta0_i,phi,Lambda,ps,k,sigmasq_y,Psi) -
         dmvnorm(eta_star,mun,Vn,log=T) + 
         dmvnorm(eta0_i,mun,Vn,log=T)
      
      
      
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
