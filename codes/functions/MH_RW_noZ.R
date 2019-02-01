MH_RW_noZ = function(phi,sigmasq_y,Lambda,ps,k,X,y,eta,acp,s,burn,Psi,delta_rw){
   #update of MH s.t. the posterior of eta_i is normal and the Covariance matrix is
   #the same for every eta_i
   
   n = nrow(X)
   p = ncol(X)
   
   a = phi%*%t(phi)/sigmasq_y + t(Lambda)%*%diag(ps)%*%Lambda + diag(k)
   
   for (i in 1:n){
      
      eta_star = as.numeric(rmvnorm(1,eta[i,],diag(k)*delta_rw))
      
      # MH update
      # they give the same result with the function or with the expression below
      # I use the expression since I compute a just one for the loop
      logr = t(eta_star)%*%(a - 2*Psi*y[i]/sigmasq_y)%*%eta_star -
         2*t(eta_star)%*%(t(Lambda)%*%diag(ps)%*%X[i,] + phi*y[i]/sigmasq_y) +
         2*t(eta_star)%*%phi*(eta_star%*%Psi%*%eta_star)/sigmasq_y + 
         (1/sigmasq_y)*(eta_star%*%Psi%*%eta_star)^2 - 
         (t(eta[i,])%*%(a - 2*Psi*y[i]/sigmasq_y)%*%eta[i,] -
         2*t(eta[i,])%*%(t(Lambda)%*%diag(ps)%*%X[i,] + phi*y[i]/sigmasq_y) +
         2*t(eta[i,])%*%phi*(eta[i,]%*%Psi%*%eta[i,])/sigmasq_y +
         (1/sigmasq_y)*(eta[i,]%*%Psi%*%eta[i,])^2)
      logr = logr*(-0.5)

      #logr = log_posterior_eta(X[i,],y[i],eta_star,phi,Lambda,ps,k,sigmasq_y,Psi) - 
      #   log_posterior_eta(X[i,],y[i],eta[i,],phi,Lambda,ps,k,sigmasq_y,Psi) 
      
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
