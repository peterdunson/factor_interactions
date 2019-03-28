MH_Langevin_noZ = function(phi,sigmasq_y,Lambda,ps,k,X,y,eta,
                           acp,s,burn,Psi,epsilon_rw){
   #update of MH s.t. the posterior of eta_i is normal and the Covariance matrix is
   #the same for every eta_i
   
   n = nrow(X)
   p = ncol(X)
   
   a = phi%*%t(phi)/sigmasq_y + t(Lambda)%*%diag(ps)%*%Lambda + diag(k)
   
   for (i in 1:n){
      
      mu_eta_i = 2*(a-2*Psi*y[i]/sigmasq_y)%*%eta[i,] -
         2*(t(Lambda)%*%diag(ps)%*%X[i,] + phi*y[i]/sigmasq_y) + 
         (2*phi%*%t(eta[i,]) + 4*as.vector(t(eta[i,])%*%phi))%*%Psi%*%eta[i,]/sigmasq_y
         4*as.vector(t(eta[i,])%*%Psi%*%eta[i,])*Psi%*%eta[i,]/sigmasq_y
      
         
      # should be \mu_eta_i / 2 ---> check
      eta_star = as.numeric(rmvnorm(1,eta[i,] - epsilon_rw*mu_eta_i/4,
                                    diag(k)*epsilon_rw))
      
      # MH update
      logr = t(eta_star)%*%(a - 2*Psi*y[i]/sigmasq_y)%*%eta_star -
         2*eta_star%*%(t(Lambda)%*%diag(ps)%*%X[i,] + phi*y[i]/sigmasq_y) +
         2*t(eta_star)%*%phi*(eta_star%*%Psi%*%eta_star)/sigmasq_y + 
         (1/sigmasq_y)*(eta_star%*%Psi%*%eta_star)^2 - 
         (t(eta[i,])%*%(a - 2*Psi*y[i]/sigmasq_y)%*%eta[i,] -
         2*eta[i,]%*%(t(Lambda)%*%diag(ps)%*%X[i,] + phi*y[i]/sigmasq_y) +
         2*t(eta[i,])%*%phi*(eta[i,]%*%Psi%*%eta[i,])/sigmasq_y +
         (1/sigmasq_y)*(eta[i,]%*%Psi%*%eta[i,])^2)
      logr = logr*(-0.5) - 
         dmvnorm(eta_star,eta[i,]-epsilon_rw*mu_eta_i/4, diag(k)*epsilon_rw, log = T) +
         dmvnorm(as.vector(epsilon_rw*mu_eta_i/4), sigma = diag(k)*epsilon_rw, log = T)
         
      
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
