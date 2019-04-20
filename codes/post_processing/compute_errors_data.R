compute_errors_data = function(hiernet,Family,PIE,RAMP,
                          y,y_test,X_test,Z_test,
                          ...){
   
   #in sample error
   err_insample = c(mean((y-mean(y))^2),mean(hiernet$err^2),mean(Family$err^2),
                    mean(PIE$err^2),mean(RAMP$err^2))
   
   #Predictive error
   err_pred = c(mean((y_test-mean(y_test))^2),mean(hiernet$err_pred^2),mean(Family$err_pred^2),
                mean(PIE$err_pred^2),mean(RAMP$err_pred^2))
   
   name = c("base case","hiernet","Family","PIE","RAMP")
   
   gibbs_list = list(...)
   t = length(gibbs_list)
   
   if(t > 0){
      for(k in 1:t){
         name_curr = paste("gibbs",k,sep="_")
         gibbs = gibbs_list[[k]]
         
         Omega_hat = apply(gibbs$Omega_bayes, c(2,3), mean)
         Phi_hat = apply(gibbs$Omega_conf, c(2,3), mean)
         alpha_hat = mean(gibbs$alpha_bayes)
         beta_hat = apply(gibbs$beta_bayes,2,mean)
         beta_Z_hat = apply(gibbs$beta_Z,2,mean)
         
         y_hat = X%*%beta_hat + 
            as.vector(diag(X%*%Omega_hat%*%t(X))) +
            alpha_hat + Z%*%beta_Z_hat +
            as.vector(diag(X%*%Phi_hat%*%t(Z)))
         y_pred = X_test%*%beta_hat + 
            as.vector(diag(X_test%*%Omega_hat%*%t(X_test))) +
            alpha_hat + Z_test%*%beta_Z_hat +
            as.vector(diag(X_test%*%Phi_hat%*%t(Z_test)))
         
         #in sample
         err_factor = y - y_hat
         err_insample = c(err_insample,mean(err_factor^2))
         
         #prediction
         err_pred_factor = y_test - y_pred
         err_pred = c(err_pred,mean(err_pred_factor^2))
         
         
         #name
         name = c(name,name_curr)
      }
   }
   
   names(err_insample) = names(err_pred) = name
   
   return(list(err_insample = err_insample,
               err_pred = err_pred))
}
