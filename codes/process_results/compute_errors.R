library(psych)

compute_errors = function(hiernet,Family,PIE,RAMP,
                          y,y_test,Omega_true,beta_true,
                          ...){
   
   #in sample error
   err_insample = c(mean((y-mean(y))^2),mean(hiernet$err^2),mean(Family$err^2),
     mean(PIE$err^2),mean(RAMP$err^2))

   #Predictive error
   err_pred = c(mean((y_test-mean(y_test))^2),mean(hiernet$err_pred^2),mean(Family$err_pred^2),
     mean(PIE$err_pred^2),mean(RAMP$err_pred^2))
   
   #Frobenius norm
   Fr_hiernet = hiernet$Omega - Omega_true; Fr_Family = Family$Omega - Omega_true
   Fr_PIE = as.matrix(PIE$Omega) - Omega_true; Fr_RAMP = RAMP$Omega - Omega_true
   FR_norm = c(sqrt(tr(t(Omega_true%*%Omega_true))),sqrt(tr(t(Fr_hiernet)%*%Fr_hiernet)),
     sqrt(tr(t(Fr_Family)%*%Fr_Family)),sqrt(tr(t(Fr_PIE)%*%Fr_PIE)),
     sqrt(tr(t(Fr_RAMP)%*%Fr_RAMP)))
   
   #Betas
   beta_MSE = c(mean(beta_true^2),mean((hiernet$beta-beta_true)^2),
     mean((Family$beta-beta_true)^2),mean((PIE$beta-beta_true)^2),
     mean((RAMP$beta-beta_true)^2))
   
   name = c("base case","hiernet","Family","PIE","RAMP")
   
   gibbs_list = list(...)
   t = length(gibbs_list)

   if(t > 0){
      for(k in 1:t){
         name_curr = paste("gibbs",k,sep="_")
         gibbs = gibbs_list[[k]]
         alpha_bayes_hat = mean(gibbs$alpha_bayes)
         betabayes_hat = apply(gibbs$beta_bayes,2,mean)
         Omegabayes_hat = apply(gibbs$Omega_bayes,c(2,3),mean)
         
         
         #in sample
         err_factor = y - X%*%betabayes_hat - 
            as.vector(diag(X%*%Omegabayes_hat%*%t(X))) - alpha_bayes_hat
         err_insample = c(err_insample,mean(err_factor^2))
         
         #prediction
         err_pred_factor = y_test - X_test%*%betabayes_hat - 
            as.vector(diag(X_test%*%Omegabayes_hat%*%t(X_test))) - alpha_bayes_hat
         err_pred = c(err_pred,mean(err_pred_factor^2))
         
         #Frob
         Fr_factor = Omegabayes_hat - Omega_true
         FR_norm = c(FR_norm,sqrt(tr(t(Fr_factor)%*%Fr_factor)))
         
         #beta MSE
         beta_MSE = c(beta_MSE,mean((betabayes_hat-beta_true)^2))
         
         #name
         name = c(name,name_curr)
      }
   }
   
   names(err_insample) = names(err_pred) = names(FR_norm) = names(beta_MSE) = name
   
   return(list(err_insample = err_insample,
               err_pred = err_pred,
               FR_norm = FR_norm,
               beta_MSE = beta_MSE))
}
