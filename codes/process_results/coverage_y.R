coverage_y = function(y_test,X_test,gibbs,alpha = 0.05){
   
   alphabayes = gibbs$alpha_bayes
   beta = gibbs$beta_bayes
   Omega = gibbs$Omega_bayes
   sigmasq = gibbs$sigmasq_st
   
   S = length(alphabayes)
   y_pred = matrix(0,nrow = S, ncol = length(y_test))
   
   for(s in 1:S){
      mu = X_test%*%beta[s,] + as.vector(diag(X_test%*%Omega[s,,]%*%t(X_test))) + alphabayes[s]
      for(i in 1:length(y_test)){
         
         y_pred[s,i] = rnorm(1, mean = mu[i], sd = sqrt(sigmasq[s]))
         
      }
   }
   
   inf = alpha/2
   sup = 1 - alpha/2
   n = length(y_test)
   y_c = apply(y_pred, 2, function(x) quantile(x,probs = c(inf,sup)))
   
   y_cov = rep(0,n)
   for (i in 1:n){
      if(y_test[i] < y_c[2,i] & y_test[i] > y_c[1,i]){
         y_cov[i] = 1
      }
   }
   
   bias = y_test - apply(y_pred,2,mean)
   
   cov_y = y_cov %>% mean()
   bias_y = bias %>% mean()
   ls_ret = list(coverage = cov_y,
                 bias = bias_y)
   return(ls_ret)
   
}




