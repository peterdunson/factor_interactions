rate_recovery_interactions = function(gibbs,gibbs2 = gibbs,alpha = 0.05,
                                      Omega_true,...){
   
   #first gibbs
   coverage = coverage_int(gibbs$beta_bayes,gibbs$Omega_bayes,
                           inf = alpha/2, sup = 1-alpha/2)
   Omega_cov = coverage[[2]]
   ind_Omega_bayes =  which(Omega_cov[lower.tri(Omega_cov,diag = T)] != 0)
   ind_Omega_bayes_zero =  which(Omega_cov[lower.tri(Omega_cov,diag=T)] == 0)
   #second gibbs
   coverage2 = coverage_int(gibbs2$beta_bayes,gibbs2$Omega_bayes,
                            inf = alpha/2, sup = 1-alpha/2)
   Omega_cov2 = coverage2[[2]]
   ind_Omega_bayes2 =  which(Omega_cov2[lower.tri(Omega_cov2,diag = T)] != 0)
   ind_Omega_bayes_zero2 =  which(Omega_cov2[lower.tri(Omega_cov2,diag=T)] == 0)
   
   #indeces of nonzero and zero values
   ind_Omega = which(Omega_true[lower.tri(Omega_true,diag = T)] != 0)
   ind_Omega_zero = which(Omega_true[lower.tri(Omega_true,diag = T)] == 0)
   
   # true positives and true negatives
   TP = sum(complete.cases(match(ind_Omega,ind_Omega_bayes)))/length(ind_Omega)
   TP2 = sum(complete.cases(match(ind_Omega,ind_Omega_bayes2)))/length(ind_Omega)
   TN = sum(complete.cases(match(ind_Omega_zero,ind_Omega_bayes_zero)))/length(ind_Omega_zero)
   TN2 = sum(complete.cases(match(ind_Omega_zero,ind_Omega_bayes_zero2)))/length(ind_Omega_zero)
   
   Other_Omega = list(...)
   t = length(Other_Omega)
   TP_vector = c(TP,TP2)
   TN_vector = c(TN,TN2)
   
   if(t > 0){
      for(k in 1:t){
         Omega_curr = as.matrix(Other_Omega[[k]])
         ind_Omega_curr =  which(Omega_curr[lower.tri(Omega_curr,diag = T)] != 0)
         ind_Omega_curr_zero =  which(Omega_curr[lower.tri(Omega_curr,diag=T)] == 0)
         
         TP_curr = sum(complete.cases(match(ind_Omega,ind_Omega_curr)))/length(ind_Omega)
         TN_curr = sum(complete.cases(match(ind_Omega_zero,ind_Omega_curr_zero)))/length(ind_Omega_zero)
         
         TP_vector = c(TP_vector,TP_curr)
         TN_vector = c(TN_vector,TN_curr)
      }
   }
   return(list(TP = TP_vector,TN = TN_vector))
}

rate_recovery_maineff = function(gibbs,gibbs2 = gibbs,alpha = 0.05,
                                         beta_true,...){
   
   #first gibbs
   coverage = coverage_int(gibbs$beta_bayes,gibbs$Omega_bayes,
                           inf = alpha/2, sup = 1-alpha/2)
   beta_cov = coverage[[1]]
   ind_beta_bayes =  which(beta_cov != 0)
   ind_beta_bayes_zero =  which(beta_cov == 0)
   #second gibbs
   coverage2 = coverage_int(gibbs$beta_bayes,gibbs$Omega_bayes,
                           inf = alpha/2, sup = 1-alpha/2)
   beta_cov2 = coverage2[[1]]
   ind_beta_bayes2 =  which(beta_cov2 != 0)
   ind_beta_bayes_zero2 =  which(beta_cov2 == 0)
   
   #indeces of nonzero and zero values
   ind_beta = which(beta_true != 0)
   ind_beta_zero = which(beta_true == 0)
   
   # true positives and true negatives
   TP = sum(complete.cases(match(ind_beta,ind_beta_bayes)))/length(ind_beta)
   TP2 = sum(complete.cases(match(ind_beta,ind_beta_bayes2)))/length(ind_beta)
   TN = sum(complete.cases(match(ind_beta_zero,ind_beta_bayes_zero)))/length(ind_beta_zero)
   TN2 = sum(complete.cases(match(ind_beta_zero,ind_beta_bayes_zero2)))/length(ind_beta_zero)
   
   Other_beta= list(...)
   t = length(Other_beta)
   TP_vector = c(TP,TP2)
   TN_vector = c(TN,TN2)
   
   if(t > 0){
      for(k in 1:t){
         beta_curr = as.matrix(Other_beta[[k]])
         ind_beta_curr =  which(beta_curr != 0)
         ind_beta_curr_zero =  which(beta_curr == 0)
         
         TP_curr = sum(complete.cases(match(ind_beta,ind_beta_curr)))/length(ind_beta)
         TN_curr = sum(complete.cases(match(ind_beta_zero,ind_beta_curr_zero)))/length(ind_beta_zero)
         
         TP_vector = c(TP_vector,TP_curr)
         TN_vector = c(TN_vector,TN_curr)
      }  
   }
   return(list(TP = TP_vector,TN = TN_vector))
   
   
}