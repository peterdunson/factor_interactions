####### Competitors Functions ######
library(hierNet)

library(RAMP)
library('glmnet')

Hiernet_fct = function(y, X, X_test = X, y_test = y, strong = T){
   p = ncol(X)
   n = nrow(X)
   fit=hierNet(X,y,lam=50, strong=strong,
               trace = 0)
   
   #extract main effects and interactions
   beta = fit$bp-fit$bn
   Omega = fit$th
   
   #prediction
   y_hat = predict(fit,X)
   err = y-y_hat
   y_pred = predict(fit,X_test)
   err_pred = y_test-y_pred
   return(list(beta = beta,Omega = Omega,
               err = err,err_pred = err_pred))
}


   




RAMP_fct = function(y, X, X_test = X, y_test = y, hier = "Strong",
                    max.iter = 100){
   
   #estimate model
   fit = RAMP(X,y,max.iter = max.iter,
              hier = hier)
   
   p = ncol(X)
   n = nrow(X)
   
   #extract coef
   #beta
   ind_main_eff = fit$mainInd
   beta = numeric(p)
   beta[ind_main_eff] = fit$beta.m
   
   #Omega
   imp_int = fit$interInd
   #int_list = sort(fit$interInd.list[[max.iter]])
   
   int_list = character(p + p*(p-1)/2)
   count = 1
   for(i in 1:p){
      for(j in i:p){
         int_list[count] = paste("X",i,"X",j,sep="")
         count = count+1
      }
   }
   
   ind = match(imp_int,int_list)
   int = numeric(p+p*(p-1)/2)
   int[ind] = fit$beta.i
   Omega = matrix(0,p,p)
   Omega[lower.tri(Omega,diag = T)] = int/2
   Omega = Omega + t(Omega)
   
   #prediction
   y_hat = predict(fit,X)
   err = y-y_hat
   y_pred = predict(fit, X_test)
   err_pred = y_test-y_pred
   
   return(list(beta = beta,Omega = Omega,
               err = err,err_pred = err_pred))
   
}


library(PIE)
PIE_fct = function(y, X, X_test = X, y_test = y){
   
   p = ncol(X)
   n = nrow(X)
   
   #estimation
   beta = as.vector(coef(cv.glmnet(X,y,nfolds = 5),
                         s="lambda.min"))[-1];  
   Omega = PIE(X,y-X%*%beta)
   
   #prediction
   err = y-X%*%beta - as.vector(diag(X%*%Omega%*%t(X)))
   err_pred = y_test-X_test%*%beta - 
      as.vector(diag(X_test%*%Omega%*%t(X_test)))
   
   return(list(beta = beta,Omega = Omega,
               err = err,err_pred = err_pred))
}