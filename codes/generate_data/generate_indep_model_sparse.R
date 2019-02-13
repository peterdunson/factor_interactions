generate_indep_model_sparse = function(p = 10,n = 100){
   
   # generate X
   X = matrix(rnorm(n*p,0,1),n,p)
   X_test = matrix(rnorm(n*p,0,1),n,p)
   X_big = scale(rbind(X,X_test))
   X = X_big[1:n,]
   X_test = X_big[(n+1):(2*n),]
   
   # true coefficients
   beta_true = numeric(p)
   beta_true[c(1,3)] = c(1,-1)
   Omega_true = matrix(0,p,p)
   Omega_true[1,2] = 1;Omega_true[3,3] = -1;
   Omega_true = Omega_true+t(Omega_true)
   
   # generate output
   y=as.vector(diag(X%*%Omega_true%*%t(X))+
                  X%*%beta_true+rnorm(n,0,0.5))
   
   y_test = as.vector(diag(X_test%*%Omega_true%*%t(X_test))+
                         X_test%*%beta_true+rnorm(n,0,0.5))
   return(list(y = y, X = X, 
               beta_true = beta_true, 
               Omega_true = Omega_true,
               y_test = y_test,
               X_test = X_test))
}