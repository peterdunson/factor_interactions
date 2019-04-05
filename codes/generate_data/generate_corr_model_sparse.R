generate_corr_model_sparse = function(p = 10,n = 100,k_true = 5){
   
   
   # generate correlated matrix of X
   L = matrix(rnorm(p*k_true,0,1),p,k_true)
   X = matrix(rnorm(n*k_true,0,1),n,k_true)%*%t(L)
   X_test = matrix(rnorm(n*k_true,0,1),n,k_true)%*%t(L)
   X_big = scale(rbind(X,X_test))8
   X = X_big[1:n,]
   X_test = X_big[(n+1):(2*n),]
   
   # true coefficients
   beta_true = numeric(p)
   beta_true[c(1,3,4)] = c(1,-1,1)
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
