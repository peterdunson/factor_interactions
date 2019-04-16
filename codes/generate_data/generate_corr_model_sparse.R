generate_corr_model_sparse = function(p = 10,n = 100,k_true = 5,
                                      ratio_Om = 0.02,ratio_beta = 0.1){
   
   
   # generate correlated matrix of X
   L = matrix(rnorm(p*k_true,0,1),p,k_true)
   X = matrix(rnorm(n*k_true,0,1),n,k_true)%*%t(L)
   X_test = matrix(rnorm(n*k_true,0,1),n,k_true)%*%t(L)
   X_big = scale(rbind(X,X_test))
   X = X_big[1:n,]
   X_test = X_big[(n+1):(2*n),]
   
   # true coefficients
   beta_true = numeric(p)
   nonzero_beta = floor(p*ratio_beta)
   nonzero_Om = floor((p^2)*ratio_Om*2/3)
   coeffs = c(-1,-0.8,-0.6,-0.4,1,0.8,0.6,0.4)
   beta_true[sample(1:p,nonzero_beta)] = sample(coeffs,nonzero_beta)
   
   Om_ind = cbind(sample(1:p,nonzero_Om,replace = T),sample(1:p,nonzero_Om,replace = T))
   Omega_true = matrix(0,p,p)
   Omega_true[Om_ind] = sample(coeffs,nonzero_Om,replace = T)
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
