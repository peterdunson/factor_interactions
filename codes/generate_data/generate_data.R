generate_data = function(p = 10,n = 100,k_true = 5,
                        ratio_Om = 0.2,ratio_beta = 0.2,
                        type = "independent", sigmasq = 1){
   
   if(type == "independent"){
      X_big = scale(mvtnorm::rmvnorm(n*2,sigma=diag(p)))
      X = X_big[1:n,]
      X_test = X_big[(n+1):(2*n),]
   }else if(type == "wishart"){
      W = stats::rWishart(1,p+1,S = diag(p))
      X_big = scale(mvtnorm::rmvnorm(n*2,sigma=W[,,1]))
      X = X_big[1:n,]
      X_test = X_big[(n+1):(2*n),]
   }else if (type == "correlated"){
      # generate correlated matrix of X
      L = matrix(rnorm(p*k_true,0,1),p,k_true)
      X = matrix(rnorm(n*k_true,0,1),n,k_true)%*%t(L)
      X_test = matrix(rnorm(n*k_true,0,1),n,k_true)%*%t(L)
      X_big = scale(rbind(X,X_test))
      X = X_big[1:n,]
      X_test = X_big[(n+1):(2*n),]
   }else if (type == "power covariance"){
      W = matrix(0,p,p)
      
      for(j in 1:p){
         for(k in 1:j){
            W[j,k] = 0.8^abs(j-k)
         }
      }
      
      W = W + t(W)
      diag(W) = diag(W)/2
      
      X_big = scale(mvtnorm::rmvnorm(n*2,sigma=W[,]))
      X = X_big[1:n,]
      X_test = X_big[(n+1):(2*n),]
      
      
      
   }else{
      return("provide proper data generating process")
   }
   
   
   # true coefficients
   beta_true = numeric(p)
   nonzero_beta = floor(p*ratio_beta)
   nonzero_Om = floor((p^2)*ratio_Om*2/3)
   coeffs = c(-1,-0.8,-0.6,-0.4,1,0.8,0.6,0.4)
   beta_true[sample(1:p,nonzero_beta,replace = T)] = sample(coeffs,nonzero_beta,replace = T)
   
   Om_ind = cbind(sample(1:p,nonzero_Om,replace = T),sample(1:p,nonzero_Om,replace = T))
   Omega_true = matrix(0,p,p)
   Omega_true[Om_ind] = sample(coeffs,nonzero_Om,replace = T)
   Omega_true = Omega_true+t(Omega_true)
   
   # generate output
   y=as.vector(diag(X%*%Omega_true%*%t(X))+
                  X%*%beta_true+rnorm(n,0,sigmasq))
   
   y_test = as.vector(diag(X_test%*%Omega_true%*%t(X_test))+
                         X_test%*%beta_true+rnorm(n,0,sigmasq))
   return(list(y = y, X = X, 
               beta_true = beta_true, 
               Omega_true = Omega_true,
               y_test = y_test,
               X_test = X_test))
}
