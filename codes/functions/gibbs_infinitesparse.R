gibbs_infinitesparse = function(y, X ,nrun, burn, thin = 1, 
                    delta_rw = 0.002, epsilon_rw = 0.5,
                    a = 1/2, k = NULL){
   n = nrow(X)
   p = ncol(X)                    # collect data attributes
   
   if(is.null(k)) k = floor(log(p)*3)
   if(length(y) != n) stop("Mismatching input lengths")
   
   VX = apply(X, 2, var)          # prepare data matrix
   X = as.matrix(scale(X))
   
   sp = floor((nrun - burn)/thin)
   
   b0 = 1            # factorization hyperparameters
   b1 = 0.0005
   epsilon = 1e-3
   prop = 1.00
   as = 1
   bs = 0.3
   df = 150
   ad1 = 3
   bd1 = 1
   ad2 = 4
   bd2 = 1
   adf = 1
   bdf = 1
   
   Omega_bayes = array(0,c(sp,p,p))  # sample storage memory allocation
   tau_st = array(0,c(sp,p))
   sigmasq_st = numeric(sp)
   alpha_bayes = numeric(sp)
   beta_bayes = matrix(0,sp,p)
   acp = numeric(n)
   k_st = numeric(sp)
   
   sigmasq_y = 1                     # initial values
   phi = numeric(k)
   ps = rgamma(p,as,bs)
   Sigma = diag(1/ps)
   Lambda = matrix(0,p,k)
   eta = matrix(rnorm(n*k),n,k)
   meta = matrix(0,n,k)
   veta = diag(rep(1,k))
   Psi = matrix(0,k,k)
   
   psijh = matrix(rgamma(p*k,df/2,df/2),p,k)
   Plam = t(t(psijh) * tauh)
   delta = c(rgamma(1,ad1,1/bd1),rgamma(k-1,ad2,1/bd2))
   tauh = cumprod(delta)
   Plam = psijh*matrix(rep(tauh,p),p,k,byrow=T)
   
   t = t0 = Sys.time()               # begin sample timing
   count = 1
   
   for(i in 1:nrun){
      
      # --- Update eta --- #
      Lambda.T = t(Lambda)
      aMH = phi%*%t(phi)/sigmasq_y + Lambda.T%*%diag(ps)%*%Lambda + diag(k)
      
      for (h in 1:n){                # Metropolis hastings step 
         
         eta_star = bayesSurv::rMVNorm(1,eta[h,],diag(k)*delta_rw)
         eta_star.T = t(eta_star)     # avoid repeated transpose calls
         eta.T = t(eta[h,])
         
         logr = eta_star.T%*%(aMH - 2*Psi*y[h]/sigmasq_y)%*%eta_star -
            2*eta_star.T%*%(Lambda.T%*%diag(ps)%*%X[h,] + phi*y[h]/sigmasq_y) +
            2*eta_star.T%*%phi*(eta_star%*%Psi%*%eta_star)/sigmasq_y + 
            (1/sigmasq_y)*(eta_star%*%Psi%*%eta_star)^2 - 
            (eta.T%*%(aMH - 2*Psi*y[h]/sigmasq_y)%*%eta[h,] -
                2*eta.T%*%(Lambda.T%*%diag(ps)%*%X[h,] + phi*y[h]/sigmasq_y) +
                2*eta.T%*%phi*(eta[h,]%*%Psi%*%eta[h,])/sigmasq_y +
                (1/sigmasq_y)*(eta[h,]%*%Psi%*%eta[h,])^2)
         logr = logr*(-0.5)
         
         logu = log(runif(1))
         
         if (logr > logu){
            eta[h,] = eta_star
            if(h>burn){
               acp[h] = acp[h] + 1
            }
         }
      }
      
      # --- Update Psi --- #
      MM = model.matrix(y~.^2 - 1,as.data.frame(eta))   # perform factorized regression
      X_reg = cbind(eta^2,MM[,(k+1):ncol(MM)])
      X_reg.T = t(X_reg)
      Lambda_n = X_reg.T%*%X_reg/sigmasq_y + diag(rep(1,ncol(X_reg)))/100
      Vcsi = solve(Lambda_n)
      Mcsi = Vcsi%*%X_reg.T%*%(y-eta%*%phi)/sigmasq_y
      csi = bayesSurv::rMVNorm(n=1,mean=Mcsi,Sigma=Vcsi)
      lambda_diag = csi[1:k]
      lambda = csi[(k+1):length(csi)]
      Psi[lower.tri(Psi)] = lambda/2
      Psi[upper.tri(Psi)] = 0
      Psi = Psi + t(Psi)
      diag(Psi) = lambda_diag
      
      # --- Update phi --- #
      eta.T = t(eta)
      Lambda_n = eta.T%*%eta/sigmasq_y + diag(rep(1,ncol(eta)))/100
      Vcsi = solve(Lambda_n)
      Mcsi = Vcsi%*%eta.T%*%(y-diag(eta%*%Psi%*%eta.T))/sigmasq_y     # using updated psi
      phi = bayesSurv::rMVNorm(n = 1, mean = Mcsi, Sigma = Vcsi)
      
      # --- Update sigmasq_y --- #
      an = 0.5 + n/2
      bn = 0.5 + 0.5*t(y-eta%*%phi-diag(eta%*%Psi%*%eta.T))%*%
         (y-eta%*%phi-diag(eta%*%Psi%*%eta.T))
      sigmasq_y = 1/rgamma(1,an,bn)
      
      # -- update Lambda (rue & held) -- #
      eta2 = t(eta) %*% eta    # prepare eta crossproduct before the loop
      zlams = rnorm(k*p)       # generate normal draws all at once 
      
      for(j in 1:p) {
         Llamt = chol(diag(Plam[j,]) + ps[j]*eta2)
         Lambda[j,] = t(solve(Llamt,
                              zlams[1:k + (j-1)*k]) + 
                           solve(Llamt,
                                 solve(t(Llamt),
                                       ps[j] * t(eta) %*% X[,j])))
      }
      
      #------Update psi_{jh}'s------#
      psijh = matrix(rgamma(p*k,
                            df/2 + 0.5,
                            df/2 + t(t(Lambda)^2 * (tauh))),
                     nrow = p, ncol = k)
      
      #------Update theta & tauh------#
      mat = psijh * Lambda^2
      ad = ad1 + 0.5*p*k
      bd = bd1 + 0.5 * theta[1] * sum(tauh*colSums(mat))
      theta[1] = 1 / rgamma(1,ad,bd)           
      tauh = 1 / cumprod(theta)
      
      
      for(h in 2:k) {
         ad = ad2 + 0.5*p*(k-h+1)
         bd = bd2 * theta[h-1] + 0.5 * theta[h] * sum(tauh[h:k]*colSums(mat[,h:k, drop = F]))
         theta[h] = 1 / rgamma(1,ad,bd)
         tauh = 1 / cumprod(theta)
      }
      
      # -- Update Sigma -- #
      Xtil = X - eta %*% t(Lambda)
      ps= rgamma(p, as + 0.5*n, bs+0.5*colSums(Xtil^2))
      Sigma=diag(1/ps)
      
      #---update precision parameters----#
      Plam = t(t(psijh) * tauh)
      
      # ----- make adaptations ----#
      prob = 1/exp(b0 + b1*i)                    # probability of adapting
      uu = runif(1)
      lind = colSums(abs(Lambda) < epsilon)/p    # proportion of elements in each column less than eps in magnitude
      vec = lind >= prop
      num = sum(vec)                             # number of redundant columns
      
      if(uu < prob) {
         if((i > 20) & (num == 0) & all(lind < 0.995)) {
            k = k + 1
            Lambda = cbind(Lambda, rep(0,p))
            eta = cbind(eta,rnorm(n))
            psijh = cbind(psijh, rgamma(p,df/2,df/2))
            theta[k] = 1 / rgamma(1, ad2,bd2)
            tauh = 1 / cumprod(theta)
            Plam = t(t(psijh) * tauh)
         } else {
            if (num > 0) {
               k = max(k - num,1)
               Lambda = Lambda[,!vec, drop = F]
               psijh = psijh[,!vec, drop = F]
               eta = eta[,!vec, drop = F]
               theta = theta[!vec]
               tauh = 1 / cumprod(theta)
               Plam = t(t(psijh) * tauh)
            }
         }
      }
      
      #store posterior samples
      if ((i %% thin == 0) & (i > burn)){
         
         # parameters of the factor model
         sigmasq_st[count] = sigmasq_y
         tau_st[count,] = tau
         
         #Bayesian estimators
         V_n = solve(Lambda.T%*%solve(Sigma)%*%Lambda+diag(rep(1,ncol(Lambda))))
         a_n = V_n%*%Lambda.T%*%solve(Sigma)
         dsVX = diag(sqrt(VX))
         Omega_bayes[count,,] = dsVX%*%t(a_n)%*%Psi%*%a_n%*%dsVX
         beta_bayes[count,] = as.vector(t(phi)%*%a_n)
         alpha_bayes[count] = tr(Psi%*%V_n)
         k_st[count] = k
         count = count + 1
      }
   }
   
   beta_bayes = beta_bayes%*%diag(sqrt(VX))
   return(list(alpha_bayes = alpha_bayes,
               beta_bayes = beta_bayes,
               Omega_bayes = Omega_bayes,
               acp = acp,
               tau = tau_st,
               sigmasq_st = sigmasq_st))
}
