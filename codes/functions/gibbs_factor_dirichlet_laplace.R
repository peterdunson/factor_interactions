gibbs_factor_dirichlet_laplace = function(y, X ,nrun = 2000, burn = 1000, thin = 1, 
                                    delta_rw = 0.002, epsilon_rw = 0.5,
                                    a = 1/2, k = floor(log(p)*3)){
   
   #tau here is the parameter of the DL prior
   
   #this is needed to sample from a Dirichlet distribution
   library(MCMCpack)
   #this is needed to sample from a GIG (generalized inverse guassian) 
   #distribution
   library(GeneralizedHyperbolic)
   
   #to sample from inverse gaussian
   library(statmod)
   
   n = nrow(X)
   p = ncol(X)
   
   
   #### scaling ####
   M = apply(X,2,mean)
   VX = apply(X, 2, var)
   X_nst = X
   C = cov(X,y)
   #trick quadratic regression
   Lambda_X = array(0,c(n,p,p))
   for(i in 1:n){
      Lambda_X[i,,] = X[i,]%*%t(X[i,])
   }
   X = as.matrix(scale(X))
   
   #### gibbs sampler ####
   sp = (nrun - burn)/thin
   kinit = k
   b0 = 1
   b1 = 0.0005
   epsilon = 1e-3
   prop = 1.00
   
   # store the posterior samples
   Omega = array(0,c(sp,p,p))
   Omega1 = array(0,c(sp,p,p))
   Omega2 = array(0,c(sp,p,p))
   Omega_bayes = array(0,c(sp,p,p))
   Lambda_st = array(0,c(sp,p,5))
   phi_st = array(0,c(sp,5))
   tau_st = array(0,c(sp,p))
   sigmasq_st = numeric(sp)
   sigmasqX_st = numeric(sp)
   alpha_bayes = numeric(sp)
   k_st = numeric(sp)
   beta_st2 = beta_bayes = matrix(0,sp,p)
   lambda_2 = matrix(0,sp,p*(p-1)/2)
   lambda_diag_2 = matrix(0,sp,p)
   Cov_Xy = Cov_bayes = matrix(0,sp,p)
   
   
   #initial values number of factors
   num = 0
   #eta_k_init = eta0_initial(y,X)
   #k = eta_k_init$k

   # prior parameters
   as = 1;bs = 0.3;
   df = 150
   ad1 = 3;bd1 = 1
   ad2 = 4;bd2 = 1
   adf = 1; bdf = 1
   
   # initial values
   sigmasq_y = 1
   phi = numeric(k)
   ps = rgamma(p,as,bs)
   Sigma = diag(1/ps)
   Lambda = matrix(0,p,k)
   eta = matrix(rnorm(n*k),n,k)
   #eta = eta_k_init$eta
   meta = matrix(0,n,k)
   veta = diag(rep(1,k))
   Psi = matrix(0,k,k)
   
   #DL prior
   tau = rgamma(p,p*a,1/2)
   psijh = matrix(rexp(p*k,1/2),p,k)
   phijh = matrix(0,p,k)
   for(j in 1:p){
      phijh[j,] = rdirichlet(1,rep(a,k))
   }
   Plam = psijh*(phijh^2)*matrix(rep(tau^2,k),p,k,byrow=F)
   
   # output
   nofout = numeric(nrun+1)
   nofout[1] = k
   nofout1 = numeric(sp)
   
   #acceptance rates of etas
   acp = numeric(n)
   
   t = t0 = Sys.time()
   
   for(s in 1:nrun){
      
      #update eta with fastest MH
      #MH = MH_fastest(phi,sigmasq_y,Lambda,ps,k,Z,alpha,X,y,eta,acp,s,burn,Psi)
      #eta = MH[[1]]
      #acp = MH[[2]]
      
      #new fast MH
      #MH = MH_fast_new_mean2_noZ(phi,sigmasq_y,Lambda,ps,k,X,y,eta,acp,s,burn,Psi)
      #eta = MH[[1]]
      #acp = MH[[2]]
      
      #new RW MH
      MH = MH_RW_noZ(phi,sigmasq_y,Lambda,ps,k,X,y,eta,acp,s,burn,Psi,delta_rw)
      eta = MH[[1]]
      acp = MH[[2]]
      
      # Langevin MH
      #MH = MH_Langevin_noZ(phi,sigmasq_y,Lambda,ps,k,X,y,eta,
      #                     acp,s,burn,Psi,epsilon_rw)
      #eta = MH[[1]]
      #acp = MH[[2]]
      
      #update eta with laplace MH, too slow
      #MH = MH_laplace(phi,sigmasq_y,Lambda,ps,k,Z,alpha,X,y,eta,acp,s,burn)
      #eta = MH[[1]]
      #acp = MH[[2]]
      
      #update eta with more precise MH
      #MH = MH_var(phi,sigmasq_y,Lambda,ps,k,Z,alpha,X,y,eta,acp,s,burn)
      #eta = MH[[1]]
      #acp = MH[[2]]
      
      #update eta with Modular Bayes approach
      #eta = MB_update_eta(phi,sigmasq_y,Lambda,ps,k,Z,alpha,X,y,eta)
      
      # update Psi
      MM = model.matrix(y~.^2 - 1,as.data.frame(eta))
      X_reg = cbind(eta^2,MM[,(k+1):ncol(MM)])
      Lambda_n = t(X_reg)%*%X_reg/sigmasq_y + diag(rep(1,ncol(X_reg)))/100
      Vcsi = solve(Lambda_n)
      Mcsi = Vcsi%*%t(X_reg)%*%(y-eta%*%phi)/sigmasq_y
      csi = as.numeric(rmvnorm(1,Mcsi,Vcsi))
      lambda_diag = csi[1:k]
      lambda = csi[(k+1):length(csi)]
      Psi[lower.tri(Psi)] = lambda/2
      Psi[upper.tri(Psi)] = 0
      Psi = Psi + t(Psi)
      diag(Psi) = lambda_diag
      
      # update phi
      Lambda_n = t(eta)%*%eta/sigmasq_y + diag(rep(1,ncol(eta)))/100
      Vcsi = solve(Lambda_n)
      Mcsi = Vcsi%*%t(eta)%*%(y-diag(eta%*%Psi%*%t(eta)))/sigmasq_y
      #Mcsi = Vcsi%*%t(eta)%*%(y-Z%*%alpha-X_reg%*%csi)/sigmasq_y
      phi = as.numeric(rmvnorm(1,Mcsi,Vcsi))
      
      #update sigmasq y when updating all together alpha phi and lambda
      an = 0.5 + n/2
      bn = 0.5 + 0.5*t(y-eta%*%phi-diag(eta%*%Psi%*%t(eta)))%*%
         (y-eta%*%phi-diag(eta%*%Psi%*%t(eta)))
      #bn = 0.5 + 0.5*t(y-Z%*%alpha-eta%*%phi-X_reg%*%csi)%*%
      #  (y-Z%*%alpha-eta%*%phi-X_reg%*%csi)
      sigmasq_y = 1/rgamma(1,an,bn)
      
      #update Lambda
      Plam = psijh*(phijh^2)*matrix(rep(tau^2,k),p,k,byrow=F)
      eta2 = t(eta)%*%eta
      for(j in 1:p){
         Qlam = diag(Plam[j,])+ps[j]*eta2
         blam = ps[j]*(t(eta)%*%X[,j])
         Llam = t(chol(Qlam))
         zlam = rnorm(k)
         vlam = solve(Llam,blam)
         mlam = solve(t(Llam),vlam)
         ylam = solve(t(Llam),zlam)
         Lambda[j,] = (ylam + mlam)
      }
      
      #update psijh
      mujh = phijh*matrix(rep(tau,k),p,k,byrow = F)/abs(Lambda)
      psijh = matrix(rinvgauss(k*p, mujh),p,k,byrow = T)
      
      #update tau
      for(j in 1:p){
         
         tau[j] = rgig(1,lambda = 1-k, psi = 1, 
                       chi = 2*sum(abs(Lambda[j,])/phijh[j,]))
         
      }
      
      #update phijh
      Tjh = numeric(k)
      for(j in 1:p){
         for(h in 1:k){
            Tjh[h] = rgig(1,lambda = a-1,psi = 1,
                          chi = 2*abs(Lambda[j,h]))
         }
         phijh[j,] = Tjh/sum(Tjh)
      }
      
      
      #update Sigma
      Ytil = X - eta%*%t(Lambda)
      ps = rgamma(p,as+0.5*n,1)
      ps = (1/(bs + 0.5*apply(Ytil^2,2,sum)))*ps
      Sigma = diag(1/ps)
      
      #no adaptations since the dirichlet variable is by row of Lambda
      # so we cannot add columns
      
      #store posterior samples
      if ((s%%thin==0)&(s>burn)){
         #Factor Model
         Omega[s-burn,,] = Lambda%*%t(Lambda)+Sigma
         Omega1[s-burn,,] = Omega[s-burn,,]*sqrt(VX%*%t(VX))
         nofout1[s-burn] = (nofout[(s-burn)/thin]-num)*(num>0)
         k_st[s-burn] = k
         sigmasq_st[s-burn] = sigmasq_y
         tau_st[s-burn,] = tau
         
         
         #Bayesian estimators
         #they work for non-standardized data
         Sigma_nst = Sigma*sqrt(VX%*%t(VX))
         V_n = solve(t(Lambda)%*%solve(Sigma)%*%Lambda+diag(rep(1,ncol(Lambda))))
         a_n = V_n%*%t(Lambda)%*%solve(Sigma)
         Omega_bayes[s-burn,,] = t(a_n)%*%Psi%*%a_n
         
         #rescale Omega_bayes
         Omega_bayes[s-burn,,] = diag(sqrt(VX))%*%Omega_bayes[s-burn,,]%*%diag(sqrt(VX))
         
         beta_bayes[s-burn,] = as.vector(t(phi)%*%a_n)
         Cov_bayes[s-burn,] = as.vector(t(phi)%*%t(Lambda))
         alpha_bayes[s-burn] = tr(Psi%*%V_n)
         
         # variance of y|X
         #an = 0.5 + n/2
         #bn = 0.5 +
         #0.5*t(y-X%*%beta_bayes[s-burn,]-diag(X%*%Omega_bayes[s-burn,,]%*%t(X)))%*%
         #(y--X%*%beta_bayes[s-burn,]-diag(X%*%Omega_bayes[s-burn,,]%*%t(X)))
         #sigmasqX_st[s-burn] = 1/rgamma(1,an,bn)
      }
      
      nofout[s+1] = k
      
      
      if (s%%100==0){
         print(paste("time for last 100 iterations:",round(as.numeric(Sys.time()-t),0),
                     "seconds",sep=" "))
         t = Sys.time()
         print(paste("iteration number",s,"out of",nrun,sep=" "))
         t_end = round(((as.numeric(difftime(Sys.time(),t0,units="secs")))/s)*(nrun-s),0)
         print(paste("estimated time to end:",t_end,"seconds",sep=" "))
      }
   }
   
   #rescale the data
   beta_bayes = beta_bayes%*%diag(sqrt(VX))
   
   
   return(list(alpha_bayes = alpha_bayes,
               beta_bayes = beta_bayes,
               Omega_bayes = Omega_bayes,
               acp = acp,
               tau = tau_st,
               sigmasq_st = sigmasq_st,
               k_st = nofout))
   
   
}