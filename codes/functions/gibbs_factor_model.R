gibbs_factor_model = function(y, X ,nrun = 2000, burn = 1000, thin = 1, 
                              delta_rw = 0.002, epsilon_rw = 0.5){
   
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
   kinit = floor(log(p)*3)
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
   k = kinit
   
   # prior parameters
   as = 1;bs = 0.3;
   df = 3
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
   
   psijh = matrix(rgamma(p*k,df/2,df/2),p,k)
   delta = c(rgamma(1,ad1,1/bd1),rgamma(k-1,ad2,1/bd2))
   tauh = cumprod(delta)
   Plam = psijh*matrix(rep(tauh,p),p,k,byrow=T)
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
      
      #update psijh's
      psijh = rgamma(p*k,df/2 + 0.5,1)
      psijh = psijh*(1/(df/2+(Lambda^2)*matrix(rep(tauh,p),p,k,byrow = T)))
      
      #update delta e tauh
      mat = psijh*(Lambda^2)
      ad = ad1 + 0.5*p*k
      bd = bd1 + 0.5*(1/delta[1])*(sum(tauh*apply(mat,2,sum)))
      delta[1] = rgamma(1,ad,bd)
      tauh = cumprod(delta)
      
      for(h in 2:(k)){
         ad = ad2 + 0.5*p*(k-h+1)
         bd = bd2 + 0.5*(1/delta[h])*
            sum(tauh[h:length(tauh)]*apply(as.matrix(mat[,h:ncol(mat)]),2,sum))
         delta[h] = rgamma(1,ad,bd)
         tauh = cumprod(delta)
      }
      
      #update Sigma
      Ytil = X - eta%*%t(Lambda)
      ps = rgamma(p,as+0.5*n,1)
      ps = (1/(bs + 0.5*apply(Ytil^2,2,sum)))*ps
      Sigma = diag(1/ps)
      
      #update precision parameters
      Plam = psijh*matrix(rep(tauh,p),p,k,byrow=T)
      
      #make adaptations
      prob = 1/(exp(b0+b1*i))
      uu = runif(1)
      lind = apply((abs(Lambda)< epsilon)/p,2,sum)
      vec = (lind >= prop)
      num = sum(vec)
      
      
      #store posterior samples
      if ((s%%thin==0)&(s>burn)){
         #Factor Model
         Omega[s-burn,,] = Lambda%*%t(Lambda)+Sigma
         Omega1[s-burn,,] = Omega[s-burn,,]*sqrt(VX%*%t(VX))
         nofout1[s-burn] = (nofout[(s-burn)/thin]-num)*(num>0)
         k_st[s-burn] = k
         sigmasq_st[s-burn] = sigmasq_y
         
         
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
      
      if (uu<prob){
         if ((i>20)&(num==0)&all(lind<0.995)){
            k = k+1
            Lambda = cbind(Lambda,0)
            phi = c(phi,0)
            eta = cbind(eta,rnorm(n))
            psijh = cbind(psijh,rgamma(p,df/2,df/2))
            delta = c(delta,rgamma(1,ad2,bd2))
            tauh = cumprod(delta)
            Plam = psijh*matrix(rep(tauh,p),p,k,byrow=T)
            Psi = rbind(Psi,0)
            Psi = cbind(Psi,0)
         }else if (num > 0){
            nonred = setdiff(1:k,which(vec !=0))
            k = max(k - num,1)
            Lambda = Lambda[,nonred]
            phi = phi[nonred]
            psijh = psijh[,nonred]
            eta = eta[,nonred]
            delta = delta[nonred]
            tauh = cumprod(delta)
            Plam = psijh*matrix(rep(tauh,p),p,k,byrow=T)
            Psi = Psi[nonred,nonred]
         }
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
               sigmasq_st = sigmasq_st,
               k_st = nofout))
   
   
}