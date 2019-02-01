### Generate distribution from the prior up to 4-th order interactions
library(MCMCpack)


S = 10000
k = 5
p = 10

# Parameters for dirichlet laplace
a_DL = 1/(p*k)
a_DL2 = 1/k
as = 1; bs = 0.3
# Storage
beta1_X_DL = beta1_X_DL2 =  matrix(0,p,S)
beta2_X_DL = beta2_X_DL2 = array(0,c(p,p,S))
beta3_X_DL = beta3_X_DL2 = array(0,c(1,2,3,S))
beta4_X_DL = beta4_X_DL2 = array(0,c(1,2,3,4,S))


for (i in 1:S){
   
   Lambda_DL = Lambda_DL2 = matrix(0,p,k)
   ps = rgamma(p,as,bs)
   for(j in 1:p){
      phi_DL = rdirichlet(1, rep(a_DL,p))
      phi_DL2 = rdirichlet(1, rep(a_DL2,p))
      for(h in 1:k){
         
         #Dirichlet-Laplace
         psi = rexp(1,0.5)
         tau_DL = rgamma(1,p*a_DL,0.5)
         tau_DL2 = rgamma(1,p*a_DL2,0.5)
         Lambda_DL[j,h] = rnorm(1,0,tau_DL*phi_DL[j]*sqrt(psi))
         Lambda_DL2[j,h] = rnorm(1,0,tau_DL2*phi_DL2[j]*sqrt(psi))
      }
   }
   
   omega1 = rnorm(k,0,10)
   omega2 = rnorm(k,0,10)
   omega3 = rnorm(k,0,10)
   omega4 = rnorm(k,0,10)
   
   A_DL = solve(t(Lambda_DL)%*%diag(ps)%*%Lambda_DL + 
                   diag(k))%*%t(Lambda_DL)%*%diag(ps)
   A_DL2 = solve(t(Lambda_DL2)%*%diag(ps)%*%Lambda_DL2 + 
                    diag(k))%*%t(Lambda_DL2)%*%diag(ps)
   
   beta1_X_DL[,i] = t(A_DL)%*%omega1 + 3*diag(1/ps)%*%t(A_DL)%*%omega3
   beta1_X_DL2[,i] = t(A_DL2)%*%omega1 + 3*diag(1/ps)%*%t(A_DL2)%*%omega3
   
   beta2_X_DL[,,i] = 2*t(A_DL)%*%diag(omega2)%*%A_DL + 
      6*diag(1/ps)%*%t(A_DL)%*%diag(omega2)%*%A_DL
   beta2_X_DL2[,,i] = 2*t(A_DL2)%*%diag(omega2)%*%A_DL2 + 
      6*diag(1/ps)%*%t(A_DL2)%*%diag(omega2)%*%A_DL2
   
   beta3_X_DL[1,2,3,i] = 6*sum(omega3*A_DL[1,]*A_DL[2,]*A_DL[3,])
   beta3_X_DL2[1,2,3,i] = 6*sum(omega3*A_DL2[1,]*A_DL2[2,]*A_DL2[3,])
   
   beta4_X_DL[1,2,3,4,i] = 24*sum(omega4*A_DL[1,]*A_DL[2,]*
                                     A_DL[3,]*A_DL[4,])
   beta4_X_DL2[1,2,3,4,i] = 24*sum(omega4*A_DL2[1,]*A_DL2[2,]*
                                      A_DL2[3,]*A_DL2[4,])
   
   
   if (i%%500==0){
      print(i)
   }
}


# DL
# > means heavier above
# tails_DL_a=1/(p*k) > tails_DL_a=1/p > tails_Sparse > tails_DL_a=1/2
xlim = c(-2,2)
par(mfrow = c(2,2))

hist(beta1_X_DL2[1,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim)
abline(v = quantile(beta1_X_DL2[1,],probs = c(0.05,0.95)),col="red",lty="dotted")
hist(beta2_X_DL2[1,2,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim)
abline(v = quantile(beta2_X_DL2[1,2,],probs = c(0.05,0.95)),col="red",lty="dotted")
hist(beta3_X_DL2[1,2,3,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim)
abline(v = quantile(beta3_X_DL2[1,2,3,],probs = c(0.05,0.95)),col="red",lty="dotted")
hist(beta4_X_DL2[1,2,3,4,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim)
abline(v = quantile(beta4_X_DL2[1,2,3,4,],probs = c(0.05,0.95)),col="red",lty="dotted")




hist(beta1_X_DL[1,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim)
abline(v = quantile(beta1_X_DL[1,],probs = c(0.05,0.95)),col="red",lty="dotted")
hist(beta2_X_DL[1,2,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim)
abline(v = quantile(beta2_X_DL[1,2,],probs = c(0.05,0.95)),col="red",lty="dotted")
hist(beta3_X_DL[1,2,3,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim)
abline(v = quantile(beta3_X_DL[1,2,3,],probs = c(0.05,0.95)),col="red",lty="dotted")
hist(beta4_X_DL[1,2,3,4,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim)
abline(v = quantile(beta4_X_DL[1,2,3,4,],probs = c(0.05,0.95)),col="red",lty="dotted")




# comparison
round(cbind(quantile(beta1_X),quantile(beta1_X_DL)),2)



# the coefficients of quadratic effects has heavier tails than the one of interactions
# both have a spike in zero

# super heavy tails if we do not consider the shrinkage parameters tau
# probabibly the distribution for beta_X and Omega_X is improper
# might be a good prior distribution

eps_seq = seq(1,0.01,by = -0.01)
probs = numeric(length(eps_seq))

for (i in 1:length(eps_seq)){
   epsilon = eps_seq[i]
   probs[i] = mean(Omega_X > -epsilon & Omega_X < epsilon)
   print(i/length(eps_seq))
}

plot(eps_seq,probs)

sigma1 = rgamma(50000,2,1)
sigma2 = 1/sigma1
par(mfrow = c(1,2))
plot(density(sigma1))
plot(density(sigma2))

