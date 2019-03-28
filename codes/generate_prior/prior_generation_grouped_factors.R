### Generate distribution from the prior up to 4-th order interactions
library(MCMCpack)
library(gridExtra)
library(grid)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(lattice)

S = 50000
k = 2
p = 6

# Parameters for dirichlet laplace
a_DL = 1/(p*k)
a_DL2 = 1/2
as = 1; bs = 0.3
v_om = 5
# Storage
# fourth order interactions in the model
beta1_X_DL = beta1_X_DL2 =  matrix(0,p,S)
beta2_X_DL = beta2_X_DL2 = array(0,c(p,p,S))
beta3_X_DL = beta3_X_DL2 = array(0,c(1,2,3,S))
beta4_X_DL = beta4_X_DL2 = array(0,c(1,2,3,4,S))

# only second order interactions
beta12_X_DL = beta12_X_DL2 =  matrix(0,p,S)
beta22_X_DL = beta22_X_DL2 = array(0,c(p,p,S))

for (i in 1:S){
   
   Lambda_DL = Lambda_DL2 = matrix(0,p,k)
   ps = rgamma(p,as,bs)
   # for(j in 1:p){
   #    phi_DL = rdirichlet(1, rep(a_DL,p))
   #    phi_DL2 = rdirichlet(1, rep(a_DL2,p))
   #    for(h in 1:k){
   # 
   #       #Dirichlet-Laplace
   #       psi = rexp(1,0.5)
   #       tau_DL = rgamma(1,k*a_DL,0.5)
   #       tau_DL2 = rgamma(1,k*a_DL2,0.5)
   #       Lambda_DL[j,h] = rnorm(1,0,tau_DL*phi_DL[j]*sqrt(psi))
   #       Lambda_DL2[j,h] = rnorm(1,0,tau_DL2*phi_DL2[j]*sqrt(psi))
   #    }
   # }
   
   Lambda_DL = matrix(c(1,1,1,0,0,0,0,0,0,1,1,1),p,k)
   
   omega1 = rnorm(k,0,v_om)
   omega2 = rnorm(k,0,v_om)
   omega3 = rnorm(k,0,v_om)
   omega4 = rnorm(k,0,v_om)
   A_DL = solve(t(Lambda_DL)%*%diag(ps)%*%Lambda_DL + 
                   diag(k))%*%t(Lambda_DL)%*%diag(ps)
   A_DL2 = solve(t(Lambda_DL2)%*%diag(ps)%*%Lambda_DL2 + 
                    diag(k))%*%t(Lambda_DL2)%*%diag(ps)
   
   
   ### only interactions up to second order
   # multiply by to so we get the actuall coefficient and not the matrix
   beta12_X_DL2[,i] = t(A_DL2)%*%omega1
   beta12_X_DL[,i] = t(A_DL)%*%omega1
   
   beta22_X_DL[,,i] = 2*t(A_DL)%*%diag(omega2)%*%A_DL
   beta22_X_DL2[,,i] = 2*t(A_DL2)%*%diag(omega2)%*%A_DL2
   
   # up to 4th order interactions
   beta1_X_DL[,i] = t(A_DL)%*%omega1 + 3*diag(ps)%*%t(A_DL)%*%omega3
   beta1_X_DL2[,i] = t(A_DL2)%*%omega1 + 3*diag(ps)%*%t(A_DL2)%*%omega3
   
   beta2_X_DL[,,i] = 2*(2*t(A_DL)%*%diag(omega2)%*%A_DL + 
                           6*diag(1/ps)%*%t(A_DL)%*%diag(omega2)%*%A_DL)
   beta2_X_DL2[,,i] = 2*(2*t(A_DL2)%*%diag(omega2)%*%A_DL2 + 
                            6*diag(1/ps)%*%t(A_DL2)%*%diag(omega2)%*%A_DL2)
   
   beta3_X_DL[1,2,3,i] = 6*sum(omega3*A_DL[,1]*A_DL[,2]*A_DL[,3])
   beta3_X_DL2[1,2,3,i] = 6*sum(omega3*A_DL2[,1]*A_DL2[,2]*A_DL2[,3])
   #beta3_X_DL[1,2,3,i] = 6*sum(omega3*A_DL[1,]*A_DL[2,]*A_DL[3,])
   
   
   
   beta4_X_DL[1,2,3,4,i] = 24*sum(omega4*A_DL[,1]*A_DL[,2]*
                                     A_DL[,3]*A_DL[,4])
   beta4_X_DL2[1,2,3,4,i] = 24*sum(omega4*A_DL2[,1]*A_DL2[,2]*
                                      A_DL2[,3]*A_DL2[,4])
   #beta4_X_DL2[1,2,3,4,i] = 24*sum(omega4*A_DL2[1,]*A_DL2[2,]*
   #                                   A_DL2[3,]*A_DL2[4,])
   
   
   if (i%%500==0){
      print(i)
   }
}


df = data.frame(main_effect = beta1_X_DL[1,], 
                main_effect2 = beta1_X_DL[2,],
                interaction_12 = beta2_X_DL[1,2,],
                interaction_13 = beta2_X_DL[1,3,],
                interaction_14 = beta2_X_DL[1,4,],
                interaction_34 = beta2_X_DL[3,4,],
                quad_11 = beta2_X_DL[1,1,],
                interaction_3rd = beta3_X_DL[1,2,3,])


X = cbind(df$interaction_12,df$interaction_34)

plot(X, xlab="X label", ylab="Y label", pch=19, cex=Cex,
     xlim = xLim,ylim =xLim)
