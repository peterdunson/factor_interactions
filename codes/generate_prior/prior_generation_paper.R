### Generate distribution from the prior up to 4-th order interactions
library(MCMCpack)
library(gridExtra)
library(grid)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(lattice)

S = 100000
k = 5
p = 20

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
   for(j in 1:p){
      phi_DL = rdirichlet(1, rep(a_DL,p))
      phi_DL2 = rdirichlet(1, rep(a_DL2,p))
      for(h in 1:k){
         
         #Dirichlet-Laplace
         psi = rexp(1,0.5)
         tau_DL = rgamma(1,k*a_DL,0.5)
         tau_DL2 = rgamma(1,k*a_DL2,0.5)
         Lambda_DL[j,h] = rnorm(1,0,tau_DL*phi_DL[j]*sqrt(psi))
         Lambda_DL2[j,h] = rnorm(1,0,tau_DL2*phi_DL2[j]*sqrt(psi))
      }
   }
   
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


# DL
# > means heavier above
# tails_DL_a=1/(p*k) > tails_DL_a=1/p > tails_Sparse > tails_DL_a=1/2

# Only two ways int in the model
xlim = c(-2,2)
par(mfrow = c(1,2))
hist(beta12_X_DL2[1,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim)
lines(density(beta12_X_DL2[1,]))
abline(v = quantile(beta12_X_DL2[1,],probs = c(0.05,0.95)),col="red",lty="dotted")
hist(beta22_X_DL2[1,2,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim)
abline(v = quantile(beta22_X_DL2[1,2,],probs = c(0.05,0.95)),col="red",lty="dotted")


xlim = c(-3,3)
par(mfrow = c(1,3))
hist(beta1_X_DL2[1,],freq = F,breaks = seq(-10000,10000,by =0.05),xlim = xlim,
     xlab = "main effect",main = "")
abline(v = quantile(beta1_X_DL2[1,],probs = c(0.05,0.95)),col="red",lty="dotted")
abline(v = quantile(beta1_X_DL2[1,],probs = c(0.25,0.75)),col="green",lty="dotted")
hist(beta2_X_DL2[1,2,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim,
     xlab = "2nd order interaction",main = "k = 5, p = 20")
abline(v = quantile(beta2_X_DL2[1,2,],probs = c(0.05,0.95)),col="red",lty="dotted")
abline(v = quantile(beta2_X_DL2[1,2,],probs = c(0.25,0.75)),col="green",lty="dotted")
hist(beta2_X_DL2[1,2,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim,
     xlab = "3rd order interaction",main = "")
abline(v = quantile(beta3_X_DL2[1,2,3,],probs = c(0.05,0.95)),col="red",lty="dotted")
abline(v = quantile(beta3_X_DL2[1,2,3,],probs = c(0.25,0.75)),col="green",lty="dotted")




df = data.frame(main_effect = beta1_X_DL2[1,], 
                main_effect2 = beta1_X_DL2[2,],
                interaction_2nd = beta2_X_DL2[1,2,],
                interaction_2nd2 = beta2_X_DL2[1,3,],
                interaction_3rd = beta3_X_DL2[1,2,3,])
p1 = ggplot(df, aes(x=main_effect)) + 
   geom_density()+
   xlim(c(-20,20))+
   xlab("")+
   ggtitle("main effect")+
   geom_vline(xintercept=quantile(df$main_effect,probs = c(0.05)), linetype="dashed", color = "red")+
   geom_vline(xintercept=quantile(df$main_effect,probs = c(0.95)), linetype="dashed", color = "red")+
   geom_vline(xintercept=quantile(df$main_effect,probs = c(0.25)), linetype="dashed", color = "green")+
   geom_vline(xintercept=quantile(df$main_effect,probs = c(0.75)), linetype="dashed", color = "green")


p2 = ggplot(df, aes(x=interaction_2nd)) + 
   geom_density()+ 
   xlim(c(-2,2))+
   xlab("")+
   ggtitle("2nd order")+
   geom_vline(xintercept=quantile(df$interaction_2nd,probs = c(0.05)), linetype="dashed", color = "red")+
   geom_vline(xintercept=quantile(df$interaction_2nd,probs = c(0.95)), linetype="dashed", color = "red")+
   geom_vline(xintercept=quantile(df$interaction_2nd,probs = c(0.25)), linetype="dashed", color = "green")+
   geom_vline(xintercept=quantile(df$interaction_2nd,probs = c(0.75)), linetype="dashed", color = "green")

p3 = ggplot(df, aes(x=interaction_3rd)) + 
   geom_density()+ 
   xlim(c(-1,1))+
   xlab("")+
   ggtitle("3rd order")+
   geom_vline(xintercept=quantile(df$interaction_3rd,probs = c(0.05)), linetype="dashed", color = "red")+
   geom_vline(xintercept=quantile(df$interaction_3rd,probs = c(0.95)), linetype="dashed", color = "red")+
   geom_vline(xintercept=quantile(df$interaction_3rd,probs = c(0.25)), linetype="dashed", color = "green")+
   geom_vline(xintercept=quantile(df$interaction_3rd,probs = c(0.75)), linetype="dashed", color = "green")


grid.arrange(p1,p2,p3,
             ncol = 3,
             top=textGrob("p = 20, k = 5",gp=gpar(fontsize=20,font=3)))

commonTheme = list(labs(color="Density",fill="Density",
                        x="RNA-seq Expression",
                        y="Microarray Expression"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))


ggplot(data=df,aes(interaction_2nd,interaction_2nd2)) + 
   geom_density2d(aes(colour=..level..)) + 
   scale_colour_gradient(low="green",high="red") + 
   commonTheme


k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))

par(mfrow = c(1,2))
## compute 2D kernel density, see MASS book, pp. 130-131
#z = kde2d(df$interaction_2nd, df$interaction_2nd2, n=500)
xLim = c(-6,6)
Cex = 0.02
X = cbind(df$interaction_2nd,df$interaction_2nd2)
plot(X, xlab="X label", ylab="Y label", pch=19, cex=Cex,
     xlim = xLim,ylim =xLim)
#contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

X = cbind(df$main_effect,df$main_effect2)
plot(X, xlab="X label", ylab="Y label", pch=19, cex=Cex,
     xlim =xLim,ylim =xLim)


par(mfrow = c(1,2))
hist(beta12_X_DL[1,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim)
abline(v = quantile(beta12_X_DL[1,],probs = c(0.05,0.95)),col="red",lty="dotted")
hist(beta22_X_DL[1,2,],freq = F,breaks = seq(-1000,1000,by =0.01),xlim = xlim)
abline(v = quantile(beta22_X_DL[1,2,],probs = c(0.05,0.95)),col="red",lty="dotted")

xlim = c(-2,2)
par(mfrow = c(2,2))

hist(beta1_X_DL2[1,],freq = F,breaks = seq(-10000,10000,by =0.05),xlim = xlim)
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

