####### CLEAN CODE for simulations ######
#lapply(.packages(all.available = TRUE), 
#       function(xx) library(xx,     character.only = TRUE))
library(mvtnorm)
library(MASS)
library(beepr)
library(psych)
library(bayesSurv)
library('PIE')
library('glmnet')
library(RAMP)
library(hierNet)
library(FAMILY)
library(RCurl)
library(R.utils)
sourceDirectory("~/Desktop/Code_factor_simulations/Functions")

# generate the data
data = generate_factor_model(p = 50,n = 200,k_true = 7)
#data = generate_indep_model(p = 10,n = 100)
y = data$y; X = data$X; beta_true = data$beta_true; 
Omega_true = data$Omega_true;
y_test = data$y_test; X_test = data$X_test
p = ncol(X);

#Factor model
nrun = 2500; burn = 1500;
gibbs = gibbs_factor_model(y, X ,nrun = nrun, burn = burn, 
                           epsilon_rw = 0.05 , delta_rw = 0.01)
plot_acp(nrun,burn,gibbs)
beep()

gibbs_CUSP = gibbs_factor_CUSP(y, X ,nrun = nrun, burn = burn,
                delta_rw = 0.0005, epsilon_rw = 0.5, k = floor(log(p)*3),
                alpha_prior = p*floor(log(p)*3)/10, theta_inf = 0.05)
plot_acp(nrun,burn,gibbs_CUSP)


gibbs_notau = gibbs_factor_model_notau(y, X ,nrun = nrun, burn = burn, 
                                       epsilon_rw = 0.1, delta_rw = 0.01,
                                       k_max =  floor(log(ncol(X))*3))
plot_acp(nrun,burn,gibbs_notau)

gibbs_DL = gibbs_factor_dirichlet_laplace(y, X ,nrun = nrun, burn = burn, thin = 1, 
                                          delta_rw = 0.0005, epsilon_rw = 0.5,
                                          a = p, k = floor(log(ncol(X))*3))
plot_acp(nrun,burn,gibbs_DL)

# for good performance max_acp < 0.3
gibbs_highint = gibbs_factor_model_highint(y, X ,nrun = nrun, burn = burn, thin = 1, 
                                      delta_rw = 0.008, epsilon_rw = 0.5,
                                      k_max = floor(log(ncol(X))*3))
plot_acp(nrun,burn,gibbs_highint)

#acp
plot_acp(nrun,burn,gibbs)

#Competitors
hiernet = Hiernet_fct(y, X, X_test, y_test)
Family = FAMILY_fct(y, X, X_test, y_test)
PIE = PIE_fct(y, X, X_test, y_test)
RAMP = RAMP_fct(y, X, X_test, y_test)

#Compare 
compute_errors(hiernet,Family,PIE,RAMP,
               y,y_test,Omega_true,beta_true,
               gibbs_DL,gibbs_CUSP)
   


gibbs_notau = gibbs_highint
#beep()
alpha_bayes_hat = mean(gibbs$alpha_bayes)
betabayes_hat = apply(gibbs$beta_bayes,2,mean)
Omegabayes_hat = apply(gibbs$Omega_bayes,c(2,3),mean)
err_factor = y - X%*%betabayes_hat - 
   as.vector(diag(X%*%Omegabayes_hat%*%t(X))) - alpha_bayes_hat
err_pred_factor = y_test - X_test%*%betabayes_hat - 
   as.vector(diag(X_test%*%Omegabayes_hat%*%t(X_test))) - alpha_bayes_hat


## Coverage
coverage = coverage_int(gibbs_DL$beta_bayes,gibbs_DL$Omega_bayes)
coverage2 = coverage_int(gibbs_notau$beta_bayes,gibbs_notau$Omega_bayes)
#betabayes_hat[!coverage[[1]]] = 0
#Omegabayes_hat[!coverage[[2]]] = 0 
#Omegabayes_hat2[!coverage2[[2]]] = 0 
#betabayes_hat2[!coverage2[[1]]] = 0

#Rate recovery
alpha = 0.2
rate_recovery_maineff(gibbs_DL,gibbs_highint,alpha = alpha,beta_true = beta_true,
                      hiernet$beta,Family$beta,PIE$beta,RAMP$beta)
rate_recovery_interactions(gibbs_DL,gibbs_highint,alpha = alpha,Omega_true=Omega_true,
                           hiernet$Omega,Family$Omega,PIE$Omega,RAMP$Omega)

alpha = seq(0.01,0.14,by = 0.01)
cov_f = cov_f2 = numeric(length(alpha))
for(i in 1:length(alpha)){
   rate = rate_recovery_interactions(gibbs,gibbs_notau,alpha = alpha[i],Omega_true=Omega_true,
                              hiernet$Omega,Family$Omega,PIE$Omega,RAMP$Omega)
   cov_f[i] = rate$TN[1]
   cov_f2[i] = rate$TN[2]
}

plot(1-alpha,cov_f, ylim = c(min(cov_f),1),pch=16)
lines(1-alpha,cov_f2,col="red",ty="p",pch=16)
abline(a = 0, b = 1)

plot(1-alpha,cov_f2,col="red",ty="p",pch=16)
abline(a = 0, b = 1)
