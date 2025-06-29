######## FACTOR MODEL FOR CLUSTER ####

#### Load Libraries #####
library(mvtnorm)
library(MASS)
library(beepr)
library(psych)
library(bayesSurv)
library('PIE')
library('glmnet')
library(RAMP)
library(hierNet)
library(RCurl)
library(stargazer)
library(R.utils)
library(GIGrvg)
library(coda)


##### Source Functions from local git repo #####
# git clone https://github.com/fedfer/factor_interactions.git
source(file.path(repo_path, "codes/functions/Gibbs_DL.R"))
source(file.path(repo_path, "codes/functions/quiet.R"))

# Source all scripts in relevant folders
R.utils::sourceDirectory(file.path(repo_path, "codes/generate_data"))
R.utils::sourceDirectory(file.path(repo_path, "codes/functions"))
#R.utils::sourceDirectory(file.path(repo_path, "codes/process_results"))
exists("generate_indep_model_notsparse")


n = 500; p = 50
# if type_model = 0 --> indep model; othw factor model with high correlation
# if sparse = 0 --> not sparse Omega; othw sparse omega
type_model = 0; sparse = 0
# number of true factors in correlated model
k_true = 17
# noise in the model
sigmasq = 1

if(type_model == 0){
   if(sparse == 0){
      ratio_Om = 0.2
      ratio_beta = 0.2
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_notsparse.rds",sep="")
      k_start = 15
      type = "independent"
      
   }else if(sparse == 1){
      ratio_Om = 0.02
      ratio_beta = 0.1
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_sparse.rds",sep="")
      k_start = 15
      type = "independent"
      
   }
}else if(type_model == 1){
   if(sparse == 0){
      ratio_Om = 0.2
      ratio_beta = 0.2
      out_name = paste("n",n,"_lp",p,"_sigmasq",sigmasq,"_corr","_notsparse.rds",sep="")
      k_start = k_true + 1
      type = "correlated"
      
   }else if (sparse == 1){
      ratio_Om = 0.01
      ratio_beta = 0.1
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_corr","_sparse.rds",sep="")
      k_start = k_true + 1
      type = "correlated"
   }
}else if(type_model == 2){
   if(sparse == 0){
      ratio_Om = 0.2
      ratio_beta = 0.2
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_wishart","_notsparse.rds",sep="")
      k_start = 23
      #type = "wishart"
      type = "power covariance"
   }else if (sparse == 1){
      ratio_Om = 0.01
      ratio_beta = 0.1
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_wishart","_sparse.rds",sep="")
      k_start = 23
      type = "power covariance"
   }
}





# generate the data
data = generate_data(
   p = p,
   n = n,
   k_true = k_true,
   ratio_Om = ratio_Om,
   ratio_beta = ratio_beta,
   type = type,
   sigmasq = sigmasq
)
y = data$y
X = data$X
beta_true = data$beta_true

cor(X) %>% abs() %>% mean()

Omega_true = data$Omega_true

y_test = data$y_test
X_test = data$X_test

#Factor models
nrun = 2000
burn = 1500
thin = 5



eig_values = eigen(cor(X))$values
plot(eig_values)

k_start = 25
# gibbs_DL_k = gibbs_DL(
#    y,
#    X ,
#    nrun,
#    burn,
#    thin = 1,
#    delta_rw = delta_k,
#    epsilon_rw = 0.5,
#    a = k_start,
#    k = k_start
# )

gibbs_DL_P = gibbs_DL_Plam(
   y,
   X ,
   nrun,
   burn,
   thin = 1,
   delta_rw = 1,
   epsilon_rw = 0.5,
   a = k_start,
   k = k_start
)
gibbs_DL_k = gibbs_DL_P


# Coverage
cov_y = coverage_y(y_test, X_test, gibbs_DL_P)
bias_pred = cov_y$bias
coverage_pred = cov_y$coverage

library(beepr)
beep()
cov_y

# apply(gibbs_DL_k$beta_bayes,2,effectiveSize)
# apply(gibbs_DL_k$Omega_bayes,c(2,3),effectiveSize)
# apply(gibbs_DL_k$beta_bayes,2,mean)
# Omega_hat = apply(gibbs_DL_k$Omega_bayes,c(2,3),mean)
# plot(gibbs_DL_P$beta_bayes[,4],ty="l")
# plot(gibbs_DL_P$Omega_bayes[,2,2],ty="l")


# Competitors
hiernet = quiet(Hiernet_fct(y, X, X_test, y_test))
#Family = quiet(FAMILY_fct(y, X, X_test, y_test))
#PIE = PIE_fct(y, X, X_test, y_test)
RAMP = RAMP_fct(y, X, X_test, y_test)
PIE = RAMP
Family = RAMP

str(hiernet)
str(RAMP)
str(Omega_true)
str(beta_true)
str(y)
str(y_test)



#Errors
errors = compute_errors(
   hiernet,
   Family,
   PIE,
   RAMP,
   y,
   y_test,
   Omega_true,
   beta_true,
   gibbs_DL_P,
   gibbs_DL_P,
   gibbs_DL_k,
   gibbs_DL_k
)
errors
library(beepr)
beep()
Sys.sleep(1)
beep()
cov_y
