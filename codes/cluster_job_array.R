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
library(FAMILY)
library(RCurl)
library(stargazer)
library(R.utils)
library(GIGrvg)
library(coda)


##### Source Functions from local git repo #####
# git clone https://github.com/fedfer/factor_interactions.git
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/functions")
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/generate_data")
#sourceDirectory("/work/sta790/ff31/factor_interactions/codes/process_results")
source("/work/sta790/ff31/factor_interactions/codes/process_results/compute_errors.R")
source("/work/sta790/ff31/factor_interactions/codes/process_results/rate_recovery.R")
source("/work/sta790/ff31/factor_interactions/codes/process_results/coverage_int.R")
#exists("generate_indep_model_notsparse")
exists("compute_errors")

##### Argument of R script ######
args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
   stop("At least one argument must be supplied", call.=FALSE)
}

n = as.numeric(args[1]); p = as.numeric(args[2])
# if type_model = 0 --> indep model; othw factor model with high correlation
# if sparse = 0 --> not sparse Omega; othw sparse omega
type_model = as.numeric(args[3]); sparse = as.numeric(args[4])
# number of true factors in correlated model
k_true = as.numeric(args[5]);
# noise in the model
sigmasq = as.numeric(args[6]);

# Get the slurm array ID
i = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

if(type_model == 0){
   if(sparse == 0){
      ratio_Om = 0.2
      ratio_beta = 0.2
      dir_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_notsparse",sep="")
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_notsparse",
                       "_iter=",i,".rds",sep="")
      k_start = 15
      type = "independent"
      
   }else if(sparse == 1){
      ratio_Om = 0.02
      ratio_beta = 0.1
      dir_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_notsparse",sep="")
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_sparse",
                       "_iter=",i,".rds",sep="")
      k_start = 15
      type = "independent"
      
   }
}else if(type_model == 1){
   if(sparse == 0){
      ratio_Om = 0.2
      ratio_beta = 0.2
      dir_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_notsparse",sep="")
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_corr","_notsparse",
                       "_iter=",i,".rds",sep="")
      k_start = k_true + 1
      type = "correlated"
      
   }else if (sparse == 1){
      ratio_Om = 0.01
      ratio_beta = 0.1
      dir_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_notsparse",sep="")
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_corr","_sparse",
                       "_iter=",i,".rds",sep="")
      k_start = k_true + 1
      type = "correlated"
   }
}else if(type_model == 2){
   if(sparse == 0){
      ratio_Om = 0.2
      ratio_beta = 0.2
      dir_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_notsparse",sep="")
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_wishart","_notsparse",
                       "_iter=",i,".rds",sep="")
      k_start = 12
      #type = "wishart"
      type = "power covariance"
   }else if (sparse == 1){
      ratio_Om = 0.01
      ratio_beta = 0.1
      dir_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_ind","_notsparse",sep="")
      out_name = paste("n",n,"_p",p,"_sigmasq",sigmasq,"_wishart","_sparse",
                       "_iter=",i,".rds",sep="")
      type = "power covariance"
   }
}

alpha = 0.05
delta_05 = delta_k = delta_10 = delta_50 = 0.2

   
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

Omega_true = data$Omega_true

y_test = data$y_test
X_test = data$X_test

#Factor models
nrun = 5000
burn = 4000
thin = 5

gibbs_DL_P = gibbs_DL_Plam(
   y,
   X ,
   nrun,
   burn,
   thin = 1,
   delta_rw = delta_k,
   epsilon_rw = 0.5,
   a = k_start,
   k = k_start
)



# Competitors
hiernet = quiet(Hiernet_fct(y, X, X_test, y_test))
Family = quiet(FAMILY_fct(y, X, X_test, y_test))
PIE = PIE_fct(y, X, X_test, y_test)
RAMP = RAMP_fct(y, X, X_test, y_test)

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
   gibbs_DL_P,
   gibbs_DL_P
)



#Error
err = errors$err_insample[-1]
err_pred = errors$err_pred[-1]
FR = errors$FR_norm[-1]
err_beta = errors$beta_MSE[-1]

#TP and TNs
rate_main2 = rate_recovery_maineff(
   gibbs_DL_P,
   gibbs_DL_P,
   alpha = alpha,
   beta_true = beta_true,
   hiernet$beta,
   Family$beta,
   PIE$beta,
   RAMP$beta
)
rate_int2 = rate_recovery_interactions(
   gibbs_DL_P,
   gibbs_DL_P,
   alpha = alpha,
   Omega_true = Omega_true,
   hiernet$Omega,
   Family$Omega,
   PIE$Omega,
   RAMP$Omega
)
rate_main1 = rate_recovery_maineff(
   gibbs_DL_P,
   gibbs_DL_P,
   alpha = alpha,
   beta_true = beta_true,
   hiernet$beta,
   Family$beta,
   PIE$beta,
   RAMP$beta
)
rate_int1 = rate_recovery_interactions(
   gibbs_DL_P,
   gibbs_DL_P,
   alpha = alpha,
   Omega_true = Omega_true,
   hiernet$Omega,
   Family$Omega,
   PIE$Omega,
   RAMP$Omega
)

TP_main = c(rate_main2$TP[1:2], rate_main1$TP)
TN_main = c(rate_main2$TN[1:2], rate_main1$TN)
TP_int = c(rate_int2$TP[1:2], rate_int1$TP)
TN_int = c(rate_int2$TN[1:2], rate_int1$TN)

TP_main = TP_main[c(5:8,1:4)]; TN_main = TN_main[c(5:8,1:4)]
TP_int = TP_int[c(5:8,1:4)]; TN_int = TN_int[c(5:8,1:4)]
#col_names = c("Hiernet","Family","Pie","RAMP","DL_P","DL_P","CUSP_10","CUSP_50")
#col_names = c("DL_05","DL_k","CUSP_10","CUSP_50","Hiernet","Family","Pie","RAMP")
#colnames(err_beta) = colnames(err_pred) = colnames(err) = colnames(FR) = col_names
#colnames(TP_main) = colnames(TN_main) = colnames(TP_int) = colnames(TN_int) = col_names


list_res = list(
   TP_main = TP_main,
   TN_main = TN_main,
   TP_int = TP_int,
   TN_int = TN_int,
   MSE_beta = err_beta,
   err_pred = err_pred,
   err_test = err,
   FR = FR
)


# save
dir_path = paste("/work/sta790/ff31/factor_interactions/results/array_jobs",
                 dir_name,sep="")

if(dir.exists(dir_path) == F){
   dir.create(dir_path)
}

saveRDS(list_res, file.path(dir_path, out_name))





