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


##### Source Functions from local git repo #####
# git clone https://github.com/fedfer/factor_interactions.git
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/functions")
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/generate_data")
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/post_processing")
#source("/Users/felpo/factor_interactions/codes/functions/Gibbs_DL.R")
#source("/Users/felpo/factor_interactions/codes/functions/quiet.R")
#source("/Users/felpo/factor_interactions/codes/functions/Gibbs_CUSP.R")
#sourceDirectory("/Users/felpo/factor_interactions/codes/generate_data")
#sourceDirectory("/Users/felpo/factor_interactions/codes/post_processing")
exists("generate_indep_model_notsparse")



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
k_true = as.numeric(args[5])

if(type_model == 0){
   if(sparse == 0){
      ratio_Om = 0.2
      ratio_beta = 0.2
      out_name = paste("n",n,"_p",p,"_ind","_notsparse.rds",sep="")
      k_start = 15
      type = "independent"
      
   }else if(sparse == 1){
      ratio_Om = 0.02
      ratio_beta = 0.1
      out_name = paste("n",n,"_p",p,"_ind","_sparse.rds",sep="")
      k_start = 15
      type = "independent"
      
   }
}else if(type_model == 1){
   if(sparse == 0){
      ratio_Om = 0.2
      ratio_beta = 0.2
      out_name = paste("n",n,"_p",p,"_corr","_notsparse.rds",sep="")
      k_start = k_true + 1
      type = "correlated"
      
   }else if (sparse == 1){
      ratio_Om = 0.01
      ratio_beta = 0.1
      out_name = paste("n",n,"_p",p,"_corr","_sparse.rds",sep="")
      k_start = k_true + 1
      type = "correlated"
   }
}else if(type_model == 2){
   if(sparse == 0){
      ratio_Om = 0.2
      ratio_beta = 0.2
      out_name = paste("n",n,"_p",p,"_power","_notsparse.rds",sep="")
      k_start = 12
      type = "power covariance"
   }else if (sparse == 1){
      ratio_Om = 0.01
      ratio_beta = 0.1
      out_name = paste("n",n,"_p",p,"_power","_sparse.rds",sep="")
      k_start = 12
      type = "power covariance"
   }
}


S = 20
err = err_pred = FR = err_beta = matrix(0,ncol = 8, nrow = S)
alpha = 0.05
TP_main = TN_main = TP_int = TN_int = matrix(0,ncol = 8, nrow = S)
acp_min = acp_max = acp_mean_st = matrix(0,ncol = 4,nrow = S)
delta_05 = delta_k = delta_10 = delta_50 = 0.2

for(s in 1:S){
   
   # generate the data
   data = generate_data(p = p,n = n,k_true = k_true,
                        ratio_Om = ratio_Om,ratio_beta = ratio_beta,
                        type = type)
   y = data$y; X = data$X; beta_true = data$beta_true; 
   Omega_true = data$Omega_true;
   y_test = data$y_test; X_test = data$X_test
   
   #Factor models
   nrun = 5000; burn = 4000; thin = 1; 
   
   
   # gibbs_DL_05 = gibbs_DL(y, X ,nrun, burn, thin = 1, 
   #                         delta_rw = delta_05, epsilon_rw = 0.5,
   #                        a = 10, k = NULL)
   
   
   gibbs_DL_k = gibbs_DL(y, X ,nrun, burn, thin = 1, 
                          delta_rw = delta_k, epsilon_rw = 0.5,
                          a = k_start, k = k_start )
   
   # apply(gibbs_DL_k$beta_bayes,2,mean)
   # Omega_hat = apply(gibbs_DL_k$Omega_bayes,c(2,3),mean)
   # plot(gibbs_DL_k$beta_bayes[,4],ty="l")
   # plot(gibbs_DL_k$Omega_bayes[,2,2],ty="l")
   
   
   # gibbs_CUSP_10 = gibbs_CUSP(y, X ,nrun, burn, thin = 1, 
   #                        delta_rw = delta_10, epsilon_rw = 0.5, 
   #                        k = NULL, alpha_prior = p*floor(log(p)*3)/10,
   #                        theta_inf = 0.05)
   # 
   # gibbs_CUSP_50 = gibbs_CUSP(y, X ,nrun, burn, thin = 1, 
   #                            delta_rw = delta_50, epsilon_rw = 0.5, 
   #                            k = NULL, alpha_prior = p*floor(log(p)*3)/2,
   #                            theta_inf = 0.05)
   
   # acp 
   # acp_1 = gibbs_DL_05$acp;acp_2 = gibbs_DL_k$acp;acp_3 = gibbs_CUSP_10$acp;acp_4 = gibbs_CUSP_50$acp;
   # acp_min[s,1] = min(acp_1/(nrun-burn));acp_max[s,1] = max(acp_1/(nrun-burn));acp_mean_st[s,1] = mean(acp_1/(nrun-burn))
   # acp_min[s,2] = min(acp_2/(nrun-burn));acp_max[s,2] = max(acp_2/(nrun-burn));acp_mean_st[s,2] = mean(acp_2/(nrun-burn))
   # acp_min[s,3] = min(acp_3/(nrun-burn));acp_max[s,3] = max(acp_3/(nrun-burn));acp_mean_st[s,3] = mean(acp_3/(nrun-burn))
   # acp_min[s,4] = min(acp_4/(nrun-burn));acp_max[s,4] = max(acp_4/(nrun-burn));acp_mean_st[s,4] = mean(acp_4/(nrun-burn))
   # 
   
   # Competitors
   hiernet = quiet(Hiernet_fct(y, X, X_test, y_test))
   Family = quiet(FAMILY_fct(y, X, X_test, y_test))
   PIE = PIE_fct(y, X, X_test, y_test)
   RAMP = RAMP_fct(y, X, X_test, y_test)
   
   #Errors
   errors = compute_errors(hiernet,Family,PIE,RAMP,
                           y,y_test,Omega_true,beta_true,
                           gibbs_DL_k,gibbs_DL_k,
                           gibbs_DL_k,gibbs_DL_k)
   
   
   # errors = compute_errors(hiernet,Family,PIE,RAMP,
   #                         y,y_test,Omega_true,beta_true,
   #                         gibbs_DL_05,gibbs_DL_k,
   #                         gibbs_CUSP_10,gibbs_CUSP_50)
   
   #Error 
   err[s,] = errors$err_insample[-1]
   err_pred[s,] = errors$err_pred[-1]
   FR[s,] = errors$FR_norm[-1]
   err_beta[s,] = errors$beta_MSE[-1]
   
   #TP and TNs
   rate_main2 = rate_recovery_maineff(gibbs_DL_k,gibbs_DL_k,alpha = alpha,beta_true = beta_true,
                                      hiernet$beta,Family$beta,PIE$beta,RAMP$beta)
   rate_int2 = rate_recovery_interactions(gibbs_DL_k,gibbs_DL_k,alpha = alpha,Omega_true=Omega_true,
                                          hiernet$Omega,Family$Omega,PIE$Omega,RAMP$Omega)
   rate_main1 = rate_recovery_maineff(gibbs_DL_k,gibbs_DL_k,alpha = alpha,beta_true = beta_true,
                                      hiernet$beta,Family$beta,PIE$beta,RAMP$beta)
   rate_int1 = rate_recovery_interactions(gibbs_DL_k,gibbs_DL_k,alpha = alpha,Omega_true=Omega_true,
                                          hiernet$Omega,Family$Omega,PIE$Omega,RAMP$Omega)
   
   # rate_main2 = rate_recovery_maineff(gibbs_DL_05,gibbs_DL_k,alpha = alpha,beta_true = beta_true,
   #                                    hiernet$beta,Family$beta,PIE$beta,RAMP$beta)
   # rate_int2 = rate_recovery_interactions(gibbs_DL_05,gibbs_DL_k,alpha = alpha,Omega_true=Omega_true,
   #                                        hiernet$Omega,Family$Omega,PIE$Omega,RAMP$Omega)
   # rate_main1 = rate_recovery_maineff(gibbs_CUSP_10,gibbs_CUSP_50,alpha = alpha,beta_true = beta_true,
   #                                    hiernet$beta,Family$beta,PIE$beta,RAMP$beta)
   # rate_int1 = rate_recovery_interactions(gibbs_CUSP_10,gibbs_CUSP_50,alpha = alpha,Omega_true=Omega_true,
   #                                        hiernet$Omega,Family$Omega,PIE$Omega,RAMP$Omega)
   
   TP_main[s,] = c(rate_main2$TP[1:2],rate_main1$TP)
   TN_main[s,] = c(rate_main2$TN[1:2],rate_main1$TN)
   TP_int[s,] = c(rate_int2$TP[1:2],rate_int1$TP)
   TN_int[s,] = c(rate_int2$TN[1:2],rate_int1$TN)
   
   print(s)
}

TP_main = TP_main[,c(5:8,1:4)]; TN_main = TN_main[,c(5:8,1:4)]
TP_int = TP_int[,c(5:8,1:4)]; TN_int = TN_int[,c(5:8,1:4)]
col_names = c("Hiernet","Family","Pie","RAMP","DL_05","DL_k","CUSP_10","CUSP_50")
#col_names = c("DL_05","DL_k","CUSP_10","CUSP_50","Hiernet","Family","Pie","RAMP")
colnames(err_beta) = colnames(err_pred) = colnames(err) = colnames(FR) = col_names
colnames(TP_main) = colnames(TN_main) = colnames(TP_int) = colnames(TN_int) = col_names


list_res = list(
   TP_main = TP_main,
   TN_main = TN_main,
   TP_int = TP_int,
   TN_int = TN_int,
   MSE_beta = err_beta,
   err_pred = err_pred,
   err_test = err,
   FR = FR
   #acp_min = acp_min,
   #acp_mean = acp_mean_st,
   #acp_max = acp_max,
   #delta = c(delta_k,delta_05,delta_10,delta_50)
)

results_dir = file.path("/work/sta790/ff31/results_fact")
#dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list_res, file.path(results_dir, out_name))
#readRDS("results/n100_p10_ind_notsparse")






