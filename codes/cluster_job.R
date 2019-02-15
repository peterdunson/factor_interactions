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
      data_gen = generate_indep_model_notsparse
      out_name = paste("n",n,"_p",p,"_ind","_notsparse",sep="")
   }else if(sparse == 1){
      data_gen = generate_indep_model_sparse
      out_name = paste("n",n,"_p",p,"_ind","_sparse",sep="")
   }
}else if(type_model == 1){
   if(sparse == 0){
      data_gen = generate_corr_model_notsparse
      out_name = paste("n",n,"_p",p,"_corr","_notsparse",sep="")
   }else if (sparse == 1){
      data_gen = generate_corr_model_sparse
      out_name = paste("n",n,"_p",p,"_corr","_sparse",sep="")
   }
}


S = 10
err = err_pred = FR = err_beta = matrix(0,ncol = 8, nrow = S)
alpha = 0.05
TP_main = TN_main = TP_int = TN_int = matrix(0,ncol = 8, nrow = S)
acp_min = acp_max = acp_mean_st = matrix(0,ncol = 4,nrow = S)
delta_05 = delta_k = delta_10 = delta_50 = 0.2

for(s in 1:S){
   
   # generate the data
   data = data_gen(p = p,n = n, k = k_true)
   y = data$y; X = data$X; beta_true = data$beta_true; 
   Omega_true = data$Omega_true;
   y_test = data$y_test; X_test = data$X_test
   p = ncol(X)
   
   #Factor models
   nrun = 15000; burn = 10000; thin = 10; 
   
   
   gibbs_DL_05 = gibbs_DL(y, X ,nrun, burn, thin = 1, 
                          delta_rw = delta_05, epsilon_rw = 0.5,
                          a = 1/2, k = NULL)
   acp_mean = mean(gibbs_DL_05$acp)
   while(acp_mean > 0.3 | acp_mean < 0.2){
      if(acp_mean > 0.3){
         delta_05 = delta_05*2
      }else if(acp_mean < 0.2){
         delta_05 = delta_05*2/3
      }
      gibbs_DL_05 = gibbs_DL(y, X ,nrun, burn, thin = 1, 
                             delta_rw = delta_05, epsilon_rw = 0.5,
                             a = 1/2, k = NULL)
      acp_mean = mean(gibbs_DL_05$acp)
      print(acp_mean)
   }
   
   if(s == 1){
      delta_k = delta_10 = delta_50 = delta_05
   }
   
   
   gibbs_DL_k = gibbs_DL(y, X ,nrun, burn, thin = 1, 
                          delta_rw = delta_k, epsilon_rw = 0.5,
                          a = floor(log(p)*3), k = NULL)
   acp_mean = mean(gibbs_DL_k$acp)
   while(acp_mean > 0.3 | acp_mean < 0.2){
      if(acp_mean > 0.3){
         delta_k = delta_k*2
      }else if(acp_mean < 0.2){
         delta_k = delta_k*2/3
      }
      gibbs_DL_k = gibbs_DL(y, X ,nrun, burn, thin = 1, 
                             delta_rw = delta_k, epsilon_rw = 0.5,
                             a = floor(log(p)*3), k = NULL)
      acp_mean = mean(gibbs_DL_k$acp)
      print(acp_mean)
   }
   
   
   
   gibbs_CUSP_10 = gibbs_CUSP(y, X ,nrun, burn, thin = 1, 
                          delta_rw = delta_10, epsilon_rw = 0.5, 
                          k = NULL, alpha_prior = p*floor(log(p)*3)/10,
                          theta_inf = 0.05)
   acp_mean = mean(gibbs_CUSP_10$acp)
   while(acp_mean > 0.3 | acp_mean < 0.2){
      if(acp_mean > 0.3){
         delta_10 = delta_10*2
      }else if(acp_mean < 0.2){
         delta_10 = delta_10*2/3
      }
      gibbs_CUSP_10 = gibbs_CUSP(y, X ,nrun, burn, thin = 1, 
                                 delta_rw = delta_10, epsilon_rw = 0.5, 
                                 k = NULL, alpha_prior = p*floor(log(p)*3)/10,
                                 theta_inf = 0.05)
      acp_mean = mean(gibbs_CUSP_10$acp)
      print(acp_mean)
   }
   
   gibbs_CUSP_50 = gibbs_CUSP(y, X ,nrun, burn, thin = 1, 
                              delta_rw = delta_50, epsilon_rw = 0.5, 
                              k = NULL, alpha_prior = p*floor(log(p)*3)/2,
                              theta_inf = 0.05)
   acp_mean = mean(gibbs_CUSP_50$acp)
   while(acp_mean > 0.3 | acp_mean < 0.2){
      if(acp_mean > 0.3){
         delta_50 = delta_50*2
      }else if(acp_mean < 0.2){
         delta_50 = delta_50*2/3
      }
      gibbs_CUSP_50 = gibbs_CUSP(y, X ,nrun, burn, thin = 1, 
                                delta_rw = delta_50, epsilon_rw = 0.5, 
                                k = NULL, alpha_prior = p*floor(log(p)*3)/2,
                                theta_inf = 0.05)
      acp_mean = mean(gibbs_CUSP_50$acp)
      print(acp_mean)
   }

   
   # acp 
   acp_1 = gibbs_DL_05$acp;acp_2 = gibbs_DL_k$acp;acp_3 = gibbs_CUSP_10$acp;acp_4 = gibbs_CUSP_50$acp;
   acp_min[s,1] = min(acp_1/(nrun-burn));acp_max[s,1] = max(acp_1/(nrun-burn));acp_mean_st[s,1] = mean(acp_1/(nrun-burn))
   acp_min[s,2] = min(acp_2/(nrun-burn));acp_max[s,2] = max(acp_2/(nrun-burn));acp_mean_st[s,2] = mean(acp_2/(nrun-burn))
   acp_min[s,3] = min(acp_3/(nrun-burn));acp_max[s,3] = max(acp_3/(nrun-burn));acp_mean_st[s,3] = mean(acp_3/(nrun-burn))
   acp_min[s,4] = min(acp_4/(nrun-burn));acp_max[s,4] = max(acp_4/(nrun-burn));acp_mean_st[s,4] = mean(acp_4/(nrun-burn))
   
   # Competitors
   hiernet = Hiernet_fct(y, X, X_test, y_test)
   Family = FAMILY_fct(y, X, X_test, y_test)
   PIE = PIE_fct(y, X, X_test, y_test)
   RAMP = RAMP_fct(y, X, X_test, y_test)
   
   #Errors
   errors = compute_errors(hiernet,Family,PIE,RAMP,
                           y,y_test,Omega_true,beta_true,
                           gibbs_DL_05,gibbs_DL_k,
                           gibbs_CUSP_10,gibbs_CUSP_50)
   
   #Error 
   err[s,] = errors$err_insample[-1]
   err_pred[s,] = errors$err_pred[-1]
   FR[s,] = errors$FR_norm[-1]
   err_beta[s,] = errors$beta_MSE[-1]
   
   #TP and TNs
   rate_main2 = rate_recovery_maineff(gibbs_DL_05,gibbs_DL_k,alpha = alpha,beta_true = beta_true,
                                      hiernet$beta,Family$beta,PIE$beta,RAMP$beta)
   rate_int2 = rate_recovery_interactions(gibbs_DL_05,gibbs_DL_k,alpha = alpha,Omega_true=Omega_true,
                                          hiernet$Omega,Family$Omega,PIE$Omega,RAMP$Omega)
   rate_main1 = rate_recovery_maineff(gibbs_CUSP_10,gibbs_CUSP_50,alpha = alpha,beta_true = beta_true,
                                      hiernet$beta,Family$beta,PIE$beta,RAMP$beta)
   rate_int1 = rate_recovery_interactions(gibbs_CUSP_10,gibbs_CUSP_50,alpha = alpha,Omega_true=Omega_true,
                                          hiernet$Omega,Family$Omega,PIE$Omega,RAMP$Omega)
   
   TP_main[s,] = c(rate_main2$TP[1:2],rate_main1$TP)
   TN_main[s,] = c(rate_main2$TN[1:2],rate_main1$TN)
   TP_int[s,] = c(rate_int2$TP[1:2],rate_int1$TP)
   TN_int[s,] = c(rate_int2$TN[1:2],rate_int1$TN)
   
   print(s)
}

TP_main = TP_main[,c(5:8,1:4)]; TN_main = TN_main[,c(5:8,1:4)]
TP_int = TP_int[,c(5:8,1:4)]; TN_int = TN_int[,c(5:8,1:4)]
col_names = c("Hiernet","Family","Pie","RAMP","DL_05","DL_k","CUSP_10","CUSP_50")
colnames(err_beta) = colnames(err_pred) = colnames(err) = colnames(FR) = 
   colnames(TP_main) = colnames(TN_main) = colnames(TP_int) = colnames(TN_int) = col_names


list_res = list(
   TP_main = TP_main,
   TN_main = TN_main,
   TP_int = TP_int,
   TN_int = TN_int,
   MSE_beta = err_beta,
   err_pred = err_pred,
   err_test = err,
   FR = FR,
   acp_min = acp_min,
   acp_mean = acp_mean_st,
   acp_max = acp_max
)

results_dir = file.path("/work/sta790/ff31/results")
#dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list_res, file.path(results_dir, out_name))
#readRDS("results/n100_p10_ind_notsparse")






