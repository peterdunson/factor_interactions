### Postprocessing results simulations
library(tidyverse)
library(stargazer)
library(plotly)

# Organize in a folder
# cd /Users/felpo/factor_interactions/results/array_jobs/
# mkdir n500_p25_sigmasq1_corr_notsparse
# mv n500_p25_sigmasq1_corr_notsparse_* n500_p25_sigmasq1_corr_notsparse/


folder = "n500_p25_sigmasq1_ind_sparse"
file_name = "n500_p25_sigmasq1_ind_sparse"


# create matrices for results
FR = matrix(0,nrow = 50, ncol = 5)
col_names = c("hiernet","Family","PIE","RAMP","FIN")
colnames(FR) = col_names
TP_main = TN_main = TP_int = TN_int = 
   MSE_beta = err_pred = err_train = FR
cov_pred = bias_pred = numeric(50)
zeros = c()
for(i in 1:50){
   out = paste("~/factor_interactions/results/array_jobs_cov/",folder,
               "/",folder,"_iter=",
               i,".rds", sep = "")
   out = paste("~/factor_interactions/results/array_jobs_cov/",folder,
               "/",file_name,"_iter=",
               i,".rds", sep = "")
   
   if(file.exists(out)){
      
      results_curr = readRDS(out) 
      TP_main[i,] = results_curr$TP_main[1:5]
      TN_main[i,] = results_curr$TN_main[1:5]
      TP_int[i,] = results_curr$TP_int[1:5]
      TN_int[i,] = results_curr$TN_int[1:5]
      MSE_beta[i,] = results_curr$MSE_beta[1:5]
      err_pred[i,] = results_curr$err_pred[1:5]
      err_train[i,] = results_curr$err_test[1:5]
      FR[i,] = results_curr$FR[1:5]
      cov_pred[i] = results_curr$coverage_pred
      bias_pred[i] = results_curr$bias_pred
      
      
   }else{
      zeros = c(zeros,i)
   }
   
}

if(is.null(zeros) == F){
   TP_main = TP_main[-zeros,];TN_main = TN_main[-zeros,]
   TP_int = TP_int[-zeros,];TN_int = TN_int[-zeros,]
   MSE_beta = MSE_beta[-zeros,];err_pred = err_pred[-zeros,]
   err_train = err_train[-zeros,];FR = FR[-zeros,]
   
}

results = list(TP_main = TP_main, TN_main = TN_main,
               TP_int = TP_int, TN_int = TN_int, 
               MSE_beta = MSE_beta,err_pred = err_pred, 
               err_test = err_train, FR = FR)


source("~/factor_interactions/codes/process_results/return_G.R")
G = return_G(results)
G

cov_pred %>% mean()
bias_pred %>% mean()






