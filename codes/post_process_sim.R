### Postprocessing results simulations
library(tidyverse)
library(plotly)
#results = readRDS("~/factor_interactions/results_fact/n100_p10_corr_notsparse.RDS")
#results = readRDS("~/factor_interactions/results_fact/n100_p10_corr_sparse.RDS")
#results = readRDS("~/factor_interactions/results_fact/n100_p10_ind_notsparse.RDS")
#results = readRDS("~/factor_interactions/results_fact/n100_p10_ind_sparse.RDS")

# --- quantile functions ---
quant = function(x){
   quantile(x,probs = c(0.05,0.95))
}
quant_inf = function(x){
   quantile(x,probs = c(0.25))
}
quant_sup = function(x){
   quantile(x,probs = c(0.75))
}

apply(results$MSE_beta,2,mean)
apply(results$FR,2,mean)
apply(results$err_pred,2,mean)
apply(results$err_test,2,mean)
apply(results$TP_main,2,mean)
apply(results$TN_main,2,mean)
apply(results$TP_int,2,mean)
apply(results$TN_int,2,mean)


apply(TP_main,2,mean)
apply(TN_main,2,mean)
apply(TP_int,2,mean)
apply(TN_int,2,mean)

apply(err,2,sd)
min_err = apply(err,1, function(x) (x == min(x)) )
min_err = apply(min_err,1,sum)

apply(err_pred,2,mean)
apply(err_pred,2,sd)
min_err_pred = apply(err_pred,1, function(x) (x == min(x)) )
min_err_pred = apply(min_err_pred,1,sum)

apply(FR,2,mean)
apply(FR,2,sd)
min_FR = apply(FR,1, function(x) (x == min(x)) )
min_FR = apply(min_FR,1,sum)

apply(err_beta,2,mean)
apply(err_beta,2,sd)
min_beta = apply(err_beta,1, function(x) (x == min(x)) )
min_beta = apply(min_beta,1,sum)

G = rbind(apply(err,2,mean),min_err/100,
          apply(err_pred,2,mean),min_err_pred/100,
          apply(FR,2,mean),min_FR/100,
          apply(err_beta,2,mean),min_beta/100)
rownames(G) = c("MSE","","MSE prediction","",
                "Frobenious norm","","MSE beta","")

stargazer(G)



