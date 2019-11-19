### Postprocessing results simulations
library(tidyverse)
library(stargazer)
#library(plotly)
#results = readRDS("~/factor_interactions/results_fact3/n500_p25_ind_sparse.RDS")
#results = readRDS("~/factor_interactions/results_fact_ind/n500_p25_ind_sparse.RDS")

#results = readRDS("~/factor_interactions/results_fact/n100_p10_corr_sparse.RDS")
#results = readRDS("~/factor_interactions/results_fact/n100_p10_ind_notsparse.RDS")
#results = readRDS("~/factor_interactions/results_fact/n100_p10_ind_sparse.RDS")

results$TP_main
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

source("~/factor_interactions/codes/post_processing/return_G.R")
G = return_G(results)

results = readRDS("~/factor_interactions/results_fact2/n500_p25_corr_sparse.RDS")
G1 = return_G(results)
results = readRDS("~/factor_interactions/results_fact3/n500_p25_corr_notsparse.RDS")
G2 = return_G(results)
results = readRDS("~/factor_interactions/results_fact_ind/n500_p25_ind_sparse.RDS")
G3 = return_G(results)
results = readRDS("~/factor_interactions/results_fact_ind/n500_p25_ind_notsparse.RDS")
G4 = return_G(results)
G = rbind(G1,G2,G3,G4)
stargazer(G[,1:5])



