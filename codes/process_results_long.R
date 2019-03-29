#### Processing of results
library(tidyverse)
library(plotly)
results = readRDS("~/factor_interactions/results/long.RDS")

y_pred_mat = results$y_pred
y_hat = apply(y_pred_mat,2,mean)

plot(results$sigmasq_st,type="l")

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

apply(results$beta_bayes,2,quant)
apply(results$beta_bayes,2,mean)
apply(results$beta_z,2,mean)
apply(results$beta_z,2,quant)
apply(results$Omega_bayes,c(2,3),mean)
apply(results$Omega_bayes,c(2,3),quant)

plot(results$beta_bayes[,7],ty="l")
plot(results$Omega_bayes[,5,10],ty="l")

