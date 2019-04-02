#### Processing of results
library(tidyverse)
library(plotly)
results = readRDS("~/factor_interactions/results/long_a133.RDS")
df_chem = readRDS("~/factor_interactions/data/df_chem.RDS")

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
apply(results$beta_Z,2,mean)
apply(results$beta_Z,2,quant)
apply(results$Omega_bayes,c(2,3),mean)
apply(results$Omega_bayes,c(2,3),quant)
apply(results$Psi,c(2,3),mean)


d_ind = 4
plot(results$Omega_bayes[,d_ind,d_ind],ty="l")
d_ind = d_ind + 1
plot(results$beta_bayes[,2],ty="l")
plot(results$Omega_bayes[,d_ind,d_ind],ty="l")
plot(results$beta_Z[,3],ty="l")
plot(results$Omega[,1,5],ty="l")

#### plot covariance matrix for paper
sourceDirectory("~/factor_interactions/codes/post_processing/coverage_int")
cover = coverage_int(results$beta_bayes,results$Omega_bayes)
beta_hat = apply(results$beta_bayes,2,mean)
Omega_hat = apply(results$Omega_bayes,c(2,3),mean)
alpha_hat = apply(results$beta_Z,2,mean)
ind_b = which(cover[[1]] == 0); ind_i = which(cover[[2]] == 0) 
beta_hat[ind_b] = 0; Omega_hat[ind_i] = 0

#### plot interactions
# library(ggcorrplot)
# colnames(Omega_hat) = rownames(Omega_hat) = colnames(X)
# max_o = max(Omega_hat); min_o = min(Omega_hat)
# Omega_hat2 = (2*Omega_hat - (max_o + min_o))/(max_o - min_o)
# ggcorrplot(Omega_hat2,
#            outline.col = "white",
#            #ggtheme = ggplot2::theme_gray,
#            colors = c("#6D9EC1", "white", "#E46726"))


### competitors
df_chem = readRDS("~/factor_interactions/data/df_chem.RDS")
sourceDirectory("~/factor_interactions/codes/post_processing/compute_errors")
y = as.numeric(scale(df_chem$y)); X = as.matrix(scale(df_chem$X)); Z = df_chem$Z
X_test = X; y_test = y; Omega_true = Omega_hat; beta_true = beta_hat
hiernet = quiet(Hiernet_fct(y, X, X_test, y_test))
Family = quiet(FAMILY_fct(y, X, X_test, y_test))
PIE = PIE_fct(y, X, X_test, y_test)
RAMP = RAMP_fct(y, X, X_test, y_test)

errors = compute_errors(hiernet,Family,PIE,RAMP,
                        y,y_test-Z%*%alpha_hat,Omega_true,beta_true,
                        results)


