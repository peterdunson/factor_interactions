# -------- Ferrari Workflow With NHANES Data --------
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
library(tidyverse)
library(dplyr)
# ... Add other necessary libraries here

# ---- Source Functions ----
repo_path <- "/Users/peterdunson/factor_interactions"
source(file.path(repo_path, "codes/functions/Gibbs_DL.R"))
source(file.path(repo_path, "codes/functions/quiet.R"))
R.utils::sourceDirectory(file.path(repo_path, "codes/functions"))
# If you need coverage_y or compute_errors, source them as well



# ---- Load Libraries ----
library(tidyverse)
library(plyr)
library(mice)

# ---- Set Data Directory ----
repo_path <- "/Users/peterdunson/factor_interactions"  # adjust as needed

# ---- Load Data ----
load(file.path(repo_path, "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_cov_1516.RData"))
load(file.path(repo_path, "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_out_1516.RData"))
load(file.path(repo_path, "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_metals_1516.RData"))
load(file.path(repo_path, "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_phalates_pfas_1516.RData"))

# ---- Log-transform Chemicals ----
df_metals_log <- df_metals %>%
   dplyr::select(-SEQN) %>%
   log(., base = 10) %>%
   cbind(df_metals$SEQN, .) %>%
   mutate(SEQN = df_metals$SEQN) %>%
   dplyr::select(-"df_metals$SEQN")

df_phalates_pfas_log <- df_phalates_pfas %>%
   dplyr::select(-SEQN) %>%
   log(., base = 10) %>%
   cbind(df_phalates_pfas$SEQN, .) %>%
   mutate(SEQN = df_phalates_pfas$SEQN) %>%
   dplyr::select(-"df_phalates_pfas$SEQN")

# ---- Remove Non-normal Metals ----
df_metals_log_norm <- df_metals_log %>%
   dplyr::select(-c(LBXBCD, LBXBCR, LBXBGE, LBXBGM, LBXIHG, LBXTHG, URXUHG,
                    URXUCD, URXUMN, URXUSB, URXUAB, URXUAC, URXUAS3, URXUAS5,
                    URXUDMA, URXUMMA, URXUUR, URXUTU))

# ---- Select Normalized Phalates/PFAS ----
df_phalates_pfas_log_norm <- df_phalates_pfas_log %>%
   dplyr::select(c(LBXMFOS, LBXNFOA, LBXNFOS, LBXPFHS,
                   URXCNP, URXCOP, URXECP, URXHIBP, URXMBP, URXMEP,
                   URXMHH, URXMIB, URXMOH, URXMZP, SEQN))

# ---- Assemble Final Dataframe ----
chem_names <- c(colnames(df_metals_log_norm), colnames(df_phalates_pfas_log_norm))
df_out_analysis <- df_out %>% dplyr::select(SEQN, BMXBMI)
df <- join_all(list(df_cov, df_out_analysis, df_metals_log_norm, df_phalates_pfas_log_norm),
               by = 'SEQN', type = 'full')

# ---- Remove Rows With Missing Outcome ----
df <- df[!is.na(df$BMXBMI), ]

# ---- Create Design Matrix (Chemicals) ----
X <- df %>% dplyr::select(all_of(chem_names)) %>%
   transmute(
      cobalt = LBXBCO, copper_serum = LBXSCU, selenium_serum = LBXSSE,
      zinc_serum = LBXSZN, lead_blood = LBXBPB, selenium_blood = LBXBSE,
      manganese_blood = LBXBMN, barium_urine = URXUBA, cobalt_urine = URXUCO,
      cesium_urine = URXUCS, molybdenum_urine = URXUMO, lead_urine = URXUPB,
      tin_urine = URXUSN, strontium_urine = URXUSR, thallium_urine = URXUTL,
      sm_pfos = LBXMFOS, n_perfluorooctanoic = LBXNFOA, n_perfluorooctane = LBXNFOS,
      perfluorohexane = LBXPFHS, mono_carboxyisononyl = URXCNP,
      mono_carboxyisoctyl = URXCOP, mono_2_ethyl_5_carboxypentyl = URXECP,
      mhibp = URXHIBP, mono_n_butyl = URXMBP, mono_ethyl = URXMEP,
      mono_2_ethyl_5_hydroxyhexyl = URXMHH, mono_isobutyl = URXMIB,
      mono_2_ethyl_5_oxohexyl = URXMOH, mono_benzyl = URXMZP
   ) %>%
   select(-selenium_serum) # Exclude for exact match to Ferrari

# ---- Remove Rows Missing All Chemicals ----
X_na <- is.na(X)
ind_na_all <- which(apply(X_na, 1, mean) == 1)
X <- X[-ind_na_all, ] %>% scale() %>% as.data.frame()
y <- df$BMXBMI[-ind_na_all] %>% scale() %>% as.numeric()

# At this point:
#   X = scaled chemicals, as data.frame or matrix
#   y = scaled BMI


# ---- Load and Prepare NHANES Data (as before) ----
# ... Use your processing code to get X and y ...
# ... After all scaling, NA handling, etc. ...
# X: n x p matrix (chemicals, scaled)
# y: vector (scaled BMI)
X <- as.matrix(X)
y <- as.numeric(y)

# ---- Split Into Train/Test Sets ----
set.seed(2024)
n <- nrow(X)
test_prop <- 0.2
test_ind <- sample(seq_len(n), size = floor(test_prop * n))
X_train <- X[-test_ind, ]
y_train <- y[-test_ind]
X_test <- X[test_ind, ]
y_test <- y[test_ind]

# ---- Bayesian Factor Model ----
nrun <- 2000
burn <- 1500
thin <- 1
k_start <- 10

gibbs_DL_P <- gibbs_DL_Plam(
   y_train,
   X_train,
   nrun,
   burn,
   thin,
   delta_rw = 1,
   epsilon_rw = 0.5,
   a = k_start,
   k = k_start
)
gibbs_DL_k <- gibbs_DL_P

# ---- Competitors (HierNet, RAMP, etc) ----
hiernet <- quiet(Hiernet_fct(y_train, X_train, X_test, y_test))
RAMP <- RAMP_fct(y_train, X_train, X_test, y_test)
PIE <- RAMP
Family <- RAMP

# ---- Compute Errors (using Ferrari's compute_errors) ----
# Note: Omega_true and beta_true are unknown for real data,
# so set to NA or skip metrics that require them.
compute_errors <- function(hiernet, Family, PIE, RAMP, y, y_test, Omega_true, beta_true, ...) {
   # In-sample error
   err_insample = c(mean((y-mean(y))^2), mean(hiernet$err^2), mean(Family$err^2),
                    mean(PIE$err^2), mean(RAMP$err^2))
   # Predictive error
   err_pred = c(mean((y_test-mean(y_test))^2), mean(hiernet$err_pred^2), mean(Family$err_pred^2),
                mean(PIE$err_pred^2), mean(RAMP$err_pred^2))
   name = c("base case","hiernet","Family","PIE","RAMP")
   gibbs_list = list(...)
   t = length(gibbs_list)
   if(t > 0){
      for(k in 1:t){
         name_curr = paste("gibbs",k,sep="_")
         gibbs = gibbs_list[[k]]
         alpha_bayes_hat = mean(gibbs$alpha_bayes)
         betabayes_hat = apply(gibbs$beta_bayes,2,mean)
         Omegabayes_hat = apply(gibbs$Omega_bayes,c(2,3),mean)
         # in-sample
         err_factor = y - X %*% betabayes_hat -
            as.vector(diag(X %*% Omegabayes_hat %*% t(X))) - alpha_bayes_hat
         err_insample = c(err_insample,mean(err_factor^2))
         # prediction
         err_pred_factor = y_test - X_test %*% betabayes_hat -
            as.vector(diag(X_test %*% Omegabayes_hat %*% t(X_test))) - alpha_bayes_hat
         err_pred = c(err_pred,mean(err_pred_factor^2))
         name = c(name,name_curr)
      }
   }
   names(err_insample) = names(err_pred) = name
   return(list(err_insample = err_insample,
               err_pred = err_pred))
}

errors <- compute_errors(
   hiernet,
   Family,
   PIE,
   RAMP,
   y_train,
   y_test,
   Omega_true = matrix(NA, ncol(X_train), ncol(X_train)),
   beta_true = rep(NA, ncol(X_train)),
   gibbs_DL_P, gibbs_DL_k
)
print(errors)

# ---- Output Example ----
cat("In-sample MSE:\n")
print(errors$err_insample)
cat("\nTest/Prediction MSE:\n")
print(errors$err_pred)

# ---- If you want coverage and bias on test set ----
coverage_y <- function(y_test, X_test, gibbs) {
   beta_hat <- apply(gibbs$beta_bayes, 2, mean)
   Omega_hat <- apply(gibbs$Omega_bayes, c(2, 3), mean)
   alpha_hat <- mean(gibbs$alpha_bayes)
   y_pred <- X_test %*% beta_hat +
      diag(X_test %*% Omega_hat %*% t(X_test)) +
      alpha_hat
   bias <- mean(y_test - y_pred)
   coverage <- mean((y_test - y_pred)^2)
   list(coverage = coverage, bias = bias)
}
cov_y <- coverage_y(y_test, X_test, gibbs_DL_P)
print(cov_y)

# ---- Done! ----
