#### Load Libraries #####
library(mvtnorm)
library(MASS)
library(beepr)
library(psych)
library(bayesSurv)
library('PIE')
library(R.utils)
library('glmnet')
library(RAMP)
library(hierNet)
library(FAMILY)
library(RCurl)
library(stargazer)
library(R.utils)
library(GIGrvg)

##### Source Functions from local git repo #####
# git clone https://github.com/fedfer/gp.git
# git pull when in the folder
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/functions")
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/generate_data")
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/post_processing")
#source("~/factor_interactions/codes/functions/gibbs_DL_confounder.R")
sourceDirectory("~/factor_interactions/codes/generate_data")
sourceDirectory("~/factor_interactions/codes/post_processing")
exists("generate_indep_model_notsparse")



##### Argument of R script ######
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
   stop("At least one argument must be supplied", call.=FALSE)
}
nrun = as.numeric(args[1]); 
burn = as.numeric(args[2]); 
thin = as.numeric(args[3]);
a = 1/as.numeric(args[4]);
if (length(a) == 0){
   a = 1/13
}
if (length(nrun) == 0){
   nrun = 1000
   burn = 500
   thin = 1
}

###### Read Data from local git repo ##### 
#cluster
df_chem = readRDS("/work/sta790/ff31/factor_interactions/data/df_chem.RDS")
#local macbook
if (!exists("df_chem")){
   df_chem = readRDS("~/factor_interactions/data/df_chem.RDS")
}
exists("df_chem")


###### Run Algorithm 
### Parameters
y = as.numeric(scale(df_chem$y)); X = as.matrix(scale(df_chem$X)); Z = df_chem$Z
delta_05 = 0.126749
res = gibbs_DL_confounder(y, X, Z, nrun, burn, thin = thin,
                          delta_rw = delta_05, epsilon_rw = 0.5,
                          a = a, k = 7)
# res = gibbs_DL(y, X, nrun, burn, thin = thin, 
#                           delta_rw = delta_05, epsilon_rw = 0.5,
#                           a = a, k = NULL)


####### Save results in cluster folder
results_dir = file.path("/work/sta790/ff31/factor_interactions/results")
#dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(res,   file.path(results_dir, "long_01_02.rds"))
