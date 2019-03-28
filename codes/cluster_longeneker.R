######## CLUSTER JOB ####### 
#### Load all the libraries
library(BayesLogit)
library(MASS)
library(fields)
library(mvtnorm)
library(microbenchmark)
library(R.utils)
library(psych)
library(hierNet)
library(bkmr)
library(FAMILY)
library(PIE)
library(tidyverse)
library(glmnet)
library(RAMP)


##### Source Functions from local git repo #####
# git clone https://github.com/fedfer/gp.git
# git pull when in the folder
sourceDirectory("/work/sta790/ff31/gp/functions")
sourceDirectory("~/gp/functions")
sourceDirectory("/work/sta790/ff31/gp/generate_data")
sourceDirectory("~/gp/generate_data")


##### Argument of R script ######
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
   stop("At least one argument must be supplied", call.=FALSE)
}
nrun = as.numeric(args[1]); 
burn = as.numeric(args[2]); 
thin = as.numeric(args[3]);

###### Read Data from local git repo ##### 
#cluster
df_chem = read.table("/work/sta790/ff31/gp/data/finaldde2.txt",
                     header = T,sep = ",",na.strings = ".")
#local macbook
if (!exists("df_chem")){
   df_chem = read.table("~/gp/data/finaldde2.txt",
                        header = T,sep = ",",na.strings = ".")
}
exists("gibbs_binary_gp_v3")
exists("df_chem")
df_chem = df_chem[-1861,]


###### create matrix y, X, Z
mylogit = glm(PRETERM ~ (DDE_A + P028_A1 + P052_A1 + P074_A1 +
                            P105_A1 + P118_A1+P153_A1+P170_A1+
                            P138_A1+ P180_A1+ P194_A1+P203_A1+ TOT_CHOL)^2 + 
                 RACEC1 + V_SMKNOW + V_SEINDX + V_MHGT + TRIGLYC +
                 V_MAGE + BMICAT, data = df_chem, family = "binomial")

y = as.numeric(df_chem$PRETERM)
X = scale(model.matrix(mylogit)[,c(2:14)])
Z = scale(model.matrix(mylogit)[,c(15:21)])


###### Run Algorithm 
### Parameters
tau = 0.08; c = 0.2; delta_sigf = 0.5
dnorm(0.03,0,(c*tau))
dnorm(0.03,0,(tau))

res = gibbs_binary_gp_v3(y, X , Z = Z,
                         nrun = 1500, burn = 1000,thin = 1,
                         tau = 0.3,c = 2,delta_sigf = 0.5, exp_C = 2,
                         heredity = "strong")


####### Save results in cluster folder
results_dir = file.path("/work/sta790/ff31/gp/results_0322")
#dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(res,   file.path(results_dir, "chain5.rds"))