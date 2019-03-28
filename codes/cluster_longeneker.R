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
# git clone https://github.com/fedfer/gp.git
# git pull when in the folder
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/functions")
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/generate_data")
sourceDirectory("/work/sta790/ff31/factor_interactions/codes/post_processing")
sourceDirectory("~/factor_interactions/codes/functions")
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

###### Read Data from local git repo ##### 
#cluster
df_chem = read.table("/work/sta790/ff31/factor_interactions/data/finaldde2.txt",
                     header = T,sep = ",",na.strings = ".")
#local macbook
if (!exists("df_chem")){
   df_chem = read.table("~/factor_interactions/data/finaldde2.txt",
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
delta_05 = 0.2
a = 1/2

res = gibbs_DL_confounder(y, X, Z, nrun, burn, thin = thin, 
                             delta_rw = delta_05, epsilon_rw = 0.5,
                             a = a, k = NULL)


####### Save results in cluster folder
results_dir = file.path("/work/sta790/ff31/results")
#dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(res,   file.path(results_dir, "long.rds"))
