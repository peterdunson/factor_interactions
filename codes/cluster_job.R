######## FACTOR MODEL FOR CLUSTER ####

#### Load Libraries #####
library(mvtnorm)
library(MASS)
library(beepr)
library(psych)
library(bayesSurv)
#library('PIE')
library('glmnet')
library(RAMP)
library(hierNet)
library(FAMILY)
library(RCurl)
library(stargazer)
library(R.utils)

##### Source Functions from local git repo #####
# git clone https://github.com/fedfer/factor_interactions.git
sourceDirectory("~/factor_interactions/codes/functions")


##### Argument of R script ######

args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
   stop("At least one argument must be supplied", call.=FALSE)
}

n = args[1]; p = args[2]
# if type_model = 0 --> indep model; othw factor model with high correlation
# if sparse = 0 --> not sparse Omega; othw sparse omega
type_model = args[3]; sparse = args[4]

if(type_model == 0){
   if(sparse == 0){
      #data_gen = generate_factor_model
      out_name = paste("n",n,"_p",p,"_ind","_sparse",sep="")
   }else if(sparse == 1){
      #data_gen = generate_factor_model
      out_name = paste("n",n,"_p",p,"_ind","_notsparse",sep="")
   }
}else if(type_model == 1){
   if(sparse == 0){
      #data_gen = generate_factor_model
      out_name = paste("n",n,"_p",p,"_cor","_sparse",sep="")
   }else if (sparse == 1){
      #data_gen = generate_factor_model
      out_name = paste("n",n,"_p",p,"_cor","_notsparse",sep="")
   }
}














