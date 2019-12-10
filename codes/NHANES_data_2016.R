#### Data analysis ####
library(tidyverse)
library(plyr)

# load data
load(file = "/work/sta790/ff31/factor_interactions/data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_cov_1516.RData")
load(file = "/work/sta790/ff31/factor_interactions/data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_out_1516.RData")
load(file = "/work/sta790/ff31/factor_interactions/data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_metals_1516.RData")
load(file = "/work/sta790/ff31/factor_interactions/data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_phalates_pfas_1516.RData")
# load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_cov_1516.RData")
# load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_out_1516.RData")
# load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_metals_1516.RData")
# load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_phalates_pfas_1516.RData")

# log trasform chemicals
# metals
df_metals_log = df_metals %>% 
   dplyr::select(-SEQN) %>% 
   log(., base = 10) %>% 
   cbind(df_metals$SEQN, .) %>%  #but the name of the column is not pretty
   mutate(SEQN = df_metals$SEQN) %>% 
   dplyr::select(- "df_metals$SEQN")#so we append new column named SEQN to the end of dataset
# phalates and pfs
df_phalates_pfas_log = df_phalates_pfas %>% 
   dplyr::select(-SEQN) %>% 
   log(., base = 10) %>% 
   cbind(df_phalates_pfas$SEQN, .) %>%  #but the name of the column is not pretty
   mutate(SEQN = df_phalates_pfas$SEQN) %>% #so we append new column named SEQN to the end of dataset
   dplyr::select(- "df_phalates_pfas$SEQN")

# get rid of the metals that do not respect the normality assumption
# function to plot
hist_chem = function(df,col_number){
   df[,col_number] %>%
      gather() %>% 
      ggplot(aes(value)) +
      facet_wrap(~ key, scales = "free") +
      geom_histogram()
}

# metals
#hist_chem(df_metals_log,1:12)
#hist_chem(df_metals_log,13:25)
#hist_chem(df_metals_log,26:34)
df_metals_log_norm = df_metals_log %>% 
   dplyr::select(-c(LBXBCD, LBXBCR, LBXBGE,
             LBXBGM, LBXIHG,LBXTHG,URXUHG,
             URXUCD,URXUMN,URXUSB,URXUAB,URXUAC,
             URXUAS3,URXUAS5,URXUDMA,URXUMMA,
             URXUUR,URXUTU
             ))

# phalates and pfs
#hist_chem(df_phalates_pfas_log,1:12)
#hist_chem(df_phalates_pfas_log,13:27)
df_phalates_pfas_log_norm = df_phalates_pfas_log %>% 
   dplyr::select(c(LBXMFOS,LBXNFOA,LBXNFOS,LBXPFHS,
            URXCNP,URXCOP,URXECP,URXHIBP,URXMBP,
            URXMEP,URXMHH,URXMIB,URXMOH,URXMZP,SEQN
   ))



# join data
# df_out_analysis = df_out %>% select(SEQN,BPXSY1,BPXDI1,
#                                     BMXWAIST,BMXBMI)

# save the names for later
chem_names = c(colnames(df_metals_log_norm),
               colnames(df_phalates_pfas_log_norm))
cov_names = colnames(df_cov)
col_Z = c("RIDAGEYR","RIDRETH1","INDFMPIR","DMDMARTL",
          "RIDEXPRG","LBXTC","URXUCR")

# join 
df_out_analysis = df_out %>% dplyr::select(SEQN,BMXBMI)
df = join_all(list(df_cov,
                   df_out_analysis,
                   df_metals_log_norm,
                   df_phalates_pfas_log_norm),
              by='SEQN', type='full')

# no NA's for the outcome
ind_na = df$BMXBMI %>% is.na %>% which(. == T)
df = df[-ind_na,]

# create the matrices for analysis
y = df$BMXBMI
X = df %>% dplyr::select(. , chem_names) %>% 
   transmute(
      cobalt = LBXBCO,copper_serum = LBXSCU,selenium_serum = LBXSSE,
      zinc_serum = LBXSZN, lead_blood = LBXBPB, selenium_blood = LBXBSE,
      manganese_blood = LBXBMN, barium_urine = URXUBA, cobalt_urine = URXUCO,
      cesium_urine = URXUCS,molybdenum_urine = URXUMO,lead_urine = URXUPB,
      tin_urine = URXUSN,strontium_urine = URXUSR,thallium_urine = URXUTL,
      sm_pfos = LBXMFOS, n_perfluorooctanoic = LBXNFOA, n_perfluorooctane = LBXNFOS,
      perfluorohexane = LBXPFHS, mono_carboxyisononyl = URXCNP,
      mono_carboxyisoctyl = URXCOP, mono_2_ethyl_5_carboxypentyl= URXECP,
      mhibp = URXHIBP, mono_n_butyl = URXMBP, mono_ethyl = URXMEP,
      mono_2_ethyl_5_hydroxyhexyl = URXMHH, mono_isobutyl = URXMIB,
      mono_2_ethyl_5_oxohexyl = URXMOH, mono_benzyl = URXMZP
   ) %>%
   select(- selenium_serum)

Z = df %>% 
   dplyr::select(. , col_Z) %>% 
   transmute(age = RIDAGEYR,gender = RIDRETH1, 
             race = RIDRETH1 ,chol = log(LBXTC),
             # ratio_income = INDFMPIR,marital_status = DMDMARTL,
             # albuminum = URXUMA,
             creatinine = log(URXUCR))

# observations that are missing all chemicals
X_na = X %>% is.na()
na_mean = apply(X_na, 1, mean)
ind_na_all = which(na_mean == 1)
X = X[-ind_na_all,] %>% scale() %>% as.data.frame()
Z = Z[-ind_na_all,] 
y = y[-ind_na_all] %>% scale()

# NAs
X_na = X %>% is.na()
Z_na = Z %>% is.na()
apply(X_na, 2, sum)
apply(Z_na, 2, sum)
# X_na %>% as.matrix() %>% image()
# Z_na %>% as.matrix() %>% image()

# impute the mean for chol, albuminum and creatinine
n = nrow(Z)
mean_Z = apply(Z, 2, function(x) mean(x,na.rm = T))
Z_pred = matrix(mean_Z, ncol = ncol(Z), nrow = n, byrow = T)
Z_imputed = Z; Z_imputed[Z_na] = Z_pred[Z_na]
Z_imputed = Z_imputed %>% scale() %>% as.matrix()

source("/work/sta790/ff31/factor_interactions/codes/functions/gibbs_DL_confounder_NA.R")

gibbs = gibbs_DL_confounder_NA(y, X, X_na, Z_imputed,
                       nrun = 100,burn = 50, k = 13)
results_dir = "/work/sta790/ff31/factor_interactions/"
saveRDS(gibbs,   file.path(results_dir, "metals_pfas_plam.rds"))
print("saved?")





