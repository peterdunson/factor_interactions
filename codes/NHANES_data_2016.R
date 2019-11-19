#### Data analysis ####
library(tidyverse)
library(plyr)

# load data
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_cov_1516.RData")
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_out_1516.RData")
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_metals_1516.RData")
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16//nhanes_phalates_pfas_1516.RData")

# log trasform chemicals
# metals
df_metals_log = df_metals %>% 
  select(-SEQN) %>% 
  log(., base = 10) %>% 
  cbind(df_metals$SEQN, .) %>%  #but the name of the column is not pretty
  mutate(SEQN = df_metals$SEQN) %>% #so we append new column named SEQN to the end of dataset
  select(- "df_metals$SEQN") #and delete column with ugly name
# phalates and pfs
df_phalates_pfas_log = df_phalates_pfas %>% 
   select(-SEQN) %>% 
   log(., base = 10) %>% 
   cbind(df_phalates_pfas$SEQN, .) %>%  #but the name of the column is not pretty
   mutate(SEQN = df_phalates_pfas$SEQN) %>% #so we append new column named SEQN to the end of dataset
   select(- "df_phalates_pfas$SEQN")

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
   select(-c(LBXBCD, LBXBCR, LBXBGE,
             LBXBGM, LBXIHG,LBXTHG,URXUHG,
             URXUCD,URXUMN,URXUSB,URXUAB,URXUAC,
             URXUAS3,URXUAS5,URXUDMA,URXUMMA,
             URXUUR,URXUTU
             ))

# phalates and pfs
#hist_chem(df_phalates_pfas_log,1:12)
#hist_chem(df_phalates_pfas_log,13:27)
df_phalates_pfas_log_norm = df_phalates_pfas_log %>% 
   select(c(LBXMFOS,LBXNFOA,LBXNFOS,LBXPFHS,
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
          "RIDEXPRG","LBXTC","URXUMA","URXUCR")

# join 
df_out_analysis = df_out %>% select(SEQN,BMXBMI)
df = join_all(list(df_cov,
                   df_out_analysis,
                   df_metals_log_norm,
                   df_phalates_pfas_log_norm),
              by='SEQN', type='full')

# no NA's for the outcome
ind_na = df$BMXBMI %>% is.na 
df = df[-ind_na,]

# create the matrices for analysis
y = df$BMXBMI
X = df %>% select(. , chem_names) %>% 
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
   )

Z = df %>% 
   select(. , col_Z) %>% 
   transmute(age = RIDAGEYR,gender = RIDRETH1, 
             race = RIDRETH1 ,chol = LBXTC,
             # ratio_income = INDFMPIR,marital_status = DMDMARTL,
             albuminum = URXUMA,creatinine = URXUCR)

# observations that are missing all chemicals
X_na = X %>% is.na()
na_mean = apply(X_na, 1, mean)
ind_na_all = which(na_mean == 1)
X = X[-ind_na_all,]
Z = Z[-ind_na_all,]
y = y[-ind_na_all]

# NAs
X_na = X %>% is.na()
Z_na = Z %>% is.na()
apply(X_na, 2, sum)
apply(Z_na, 2, sum)
# X_na %>% as.matrix() %>% image()
# Z_na %>% as.matrix() %>% image()

# impute the mean for chol, albuminum and creatinine




