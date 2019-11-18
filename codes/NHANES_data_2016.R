#### Data analysis ####
library(tidyverse)
library(plyr)

# load data
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_cov_1516.RData")
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_out_1516.RData")
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_metals_1516.RData")
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16//nhanes_phalates_pfas_1516.RData")

# log trasform chemicals
df_metals_log = df_metals %>% 
  select(-SEQN) %>% 
  log(., base = 10) %>% 
  cbind(df_metals$SEQN, .) %>%  #but the name of the column is not pretty
  mutate(SEQN = df_metals$SEQN) %>% #so we append new column named SEQN to the end of dataset
  select(- "df_metals$SEQN") #and delete column with ugly name

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
df_out_analysis = df_out %>% select(SEQN,BMXBMI)
df = join_all(list(df_cov,
                   df_out_analysis,
                   df_metals_log_norm,
                   df_phalates_pfas_log_norm), 
              by='SEQN', type='full')
df = df %>% select(-SEQN)




  




