library(tidyverse)
library(reshape2)
library(latex2exp)
library(glmnet)
library(RAMP)
library(hierNet)
library(FAMILY)
library(GIGrvg)
library(R.utils)
library(naniar)
library(plyr)

##### Source Functions from local git repo #####
# git clone https://github.com/fedfer/factor_interactions.git
source("/Users/felpo/factor_interactions/codes/functions/quiet.R")
sourceDirectory("/Users/felpo/factor_interactions/codes/functions")
source("/Users/felpo/factor_interactions/codes/process_results/compute_errors.R")
source("/Users/felpo/factor_interactions/codes/process_results/coverage_int.R")


#### Read the data in ####
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_cov_1516.RData")
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_out_1516.RData")
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_metals_1516.RData")
load(file = "data/data_nhanes_15-16/Rdata_nhanes_15-16/nhanes_phalates_pfas_1516.RData")

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
# pfas start with L and Phalathes start with U


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

pfas = X %>% select(sm_pfos,n_perfluorooctanoic,n_perfluorooctane,perfluorohexane)
C_pfas = cor(pfas, use = "pairwise.complete.obs")
C_pfas[lower.tri(C_pfas)] %>% mean()
complete.cases(pfas) %>% sum()

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
X_visNA = X[,c(19:28,15:18,1:14)]
visdat::vis_miss(X_visNA, warn_large_data = F)+
   theme(plot.margin=unit(c(1,3.5,2,2),"cm"))

# LOD
lod =  X %>% apply(., 2, function(x) min(x, na.rm = T))
chem_lod = X %>% apply(., 2, function(x) min(x, na.rm = T)) %>%
   rep(., times = nrow(X)) %>%
   matrix(., nrow = nrow(X), byrow = T)
X_lod = X
X_lod[X_na] = max(X, na.rm = T)
df_lod = (X_lod == chem_lod) %>% as.matrix() 
# df_lod %>% as.matrix() %>% image()
mean(df_lod)

# impute the mean for chol, albuminum and creatinine
n = nrow(Z)
mean_Z = apply(Z, 2, function(x) mean(x,na.rm = T))
Z_pred = matrix(mean_Z, ncol = ncol(Z), nrow = n, byrow = T)
Z_imputed = Z; Z_imputed[Z_na] = Z_pred[Z_na]
Z_imputed = Z_imputed %>% scale() %>% as.matrix()


## Correlation plot
Cor_plot = cor(X %>% as.matrix(),use = "pairwise.complete.obs")
colnames(Cor_plot) = rownames(Cor_plot) = colnames(X)
Cor_plot = melt(Cor_plot)
ggplot(Cor_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "#191970", high = "#006400", mid = "white") +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.ticks = element_blank(),
         #axis.text = element_blank(),
         axis.text.x = element_text(angle = 90, hjust = 1),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") + 
   ggtitle("Correlation Chemicals") + 
   scale_y_discrete(limits = rev(levels(Cor_plot$Var2))) 


## Correlation plot BRACES

bracketsGrob <- function(...){
   l <- list(...)
   e <- new.env()
   e$l <- l
   grid:::recordGrob(  {
      do.call(grid.brackets, l)
   }, e)
}
#?grid.brackets
b1 <- bracketsGrob(x1 = 0.03, y1 = 0.035, x2 = 0.03,y2 = 0.355, h=0.01, lwd=2, col="black")
b2 <- bracketsGrob(x1 = 0.03, y1 = 0.355, x2 = 0.03,y2 = 0.483, h=0.01, lwd=2, col="black")
b3 <- bracketsGrob(x1 = - 0.05, y1 = 0.483, x2 = - 0.05,y2 = 0.93, h=0.01, lwd=2, col="black")
b4 <- bracketsGrob(x1 =  0.03, y1 = 0.483, x2 =  0.03,y2 = 0.74, h=0.01, lwd=2, col="black")

Cor_plot = cor(X %>% as.matrix(),use = "pairwise.complete.obs")
colnames(Cor_plot) = rownames(Cor_plot) = colnames(X)
Cor_plot = melt(Cor_plot)
ggplot(Cor_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "#191970", high = "#006400", mid = "white") +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.text = element_blank(),
         #axis.text.x = element_text(angle = 90, hjust = 1, colour = "white"),
         #axis.text.y = element_text(hjust = 1, colour = "white"),
         axis.ticks = element_blank(),
         axis.ticks.length = unit(.85, "cm"),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") + 
   ggtitle("Correlation Chemicals") + 
   scale_y_discrete(limits = rev(levels(Cor_plot$Var2))) + 
   annotation_custom(b1)+
   annotation_custom(b2)+
   annotation_custom(b3)+
   annotation_custom(b4)+
   scale_x_discrete("",breaks=c(1,4),labels=c("Low types","High types") ) + 
   annotate("text", x=-1, y = 5, label= "phthalates",angle = 90)+
   annotate("text", x=-1, y = 12.5, label= "pfas",angle = 90)+
   annotate("text", x=-3.5, y = 21.5, label= "metals",angle = 90)+
   annotate("text", x=-1, y = 19, label= "metals urine",angle = 90)+
   coord_cartesian(ylim=c(0,30),xlim=c(0,30), clip="off")


# Eigen values
Cor_noNA = cor(X %>% as.matrix(),use = "pairwise.complete.obs")
ind = is.na(Cor_noNA)
Cor_noNA[ind] = mean(Cor_noNA, na.rm = T)
eig = eigen(Cor_noNA)
plot(eig$values)

# Eigen values
Cor_imp = cor(X_imputed %>% as.matrix(),use = "pairwise.complete.obs")
ind = is.na(Cor_imp)
Cor_imp[ind] = mean(Cor_imp, na.rm = T)
eig = eigen(Cor_imp)
plot(eig$values)
cumsum(eig$values)/sum(eig$values)

#### read simulation from cluster ####
# gibbs = readRDS("~/factor_interactions/metals_pfas.rds")
gibbs = readRDS("~/factor_interactions/metals_pfas.rds")
plot(gibbs$alpha_bayes)
plot(gibbs$beta_bayes[,1])

Omega_hat = apply(gibbs$Omega_bayes, c(2,3), mean)
Phi_hat = apply(gibbs$Omega_conf, c(2,3), mean)
alpha_hat = mean(gibbs$alpha_bayes)
beta_hat = apply(gibbs$beta_bayes,2,mean)
beta_Z_hat = apply(gibbs$beta_Z,2,mean)

colnames(Omega_hat) = rownames(Omega_hat) = colnames(X)
apply(gibbs$beta_bayes,2,quantile)
plot(gibbs$Omega_bayes[,3,3],ty="l")
plot(gibbs$Omega_conf[,1,1],ty="l")
plot(gibbs$Omega_conf[,1,2],ty="l")
plot(gibbs$Omega_conf[,1,4],ty="l")

# --- plot omega_hat --- #
cov = coverage_int(gibbs$beta_bayes,gibbs$Omega_bayes,inf = 0.005,sup = 0.995)
p = ncol(X)
Omega_plot = matrix(0,p,p);
ind = which(cov[[2]] ==1)
Omega_plot[ind] = Omega_hat[ind]
colnames(Omega_plot) = rownames(Omega_plot) = colnames(X)
Omega_plot = Omega_plot[1:(p-2),1:(p-2)]
label = "Interactions Chemicals"
Omega_plot = melt(Omega_plot)
Omega_plot = cbind(Omega_plot,label)
ggplot(Omega_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "#191970", high = "#006400", mid = "white") +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.ticks = element_blank(),
         #axis.text = element_blank(),
         axis.text.x = element_text(angle = 45, hjust = 1),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") + 
   ggtitle(TeX("Interactions Chemicals")) + 
   scale_y_discrete(limits = rev(levels(Omega_plot$Var2)))



# --- plot Phi_hat --- #
cov = coverage_int(gibbs$beta_bayes,gibbs$Omega_conf,inf = 0.005,sup = 0.995)
Phi_plot = matrix(0,p,ncol(Phi_hat)); 
ind = which(cov[[2]] ==1)
Phi_plot[ind] = Phi_hat[ind]
rownames(Phi_plot) = colnames(X)
colnames(Phi_plot) = colnames(Z)
#Phi_plot[cov[[2]]] = Phi_plot[cov[[2]]]
Phi_plot = melt(Phi_plot)
ggplot(Phi_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "#191970", high = "#006400", mid = "white") +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.ticks = element_blank(),
         #axis.text = element_blank(),
         axis.text.x = element_text(angle = 90, hjust = 1),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") +
   ggtitle(TeX("Interactions Chemicals-Covariates"))+ 
   scale_y_discrete(limits = rev(levels(Phi_plot$Var1)))

# --- Competitors --- #
library(mice)
source("/Users/felpo/factor_interactions/codes/functions/quiet.R")
source("/Users/felpo/factor_interactions/codes/functions/Competitors_fcts.R")
W = cbind(X,Z_imputed)
# W_imputed = mice(W,m=1,maxit=50,meth='pmm',seed=500)
# W_imputed = complete(W_imputed,1) 
# W_imputed = W_imputed %>% as.matrix()
# saveRDS(W_imputed,
#         "/Users/felpo/factor_interactions/data/data_nhanes_15-16/data_imputed.rds")

# hiernet = quiet(Hiernet_fct(as.numeric(y),W_imputed))
# Family = quiet(FAMILY_fct(as.vector(y),W_imputed))
# RAMP = quiet(RAMP_fct(y,W_imputed))
# PIE = quiet(PIE_fct(y,W_imputed))
#PIE = RAMP
#RAMP = hiernet
# hiernet = quiet(Hiernet_fct(y,W,W_test,y_test))
# Family = quiet(FAMILY_fct(y,W,W_test,y_test))
# RAMP = quiet(RAMP_fct(y,X,X_test,y_test))
# PIE = quiet(PIE_fct(y,X,X_test,y_test))
X_imputed = mice(X,m=1,maxit=50,meth='pmm',seed=500)
X_imputed = complete(X_imputed,1) 
X_imputed = X_imputed %>% as.matrix()
saveRDS(W_imputed,
        "/Users/felpo/factor_interactions/data/data_nhanes_15-16/chem_imputed.rds")

competitors = readRDS("~/factor_interactions/competitors.rds")
PIE = competitors$PIE
hiernet = competitors$hiernet
Family = competitors$Family
RAMP = competitors$RAMP


# --- Compute errors --- #
errors = compute_errors_data(hiernet,Family,RAMP,PIE,
                             y = y,y_test = y,X_test = W_imputed,Z_test = Z_imputed)

y_hat = gibbs$y_pred %>% apply(. , 2, mean)
(y - y_hat) %>% .^2 %>% mean()


# --- non zero interactions --- #
cov = coverage_int(gibbs$beta_bayes,gibbs$Omega_bayes)
Omega_hat01 = cov[[2]]
colnames(Omega_hat01) = rownames(Omega_hat01) = colnames(X)
Omega_plot = melt(Omega_hat01)
ggplot(Omega_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "white", high = "black") +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.ticks = element_blank(),
         #axis.text = element_blank(),
         axis.text.x = element_text(angle = 90, hjust = 1),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") +
   ggtitle(TeX("$\\Omega$"))


cov = coverage_int(gibbs$beta_bayes,gibbs$Omega_conf)
Phi01 = cov[[2]]
colnames(Phi01) = colnames(Z) 
rownames(Phi01) = colnames(X)
Phi_plot = melt(Phi01)
ggplot(Phi_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "white", high = "black") +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.ticks = element_blank(),
         axis.text = element_blank(),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") +
   ggtitle(TeX("$\\Omega$"))


# --- Postprocessing Lambda --- #
source("/Users/felpo/factor_interactions/codes/process_results/signer.R")
source("/Users/felpo/factor_interactions/codes/process_results/permuter.R")
source("/Users/felpo/factor_interactions/codes/process_results/mcrotfact.R")
source("/Users/felpo/factor_interactions/codes/process_results/permsignfact.R")
source("/Users/felpo/factor_interactions/codes/process_results/spclone.R")

lambda_sample = gibbs$Lambda[,1:p,]
lambda_sample = lapply(1:500, function(ind) lambda_sample[ind,,])
sample_mean = reduce(lambda_sample, `+`)/length(lambda_sample)
rownames(sample_mean) = colnames(X)
rotated = mcrotfact(lambda_sample, method = "varimax", file = FALSE)

# ClustAlign
#aligned = clustalignplus(rotated$samples, itermax = 500)

# MatchAlign
saveRDS(rotated$samples,"Lambda_sample.rds")
sprotated = spclone("Lambda_sample.rds", method = "BADFM", maxiter = 100, ncores = 4, tol = 1e-5)


label = "Sample Mean"
SampleMean = cbind(melt(sample_mean), label)
label = "Aligned Sample Mean"
#ProcessMean = Reduce("+", aligned)/length(aligned)
ProcessMean = Reduce("+", sprotated$samples)/length(sprotated$samples)
rownames(ProcessMean) = colnames(X)
ProcessMean = cbind(melt(ProcessMean), label)
ggdf = rbind(SampleMean, ProcessMean)
ggplot(ggdf, aes(x = Var2, y = Var1)) + 
   facet_grid(cols = vars(label)) +
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "#800000", high = "#006400", mid = "white") +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.ticks = element_blank(),
         #axis.text = element_blank(),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ")+ 
   scale_y_discrete(limits = rev(levels(ProcessMean$Var1)))


# --- GIF lambda --- #
#lambda_sample
#aligned
factl = array(unlist(lambda_sample), dim = c(11,7,500))
# factl = array(unlist(aligned), dim = c(11,7,500))
img = image_graph(600, 340, res = 96)
for(i in seq(1,500,by = 5)){
   curr = factl[,,i,  drop = F] / max(factl)
   rownames(curr) = colnames(X[,1:(p-2)])
   SampleMean = melt(curr)
   plot = ggplot(SampleMean, aes(x = Var2, y = Var1)) + 
      geom_tile(aes(fill=value), colour="grey20") + 
      scale_fill_gradient2(low = "#800000", high = "#006400", mid = "white") +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            #axis.text = element_blank(),
            legend.title = element_text(),
            plot.title = element_text(hjust = 0.5)) + 
      labs(fill = " ")
   print(plot)
   
   #lambda.add =  image_raster(plot)
   #lambda.image =  c(lambda.image, lambda.add)
}
lambda.animated = image_animate(img, fps = 20, dispose = "background")
image_write(lambda.animated, "Lambda_aligned.gif")

# For latex
# brew install ImageMagick
# magick convert -coalesce Lambda_unrotated.gif Lambda_unrotated.png


# Correlation plot
Cor_plot = cor(X[,1:(p-2)])
colnames(Cor_plot) = rownames(Cor_plot) = colnames(X[,1:(p-2)])
Cor_plot = melt(Cor_plot)
ggplot(Cor_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "#191970", high = "#006400", mid = "white") +
   theme(axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.ticks = element_blank(),
         #axis.text = element_blank(),
         axis.text.x = element_text(angle = 90, hjust = 1),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") + 
   ggtitle(TeX("Correlation matrix"))


# Plot main effects
#
names = as.character(colnames(X))
label = "FIN"
q_sup = apply(gibbs$beta_bayes,2,function(x) quantile(x,probs = 0.95))
q_inf = apply(gibbs$beta_bayes,2,function(x) quantile(x,probs = 0.05))
gibbs_plot = data.frame(Values = beta_hat, 
                        Variables = names,
                        model = label,
                        q_inf = q_inf, q_sup = q_sup)
label = "PIE"
PIE_plot = data.frame(Values = PIE$beta[1:p], Variables = names,model = label, q_inf = PIE$beta[1:p],
                      q_sup = PIE$beta[1:p])
label = "RAMP"
RAMP_plot = data.frame(Values = RAMP$beta[1:p], Variables = names,model = label, q_inf = RAMP$beta[1:p],
                       q_sup = RAMP$beta[1:p])
label = "Family"
Family_plot = data.frame(Values = Family$beta[1:p], Variables = names,model = label,q_inf = Family$beta[1:p],
                         q_sup = Family$beta[1:p])
label = "HierNet"
hiernet_plot = data.frame(Values = hiernet$beta[1:p], Variables = names,model = label,q_inf = hiernet$beta[1:p],
                          q_sup = hiernet$beta[1:p])

beta_plot = rbind(gibbs_plot,PIE_plot,RAMP_plot,Family_plot,hiernet_plot)
ggplot(beta_plot, aes(x = Variables, y = Values, color = model,
                      shape = model))+
   geom_point(size = 2.3)+
   theme(axis.title.x = element_blank(),
         #axis.title.y = element_blank(),
         #panel.grid.major = element_blank(),
         panel.border = element_blank(),
         #panel.background = element_blank(),
         #axis.ticks = element_blank(),
         #axis.text = element_blank(),
         text = element_text(size=17),
         axis.text.x = element_text(angle = 90, hjust = 1),
         #legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") + 
   ylab("Estimated Main Effects") +
   # ggtitle("Estimated Main Effects")+
   geom_errorbar(aes(ymin=q_inf, ymax=q_sup), width=.1)



