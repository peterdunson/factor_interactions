#### Processing of results
library(tidyverse)
library(plotly)
library(R.utils)
results = readRDS("~/factor_interactions/results/long_a133.rds")
df_chem = readRDS("~/factor_interactions/data/df_chem.RDS")

plot(results$sigmasq_st,type="l")

# --- quantile functions ---
quant = function(x){
   quantile(x,probs = c(0.05,0.95))
}
quant_inf = function(x){
   quantile(x,probs = c(0.05))
}
quant_sup = function(x){
   quantile(x,probs = c(0.95))
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
source("~/factor_interactions/codes/post_processing/coverage_int.R")
cover = coverage_int(results$beta_bayes,results$Omega_bayes)
beta_hat = apply(results$beta_bayes,2,mean)
Omega_hat = apply(results$Omega_bayes,c(2,3),mean)
alpha_hat = apply(results$beta_Z,2,mean)
ind_b = which(cover[[1]] == 0); ind_i = which(cover[[2]] == 0) 
beta_hat[ind_b] = 0; Omega_hat[ind_i] = 0


library(RColorBrewer)
colors1 = brewer.pal(n = 8, name = "RdBu")
colors1[4:5] = c("#FFFFFF","#FFFFFF")
colors = c("#B2182B","#F4A582","#FFFFFF","#92C5DE","#2166AC")
image(Omega_hat,col  = colors1)


library(plot.matrix)
plot(Omega_hat,col=heat.colors(10))
#### plot interactions
library(ggcorrplot)
colnames(Omega_hat) = rownames(Omega_hat) = colnames(X)
max_o = max(Omega_hat); min_o = min(Omega_hat)
Omega_hat2 = (2*Omega_hat - (max_o + min_o))/(max_o - min_o)
ggcorrplot(Omega_hat2,
            outline.col = "white",
            #ggtheme = ggplot2::theme_gray,
            colors = c("#6D9EC1", "white", "#E46726"))


library(seriation)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colors1 = brewer.pal(n = 8, name = "RdBu")
colors1[4:5] = c("#FFFFFF","#FFFFFF")

colnames(Omega_hat) = rownames(Omega_hat) = colnames(X)
VX = sqrt(attributes(df_chem$X)$`scaled:scale`); sd_y = sd(df_chem$y)
Omega_plot = diag(VX)%*%Omega_hat%*%diag(VX)*sd_y
#data = melt(Omega_hat)
data = melt(Omega_plot)
#colors2 = brewer.pal(n = 10, name = "RdBu")
#colors2 = c(colors1[1:5],rep("#FFFFFF",11),colors1[6:10])

breaks = round(quantile(Omega_plot,probs = c(0,0.1,0.15,0.5,0.6,0.85,0.9,1)),4)
ggplot(data, aes(x = Var2, y = Var1)) + 
   geom_raster(aes(fill=value)) + 
   scale_fill_gradientn(colours= colors1,na.value = "transparent",
                        breaks=breaks,
                        labels=breaks
                        #limits=c(min(Omega_hat),max(Omega_hat))
                        )+
   # scale_color_manual(values=colors1, 
   #                    breaks =breaks) +
   #theme_bw() + 
   theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                      axis.text.y=element_text(size=9),
                      plot.title=element_text(size=11))

### competitors
df_chem = readRDS("~/factor_interactions/data/df_chem.RDS")
sourceDirectory("~/factor_interactions/codes/post_processing")
source("~/factor_interactions/codes/functions/Competitors_fcts.R")
source("~/factor_interactions/codes/functions/quiet.R")
source("~/factor_interactions/codes/post_processing/compute_errors.R")


y = as.numeric(scale(df_chem$y)); X = as.matrix(scale(df_chem$X)); Z = df_chem$Z
y_test = as.numeric(scale(df_chem$y_test)); X_test = as.matrix(scale(df_chem$X_test));
Omega_true = Omega_hat; beta_true = beta_hat
hiernet = quiet(Hiernet_fct(y, X, X_test, y_test))
Family = quiet(FAMILY_fct(y, X, X_test, y_test))
PIE = PIE_fct(y, X, X_test, y_test)
RAMP = RAMP_fct(y, X, X_test, y_test)

errors = compute_errors(hiernet,Family,PIE,RAMP,
                        y,y_test,Omega_true,beta_true,
                        results)



## plot main effects
quant_inf = function(x){
   quantile(x,probs = c(0.05))
}
quant_sup = function(x){
   quantile(x,probs = c(0.95))
}

VX = df_chem$VX; sd_y = df_chem$sd_y
Scale = diag(VX*sd_y)
#beta = c(beta_hat,hiernet$beta,Family$beta,PIE$beta,RAMP$beta)
beta = c(#apply(results$beta_bayes,2,mean)%*%Scale,
         beta_hat%*%Scale,
         hiernet$beta%*%Scale,
         Family$beta%*%Scale,
         PIE$beta%*%Scale,
         RAMP$beta%*%Scale)
model = rep(c("BFM","hiernet","Family","PIE","RAMP"),each = 13)
names = rep(colnames(X),5)
q_inf = c(apply(results$beta_bayes,2,quant_inf)%*%Scale,rep(0,13*4))
q_sup = c(apply(results$beta_bayes,2,quant_sup)%*%Scale,rep(0,13*4))
df_beta = data_frame(beta = beta,
                     q_inf = q_inf,
                     q_sup = q_sup,
                     model = model,
                     covariates = names)


ggplot(df_beta,aes(covariates,beta))+
   geom_point(aes(shape = factor(model),colour =  factor(model)))+ 
   theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
         axis.text.y=element_text(size=9),
         plot.title=element_text(size=11))

ind = which(df_beta$covariates != "DDE_A" & df_beta$covariates != "TOT_CHOL")
df_beta_pcbs = df_beta[ind,]
ggplot(df_beta_pcbs,aes(covariates,beta,shape = factor(model),colour =  factor(model)))+
   geom_point()+
   #geom_errorbar(aes(ymin=q_inf, ymax=q_sup), width=.25)+
   theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
         axis.text.y=element_text(size=9),
         plot.title=element_text(size=11))

