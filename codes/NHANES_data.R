#### NHANES data #####
library(reshape2)
library(tidyverse)
library(latex2exp)

# X[,5], log transform creatinine and cholesterol
# CHANGE DISTRIBUTION OF Y TO HAVE A T DIST
# ADD INTERACTIONS between Gender and Chemicals


# read data
load("~/factor_interactions/data/data_nhanes/nhanes_complete.RData")
head(data_nhanes)

# Chemicals in the study
chemicals = subset(data_nhanes, select = c(URXMBP,URXMIB,URXMEP,URXMZP,URXMCP,
                                          URXECP,URXMHH,URXMOH,URXMHP,
                                          URXDCB,URX14D))
C = cor(chemicals)
eig = eigen(C)

# considere only complete data
data_complete = data_nhanes[complete.cases(data_nhanes),]

# factorize Race 
data_complete$Race = as.factor(data_complete$Race)
Race = model.matrix(WAIST~Race - 1,data = data_complete)
data_complete$Race = Race

# create X and y for FIN
X = subset(data_complete,select = 
              c(URXMBP,URXMIB,URXMEP,URXMZP,URXMCP,
                URXECP,URXMHH,URXMOH,URXMHP,
                URXDCB,URX14D,
                TotChol,Creatinine))%>%
   as.matrix
X = scale(X,scale = F)
Z = as.matrix(cbind(data_complete$Race,data_complete$Gender,
                    data_complete$age))
plot(eigen(cor(X))$values)
# log of instogram looks much nicer
y = log(as.numeric(data_complete$BMI))
y = as.numeric(scale(y))

# Run model
source("~/factor_interactions/codes/functions/gibbs_DL.R")
source("~/factor_interactions/codes/functions/gibbs_DL_confounder_int.R")

nrun = 500
burn = 400
k = 7
gibbs = gibbs_DL_confounder_int(y, X, Z ,nrun, burn, thin = 1, 
                            delta_rw = 0.2, epsilon_rw = 0.5,
                            a = 1/k, k = k)

plot(gibbs$alpha_bayes)
plot(gibbs$beta_bayes[,1])

Omega_hat = apply(gibbs$Omega_bayes, c(2,3), mean)
Phi_hat = apply(gibbs$Omega_conf, c(2,3), mean)
alpha_hat = mean(gibbs$alpha_bayes)
beta_hat = apply(gibbs$beta_bayes,2,mean)
beta_Z_hat = apply(gibbs$beta_Z,2,mean)

colnames(Omega_hat) = rownames(Omega_hat) = colnames(X)
apply(gibbs$beta_bayes,2,quantile)
plot(gibbs$Omega_bayes[,1,1],ty="l")

# plot omega_hat
Omega_plot = melt(Omega_hat)
ggplot(Omega_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "#3d52bf", high = "#33961b", mid = "white") +
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


# plot Phi_hat
Phi_plot = melt(Phi_hat)
ggplot(Phi_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "#3d52bf", high = "#33961b", mid = "white") +
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
   ggtitle(TeX("$\\Phi$"))

# Competitors 
hiernet = quiet(Hiernet_fct(y,X))
Family = quiet(FAMILY_fct(y,X))
RAMP = quiet(RAMP_fct(y,X))
# PIE = RAMP

hist(RAMP$err - y)
hist(y_hat)
hist(y)

# Compute errors 
X_test = X; y_test = y - mean(gibbs$alpha_bayes) - Z%*%apply(gibbs$beta_Z,2,mean)
errors = compute_errors(hiernet,Family,RAMP,RAMP,
                        y,y_test,Omega_hat,apply(gibbs$beta_bayes, 2, mean),
                        gibbs)

y_hat = X%*%beta_hat + 
   as.vector(diag(X%*%Omega_hat%*%t(X))) +
   alpha_hat + Z%*%beta_Z_hat


# non zero interactions
source("~/factor_interactions/codes/post_processing/coverage_int.R")
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
         axis.text = element_blank(),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") +
   ggtitle(TeX("$\\Omega$"))


