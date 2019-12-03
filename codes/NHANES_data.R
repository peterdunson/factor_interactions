#### NHANES data #####

# Load libraries and functions
library(reshape2)
library(tidyverse)
library(latex2exp)
library(magick)
source("~/factor_interactions/codes/process_results/coverage_int.R")
source("~/factor_interactions/codes/functions/gibbs_DL.R")
source("~/factor_interactions/codes/functions/gibbs_DL_confounder_int.R")
source("~/factor_interactions/codes/functions/Competitors_fcts.R")
source("~/factor_interactions/codes/functions/quiet.R")
source("~/factor_interactions/codes/process_results/compute_errors_data.R")
source("~/factor_interactions/codes/process_results/mcrotfact.R")
source("~/factor_interactions/codes/process_results/clustalignplus.R")



# --- Read data --- #
load("~/factor_interactions/data/data_nhanes/nhanes_complete.RData")
head(data_nhanes)

# --- Chemicals in the study --- #
chemicals = subset(data_nhanes, select = c(URXMBP,URXMIB,URXMEP,URXMZP,URXMCP,
                                          URXECP,URXMHH,URXMOH,URXMHP,
                                          URXDCB,URX14D))
C = cor(chemicals)
eig = eigen(C)

#considere only complete data
data_complete = data_nhanes[complete.cases(data_nhanes),]

# factorize Race
data_complete$Race = as.factor(data_complete$Race)
Race = model.matrix(WAIST~Race - 1,data = data_complete)
data_complete$Race = Race

# --- create X and y for FIN --- #
X = subset(data_complete,select = 
              c(URXMBP,URXMIB,URXMEP,URXMZP,URXMCP,
                URXECP,URXMHH,URXMOH,URXMHP,
                URXDCB,URX14D,
                TotChol,Creatinine))%>%
   as.matrix
colnames(X) = c("Mono-n-butyl","Mono-isobutyl","Mono-ethyl","Mono-benzyl",
                "Mono-cyclohexyl","Mono-carbox","Mono-hydrox","Mono-oxohexyl",
                "Mono-hexyl","2,5-dich","2,4-dich","TotChol","Creatinine")
X = scale(X,scale = F)
Z = as.matrix(cbind(data_complete$Race,data_complete$Gender,
                    data_complete$age))
plot(eigen(cor(X))$values)
y = log10(as.numeric(data_complete$BMI))
y = as.numeric(scale(y))

set.seed(1)
ind = sample(1:nrow(X),200)
X_test = X[ind,]; X = X[-ind,]; 
Z_test = Z[ind,]; Z = Z[-ind,]; 
y_test = y[ind]; y = y[-ind]; 

# --- Run model --- #
nrun = 7500
burn = 5000
k = 7
gibbs = gibbs_DL_confounder_int(y, X, Z ,nrun, burn, thin = 5, 
                            delta_rw = 0.2, epsilon_rw = 0.5,
                            a = 1/2, k = k)

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
plot(gibbs$Omega_conf[,1,1],ty="l")
plot(gibbs$Omega_conf[,1,2],ty="l")

# --- plot omega_hat --- #
cov = coverage_int(gibbs$beta_bayes,gibbs$Omega_bayes)
p = ncol(X)
Omega_plot = matrix(0,p,p);
ind = which(cov[[2]] ==1 )
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
         axis.text.x = element_text(angle = 90, hjust = 1),
         legend.title = element_text(),
         plot.title = element_text(hjust = 0.5)) + 
   labs(fill = " ") + 
   ggtitle(TeX("Interactions Chemicals"))


# --- plot Phi_hat --- #
cov = coverage_int(gibbs$beta_bayes,gibbs$Omega_conf)
#Phi_plot = matrix(0,p,ncol(Phi_hat)); 
Phi_plot = Phi_hat
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
   ggtitle(TeX("Interactions X-Z"))

# --- Competitors --- #
W = cbind(X,Z)
W_test = cbind(X_test,Z_test)
hiernet = quiet(Hiernet_fct(y,W,W_test,y_test))
Family = quiet(FAMILY_fct(y,W,W_test,y_test))
RAMP = quiet(RAMP_fct(y,X,X_test,y_test))
PIE = quiet(PIE_fct(y,X,X_test,y_test))


# --- Compute errors --- #
errors = compute_errors_data(hiernet,Family,RAMP,RAMP,
                        y = y,y_test = y_test,X_test = X_test,Z_test = Z_test,
                        gibbs)

y_hat = X%*%beta_hat + 
   as.vector(diag(X%*%Omega_hat%*%t(X))) +
   alpha_hat + Z%*%beta_Z_hat +
   as.vector(diag(X%*%Phi_hat%*%t(Z)))


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
lambda_sample = gibbs$Lambda[,1:(p-2),]
lambda_sample = lapply(1:500, function(ind) lambda_sample[ind,,])
sample_mean = reduce(lambda_sample, `+`)/length(lambda_sample)
rownames(sample_mean) = colnames(X[,1:(p-2)])
rotated = mcrotfact(lambda_sample, method = "varimax", file = FALSE)
aligned = clustalignplus(rotated$samples, itermax = 500)


label = "Sample Mean"
SampleMean = cbind(melt(sample_mean), label)
label = "Aligned Sample Mean"
ProcessMean = Reduce("+", aligned)/length(aligned)
rownames(ProcessMean) = colnames(X[,1:(p-2)])
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
   labs(fill = " ")


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
names = as.character(colnames(X)[1:(p-2)])
label = "FIN"
q_sup = apply(gibbs$beta_bayes,2,function(x) quantile(x,probs = 0.95))
q_inf = apply(gibbs$beta_bayes,2,function(x) quantile(x,probs = 0.05))
gibbs_plot = data.frame(Values = beta_hat[1:(p-2)], 
                       Variables = names,
                       model = label,
                       q_inf = q_inf[1:(p-2)], q_sup = q_sup[1:(p-2)])
label = "PIE"
PIE_plot = data.frame(Values = PIE$beta[1:(p-2)], Variables = names,model = label, q_inf = PIE$beta[1:(p-2)]
                      ,q_sup = PIE$beta[1:(p-2)])
label = "RAMP"
RAMP_plot = data.frame(Values = RAMP$beta[1:(p-2)], Variables = names,model = label, q_inf = RAMP$beta[1:(p-2)]
                       ,q_sup = RAMP$beta[1:(p-2)])
label = "Family"
Family_plot = data.frame(Values = Family$beta[1:(p-2)], Variables = names,model = label,q_inf = Family$beta[1:(p-2)]
                         ,q_sup = Family$beta[1:(p-2)])
label = "HierNet"
hiernet_plot = data.frame(Values = hiernet$beta[1:(p-2)], Variables = names,model = label,q_inf = hiernet$beta[1:(p-2)]
                          ,q_sup = hiernet$beta[1:(p-2)])

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
   ggtitle(TeX("Estimated Main Effects"))+
   geom_errorbar(aes(ymin=q_inf, ymax=q_sup), width=.1)



