# Correlation NHANES
library(tidyverse)
library(reshape2)
library(lessR)
library(latex2exp)

# Metals 
load("~/factor_interactions/data/data_nhanes/nhanes_complete.RData")
Chem = subset(data_nhanes, select = -c(ID,TotChol,Creatinine,
                                         BMI,WAIST,Age,Gender,Race,Ratio_income_poverty))%>%
   data.matrix()

# Phthalates
load("~/gp/data/NHANES 2015/PhthalatesMetals_nhanes.RData")
Chem = subset(PhthalatesMetals_nhanes, select = -c(SEQN,Cholesterol,Creatinine,
                                                   BMI,Waist,Age,Sex,Race,Ratio_income))%>%
   data.matrix()

# Correlation matrix 
Chem = Chem[complete.cases(Chem),]
C = stats::cor(Chem)
colnames(C) = rownames(C) = c("Mono-n-butyl","Mono-isobutyl","Mono-ethyl","Mono-benzyl",
                              "Mono-cyclohexyl","Mono-carbox","Mono-hydrox","Mono-oxohexyl",
                              "Mono-hexyl","2,5-dich","2,4-dich")
C = corReorder(C)

C_plot = melt(C)
ggplot(C_plot, aes(x = Var2, y = Var1)) + 
   geom_tile(aes(fill=value), colour="grey20") + 
   scale_fill_gradient2(low = "#191970", high = "#800000", mid = "white") +
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
   ggtitle(TeX(""))

