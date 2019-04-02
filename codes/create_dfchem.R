### Create Dataset from Longeneker

###### Read Data from local git repo ##### 
#cluster
df_chem = read.table("/work/sta790/ff31/factor_interactions/data/finaldde2.txt",
                     header = T,sep = ",",na.strings = ".")
#local macbook
if (!exists("df_chem")){
   df_chem = read.table("~/factor_interactions/data/finaldde2.txt",
                        header = T,sep = ",",na.strings = ".")
}
exists("df_chem")


##### Remove Outliers
df_chem = df_chem[-1861,]
gest_day = as.numeric(df_chem$GESTDAY)
ind = which(gest_day > 340)
df_chem = df_chem[-ind,]
pcbs = subset(df_chem,select = c( P028_A1,P052_A1,P074_A1,
                                  P105_A1,P118_A1,P153_A1,P170_A1,
                                  P138_A1,P180_A1,P194_A1,P203_A1))
totpcb = apply(pcbs,1,sum)
ind = which(totpcb > 10)
df_chem = df_chem[-ind,]

###### create matrix y, X, Z
mylogit = glm(PRETERM ~ (DDE_A + P028_A1 + P052_A1 + P074_A1 +
                            P105_A1 + P118_A1+P153_A1+P170_A1+
                            P138_A1+ P180_A1+ P194_A1+P203_A1+ TOT_CHOL)^2 + 
                 RACEC1 + V_SMKNOW + V_SEINDX + V_MHGT + TRIGLYC +
                 V_MAGE + BMICAT, data = df_chem, family = "binomial")


y = scale(as.numeric(df_chem$GESTDAY))
X = scale(model.matrix(mylogit)[,c(2:14)])
Z = scale(model.matrix(mylogit)[,c(15:21)])

VX = attributes(X)$`scaled:scale`
sd_y = attributes(y)$`scaled:scale`

set.seed(1)
ind = sample(1:length(y),100)
y_test = y[ind]; y = y[-ind,]
X_test = X[ind,]; X = X[-ind,]
Z_test = Z[ind,]; Z = Z[-ind,]

df_list = list(y = y, y_test = y_test,
               X = X, X_test = X_test,
               Z = Z, Z_test = Z_test,
               VX = VX, sd_y = sd_y)
#saveRDS(df_list, file.path("~/factor_interactions/data/df_chem.RDS"))

XX = rbind(df_list$X,df_list$X_test) 
XX = XX[-1,-1]
XX = XX[-12,-12]
mean(abs(cor(XX))[lower.tri(cor(XX))])
