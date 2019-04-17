#### NHANES data #####

# read data
head(data_nhanes)
data_nhanes$RIDRETH1 = as.factor(data_nhanes$RIDRETH1)

# how many factors? 6 seems to be a good choice
phalates = subset(data_nhanes, select = c(URXMBP,URXMIB,URXMEP,URXMZP,URXMCP,
                                          URXECP,URXMHH,URXMOH,URXMHP))
C = cor(phalates)
#image(C)
eig = eigen(C)

# considere only complete data
data_complete = data_nhanes[complete.cases(data_nhanes),]

# factorize education
educ = model.matrix(BMXWAIST~RIDRETH1 - 1,data = data_complete)
data_complete$educ = educ







