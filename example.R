# Set working directory----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Install packages------------
# install.packages("bnlearn")
# install.packages("MASS")
# install.packages("infotheo")
library(bnlearn)
# Generate data------------------

# Class variable
set.seed(1)
C = factor(sample(c("0","1"),100,T,prob = c(0.65,0.35)))
data = data.frame(C = C)
nC = table(C)
# Features variables
set.seed(12)
means = c(0,2,-3,1)
cor_matrix0 = matrix(c(1,0.9,0.4,-0.7,
                      0.9,1,0.5,-0.6,
                      0.4,0.5,1,-0.3,
                      -0.7,-0.6,-0.3,1),
                    nrow = 4,byrow = T)
data_0 = MASS::mvrnorm(n = nC[1],means,cor_matrix0,empirical = T)
set.seed(13)
means = c(0,3,-3,1)
cor_matrix1 = matrix(c(1,0.8,0.2,-0.75,
                       0.8,1,0.6,-0.7,
                       0.2,0.6,1,-0.35,
                       -0.75,-0.7,-0.35,1),
                     nrow = 4,byrow = T)
data_1 = MASS::mvrnorm(n = nC[2],means,cor_matrix1,empirical = T)

#continuous feature variables
data$Z1[C=="0"] = data_0[,1]
data$Z1[C=="1"] = data_1[,1]
data$Z2[C=="0"] = data_0[,2]
data$Z2[C=="1"] = data_1[,2]
# Discrete features variables
data$Y1[C=="0"] = ifelse(data_0[,3]>-3,"1","0")
data$Y2[C=="0"] = ifelse(data_0[,4]>1,"1","0")
data$Y1 = as.factor(data$Y1)
data$Y1[C=="1"] = ifelse(data_1[,3]>-3,"1","0")
data$Y2[C=="1"] = ifelse(data_1[,4]>1,"1","0")
data$Y2 = as.factor(data$Y2)

# Save the dataset
save(data,file = "data.Rda")
# Conditional mutual information-------------
source("fitTANGauss.R")
MI = MI_tan_gauss(data = data,target = "C")
MI

tan = fit_tan_g(target = "C",data = data)
graphviz.plot(tan)

# TAN with Y1 as root
tan2 = fit_tan_g(target = "C",data = data,root = "Y1",mutualInfoCond = MI)
graphviz.plot(tan2)

# TAN model with only continuous feature variables
tan4 = fit_tan_g(target = "C",data = data[,c("C","Z1","Z2")],root = "Z1")
graphviz.plot(tan4)