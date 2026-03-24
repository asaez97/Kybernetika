setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Sin balancear---------------------
source("fitTANGauss.R")
load("../DataSets/foldC.RData")

library(doParallel)
cores = detectCores()
cl <- makeCluster(10)
registerDoParallel(cl)
i = 1
MIGauss = foreach(i = 1:10,.packages = c("infotheo"))%dopar%{
  MI_tan_gauss(data = folds[[i]]$Training,target = "DUP.")
}
stopCluster(cl)
save(MIGauss,file="../DataSets/MIGauss.RData")
## Balanceados---------------------------------
source("fitTANGauss.R")
load("../DataSets/balancedMoP40.RData")
library(doParallel)
cores = detectCores()
cl <- makeCluster(10)
registerDoParallel(cl)
i = 1
MIGaussB = foreach(i = 1:10,.packages = c("infotheo"))%dopar%{
  MI_tan_gauss(data = balancedFold[[i]]$Training,target = "DUP.")
}
stopCluster(cl)
save(MIGaussB,file="../DataSets/MIGaussB40.RData")

