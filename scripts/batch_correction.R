library(statTarget)

setwd("C:/Users/XXX/ForStatTarget") # set your working directory

########### preprocessing #################
## 1) batch correction and imputation using stattarget ##
## files for statTarget for drift correction
pheno <- "ForStatTarget/meta4samples.csv" 
dfile <- "ForStatTarget/data4samples.csv"
labels <- read.csv("ForStatTarget/labels4samples.csv", row.names = "SampleName")

## MLmethod = 'QCRLSC' or MLmethod = 'QCRFSC' (default)
statTarget::shiftCor(pheno, dfile,  QCspan = 0.25, Frule = 0.5, 
                     degree = 2,imputeM = "KNN", ntree=500, coCV = 50)
### conclusion of batch correction and imputation through statTarget