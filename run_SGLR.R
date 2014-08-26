# run SGLR
# example: NCI database with Mh dataset in CCLE
library(synapseClient)
synapseLogin()

# Step 01 : stepwise prior selection
source("~/SGLR/SGLR_prior.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")#,"E","C")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    SGLR_CCLE(k1,k2)
  }
}


# Step 01-2 : if you need only feature selection : bootstrapping 
KK1=c(1,2,4,6,9,11,12,13,14,17,39,41,75,82,84,87,90,94,95,103,108,110,122,123,124,126,127,129)
source_url("https://raw.githubusercontent.com/insockjang/SGLR/master/bsSGLR_prior_synapse.R")
PathwayName<-c("KEGG","BIOCARTA","NCI","GO_BP","GO_MF")
DataCombine<-c("E","C")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    bsSGLR_prior_Sanger(k1,k2,KK=KK1,bsNum = 100)
  }
  print(k1)
}


library(synapseClient)
library(devtools)
source_url("https://raw.githubusercontent.com/insockjang/SGLR/master/bsSGLR_prior_synapse.R")
KK1=c(1,2,4,6,9,11,12,13,14,17,39,41,75,82,84,87,90,94,95,103,108,110,122,123,124,126,127,129)
bsSGLR_prior_Sanger("GO_MF","E",KK=KK1,bsNum = 100,mcCoreNum=1)

# Step 02: prediction is separately run with restoreSGLR.R

source("~/SGLR/restoreSGLR.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")#,"E","C")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    restoreSGLR_CCLE(k1,k2)
  }
}

# Step 03 : Elastic Net with fixed sparsity with SGLR

source("~/SGLR/ENetSameNumber.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")#,"E","C")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    ENetSameNumber_CCLE(k1,k2)
  }
}

# Step 04 : Random Genes approach : null distribution for SGLR

source("~/SGLR/randomGene.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")#,"E","C")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    randomGene_CCLE(k1,k2)
  }
}

# Step 05 : Random pathway structure approach : null distribution for SGLR

source("~/SGLR/randomPathways.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")#,"E","C")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    randomPathway_CCLE(k1,k2)
  }
}
