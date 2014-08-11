# run SGLR
# example: NCI database with Mh dataset in CCLE

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
source("~/SGLR/bsSGLR_prior.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")#,"E","C")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    bsSGLR_prior_CCLE(k1,k2,bsNum = 100)
  }
}


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
