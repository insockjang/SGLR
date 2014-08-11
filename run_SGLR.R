# run SGLR
source("~/SGLR/PriorIncorporatedLasso_CCLE.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")#,"E","C")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    PriorIncorporatedLasso_CCLE(k1,k2)
  }
}


