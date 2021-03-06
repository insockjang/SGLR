bsSGLR_prior_CCLE<-function(pathwayName,dataCombine,KK=c(1:24),bsNum = 100,mcCoreNum = 32){
  ### DEMO Stepwise grouping Lasso
  require(predictiveModeling)
  require(synapseClient)
  require(devtools)
    
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/R5/myEnetModel1.R")
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/myData_CCLE_new.R")
  source_url("https://raw.githubusercontent.com/insockjang/SGLR/master/parallel_stepwiseDecision.R")
  
  ###################################################
  #### Load Pathways                             ####
  ###################################################
  GRAPHITE<-synGet("syn2135029")
  load(GRAPHITE@filePath)
  pathwayName<-toupper(pathwayName)
  if(is.element(pathwayName,"BIOCARTA")){
    allPathways <- BIOCARTA
  }
  if(is.element(pathwayName,"KEGG")){
    allPathways <- KEGG
  }
  if(is.element(pathwayName,"REACTOME")){
    allPathways <- REACTOME
  }
  if(is.element(pathwayName,"NCI")){
    allPathways <- NCI
  }    
  MSigDB<-synGet("syn2227979")
  load(MSigDB@filePath)
  if(is.element(pathwayName,"GO_BP")){
    allPathways <- MSigDB$C5.GO_BP
  }    
  if(is.element(pathwayName,"GO_MF")){
    allPathways <- MSigDB$C5.GO_MF
  }    
  
  ###################################################
  #### Load CCLE Molecular Feature Data from Synapse ####
  ###################################################
  dataSets<-myData_CCLE_new(dataCombine,"ActArea")
  
  #   require(graphite)
  groups=list()
  for(k in 1:length(allPathways)){
    a=allPathways[[k]]
    a1<-paste(a,"_expr",sep="")
    a2<-paste(a,"_copy",sep="")
    a3<-paste(a,"_mut",sep="")
    aa<-union(a1,union(a2,a3))
    groups[[k]]<-aa
  }
  
  
  for(kk in KK){
    filename = paste("~/SGLR_bs100_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/PriorIncorporated_bsDrug_",kk,".Rdata",sep = "")
    if(!file.exists(filename)){
      
      #########################################################################################################
      ######## Training and Testing data are scaled(normalized) vs. raw(unnormalized) #######################
      #########################################################################################################
      
      # data preprocessing for preselecting features
      filteredData<-filterPredictiveModelData(dataSets$featureData,dataSets$responseData[,kk,drop=FALSE],featureVarianceThreshold = 0.2)
      
      # filtered feature and response data
      filteredFeatureData <- filteredData$featureData
      filteredResponseData <- filteredData$responseData
      
      ## scale these data
      filteredFeatureDataScaled <- scale(filteredFeatureData)
      filteredResponseDataScaled <- scale(filteredResponseData)
      
      set.seed(2)
      # bootstrapping resampling
      bootIndices <- createResample(filteredFeatureDataScaled[,1],times = bsNum,list = TRUE)
      
      resultSTEP<-foreach(kkk = 1:length(bootIndices)) %dopar% {
        STEP<-parallel_stepwiseDecision(filteredFeatureDataScaled[bootIndices[[kkk]],], filteredResponseDataScaled[bootIndices[[kkk]]],groups,coreNum = mcCoreNum,100)
        return(STEP)
      }
      save(resultSTEP,file = filename)
    }
    
  }
}

bsSGLR_prior_Sanger<-function(pathwayName,dataCombine,KK=NA,bsNum = 100,mcCoreNum = 32){
  ### DEMO Stepwise grouping Lasso
  require(predictiveModeling)
  require(synapseClient)
  require(devtools)
  
  
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/R5/myEnetModel1.R")
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/myData_Sanger.R")
  source_url("https://raw.githubusercontent.com/insockjang/SGLR/master/parallel_stepwiseDecision.R")
  
  ###################################################
  #### Load Pathways                             ####
  ###################################################
  GRAPHITE<-synGet("syn2135029")
  load(GRAPHITE@filePath)
  pathwayName<-toupper(pathwayName)
  if(is.element(pathwayName,"BIOCARTA")){
    allPathways <- BIOCARTA
  }
  if(is.element(pathwayName,"KEGG")){
    allPathways <- KEGG
  }
  if(is.element(pathwayName,"REACTOME")){
    allPathways <- REACTOME
  }
  if(is.element(pathwayName,"NCI")){
    allPathways <- NCI
  }    
  MSigDB<-synGet("syn2227979")
  load(MSigDB@filePath)
  if(is.element(pathwayName,"GO_BP")){
    allPathways <- MSigDB$C5.GO_BP
  }    
  if(is.element(pathwayName,"GO_MF")){
    allPathways <- MSigDB$C5.GO_MF
  }    
  
  ###################################################
  #### Load CCLE Molecular Feature Data from Synapse ####
  ###################################################
  dataSets<-myData_Sanger(dataCombine,"IC50")
  
  #   require(graphite)
  groups=list()
  for(k in 1:length(allPathways)){
    a=allPathways[[k]]
    a1<-paste(a,"_expr",sep="")
    a2<-paste(a,"_copy",sep="")
    a3<-paste(a,"_mut",sep="")
    aa<-union(a1,union(a2,a3))
    groups[[k]]<-aa
  }
  
  for(kk in KK){
    filename = paste("~/SGLR_bs100_filterVar02/",dataCombine,"/Sanger/",pathwayName,"/PriorIncorporated_bsDrug_",kk,".Rdata",sep = "")
    if(!file.exists(filename)){
      
      #########################################################################################################
      ######## Training and Testing data are scaled(normalized) vs. raw(unnormalized) #######################
      #########################################################################################################
      
      # data preprocessing for preselecting features
      filteredData<-filterPredictiveModelData(dataSets$featureData,dataSets$responseData[,kk,drop=FALSE],featureVarianceThreshold = 0.2)
      
      # filtered feature and response data
      filteredFeatureData <- filteredData$featureData
      filteredResponseData <- filteredData$responseData
      
      ## scale these data
      filteredFeatureDataScaled <- scale(filteredFeatureData)
      filteredResponseDataScaled <- scale(filteredResponseData)
      
      set.seed(2)
      # bootstrapping resampling
      bootIndices <- createResample(filteredFeatureDataScaled[,1],times = bsNum,list = TRUE)
      
      resultSTEP<-foreach(kkk = 1:length(bootIndices)) %dopar% {
        STEP<-parallel_stepwiseDecision(filteredFeatureDataScaled[bootIndices[[kkk]],], filteredResponseDataScaled[bootIndices[[kkk]]],groups,coreNum = mcCoreNum,100)
        return(STEP)
      }
      save(resultSTEP,file = filename)
    }
    
  }
}

