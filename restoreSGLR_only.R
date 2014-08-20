restoreSGLR_only_CCLE<-function(pathwayName,dataCombine,KK = c(1:24)){
  ### DEMO Stepwise grouping Lasso
  require(predictiveModeling)
  require(synapseClient)
  require(devtools)
  
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/R5/myEnetModel1.R")
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/R5/myLMModel1.R")
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/myData_CCLE_new.R")
  
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
  
  #   KK<-c(4,12,6,23)
  
  for(kk in KK){
    filename1 = paste("~/Result_priorIncorporateStepwiseRegression_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/PriorIncorporated_cvDrug_",kk,".Rdata",sep = "")
    load(filename1)  
    filename = paste("~/Result_priorIncorporateStepwiseRegression_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/restoredPriorIncorporated_cvDrug_",kk,".Rdata",sep = "")
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
      foldIndices <- createFolds(filteredFeatureDataScaled[,1], k = 5, list = TRUE)
      
      resultsScale <- foreach(k = 1:length(foldIndices)) %dopar% {      
        groupNum1<-which(resultSTEP[[k]]$penalty == 0)
        if(length(groupNum1)==0){
          foldModel1 <- myEnetModel1$new()
          
          foldModel1$customTrain(filteredFeatureDataScaled[-foldIndices[[k]],], filteredResponseDataScaled[-foldIndices[[k]]], alpha = 1, nfolds = 5)
          
          res <- list(trainPredictions = foldModel1$customPredict(filteredFeatureDataScaled[-foldIndices[[k]],]), 
                      trainObservations = filteredResponseDataScaled[-foldIndices[[k]]],
                      testPredictions = foldModel1$customPredict(filteredFeatureDataScaled[foldIndices[[k]],]),
                      testObservations = filteredResponseDataScaled[foldIndices[[k]]])
          
          return(res)                   
        }
        if(length(groupNum1)==1){
          
          foldModel1 <- myLMModel1$new()
          
          foldModel1$customTrain(filteredFeatureDataScaled[-foldIndices[[k]],], filteredResponseDataScaled[-foldIndices[[k]]])
          
          res <- list(trainPredictions = foldModel1$customPredict(filteredFeatureDataScaled[-foldIndices[[k]],]), 
                      trainObservations = filteredResponseDataScaled[-foldIndices[[k]]],
                      testPredictions = foldModel1$customPredict(filteredFeatureDataScaled[foldIndices[[k]],]),
                      testObservations = filteredResponseDataScaled[foldIndices[[k]]])
          
          return(res)                   
        }
        if(length(groupNum1)>1){
          foldModel1 <- myEnetModel1$new()
          
          foldModel1$customTrain(filteredFeatureDataScaled[-foldIndices[[k]],groupNum1], filteredResponseDataScaled[-foldIndices[[k]]], alpha = 0, nfolds = 5)
          
          
          res <- list(trainPredictions = foldModel1$customPredict(filteredFeatureDataScaled[-foldIndices[[k]],groupNum1]), 
                      trainObservations = filteredResponseDataScaled[-foldIndices[[k]]],
                      testPredictions = foldModel1$customPredict(filteredFeatureDataScaled[foldIndices[[k]],groupNum1]),
                      testObservations = filteredResponseDataScaled[foldIndices[[k]]])
          
          return(res)                   
        }    
      }
      
      save(resultsScale,file = filename)
      
    }
  }
}


restoreSGLR_only_Sanger<-function(pathwayName,dataCombine,KK = sort(c(103,14,129,1,39,90,127,123,13,4,110,6,9,12,95,94,108,11,17,2,122,124,126,87,84,75,41,82))){
  ### DEMO Stepwise grouping Lasso
  require(predictiveModeling)
  require(synapseClient)
  require(devtools)
  
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/R5/myEnetModel1.R")
#   source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/R5/myLMModel1.R")  
  source("~/PredictiveModel_pipeline/R5/myLMModel1.R")
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/myData_Sanger.R")
  
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
  #### Load Sanger Molecular Feature Data from Synapse ####
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
    filename1 = paste("~/Result_priorIncorporateStepwiseRegression_filterVar02/",dataCombine,"/Sanger/",pathwayName,"/PriorIncorporated_cvDrug_",kk,".Rdata",sep = "")
    load(filename1)  
    filename = paste("~/Result_priorIncorporateStepwiseRegression_filterVar02/",dataCombine,"/Sanger/",pathwayName,"/restoredPriorIncorporated_cvDrug_",kk,".Rdata",sep = "")
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
      foldIndices <- createFolds(filteredFeatureDataScaled[,1], k = 5, list = TRUE)
      
      resultsScale <- foreach(k = 1:length(foldIndices)) %dopar% {      
        groupNum1<-which(resultSTEP[[k]]$penalty == 0)
        if(length(groupNum1)==0){
          foldModel1 <- myEnetModel1$new()
          
          foldModel1$customTrain(filteredFeatureDataScaled[-foldIndices[[k]],], filteredResponseDataScaled[-foldIndices[[k]]], alpha = 1, nfolds = 5)
          
          res <- list(trainPredictions = foldModel1$customPredict(filteredFeatureDataScaled[-foldIndices[[k]],]), 
                      trainObservations = filteredResponseDataScaled[-foldIndices[[k]]],
                      testPredictions = foldModel1$customPredict(filteredFeatureDataScaled[foldIndices[[k]],]),
                      testObservations = filteredResponseDataScaled[foldIndices[[k]]])
          
          return(res)                   
        }
        if(length(groupNum1)==1){
          
          foldModel1 <- myLMModel1$new()
          
          foldModel1$customTrain(filteredFeatureDataScaled[-foldIndices[[k]],groupNum1], filteredResponseDataScaled[-foldIndices[[k]]])
          
          res <- list(trainPredictions = foldModel1$customPredict(filteredFeatureDataScaled[-foldIndices[[k]],]), 
                      trainObservations = filteredResponseDataScaled[-foldIndices[[k]]],
                      testPredictions = foldModel1$customPredict(filteredFeatureDataScaled[foldIndices[[k]],groupNum1]),
                      testObservations = filteredResponseDataScaled[foldIndices[[k]]])
          
          return(res)                   
        }
        if(length(groupNum1)>1){
          foldModel1 <- myEnetModel1$new()
          
          foldModel1$customTrain(filteredFeatureDataScaled[-foldIndices[[k]],groupNum1], filteredResponseDataScaled[-foldIndices[[k]]], alpha = 0, nfolds = 5)
          
          
          res <- list(trainPredictions = foldModel1$customPredict(filteredFeatureDataScaled[-foldIndices[[k]],groupNum1]), 
                      trainObservations = filteredResponseDataScaled[-foldIndices[[k]]],
                      testPredictions = foldModel1$customPredict(filteredFeatureDataScaled[foldIndices[[k]],groupNum1]),
                      testObservations = filteredResponseDataScaled[foldIndices[[k]]])
          
          return(res)                   
        }    
      }
      
      save(resultsScale,file = filename)
      
    }
  }
}
