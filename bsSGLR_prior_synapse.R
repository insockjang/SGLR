bsSGLR_prior_CCLE<-function(pathwayName,dataCombine,KK=c(1:24),bsNum = 100,mcCoreNum = 32){
  ### DEMO Stepwise grouping Lasso
  
  require(predictiveModeling)
  require(synapseClient)
  require(devtools)
  
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/R5/myEnetModel1.R")
  source_url("https://raw.githubusercontent.com/insockjang/PredictiveModel_pipeline/master/myData_CCLE_new.R")
  source_url("https://raw.githubusercontent.com/insockjang/SGLR/master/parallel_stepwiseDecision.R")
  source_url("https://raw.githubusercontent.com/insockjang/Supplements/master/ListMake_CCLE.R")
  # drug annotations
  drugs<-synGet("syn2631135")
  load(drugs@filePath)
  
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
  
#   myFolder  <- Folder(name = "CCLE", parentId = "syn2575943")
#   myFolder  <- synStore(myFolder)
#   
#   myFolder.1  <- Folder(name = dataCombine, parentId = myFolder$properties$id)
#   myFolder.1  <- synStore(myFolder.1)    
#   
#   myFolder.2  <- Folder(name = pathwayName, parentId = myFolder.1$properties$id)
#   myFolder.2  <- synStore(myFolder.2)
#   
#   myFolder.3  <- Folder(name = "SGLR_prior_bootstrap", parentId = myFolder.2$properties$id)
#   myFolder.3  <- synStore(myFolder.3)
#   
  qry0<-synapseQuery(paste("select id, name from entity where entity.parentId == '","syn2575943", "'"))  
  qry1<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry0$entity.id[which(qry0$entity.name == "CCLE")], "'"))  
  qry2<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry1$entity.id[which(qry1$entity.name == dataCombine)], "'"))  
  qry3<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry2$entity.id[which(qry2$entity.name == pathwayName)], "'"))  
  qry4<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap")], "'"))  
    
  
  if(is.null(qry4)){
    
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
      name1<-drugNameCCLE[kk]
      KKK<-ListMake2("ActArea",dataCombine,pathwayName)
      plotFile  <- synStore(File(path=filename, parentId= qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap")],name = name1),
                            used=KKK,                              
                            activityName="Incoporated Priors from Stepwise forward selection : bootstrapping for features",
                            activityDescription="To execute run: bsSGLR_prior_CCLE(pathwayName,dataCombine,KK=c(1:24),bsNum = 100)")            
    }
  }else{
    if(nrow(qry4)!=length(KK)){
      
      for(kk in KK[(nrow(qry4)+1):length(KK)]){
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
        name1<-drugNameCCLE[kk]
        KKK<-ListMake2("ActArea",dataCombine,pathwayName)
        plotFile  <- synStore(File(path=filename, parentId=qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap")],name = name1),
                              used=KKK,                              
                              activityName="Incoporated Priors from Stepwise forward selection : bootstrapping for features",
                              activityDescription="To execute run: bsSGLR_prior_CCLE(pathwayName,dataCombine,KK=c(1:24),bsNum = 100)")            
      }
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
  source_url("https://raw.githubusercontent.com/insockjang/Supplements/master/ListMake_Sanger.R")
  
  # drug annotations
  drugs<-synGet("syn2631135")
  load(drugs@filePath)
  
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
  
  
  qry0<-synapseQuery(paste("select id, name from entity where entity.parentId == '","syn2575943", "'"))  
  qry1<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry0$entity.id[which(qry0$entity.name == "Sanger")], "'"))  
  qry2<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry1$entity.id[which(qry1$entity.name == dataCombine)], "'"))  
  qry3<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry2$entity.id[which(qry2$entity.name == pathwayName)], "'"))  
  qry4<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap")], "'"))  
  
  if(is.null(qry4)){
    
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
      name1<-drugNameSangerIC[kk]
      KKK<-ListMake2("IC50",dataCombine,pathwayName)
      plotFile  <- synStore(File(path=filename, parentId=qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap")],name = name1),
                            used=KKK,                              
                            activityName="Incoporated Priors from Stepwise forward selection : bootstrapping for features",
                            activityDescription="To execute run: bsSGLR_prior_CCLE(pathwayName,dataCombine,KK=c(1:24),bsNum = 100)")            
    }
  }else{
    if(nrow(qry4)!=length(KK)){
      
      for(kk in KK[(nrow(qry4)+1):length(KK)]){
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
        name1<-drugNameSangerIC[kk]
        KKK<-ListMake2("IC50",dataCombine,pathwayName)
        plotFile  <- synStore(File(path=filename, parentId = qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap")],name = name1),
                              used=KKK,                              
                              activityName="Incoporated Priors from Stepwise forward selection : bootstrapping for features",
                              activityDescription="To execute run: bsSGLR_prior_CCLE(pathwayName,dataCombine,KK=c(1:24),bsNum = 100)")            
      }
    }
  }
  
}
