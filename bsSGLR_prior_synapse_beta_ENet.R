bsSGLR_prior_CCLE_beta_ENet<-function(pathwayName,dataCombine,KK=c(1:24),bsNum = 100,mcCoreNum = 32){
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
  
    
  for(kk in KK){
    qry0<-synapseQuery(paste("select id, name from entity where entity.parentId == '","syn2575943", "'"))  
    qry1<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry0$entity.id[which(qry0$entity.name == "CCLE")], "'"))  
    qry2<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry1$entity.id[which(qry1$entity.name == dataCombine)], "'"))  
    qry3<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry2$entity.id[which(qry2$entity.name == pathwayName)], "'"))  
    qry4<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap")], "'"))  
    qry5<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap_beta_ENet")], "'"))  
    
    qq<-match(drugNameCCLE[kk],qry4$entity.name)
    pp<-match(drugNameCCLE[kk],qry5$entity.name)
    
    if(!is.na(qq) & is.na(pp)){
      filename = paste("~/SGLR_bs100_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/PriorIncorporated_bsDrug_beta_ENet_",kk,".Rdata",sep = "")
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
      
      aa1<-synGet(qry4$entity.id[qq])
      load(aa1@filePath)
            
      resultWeightFun<-function(kkk){
        bsModel1 <- myEnetModel1$new()        
        bsModel1$customTrain(filteredFeatureDataScaled[bootIndices[[kkk]],], filteredResponseDataScaled[bootIndices[[kkk]]], alpha = 0.01, nfolds = 5,penalty.factor = resultSTEP[[kkk]]$penalty)
        beta<-bsModel1$getCoefficients()      
        return(beta[-1,])
      }
      resultWeight <- mclapply(1:bsNum,function(x)resultWeightFun(x),mc.cores=mcCoreNum)
      
      save(resultWeight,file = filename)
      name1<-drugNameCCLE[kk]
      KKK<-ListMake3("ActArea",dataCombine,pathwayName,qry4$entity.name[qq],qry4$entity.id[qq])
      plotFile  <- synStore(File(path=filename, parentId= qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap_beta_ENet")],name = name1),
                            used=KKK,                              
                            activityName="Incoporated Priors from Stepwise forward selection : bootstrapping beta coefficients for features",
                            activityDescription="To execute run: bsSGLR_prior_CCLE_beta(pathwayName,dataCombine,KK=c(1:24),bsNum = 100)")            
    }
  }
}

bsSGLR_prior_Sanger_beta_ENet<-function(pathwayName,dataCombine,KK=NA,bsNum = 100,mcCoreNum = 32){
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
    
  for(kk in KK){
    qry0<-synapseQuery(paste("select id, name from entity where entity.parentId == '","syn2575943", "'"))  
    qry1<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry0$entity.id[which(qry0$entity.name == "Sanger")], "'"))  
    qry2<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry1$entity.id[which(qry1$entity.name == dataCombine)], "'"))  
    qry3<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry2$entity.id[which(qry2$entity.name == pathwayName)], "'"))  
    qry4<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap")], "'"))  
    qry5<-synapseQuery(paste("select id, name from entity where entity.parentId == '",qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap_beta_ENet")], "'"))  
    
    qq<-match(drugNameSangerIC[kk],qry4$entity.name)
    pp<-match(drugNameSangerIC[kk],qry5$entity.name)
    
    if(!is.na(qq) & is.na(pp)){
      
      filename = paste("~/SGLR_bs100_filterVar02/",dataCombine,"/Sanger/",pathwayName,"/PriorIncorporated_bsDrug_beta_ENet_",kk,".Rdata",sep = "")
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
      
      aa1<-synGet(qry4$entity.id[qq])
      load(aa1@filePath)
      
      resultWeightFun<-function(kkk){
        bsModel1 <- myEnetModel1$new()        
        bsModel1$customTrain(filteredFeatureDataScaled[bootIndices[[kkk]],], filteredResponseDataScaled[bootIndices[[kkk]]], alpha = 0.01, nfolds = 5,penalty.factor = resultSTEP[[kkk]]$penalty)
        beta<-bsModel1$getCoefficients()      
        return(beta[-1,])
      }
      
      resultWeight <- mclapply(1:bsNum,function(x)resultWeightFun(x),mc.cores=mcCoreNum)
      
      save(resultWeight,file = filename)
      
      name1<-drugNameSangerIC[kk]
      KKK<-ListMake3("IC50",dataCombine,pathwayName,qry4$entity.name[qq],qry4$entity.id[qq])
      plotFile  <- synStore(File(path=filename, parentId=qry3$entity.id[which(qry3$entity.name == "SGLR_prior_bootstrap_beta_ENet")],name = name1),
                            used=KKK,                              
                            activityName="Incoporated Priors from Stepwise forward selection : bootstrapping beta coefficients for features",
                            activityDescription="To execute run: bsSGLR_prior_Sanger_beta(pathwayName,dataCombine,KK= 1:138,bsNum = 100)")            
    }
  }
}