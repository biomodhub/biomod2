# `.Rescaler4` <-
# function(dataToRescale, ref=NULL, run, original=FALSE)
# {
#     #preparing data
#     #homogenize the format accross original predictions and new projections 
#     if(!class(dataToRescale)[1]=='RasterLayer'){
#         DataF <- as.data.frame(dataToRescale)         
#         colnames(DataF) <- "DataF"                  
#     } else{
#         names(dataToRescale) <-"DataF"
# 	      DataF <- stack(dataToRescale) 
#     }
#     
#     #Creating or loading the scaling model
#     if(original){ 
#         Rescaling_GLM = glm(ref~DataF, data=DataF, family="binomial")
#         eval(parse(text=paste("save(Rescaling_GLM, file='", getwd(), "/models/scaling_models/Rmod_", run, "', compress='xz')", sep=""))) 
#     } else
#         eval(parse(text=paste("load('", getwd(), "/models/scaling_models/Rmod_", run, "')", sep="")))
# 	 	
#     #make the scaling prediction
#     if(!class(dataToRescale)[1]=='RasterLayer') RescaledData <- predict(Rescaling_GLM, DataF, type="response") 
# 	if(class(dataToRescale)[1]=='RasterLayer')  RescaledData <- predict(model=Rescaling_GLM, DataF, type="response")    #rasters
# 	   
#    
#     return(RescaledData)
# }



.Rescaler5 <-
function(dataToRescale, ref=NULL, name, original=FALSE, weights=NULL)
{
  compress.arg = ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
#     #preparing data
#     #homogenize the format accross original predictions and new projections 
#     if(!class(dataToRescale)[1]=='RasterLayer'){
#         DataF <- as.data.frame(dataToRescale)         
#         colnames(DataF) <- "DataF"                  
#     } else{
#         names(dataToRescale) <-"DataF"
#         DataF <- stack(dataToRescale) 
#     }
#     
    #Creating or loading the scaling model
  if(original){
      if(! file.exists(paste(getwd(),"/", unlist(strsplit(name,'_'))[1], "/models/scaling_models/", sep=""))){
        dir.create(paste(getwd(),"/", unlist(strsplit(name,'_'))[1], "/models/scaling_models/", sep=""), showWarnings=F)
      }
#       Rescaling_GLM = glm(ref~DataF, data=DataF, family="binomial", mustart = rep(0.5,length(ref)))
      ## customised wgts
      if(is.null(weights)){
        Prevalence=0.5
        nbPres <- sum(ref, na.rm=TRUE)
        nbAbs <- length(ref) - nbPres
        weights <- rep(1,length(ref))
        
        if(nbAbs > nbPres){ # code absences as 1
          weights[which(ref>0)] <- (Prevalence * nbAbs) / (nbPres * (1-Prevalence))
        } else{ # code presences as 1
          weights[which(ref==0 | is.na(ref))] <- (nbPres * (1-Prevalence)) / (Prevalence * nbAbs)
        }
        
         weights = round(weights[]) # test to remove glm & gam warnings
      }

      

      Rescaling_GLM = glm(ref~pred, data=data.frame(ref=as.numeric(ref), pred = as.numeric(dataToRescale)) , family=binomial(link=probit), x=TRUE, weights=weights)
#       Rescaling_GLM = glm(ref~DataF, data=DataF, family=binomial, weights=wgts)
#       Rescaling_GLM = glm(ref~DataF, data=DataF, family=binomial, weights=wgts)
      eval(parse(text=paste("save(Rescaling_GLM, file='", getwd(),"/",
                            unlist(strsplit(name,'_'))[1], "/models/scaling_models/",
                            name, "_scaled' ,  compress=", compress.arg, ")", sep=""))) 
    } else{
      eval(parse(text=paste("load('", getwd(),"/",unlist(strsplit(name,'_'))[1],
                            "/models/scaling_models/",name,"_scaled')", sep="")))
    }
    #make the scaling prediction
    if(! inherits(dataToRescale, "Raster")){
      RescaledData <- predict(Rescaling_GLM, data.frame(pred=as.numeric(dataToRescale)), type="response")
    } else{
      RescaledData <- predict(dataToRescale, model=Rescaling_GLM, type='response')
    }
    
  
#     cat("\n\t\t original range = ", min(DataF) ," - ", max(DataF), "\t scal ranged = ", min(RescaledData), " - ", max(RescaledData) )

    
# 	  if(class(dataToRescale)[1]=='RasterLayer')  RescaledData <- predict(model=Rescaling_GLM, DataF, type="response")    #rasters
	   
    return(RescaledData)
}


.scaling_model <-
  function(dataToRescale, ref=NULL, ...)
  {
    args <- list(...)
    prevalence <- args$prevalence
    weights <- args$weights
    
#     # if no prevalence define a 0.5 is set.
#     if(is.null(prevalence)) prevalence <- 0.5
    
    # if no weights given, some are created to rise the define prevalence
    if(is.null(weights) & ! is.null(prevalence)){
      nbPres <- sum(ref, na.rm=TRUE)
      nbAbs <- length(ref) - nbPres
      weights <- rep(1,length(ref))
      
      if(nbAbs > nbPres){ # code absences as 1
        weights[which(ref>0)] <- (prevalence * nbAbs) / (nbPres * (1-prevalence))
      } else{ # code presences as 1
        weights[which(ref==0 | is.na(ref))] <- (nbPres * (1-prevalence)) / (prevalence * nbAbs)
      }
      
      weights = round(weights[]) # to remove glm & gam warnings
    } else if(is.null(weights)){
      # only 1 weights vector
      weights <- rep(1,length(ref))
    }
     
    # define a glm to scal prediction from 0 to1 
    scaling_model <- glm(ref~pred, data=data.frame(ref=as.numeric(ref), pred = as.numeric(dataToRescale)) , family=binomial(link=probit), x=TRUE, weights=weights)
    
    return(scaling_model)
  }

###################################################################
###################################################################
# .Rescaler5 <-
#   function(dataToRescale, ref=NULL, name, original=FALSE, weights=NULL)
#   {
#     #     #preparing data
#     #     #homogenize the format accross original predictions and new projections 
#     #     if(!class(dataToRescale)[1]=='RasterLayer'){
#     #         DataF <- as.data.frame(dataToRescale)         
#     #         colnames(DataF) <- "DataF"                  
#     #     } else{
#     #         names(dataToRescale) <-"DataF"
#     #         DataF <- stack(dataToRescale) 
#     #     }
#     #     
#     #Creating or loading the scaling model
#     if(original){
#       if(! file.exists(paste(getwd(),"/", unlist(strsplit(name,'_'))[1], "/models/scaling_models/", sep=""))){
#         dir.create(paste(getwd(),"/", unlist(strsplit(name,'_'))[1], "/models/scaling_models/", sep=""), showWarnings=F)
#       }
#       #       Rescaling_GLM = glm(ref~DataF, data=DataF, family="binomial", mustart = rep(0.5,length(ref)))
#       ## customised wgts
#       if(is.null(weights)){
#         weights = ref
#         weights[ref==1] <- (length(ref)-sum(ref,na.rm=TRUE))  / sum(ref,na.rm=TRUE) 
#         weights[ref!=1] <- 1
#         ## transform to integer for warning prevent..
#         weights <- as.integer(weights * 100)
#       }
#       
#       
#       
#       Rescaling_GLM = glm(ref~pred, data=data.frame(ref=as.numeric(ref), pred = as.numeric(dataToRescale)) , family=binomial(link=probit), x=TRUE, weights=weights)
#       #       Rescaling_GLM = glm(ref~DataF, data=DataF, family=binomial, weights=wgts)
#       #       Rescaling_GLM = glm(ref~DataF, data=DataF, family=binomial, weights=wgts)
#       eval(parse(text=paste("save(Rescaling_GLM, file='", getwd(),"/",
#                             unlist(strsplit(name,'_'))[1], "/models/scaling_models/",
#                             name, "_scaled' , compress='",ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
#                             ,"')", sep=""))) 
#     } else{
#       eval(parse(text=paste("load('", getwd(),"/",unlist(strsplit(name,'_'))[1],
#                             "/models/scaling_models/",name,"_scaled')", sep="")))
#     }
#     #make the scaling prediction
#     if(! inherits(dataToRescale, "Raster")){
#       RescaledData <- predict(Rescaling_GLM, data.frame(pred=as.numeric(dataToRescale)), type="response")
#     } else{
#       cat("\n*** scaler5 raster scaling")
#       RescaledData <- predict(dataToRescale, model=Rescaling_GLM, type='response')
#     }
#     
#     
#     #     cat("\n\t\t original range = ", min(DataF) ," - ", max(DataF), "\t scal ranged = ", min(RescaledData), " - ", max(RescaledData) )
#     
#     
#     #     if(class(dataToRescale)[1]=='RasterLayer')  RescaledData <- predict(model=Rescaling_GLM, DataF, type="response")    #rasters
#     
#     return(RescaledData)
#   }
####################################################################
####################################################################
