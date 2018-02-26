##' @name BIOMOD_presenceonly
##' @aliases BIOMOD_presenceonly
##' 
##' @title evaluate models with presences only metrics
##' 
##' @description This function enables to evaluate BIOMOD.models.out and 
##' BIOMOD.EnsembleModeling.out object with presence-only evaluation methods 
##' (Boyce index and Minimal Predicted Area MPA)
##' 
##' @param modeling.output  "BIOMOD.models.out" object produced by a BIOMOD_Modeling run
##' @param EM.output        a "BIOMOD.EnsembleModeling.out" returned by BIOMOD_EnsembleModeling
##' @param save.output      logical. If TRUE (Default) the output is saved to the ".BIOMOD_DATA" folder
##' 
##' @details
##' 'em.by' of 'BIOMOD.EnsembleModeling' must be 'PA_dataset+repet' to have an 
##' ensemble for each RUN of the 'NbRunEval' argument (BIOMOD_Modeling funtion)
##' for evaluation.
##' The Boyce index returns NA values for 'SRE' models because it is not possible
##' to be calculated with binary predictions. This is also the reason why there 
##' are sometimes NA values for 'GLM' models if they don not converge.
##' 
##' @return
##' data.frame containing evaluation scores for the evaluation metrics used for 
##' the BIOMOD_Modeling function and additional Boyce index and MPA
##' 
##' @references
##' Engler, R., A. Guisan, and L. Rechsteiner. 2004. An improved approach for predicting the distribution of rare and endangered species from occurrence and pseudo-absence data. Journal of Applied Ecology.
##' 
##' @author
##' Frank Breiner
##' 
##' @seealso
##' ecospat.boyce, ecospat.mpa, BIOMOD_Modeling, BIOMOD_EnsembleModeling 
##' 
##' @examples
##' \dontrun{
##' require(PresenceAbsence)
##' 
##' # species occurrences
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##'                                     package="biomod2"), row.names = 1)
##' head(DataSpecies)
##' 
##' # the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##' 
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio4.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio7.grd", 
##'                              package="biomod2"),  
##'                 system.file( "external/bioclim/current/bio11.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio12.grd", 
##'                              package="biomod2"))
##' 
##' # 1. Formatting Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' # 3. Doing Modelisation
##' 
##' myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##'                                      models = c('SRE','CTA','RF'), 
##'                                      models.options = myBiomodOption, 
##'                                      NbRunEval=1, 
##'                                      DataSplit=80, 
##'                                      Yweights=NULL, 
##'                                      VarImport=3, 
##'                                      models.eval.meth = c('TSS','ROC'),
##'                                      SaveObj = TRUE,
##'                                      rescal.all.models = FALSE,
##'                                      do.full.models = FALSE)
##' 
##' # 4. Doing Ensemble Modelling
##' myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
##'                                        chosen.models = 'all',
##'                                        em.by = 'PA_dataset+repet',
##'                                        eval.metric = c('TSS'),
##'                                        eval.metric.quality.threshold = c(0.7),
##'                                        models.eval.meth = c('TSS','ROC'),
##'                                        prob.mean = TRUE,
##'                                        prob.cv = FALSE,
##'                                        prob.ci = FALSE,
##'                                        prob.ci.alpha = 0.05,
##'                                        prob.median = FALSE,
##'                                        committee.averaging = FALSE,
##'                                        prob.mean.weight = TRUE,
##'                                        prob.mean.weight.decay = 'proportional' )   
##' 
##' # evaluate Biomod models with the Boyce index and MPA
##' pres.only.eval <- BIOMOD_presenceonly(myBiomodModelOut, myBiomodEM)
##' pres.only.eval$eval
##' }

BIOMOD_presenceonly <- function(modeling.output = NULL, EM.output = NULL, save.output = T){
  
  ## note : this function need ecospat.mpa function not yet available in ecospat package
  ## TODO(damien) remove this function and make the link with the ecospat.mpa package
  ecospat.mpa <- function(Pred, Sp.occ.xy, perc = 0.9){
    perc <- 1-perc
    if(class(Pred)=="RasterLayer"){
      Pred <- extract(Pred, Sp.occ.xy)}
    round(quantile(Pred, probs = perc, na.rm = TRUE ), 3)}
  
  ecospat.boyce <- function (fit, obs, nclass = 0, window.w = "default", res = 100, PEplot = T){
    boycei<-function(interval,obs,fit){
      
      fit.bin<-fit
      obs.bin<-obs
      fit.bin[fit[]>=interval[1]&fit[]<=interval[2]]<-"i";fit.bin[fit.bin!="i"]<-0
      obs.bin[obs[]>=interval[1]&obs[]<=interval[2]]<-"i";obs.bin[obs.bin!="i"]<-0
      
      pi<-length(which(obs.bin=="i"))/length(obs)
      ei<-length(which(fit.bin=="i"))/length(fit.bin)
      fi<-pi/ei
      
      return(fi)
    }
    
    if (window.w == "default") {
      window.w <- (max(fit) - min(fit))/10
    }
    interval <- c(min(fit), max(fit))
    mini <- interval[1]
    maxi <- interval[2]
    if (nclass == 0) {
      vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - 
                                                                mini - window.w)/res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1
      interval <- cbind(vec.mov, vec.mov + window.w)
    }
    else if (length(nclass) > 1) {
      vec.mov <- c(mini, nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
    else if (nclass > 0 & length(nclass) < 2) {
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
    }
    f <- apply(interval, 1, boycei, obs, fit)
    if (length(f[which(f != "NaN")]) <= 2) {
      b <- NA
    }
    else {
      b <- cor(f[which(f != "NaN")], vec.mov[which(f != "NaN")], 
               method = "spearman")
    }
    ID <- seq(1:(length(vec.mov)))
    HS <- apply(interval, 1, sum)/2
    if (PEplot == T) 
      plot((apply(interval[which(f != "NaN"), ], 1, sum)/2), 
           f[which(f != "NaN")], xlab = "Habitat suitability", 
           ylab = "Predicted/Expected ratio")
    results <- list(F.ratio = f, Spearman.cor = round(b, 3), 
                    HS = HS, ID = ID)
    return(results)
  }
  
  
  #   if(!require(PresenceAbsence)){stop("PresenceAbsence package required!")}
  
  if(!is.null(modeling.output)){  
    myModelEval <- get_evaluations(modeling.output,as.data.frame=T)
    myModelEval[,1] <- as.character(myModelEval[,1])
    for(i in 1:nrow(myModelEval)){myModelEval[i,1] <- paste(c(modeling.output@sp.name,strsplit(as.character(myModelEval[i,1]),split="_")[[1]][3:1]),collapse="_")  } 
    myModelPred <- get_predictions(modeling.output,as.data.frame=T)
    if(modeling.output@has.evaluation.data == T){
      myModelPred.eval  <- as.data.frame(get(load(paste(modeling.output@"sp.name","/.BIOMOD_DATA/",modeling.output@modeling.id,"/models.prediction.eval", sep=""))))
      for(i in 1:ncol(myModelPred.eval)){colnames(myModelPred.eval)[i] <- paste(c(modeling.output@sp.name,strsplit(colnames(myModelPred.eval)[i],split="[.]")[[1]][3:1]),collapse="_")  }       
    }
  }
  
  if(!is.null(EM.output)){
    if(EM.output@em.by!='PA_dataset+repet'){stop("em.by of 'BIOMOD.EnsembleModeling' must be 'PA_dataset+repet'")} 
    myModelEvalEF <- get_evaluations(EM.output,as.data.frame=T)
    myModelEvalEF[,1] <- paste(modeling.output@sp.name,as.character(myModelEvalEF[,1]),sep="_")
    #for(i in 1:nrow(myModelEvalEF)){myModelEvalEF[i,1] <- paste("EF",strsplit(as.character(myModelEvalEF[,1]),split="_")[[i]][2],"AllData",sep="_")}
    if(!is.null(modeling.output)){
      myModelEval <- rbind(myModelEval, myModelEvalEF)
    }
    # DOES NOT EXIST YET BUT IN NEXT BIOMOD2 VERSION:
    #myModelPredEF <- get_predictions(EM.output,as.data.frame=T)    
    ###### than this part gets redundand:
    myBiomodProjFF <- BIOMOD_EnsembleForecasting(
      new.env = get(load(modeling.output@formated.input.data@link))@data.env.var,
      proj.name = paste(modeling.output@modeling.id,"cv_EF",sep="_"),   
      EM.output = EM.output)
    
    myBiomodProjFF <- as.data.frame(myBiomodProjFF@proj@val)
    #colnames(myBiomodProjFF) <- gsub("AllAlgos_ROC_EMwmean","EF",  colnames(myBiomodProjFF))
    myModelPred <- cbind(myModelPred, myBiomodProjFF)
    ###### end of redundancy
    
    if(modeling.output@has.evaluation.data == T){
      myBiomodProjFF.eval <- BIOMOD_EnsembleForecasting(
        new.env = get(load(modeling.output@formated.input.data@link))@eval.data.env.var,
        proj.name = paste(modeling.output@modeling.id,"cv_EF_eval",sep="_"),
        EM.output = EM.output)
      
      myBiomodProjFF.eval <- as.data.frame(myBiomodProjFF.eval@proj@val)
      #colnames(myBiomodProjFF.eval) <- gsub("AllAlgos_ROC_EMwmean","EF",  colnames(myBiomodProjFF.eval))
      myModelPred.eval <- cbind(myModelPred.eval, myBiomodProjFF.eval)      
    }  
  }
  
  mpa.eval <- boyce.eval <- myModelEval[!duplicated(myModelEval[,1]),]
  boyce.eval$Eval.metric <- "boyce"
  mpa.eval$Eval.metric <- "mpa"
  boyce.eval[,3:7]<-mpa.eval[,3:7]<-NA
  

  calib.lines <- get(load(modeling.output@calib.lines@link))
#   ## keep the names of calib lines dimensions
#   calib.lines.dimnames <- dimnames(calib.lines)
#   ## extarct the first 2 dimensions
#   calib.lines <- calib.lines[,,1, drop = FALSE]
#   dimnames(calib.lines) <- calib.lines.dimnames[1:2]
#   cat("\n*** start test \n")
#   cat(modeling.output@calib.lines@link)
#   print(dim(calib.lines))
#   print(dimnames(calib.lines))
#   cat("\n*** end test \n")

  myResp <- get(load(modeling.output@formated.input.data@link))@data.species
  ###MPA & BOYCE     
  for(i in 1:nrow(boyce.eval)){
    n <- length(strsplit(as.character(boyce.eval[i,1]),split="_")[[1]])
    tec <- paste(strsplit(as.character(boyce.eval[i,1]),split="_")[[1]][3:n],collapse="_") 
    Model.name <- boyce.eval[i,1]
    run <- strsplit(Model.name,split="_")[[1]][c(grep("RUN",strsplit(Model.name,split="_")[[1]]),grep("Full",strsplit(Model.name,split="_")[[1]]))]
    
    if(sum(!calib.lines[,paste("_",run, sep=""), 1])==0){      #this is the full model
      test <- myResp[calib.lines[, paste("_",run, sep=""), 1]]      
      Pred<-myModelPred[calib.lines[, paste("_",run, sep=""), 1],Model.name]                         
    }else{
      test <- myResp[!calib.lines[,paste("_",run, sep=""), 1]]
      Pred<-myModelPred[!calib.lines[,paste("_",run, sep=""), 1],Model.name]                         
    }
    
    if(sum(!is.na(Pred))>1){
#       boy <- ecospat::ecospat.boyce(fit=Pred,obs=Pred[test==1], PEplot=F)
      boy <- ecospat.boyce(fit=Pred,obs=Pred[test==1], PEplot=F)
      boyce.eval[boyce.eval[,1] == Model.name, 3] <- boy$Spearman.cor
      if( sum(boy$F.ratio<1,na.rm=T)>0){
        boyce.eval[boyce.eval[,1]==Model.name,5] <-  round(boy$HS[max(which(boy$F.ratio<1))],0)
        DATA<-cbind(1:length(Pred), test, Pred/1000)
        EVAL<- PresenceAbsence::presence.absence.accuracy(DATA, threshold=round(boy$HS[max(which(boy$F.ratio<1))],0)/1000) 
        boyce.eval[boyce.eval[,1]==Model.name,6] <-  EVAL$sensitivity 
        boyce.eval[boyce.eval[,1]==Model.name,7] <-  EVAL$specificity 
      }else{
        boyce.eval[boyce.eval[,1]==Model.name,7] <-  boyce.eval[boyce.eval[,1]==Model.name,6] <-  boyce.eval[boyce.eval[,1]==Model.name,5] <- NA 	
      }
      
      mpa.eval[mpa.eval[,1]==Model.name,5] <- ecospat.mpa(Pred[test==1])
      EVAL <- PresenceAbsence::presence.absence.accuracy(DATA, threshold=ecospat.mpa(Pred[test])/1000) 
      mpa.eval[mpa.eval[,1]==Model.name,6] <-  EVAL$sensitivity 
      mpa.eval[mpa.eval[,1]==Model.name,7] <-  EVAL$specificity  
    }
    
    if(modeling.output@has.evaluation.data == T){
      myResp.eval <- get(load(modeling.output@formated.input.data@link))@eval.data.species
      Pred.eval<-myModelPred.eval[,Model.name]   
      boy <- ecospat.boyce(fit=Pred.eval,obs=Pred.eval[myResp.eval==1], PEplot=F)
      boyce.eval[boyce.eval[,1]==Model.name,"Evaluating.data"] <- boy$Pearson.cor
      
      mpa.eval[mpa.eval[,1]==Model.name,"Evaluating.data"] <- ecospat.mpa(Pred.eval[myResp.eval==1]) 
    }
  }
  myModelEval[,6:7] <- round(myModelEval[,6:7],1)
  boyce.eval[,6:7] <- round(boyce.eval[,6:7]*100,1)
  mpa.eval[,6:7] <- round(mpa.eval[,6:7]*100,1)  
  
  if(modeling.output@has.evaluation.data == T){
    output <- list(eval=rbind(myModelEval,boyce.eval,mpa.eval),myBiomodProjFF=myBiomodProjFF,myBiomodProjEF.eval=myBiomodProjFF.eval) 
    if(save.output){
      if(!is.null(modeling.output)){
      sp <- modeling.output@sp.name
      modeling.id <- modeling.output@modeling.id
      }
      if(!is.null(EM.output)){
        sp <- EM.output@sp.name
        modeling.id <- EM.output@modeling.id
      }      
      save(output, paste(sp, "/.BIOMOD_DATA/", modeling.id, "/presenceonlyevaluation_",sp, sep=""))
    }
    return(output)
  }else{
    output <- list(eval=rbind(myModelEval,boyce.eval,mpa.eval),myBiomodProjFF=myBiomodProjFF)
    if(save.output){
      if(!is.null(modeling.output)){
        sp <- modeling.output@sp.name
        modeling.id <- modeling.output@modeling.id
      }
      if(!is.null(EM.output)){
        sp<-EM.output@sp.name
        modeling.id <- EM.output@modeling.id
      }      
      save(output, file=paste(sp, "/.BIOMOD_DATA/", modeling.id,"/presenceonlyevaluation_",sp, sep=""))
      }
    return(output)
    }
}


