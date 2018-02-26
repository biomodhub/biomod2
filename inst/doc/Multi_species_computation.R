### R code from vignette source 'Multi_species_computation.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: options
###################################################
options(prompt = " ", continue = "  ", width = 60, digits=4)
.CurFileName <- "simple_species"
# .PrefixName <- strsplit(.CurFileName, "\\.")[[1]][1]
.PrefixName <- file.path("figs",.CurFileName)
.RversionName <- R.version.string
.PkgVersion <- packageDescription("biomod2")$Version
.SupportedDataVignette <- paste("run:",system.file('doc/Simple_species_modelling.pdf',package='biomod2'),sep="")


###################################################
### code chunk number 2: LoadSp_1
###################################################
# 1. loading species occurrences data
# 1. loading species occurrences data
library(biomod2)

DataSpecies <- read.csv( system.file( 
                          "external/species/mammals_table.csv", 
                          package="biomod2"))
                            
head(DataSpecies)



###################################################
### code chunk number 3: LoadEnv_1
###################################################
# 2. loading environmental data

# Environmental variables extracted from Worldclim (bio_3, bio_4, 
# bio_7, bio_11 & bio_12)
require(raster)
myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
                             package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))



###################################################
### code chunk number 4: Loop_1
###################################################

# define the species of interest
sp.names <- c("ConnochaetesGnou", "GuloGulo", "PantheraOnca", "PteropusGiganteus", 
              "TenrecEcaudatus", "VulpesVulpes")

# loop on species == applying the same functions to each species
for(sp.n in sp.names){
  myRespName = sp.n
  
  cat('\n',myRespName,'modeling...')  
  ### definition of data 
  ## i.e keep only the column of our species
  myResp <- as.numeric(DataSpecies[,myRespName])
  
  myRespCoord = DataSpecies[c('X_WGS84','Y_WGS84')]  
  
  
  
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName)
    
  
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  ### Modelling 
  myBiomodModelOut <- BIOMOD_Modeling( 
                        myBiomodData, 
                        models = c('SRE','CTA','RF','MARS','FDA'), 
                        models.options = myBiomodOption, 
                        NbRunEval=3, 
                        DataSplit=80, 
                        Prevalence=0.5, 
                        VarImport=3,
                        models.eval.meth = c('TSS','ROC'),
                        SaveObj = TRUE,
                        rescal.all.models = TRUE,
                        do.full.models = FALSE,
                        modeling.id = paste(myRespName,"FirstModeling",sep=""))
  
  ### save models evaluation scores and variables importance on hard drive
  capture.output(get_evaluations(myBiomodModelOut),
      file=file.path(myRespName, 
                     paste(myRespName,"_formal_models_evaluation.txt", sep="")))
  
  capture.output(get_variables_importance(myBiomodModelOut),
      file=file.path(myRespName, 
                     paste(myRespName,"_formal_models_variables_importance.txt", sep="")))               
                 
  
  ### Building ensemble-models
  myBiomodEM <- BIOMOD_EnsembleModeling( 
                  modeling.output = myBiomodModelOut,
                  chosen.models = 'all',
                  em.by='all',
                  eval.metric = c('TSS'),
                  eval.metric.quality.threshold = c(0.7),
                  prob.mean = T,
                  prob.cv = T,
                  prob.ci = T,
                  prob.ci.alpha = 0.05,
                  prob.median = T,
                  committee.averaging = T,
                  prob.mean.weight = T,
                  prob.mean.weight.decay = 'proportional' )
  
  ### Make projections on current variable
  myBiomodProj <- BIOMOD_Projection(
                      modeling.output = myBiomodModelOut,
                      new.env = myExpl,
                      proj.name = 'current',
                      selected.models = 'all',
                      binary.meth = 'TSS',
                      compress = 'xz',
                      clamping.mask = F,
                      output.format = '.grd')
  
  ### Make ensemble-models projections on current variable
  myBiomodEF <- BIOMOD_EnsembleForecasting( 
                  projection.output = myBiomodProj,
                  EM.output = myBiomodEM)
}



###################################################
### code chunk number 5: alpha1
###################################################
# define a mask of studied
alphaMap <- reclassify(subset(myExpl,1), c(-Inf,Inf,0))

# # add all other species map
for(sp.n in sp.names){
  # add layer
  alphaMap <- 
    alphaMap + 
    subset(stack(file.path(sp.n,
                           "proj_current", 
                           paste("proj_current_",
                                 sp.n, 
                                 "_TotalConsensus_EMbyTSS_TSSbin.grd", sep=""))), 1)
}

# summary of created raster
alphaMap



###################################################
### code chunk number 6: alpha_2
###################################################
plot(alphaMap, main = expression( paste(alpha, "-diversity based on",
                                " TotalConsensus_EMbyTSS_TSSbin outputs")))



###################################################
### code chunk number 7: Lapply_1
###################################################
MyBiomodSF <- function(sp.n){
  
  myRespName = sp.n
  
  cat('\n',myRespName,'modeling...')  
  ### definition of data 
  ## i.e keep only the column of our species
  myResp <- as.numeric(DataSpecies[,myRespName])
  
  myRespCoord = DataSpecies[c('X_WGS84','Y_WGS84')]
  
  
  
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName)
  
  
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  ### Modelling 
  myBiomodModelOut <- BIOMOD_Modeling( 
    myBiomodData, 
    models = c('SRE','CTA','RF','MARS','FDA'), 
    models.options = myBiomodOption, 
    NbRunEval=3, 
    DataSplit=80, 
    Prevalence=0.5, 
    VarImport=3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))
  
  ### save models evaluation scores and variables importance on hard drive
  capture.output(get_evaluations(myBiomodModelOut),
     file=file.path(myRespName, 
                    paste(myRespName,"_formal_models_evaluation.txt", sep="")))
  
  capture.output(get_variables_importance(myBiomodModelOut),
     file=file.path(myRespName, 
                    paste(myRespName,"_formal_models_variables_importance.txt", sep="")))               
  
  
  ### Building ensemble-models
  myBiomodEM <- BIOMOD_EnsembleModeling( 
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by='all',
    eval.metric = c('TSS'),
    eval.metric.quality.threshold = c(0.7),
    prob.mean = T,
    prob.cv = T,
    prob.ci = T,
    prob.ci.alpha = 0.05,
    prob.median = T,
    committee.averaging = T,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional' )
  
  ### Make projections on current variable
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExpl,
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = 'TSS',
    compress = 'xz',
    clamping.mask = F,
    output.format = '.grd')
  
  ### Make ensemble-models projections on current variable
  myBiomodEF <- BIOMOD_EnsembleForecasting( 
    projection.output = myBiomodProj,
    EM.output = myBiomodEM)
                        
}



###################################################
### code chunk number 8: Lapply_2
###################################################
myLapply_SFModelsOut <- lapply( sp.names, MyBiomodSF)


###################################################
### code chunk number 9: SnowFold_1 (eval = FALSE)
###################################################
## install.packages('snowfall', dependencies=TRUE)


###################################################
### code chunk number 10: SnowFold_2
###################################################
library(snowfall)


###################################################
### code chunk number 11: SnowFold_3 (eval = FALSE)
###################################################
## 
## ## Init snowfall
## library(snowfall)
## sfInit(parallel=TRUE, cpus=2 )  ## we select 2 CPUs. If you have 8 CPUs, put 8. 
## 
## ## Export packages
## sfLibrary('biomod2', character.only=TRUE)
## 
## ## Export variables
## sfExport('mySpeciesOcc')
## sfExport('myExpl')
## sfExport('sp.names')
## 
## # you may also use sfExportAll() to export all your workspace variables
## 
## ## Do the run
## mySFModelsOut <- sfLapply( sp.names, MyBiomodSF)
## 
## ## stop snowfall
## sfStop( nostop=FALSE )
## 


