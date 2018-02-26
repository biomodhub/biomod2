### R code from vignette source 'Simple_species_modelling.Rnw'
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
### code chunk number 2: loading_data
###################################################
# load the library
library(biomod2)

# load our species data
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"))

head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# load the environmental raster layers (could be .img, ArcGIS 
# rasters or any supported format by the raster package)

# Environmental variables extracted from Worldclim (bio_3, bio_4, 
# bio_7, bio_11 & bio_12)
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
### code chunk number 3: formating_data
###################################################
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)


###################################################
### code chunk number 4: print_formating_data
###################################################
myBiomodData


###################################################
### code chunk number 5: plot_formating_data
###################################################
plot(myBiomodData)


###################################################
### code chunk number 6: modeling_options
###################################################
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()


###################################################
### code chunk number 7: modeling
###################################################
# 3. Computing the models 

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



###################################################
### code chunk number 8: modeling_summary
###################################################
myBiomodModelOut 


###################################################
### code chunk number 9: modeling_model_evaluation
###################################################
# get all models evaluation                                     
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
                                     
# print the dimnames of this object
dimnames(myBiomodModelEval)
                                     
# let's print the TSS scores of Random Forest
myBiomodModelEval["TSS","Testing.data","RF",,]

# let's print the ROC scores of all selected models
myBiomodModelEval["ROC","Testing.data",,,]



###################################################
### code chunk number 10: modeling_variable_importance
###################################################
# print variable importances                                    
get_variables_importance(myBiomodModelOut)


###################################################
### code chunk number 11: ensemble_modeling
###################################################
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


###################################################
### code chunk number 12: ensemble_modeling_outputs
###################################################
# print summary                     
myBiomodEM
                     
# get evaluation scores
get_evaluations(myBiomodEM)


###################################################
### code chunk number 13: projection_curent
###################################################
# projection over the globe under current conditions
myBiomodProj <- BIOMOD_Projection(
                         modeling.output = myBiomodModelOut,
                         new.env = myExpl,
                         proj.name = 'current',
                         selected.models = 'all',
                         binary.meth = 'TSS',
                         compress = 'xz',
                         clamping.mask = F,
                         output.format = '.grd')

# summary of crated oject
myBiomodProj

# files created on hard drive
list.files("GuloGulo/proj_current/")



###################################################
### code chunk number 14: projection_curent_plot
###################################################
# make some plots sub-selected by str.grep argument
plot(myBiomodProj, str.grep = 'MARS')


###################################################
### code chunk number 15: projection_curent_getProj
###################################################
# if you want to make custom plots, you can also get the projected map
myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj


###################################################
### code chunk number 16: projection_future
###################################################
# load environmental variables for the future. 
myExplFuture = stack( system.file( "external/bioclim/future/bio3.grd",
                                 package="biomod2"),
                    system.file( "external/bioclim/future/bio4.grd",
                                 package="biomod2"),
                    system.file( "external/bioclim/future/bio7.grd",
                                 package="biomod2"),
                    system.file( "external/bioclim/future/bio11.grd",
                                 package="biomod2"),
                    system.file( "external/bioclim/future/bio12.grd",
                                 package="biomod2"))

myBiomodProjFuture <- BIOMOD_Projection(
                              modeling.output = myBiomodModelOut,
                              new.env = myExplFuture,
                              proj.name = 'future',
                              selected.models = 'all',
                              binary.meth = 'TSS',
                              compress = 'xz',
                              clamping.mask = T,
                              output.format = '.grd')
                              



###################################################
### code chunk number 17: projection_current_plot
###################################################
# make some plots, sub-selected by str.grep argument
plot(myBiomodProjFuture, str.grep = 'MARS')


###################################################
### code chunk number 18: EnsembleForecasting_current
###################################################
myBiomodEF <- BIOMOD_EnsembleForecasting( 
                      EM.output = myBiomodEM,
                      projection.output = myBiomodProj)


###################################################
### code chunk number 19: EnsembleForecasting_loading_res
###################################################
myBiomodEF


###################################################
### code chunk number 20: EnsembleForecasting_plotting_res
###################################################
# reduce layer names for plotting convegences
plot(myBiomodEF)


