## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE------------------------------------------------------------
library(rangeModelMetadata)
library(sf)
library(spocc)
library(dplyr)

## -----------------------------------------------------------------------------
rmm1 <- rmmTemplate(family=c('base')) 
str(rmm1)

## -----------------------------------------------------------------------------
rmm2 <- rmmTemplate(family=NULL)
str(rmm2)

## -----------------------------------------------------------------------------
rmmSuggest('dataPrep',fullFieldDepth=FALSE)
rmmSuggest('dataPrep',fullFieldDepth=TRUE) # for all fields below the specified one
rmmSuggest('dataPrep$biological$duplicateRemoval')
rmmSuggest('dataPrep$biological$duplicateRemoval$rule')

## -----------------------------------------------------------------------------
rmmSuggest('model')
rmmSuggest('model$algorithm$maxent')
rmmSuggest('$model$algorithm$maxent$featureSet')

## -----------------------------------------------------------------------------
rmm <- rmmTemplate()
rmm <- rmmAutofillPackageCitation(rmm,c('terra','sf'))
# search GBIF for occurrence data to demonstrate the autofill function
bv <- spocc::occ('Bradypus variegatus', 'gbif', limit=50, has_coords=TRUE)
rmm <- rmmAutofillspocc(rmm,bv$gbif)
# get some env layers to demonstrate the autofill function
rasterFiles <- list.files(path=paste(system.file(package='dismo'), '/ex', sep=''),
                       pattern='grd', full.names=TRUE)
# make a stack of the rasters
env <- terra::rast(rasterFiles)
rmm <- rmmAutofillEnvironment(rmm,env,transfer=0) # for fitting environment
# just using the same rasters for demonstration; in practice these are different
rmm <- rmmAutofillEnvironment(rmm,env,transfer=1) # for transfer environment 1
rmm <- rmmAutofillEnvironment(rmm,env,transfer=2) # for transfer environment 2 

## -----------------------------------------------------------------------------
empties <- rmmCheckEmpty(rmm)

## -----------------------------------------------------------------------------
# Make an empty template
rmm1 <- rmmTemplate() 
# Add a new, non-standard field
rmm1$dataPrep$biological$taxonomicHarmonization$taxonomy_source <- "The Plant List" # # Checking the names identifies the new, non-standard field we've added ("taxonomy_source")
rmm1 <- rmmCheckName(rmm1) 

## -----------------------------------------------------------------------------
#First, we create an empty rmm template
rmm1 <- rmmTemplate() 
#We add 3 of the bioclim layers, including a spelling error (an extra space) in bio2, and a word that is clearly not a climate layer, 'cromulent'.
rmm1$data$environment$variableNames <- c("bio1", "bio 2", "bio3", "cromulent") 
#Now, when we check the values, we see that bio1 and bio2 are reported as exact matches, while 'bio 2' is flagged as a partial match with a suggested value of 'bio2', and 'cromulent' is flagged as not matched at all.
rmmCheckValue(rmm = rmm1) 
#If we'd like to return a dataframe containing this information in a perhaps more useful format:
rmmCheckValueOutput <- rmmCheckValue(rmm = rmm1, returnData = TRUE)

## -----------------------------------------------------------------------------
rmmCheckFinalize(rmm, family='base')

## ----eval=F-------------------------------------------------------------------
#  outFile <- '~/Desktop/demo_rmmToCSV.csv'
#  rmmObj <- rmmTemplate()
#  rmmToCSV(rmmObj, filename=outFile)
#  system(paste0('open ', outFile, ' -a "Microsoft Excel"'))

## -----------------------------------------------------------------------------
dd <- rmmDataDictionary()
str(dd)
# rmmDataDictionary(excel=TRUE) # try this if you have excel

