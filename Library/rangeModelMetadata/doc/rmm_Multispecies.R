## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache= FALSE)

## ----message=FALSE------------------------------------------------------------
library(rangeModelMetadata)
library(terra)

## -----------------------------------------------------------------------------
rmm <- rmmTemplate(family=c('base'))

rmm$authorship$rmmName <- 'Owens_2015_Gadids'
rmm$authorship$names <- 'Owens, Hannah'
rmm$authorship$contact <- 'hannah.owens@gmail.com'
rmm$authorship$relatedReferences <- '@article{, title={Evolution of codﬁshes (Teleostei: Gadinae) in geographical and ecological space: evidence that physiological limits drove diversiﬁcation of subarctic ﬁshes},author={Owens, Hannah},journal={Journal of Biogeography}, year={2015}, publisher={Wiley-Blackwell}}'

## -----------------------------------------------------------------------------
rmm$data$occurrence$taxon <- c("Arctogadus borisovi", "Arcogadus glacialis", "Boreogadus saida", "Eleginus gracilis", "Gadus macrocephalus", "Gadus morhua", "Gadus ogac", "Melanogrammus aeglefinus", "Merlangius merlangus", "Microgadus proximus", "Microgadus tomcod", "Pollachius pollachius", "Pollachius virens", "Gadus chalcogrammus", "Trisopterus esmarkii", "Trisopterus luscus", "Trisopterus minutus")
rmm$data$occurrence$dataType <- c("Presence only")
rmm$data$occurrence$sources <- c("Global Biodiversity
Information Facility (GBIF)", "Ocean Biogeographic Information System (OBIS)", "Zoological Institute at the Russian Academy of Sciences (ZIN)", "Natural History Museum in London (NHMUK)")
rmm$data$occurrence$presenceSampleSize <- c(16, 41, 534, 121, 249, 448, 116, 273, 151, 46, 25, 68, 200, 282, 4, 5, 128)
rmm$data$occurrence$backgroundSampleSizeSet <- 10000

rmm$data$environment$variableNames <- c("Minimum Sea Ice Concentration", "Maximum Sea Ice Concentration", "Mean Mixed Layer Depth", "Minimum Mixed Layer Depth", "Maximum Mixed Layer Depth", "Mean Bottom Salinity", "Minimum Bottom Salinity", "Maximum Bottom Salinity", "Mean Surface Salinity", "Minimum Surface Salinity", "Maximum Surface Salinity", "Mean Bottom Temperature", "Minimum Bottom Temperature", "Maximum Bottom Temperature", "Mean Surface Temperature", "Minimum Surface Temperature", "Maximum Surface Temperature")

## -----------------------------------------------------------------------------
#Getting and documenting extents
# This code is provided as an example of how to extract extent information from rasters; below we manaully input the results so that we don't need to distribute these layers with this pacakge. 
# setwd("~/Dropbox/rmm/DataForGadidDocumentation/ProjectDirectory/MTrimmedLayers/")
# trainingRegionList <- list.files(pattern = "etopo", recursive = T)
# extVect <- list()
# for(n in trainingRegionList){
#     tmp <- rast(n)
#     extVect[[length(extVect)+1]] <- ext(tmp)
# }
# rmm$data$environment$extentSet <- extVect
# rm(trainingRegionList, tmp, extVect, n)
# Manual alternative
rmm$data$environment$extentSet <- list(
  ext(-180,180,50.06732,90.06732),
  ext(-180,180,54.77051,89.77051), 
  ext(-180,180,36.77808,89.77808), 
  ext(-180,180,39.46948,81.46948), 
  ext(-180,180,23.25928 ,74.25928), 
  ext(-80.8255,52.1745,32.15668,84.15668), 
  ext( -169.2665,-31.26648,35.74249,82.74249), 
  ext(-77.09869,56.90131,35.9895,82.9895), 
  ext(-27.11713,42.88287,29.41553 ,72.41553), 
  ext(-180,-115,31.06689 ,67.06689), 
  ext(-74.711,-54.711,38.268,52.268), 
  ext(-25.20459,34.79541,41.92908,75.92908), 
  ext(-79.63971,33.36029,32.28247,84.28247), 
  ext(-180,180,28.88068,75.88068),
  ext(-26.31464,28.68536,42.62609,79.62609),
  ext(-20.26068,26.73932,23.31372,64.31372),
  ext(-20.40668,39.593,25.8183,71.8183))

## -----------------------------------------------------------------------------
rmm$data$environment$resolution <- "1 X 1 degree"
rmm$data$environment$sources <- c("NOAA National Geophysical Data Center", "NOAA World Ocean Atlas", "NOAA National Snow and Ice Data Center")

# omitted for the vignette to avoid distributing the data
# tmp <- rast("../EnvProjectionData/etopo.asc")
# rmm$data$transfer$environment1$extentSet <- tmp |> ext() |> as.character()
rmm$data$transfer$environment1$resolution <- "1 X 1 degree"
rmm$data$transfer$environment1$sources <- c("NOAA National Geophysical Data Center", "NOAA World Ocean Atlas", "NOAA National Snow and Ice Data Center")

## -----------------------------------------------------------------------------
rmm$dataPrep$duplicateRemoval$rule <- "coordinate"
rmm$dataPrep$questionablePointRemoval$notes <- "Points outside known distribution of species removed."
rmm$dataPrep$pointInPolygon$rule <- "Remove points outside training region of species."
rmm$dataPrep$spatialThin$rule <- "Reduced spatial resolution of points to match resolution of environmental data (1 X 1 resolution)."

## -----------------------------------------------------------------------------
rmm$modelFit$algorithm <- "maxent"
rmm$modelFit$algorithmCitation <- '@inproceedings{phillips2004maximum, title={A maximum entropy approach to species distribution modeling}, author={Phillips, Steven J and Dudik, Miroslav and Schapire, Robert E}, booktitle={Proceedings of the twenty-first international conference on Machine learning},pages={83},year={2004},organization={ACM}}'
rmm$modelFit$maxent$featureSet <- "LQP"
rmm$modelFit$maxent$notes <- "Ten bootstrap replicates trained with 50% of occurrence points chosen using random seed, maximum of 10000 iterations"

## -----------------------------------------------------------------------------
rmm$prediction$continuous$units <- "absolute probability"
rmm$prediction$transfer$environment1$units <- "absolute probability"
rmm$prediction$transfer$environment1$extrapolation <- "No clamping or extrapolation"

## -----------------------------------------------------------------------------
rmm$evaluation$notes <- "Inferred distribution congruent with known ranges for all species."

## -----------------------------------------------------------------------------
rmmClean <- rmmCleanNULLs(rmm)
rmmCheckFinalize(rmmClean)

