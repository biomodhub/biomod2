#' Single models package and functions
#'
#' A \code{data.frame} containing for each single model available in \pkg{biomod2} 
#' the package and functions to be called.
#'
#' @format A \code{data.frame} object with 12 rows and 5 variables:
#' \describe{
#'   \item{model}{all single models that can be computed in \pkg{biomod2}}
#'   \item{type}{data type associated to the models}
#'   \item{package}{\code{R} package used}
#'   \item{func}{function used in the \code{R} package}
#'   \item{train}{function called by \pkg{caret} for the tuning}
#' }
#' 
#' All single models available are the following : 
#' 
#' \itemize{
#'   \item ANN (\code{\link[nnet]{nnet}})
#'   \item CTA (\code{\link[rpart]{rpart}})
#'   \item FDA (\code{\link[mda]{fda}})
#'   \item GAM (\code{\link[gam]{gam}}, \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}})
#'   \item GBM (\code{\link[gbm]{gbm}})
#'   \item GLM (\code{\link[stats]{glm}})
#'   \item MARS (\code{\link[earth]{earth}})
#'   \item MAXENT (\url{https://biodiversityinformatics.amnh.org/open_source/maxent/})
#'   \item MAXNET (\code{\link[maxnet]{maxnet}})
#'   \item RF (\code{\link[randomForest]{randomForest}})
#'   \item SRE (\code{\link{bm_SRE}})
#'   \item XGBOOST (\code{\link[xgboost]{xgboost}})
#' }

"ModelsTable"

# ModelsTable <- data.frame(model = c('ANN', 'CTA', 'FDA', 'GAM', 'GAM', 'GAM', 'GBM', 'GLM'
#                                      , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
#                            , type = 'binary'
#                            , package = c('nnet', 'rpart', 'mda', 'gam', 'mgcv', 'mgcv', 'gbm', 'stats'
#                                          , 'earth', 'MAXENT', 'maxnet', 'randomForest', 'biomod2', 'xgboost')
#                            , func = c('nnet', 'rpart', 'fda', 'gam', 'bam', 'gam', 'gbm', 'glm'
#                                       , 'earth', 'MAXENT', 'maxnet', 'randomForest', 'bm_SRE', 'xgboost')
#                            , train = c('avNNet', 'rpart', 'fda', 'gamLoess', 'bam', 'gam', 'gbm', 'glm'
#                                        , 'earth', 'ENMevaluate', 'maxnet', 'rf', 'bm_SRE', 'xgbTree'))

# usethis::use_data(ModelsTable, overwrite = TRUE)
# usethis::use_data(ModelsTable, overwrite = TRUE, internal = TRUE)

#' Bigboss pre-defined parameter values for single models
#'
#' A \code{\link{BIOMOD.models.options}} object containing for each single model available in 
#' \pkg{biomod2} the parameter values pre-defined by \pkg{biomod2} team.
#'
#' @format A \code{\link{BIOMOD.models.options}} object with some changed values :
#' 
#' \describe{
#'    \item{\code{ANN.binary.nnet.nnet}}{
#'      \itemize{
#'        \item \code{size = 5}
#'        \item \code{decay = 5}
#'        \item \code{trace = FALSE}
#'        \item \code{rang = 0.1}
#'        \item \code{maxit = 200}
#'      }
#'    }
#'    \item{\code{CTA.binary.rpart.rpart}}{
#'      \itemize{
#'        \item \code{method = 'class'}
#'        \item \code{control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25)}
#'        \item \code{cost = NULL}
#'      }
#'    }
#'    \item{\code{FDA.binary.mda.fda}}{
#'      \itemize{
#'        \item \code{method = 'mars'}
#'      }
#'    }
#'    \item{\code{GAM.binary.gam.gam}}{
# \itemize{
#   \item ...
# }
#'    }
#'    \item{\code{GAM.binary.mgcv.bam}}{
     # \itemize{
     #   \item ...
     # }
#'    }
#'    \item{\code{GAM.binary.mgcv.gam}}{
#'      \itemize{
#'        \item \code{family = binomial(link = 'logit')}
#'        \item \code{method = 'GCV.Cp'}
#'        \item \code{control = list(epsilon = 1e-06, trace = FALSE, maxit = 100)}
#'      }
#'    }
#'    \item{\code{GBM.binary.gbm.gbm}}{
#'      \itemize{
#'        \item \code{n.trees = 2500}
#'        \item \code{interaction.depth = 7}
#'        \item \code{n.minobsinnode = 5}
#'        \item \code{shrinkage = 0.001}
#'        \item \code{cv.folds = 3}
#'        \item \code{keep.data = FALSE}
#'        \item \code{n.cores = 1}
#'      }
#'    }
#'    \item{\code{GLM.binary.stats.glm}}{
#'      \itemize{
#'        \item \code{family = binomial(link = 'logit')}
#'        \item \code{mustart = 0.5}
#'        \item \code{control = glm.control(maxit = 50)}
#'      }
#'    }
#'    \item{\code{MARS.binary.earth.earth}}{
#'      \itemize{
#'        \item \code{glm = list(family = binomial(link = 'logit'))}
#'        \item \code{ncross = 0}
#'        \item \code{nk = NULL}
#'        \item \code{penalty = 2}
#'        \item \code{thresh = 0.001}
#'        \item \code{nprune = NULL}
#'        \item \code{pmethod = 'backward'}
#'      }
#'    }
#'    \item{\code{MAXENT.binary.MAXENT.MAXENT}}{
#'      \itemize{
#'        \item \code{path_to_maxent.jar = '.'}
#'      }
#'    }
#'    \item{\code{RF.binary.randomForest.randomForest}}{
#'      \itemize{
#'        \item \code{type = 'classification'}
#'        \item \code{ntree = 500}
#'        \item \code{mtry = NULL}
#'        \item \code{strata = factor(c(0, 1))}
#'        \item \code{sampsize = NULL}
#'        \item \code{nodesize = 5}
#'        \item \code{maxnodes = NULL}
#'      }
#'    }
#'    \item{\code{SRE.binary.biomod2.bm_SRE}}{
#'      \itemize{
#'        \item \code{do.extrem = TRUE}
#'      }
#'    }
#'    \item{\code{XGBOOST.binary.xgboost.xgboost}}{
#'      \itemize{
#'        \item \code{params = list(max_depth = 2, eta = 1)}
#'        \item \code{nthread = 2}
#'        \item \code{nrounds = 4}
#'        \item \code{objective = 'binary:logistic'}
#'      }
#'    }
#' }
#' 
#' 

"OptionsBigboss"

# bm.opt <- bm_ModelingOptions(data.type = "binary", strategy = "default")
# bm.opt_gam2 <- bm_ModelingOptions(data.type = "binary", strategy = "default", models = "GAM.gam.gam")
# bm.opt_gam3 <- bm_ModelingOptions(data.type = "binary", strategy = "default", models = "GAM.mgcv.bam")
# bm.opt@models <- sort(unique(c(bm.opt@models, bm.opt_gam2@models, bm.opt_gam3@models)))
# bm.opt@options <- c(bm.opt@options, bm.opt_gam2@options, bm.opt_gam3@options)
# bm.opt@options <- bm.opt@options[bm.opt@models]
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$size = 5 #NULL
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$decay = 5
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$trace = FALSE
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$rang = 0.1
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$maxit = 200
# bm.opt@options$CTA.binary.rpart.rpart@args.values[['_allData_allRun']]$method = "class"
# bm.opt@options$CTA.binary.rpart.rpart@args.values[['_allData_allRun']]$control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25)
# bm.opt@options$CTA.binary.rpart.rpart@args.values[['_allData_allRun']]$cost = NULL
# bm.opt@options$FDA.binary.mda.fda@args.values[['_allData_allRun']]$method = "mars"
# bm.opt@options$GAM.binary.mgcv.gam@args.values[['_allData_allRun']]$family = binomial(link = 'logit')
# bm.opt@options$GAM.binary.mgcv.gam@args.values[['_allData_allRun']]$method = "GCV.Cp"
# bm.opt@options$GAM.binary.mgcv.gam@args.values[['_allData_allRun']]$control = list(epsilon = 1e-06, trace = FALSE, maxit = 100)
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$n.trees = 2500
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$interaction.depth = 7
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$n.minobsinnode = 5
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$shrinkage = 0.001
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$cv.folds = 3
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$keep.data = FALSE
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$n.cores = 1
# bm.opt@options$GLM.binary.stats.glm@args.values[['_allData_allRun']]$family = binomial(link = 'logit')
# bm.opt@options$GLM.binary.stats.glm@args.values[['_allData_allRun']]$mustart = 0.5
# bm.opt@options$GLM.binary.stats.glm@args.values[['_allData_allRun']]$control = glm.control(maxit = 50)
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$glm = list(family = binomial(link = 'logit'))
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$ncross = 0
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$nk = NULL
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$penalty = 2
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$thresh = 0.001
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$nprune = NULL
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$pmethod = 'backward'
# bm.opt@options$MAXENT.binary.MAXENT.MAXENT@args.default$path_to_maxent.jar = '.'
# bm.opt@options$MAXENT.binary.MAXENT.MAXENT@args.values[['_allData_allRun']]$path_to_maxent.jar = '.'
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$type = 'classification'
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$ntree = 500
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$mtry = NULL
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$strata = factor(c(0, 1))
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$sampsize = NULL
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$nodesize = 5
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$maxnodes = NULL
# bm.opt@options$SRE.binary.biomod2.bm_SRE@args.values[['_allData_allRun']]$do.extrem = TRUE
# bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']]$params = list(max_depth = 2, eta = 1)
# bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']]$nthread = 2
# bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']]$nrounds = 4
# bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']]$objective = "binary:logistic"
# OptionsBigboss <- bm.opt

# usethis::use_data(OptionsBigboss, overwrite = TRUE, internal = TRUE)
# usethis::use_data(OptionsBigboss, ModelsTable, overwrite = TRUE, internal = TRUE)


#' Presence-Absence data to build test SDM
#'
#' A dataset covering all the continent with presence/absence data for 6 mammal
#' species. Presence/absence were derived from range maps downloaded at 
#' \href{https://www.iucnredlist.org/}{IUCN}.
#'
#' @format A \code{data.frame} object with 2488 rows and 10 variables:
#' \describe{
#'   \item{X_WGS84}{Longitude}
#'   \item{Y_WGS84}{Latitude}
#'   \item{ConnochaetesGnou}{Presence (1) or Absence (0) for black wildebeest}
#'   \item{GuloGulo}{Presence (1) or Absence (0) for wolverine}
#'   \item{PantheraOnca}{Presence (1) or Absence (0) for jaguar}
#'   \item{PteropusGiganteus}{Presence (1) or Absence (0) for indian flying fox}
#'   \item{TenrecEcaudatus}{Presence (1) or Absence (0) for tailless tenrec}
#'   \item{VulpesVulpes}{Presence (1) or Absence (0) for red fox}
#' }

"DataSpecies"

# DataSpecies  <-
#   read.csv(
#     "../biomod2_old_inst_folder/external/species/mammals_table.csv",
#     row.names = 1)
# 
# usethis::use_data(DataSpecies, overwrite = TRUE)


#' Bioclimatic variables for SDM based on current condition
#'
#' A \code{\link[terra:rast]{SpatRaster}} with 5 bioclimatic variables commonly
#' used for SDM and describing current climate. Additional information available
#' at \href{https://www.worldclim.org/data/bioclim.html}{worldclim}
#'
#' @format A \code{\link[terra:rast]{SpatRaster}} with 5 layers:
#' \describe{
#'   \item{bio3}{Isothermality}
#'   \item{bio4}{Temperature Seasonality}
#'   \item{bio7}{Temperature Annual Range}
#'   \item{bio11}{Mean Temperature of Coldest Quarter}
#'   \item{bio12}{Annual Precipitation}
#' }

"bioclim_current"

# myFiles <- paste0('../biomod2_old_inst_folder/external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
# need to go through the use of raster because grd files are stored with 
# FLT4S datatype which throws an error when directly loaded with terra::rast
# for an unknown reason
# bioclim_current <-
#   wrap( 
#     terra::rast(
#       raster::stack(myFiles)
#     )
#   )
# usethis::use_data(bioclim_current, overwrite = TRUE )

#' Bioclimatic variables for SDM based on future condition
#'
#' A \code{\link[terra:rast]{SpatRaster}} with 5 bioclimatic variables commonly
#' used for SDM and describing future climate based on old RCP scenarios at the
#' horizon 2080.
#' 
#' @format A \code{\link[terra:rast]{SpatRaster}} with 5 layers:
#' \describe{
#'   \item{bio3}{Isothermality}
#'   \item{bio4}{Temperature Seasonality}
#'   \item{bio7}{Temperature Annual Range}
#'   \item{bio11}{Mean Temperature of Coldest Quarter}
#'   \item{bio12}{Annual Precipitation}
#' }

"bioclim_future"

# myFiles <- paste0('../biomod2_old_inst_folder/external/bioclim/future/bio', c(3, 4, 7, 11, 12), '.grd')
# # need to go through the use of raster because grd files are stored with
# # FLT4S datatype which throws an error when directly loaded with terra::rast
# # for an unknown reason
# bioclim_future <-
#   terra::wrap(
#     terra::rast(
#       raster::stack(myFiles)
#     )
#   )
# usethis::use_data(bioclim_future, overwrite = TRUE)
