###################################################################################################
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
#'   \item DNN (\code{\link[cito]{cito}})
#'   \item FDA (\code{\link[mda]{fda}})
#'   \item GAM (\code{\link[gam]{gam}}, \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}})
#'   \item GBM (\code{\link[gbm]{gbm}})
#'   \item GLM (\code{\link[stats]{glm}})
#'   \item MARS (\code{\link[earth]{earth}})
#'   \item MAXENT (\href{https://biodiversityinformatics.amnh.org/open_source/maxent/}{see Maxent website})
#'   \item MAXNET (\code{\link[maxnet]{maxnet}})
#'   \item RF (\code{\link[randomForest]{randomForest}})
#'   \item RFd (\code{\link[randomForest]{randomForest} downsampled})
#'   \item SRE (\code{\link{bm_SRE}})
#'   \item XGBOOST (\code{\link[xgboost]{xgboost}})
#' }
#' 
###################################################################################################

"ModelsTable"

# ModelsTable <- data.frame(model = c('ANN', 'CTA', 'DNN', 'FDA', 'GAM', 'GAM', 'GAM', 'GBM', 'GLM'
#                                     , 'MARS', 'MAXENT', 'MAXNET', 'RF','RFd', 'SRE', 'XGBOOST', 
#                                     'CTA', 'DNN', 'FDA', 'GAM', 'GAM', 'GAM', 'GBM', 'GLM'
#                                     , 'MARS', 'RF', 'XGBOOST')
#                           , type = c(rep('binary',16), rep('nonbinary',11))
#                           , package = c('nnet', 'rpart', 'cito', 'mda', 'gam', 'mgcv', 'mgcv', 'gbm', 'stats'
#                                         , 'earth', 'MAXENT', 'maxnet', 'randomForest','randomForest', 'biomod2', 'xgboost', 
#                                         'rpart', 'cito', 'mda', 'gam', 'mgcv', 'mgcv', 'gbm', 'stats'
#                                         , 'earth', 'randomForest', 'xgboost')
#                           , func = c('nnet', 'rpart', 'dnn', 'fda', 'gam', 'bam', 'gam', 'gbm', 'glm'
#                                      , 'earth', 'MAXENT', 'maxnet', 'randomForest','randomForest', 'bm_SRE', 'xgboost',
#                                      'rpart', 'dnn', 'fda', 'gam', 'bam', 'gam', 'gbm', 'glm'
#                                      , 'earth', 'randomForest', 'xgboost')
#                           , train = c('avNNet', 'rpart', 'tune', 'fda', 'gamLoess', 'bam', 'gam', 'gbm', 'glm'
#                                       , 'earth', 'ENMevaluate', 'maxnet', 'rf','rf', 'bm_SRE', 'xgbTree',
#                                       'rpart', 'tune', 'fda', 'gamLoess', 'bam', 'gam', 'gbm', 'glm'
#                                       , 'earth', 'rf', 'xgbTree'))

# usethis::use_data(ModelsTable, overwrite = TRUE)
# usethis::use_data(ModelsTable, overwrite = TRUE, internal = TRUE)


###################################################################################################
#' Bigboss pre-defined parameter values for single models
#'
#' A \code{\link{BIOMOD.models.options}} object containing for each single model available in 
#' \pkg{biomod2} the parameter values pre-defined by \pkg{biomod2} team.
#'
#' @format A \code{\link{BIOMOD.models.options}} object with some changed values :
#' 
#' \describe{
#'    \item{\code{ANN.nnet.nnet}}{
#'      \itemize{
#'        \item \code{size = 5}
#'        \item \code{decay = 0.1}
#'        \item \code{trace = FALSE}
#'        \item \code{rang = 0.1}
#'        \item \code{maxit = 200}
#'      }
#'    }
#'    \item{\code{CTA.rpart.rpart}}{
#'      \itemize{
#'        \item \code{control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 10)}
#'        \item \code{cost = NULL}
#'      }
#'    }
#'    \item{\code{DNN.cito.dnn}}{
#'      \itemize{
#'        \item \code{batchsize = 100L}
#'        \item \code{epochs = 150L}
#'        \item \code{hidden = c(100L, 100L) }
#'        \item \code{lr = 0.05}
#'        \item \code{optimizer = "adam"}
#'        \item \code{lambda = 0.001}
#'        \item \code{alpha = 1.0}
#'        \item \code{validation = 0.2}
#'        \item \code{lr_scheduler = config_lr_scheduler("reduce_on_plateau", patience = 7)}
#'        \item \code{early_stopping = 14}
#'      }
#'    }
#'    \item{\code{FDA.mda.fda}}{
#     \itemize{
#        \item ...
#     }
#'    }
#'    \item{\code{GAM.gam.gam}}{
# \itemize{
#   \item ...
# }
#'    }
#'    \item{\code{GAM.mgcv.bam}}{
# \itemize{
#   \item ...
# }
#'    }
#'    \item{\code{GAM.mgcv.gam}}{
#'      \itemize{
#'        \item \code{method = 'GCV.Cp'}
#'        \item \code{control = list(epsilon = 1e-06, trace = FALSE, maxit = 100)}
#'      }
#'    }
#'    \item{\code{GBM.gbm.gbm}}{
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
#'    \item{\code{GLM.stats.glm}}{
#'      \itemize{
#'        \item \code{mustart = 0.5}
#'        \item \code{control = glm.control(maxit = 50)}
#'      }
#'    }
#'    \item{\code{MARS.earth.earth}}{
#'      \itemize{
#'        \item \code{ncross = 0}
#'        \item \code{nk = NULL}
#'        \item \code{penalty = 2}
#'        \item \code{thresh = 0.001}
#'        \item \code{nprune = NULL}
#'        \item \code{pmethod = 'backward'}
#'      }
#'    }
#'    \item{\code{MAXENT.MAXENT.MAXENT}}{
#'      \itemize{
#'        \item \code{path_to_maxent.jar = '.'}
#'      }
#'    }
#'    \item{\code{RF.randomForest.randomForest}}{
#'      \itemize{
#'        \item \code{ntree = 500}
#'        \item \code{mtry = 2}
#'        \item \code{sampsize = NULL}
#'        \item \code{nodesize = 5}
#'        \item \code{maxnodes = NULL}
#'      }
#'    }
#'    \item{\code{RFd.randomForest.randomForest}}{
#'      \itemize{
#'        \item \code{type = 'classification'}
#'        \item \code{ntree = 500}
#'        \item \code{mtry = 2}
#'        \item \code{strata = factor(c(0, 1))}
#'        \item \code{sampsize = NULL}
#'        \item \code{nodesize = 5}
#'        \item \code{maxnodes = NULL}
#'      }
#'    }
#'    \item{\code{SRE.biomod2.bm_SRE}}{
#'      \itemize{
#'        \item \code{do.extrem = TRUE}
#'      }
#'    }
#'    \item{\code{XGBOOST.xgboost.xgboost}}{
#'      \itemize{
#'        \item \code{params = list(max_depth = 2, eta = 1)}
#'        \item \code{nthread = 2}
#'        \item \code{nrounds = 4}
#'      }
#'    }
#' }
#' 
###################################################################################################

"OptionsBigboss"

# bm.opt <- bm_ModelingOptions(data.type = "binary", strategy = "default")
# bm.opt_gam2 <- bm_ModelingOptions(data.type = "binary", strategy = "default", models = "GAM.gam.gam")
# bm.opt_gam3 <- bm_ModelingOptions(data.type = "binary", strategy = "default", models = "GAM.mgcv.bam")
# bm.opt@models <- sort(unique(c(bm.opt@models, bm.opt_gam2@models, bm.opt_gam3@models)))
# bm.opt@options <- c(bm.opt@options, bm.opt_gam2@options, bm.opt_gam3@options)
# bm.opt@options <- bm.opt@options[bm.opt@models]
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$size = 5 #NULL
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$decay = 0.1
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$trace = FALSE
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$rang = 0.1
# bm.opt@options$ANN.binary.nnet.nnet@args.values[['_allData_allRun']]$maxit = 200
# # bm.opt@options$CTA.binary.rpart.rpart@args.values[['_allData_allRun']]$method = "class"
# bm.opt@options$CTA.binary.rpart.rpart@args.values[['_allData_allRun']]$control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 10)
# bm.opt@options$CTA.binary.rpart.rpart@args.values[['_allData_allRun']]$cost = NULL
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$batchsize = 100L
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$epochs = 150L
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$hidden = c(100L, 100L) ## To wide ?
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$lr = 0.05
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$optimizer = "adam"
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$lambda = 0.001
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$alpha = 1.0 #L2 regularization
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$validation = 0.2
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$lr_scheduler = config_lr_scheduler("reduce_on_plateau", patience = 7)
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]$early_stopping = 14
# # bm.opt@options$FDA.binary.mda.fda@args.values[['_allData_allRun']]$method = "mars"
# # bm.opt@options$GAM.binary.mgcv.gam@args.values[['_allData_allRun']]$family = binomial(link = 'logit')
# bm.opt@options$GAM.binary.mgcv.gam@args.values[['_allData_allRun']]$method = "GCV.Cp"
# bm.opt@options$GAM.binary.mgcv.gam@args.values[['_allData_allRun']]$control = list(epsilon = 1e-06, trace = FALSE, maxit = 100)
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$n.trees = 2500
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$interaction.depth = 7
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$n.minobsinnode = 5
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$shrinkage = 0.001
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$cv.folds = 3
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$keep.data = FALSE
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]$n.cores = 1
# # bm.opt@options$GLM.binary.stats.glm@args.values[['_allData_allRun']]$family = binomial(link = 'logit')
# bm.opt@options$GLM.binary.stats.glm@args.values[['_allData_allRun']]$mustart = 0.5
# bm.opt@options$GLM.binary.stats.glm@args.values[['_allData_allRun']]$control = glm.control(maxit = 50)
# # bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$glm = list(family = binomial(link = 'logit'))
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$ncross = 0
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$nk = NULL
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$penalty = 2
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$thresh = 0.001
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$nprune = NULL
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]$pmethod = 'backward'
# bm.opt@options$MAXENT.binary.MAXENT.MAXENT@args.default$path_to_maxent.jar = '.'
# bm.opt@options$MAXENT.binary.MAXENT.MAXENT@args.values[['_allData_allRun']]$path_to_maxent.jar = '.'
# # bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$type = 'classification'
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$ntree = 500
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$mtry = 2
# # bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$strata = factor(c(0, 1))
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$sampsize = NULL
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$nodesize = 5
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]$maxnodes = NULL
# bm.opt@options$RFd.binary.randomForest.randomForest@args.values[['_allData_allRun']]$type = 'classification'
# bm.opt@options$RFd.binary.randomForest.randomForest@args.values[['_allData_allRun']]$ntree = 500
# bm.opt@options$RFd.binary.randomForest.randomForest@args.values[['_allData_allRun']]$mtry = 2
# bm.opt@options$RFd.binary.randomForest.randomForest@args.values[['_allData_allRun']]$strata = factor(c(0, 1))
# bm.opt@options$RFd.binary.randomForest.randomForest@args.values[['_allData_allRun']]$sampsize = NULL
# bm.opt@options$RFd.binary.randomForest.randomForest@args.values[['_allData_allRun']]$nodesize = 5
# bm.opt@options$RFd.binary.randomForest.randomForest@args.values[['_allData_allRun']]$maxnodes = NULL
# bm.opt@options$SRE.binary.biomod2.bm_SRE@args.values[['_allData_allRun']]$do.extrem = TRUE
# bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']]$params = list(max_depth = 2, eta = 1)
# bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']]$nthread = 2
# bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']]$nrounds = 4
# # bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']]$objective = "binary:logistic"

# # Remove things to be adapted to "count", "abundance" and "compositional"
# bm.opt@options$CTA.binary.rpart.rpart@args.values[['_allData_allRun']] <- bm.opt@options$CTA.binary.rpart.rpart@args.values[['_allData_allRun']][-which(names(bm.opt@options$CTA.binary.rpart.rpart@args.values[['_allData_allRun']]) == "method")]
# bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']] <- bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']][-which(names(bm.opt@options$DNN.binary.cito.dnn@args.values[['_allData_allRun']]) == "loss")]
# bm.opt@options$FDA.binary.mda.fda@args.values[['_allData_allRun']] <- bm.opt@options$FDA.binary.mda.fda@args.values[['_allData_allRun']][-which(names(bm.opt@options$FDA.binary.mda.fda@args.values[['_allData_allRun']]) == "method")]
# bm.opt@options$GAM.binary.mgcv.gam@args.values[['_allData_allRun']] <- bm.opt@options$GAM.binary.mgcv.gam@args.values[['_allData_allRun']][-which(names(bm.opt@options$GAM.binary.mgcv.gam@args.values[['_allData_allRun']]) == "family")]
# bm.opt@options$GAM.binary.mgcv.bam@args.values[['_allData_allRun']] <- bm.opt@options$GAM.binary.mgcv.bam@args.values[['_allData_allRun']][-which(names(bm.opt@options$GAM.binary.mgcv.bam@args.values[['_allData_allRun']]) == "family")]
# bm.opt@options$GAM.binary.gam.gam@args.values[['_allData_allRun']] <- bm.opt@options$GAM.binary.gam.gam@args.values[['_allData_allRun']][-which(names(bm.opt@options$GAM.binary.gam.gam@args.values[['_allData_allRun']]) == "family")]
# bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']] <- bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']][-which(names(bm.opt@options$GBM.binary.gbm.gbm@args.values[['_allData_allRun']]) == "distribution")]
# bm.opt@options$GLM.binary.stats.glm@args.values[['_allData_allRun']] <- bm.opt@options$GLM.binary.stats.glm@args.values[['_allData_allRun']][-which(names(bm.opt@options$GLM.binary.stats.glm@args.values[['_allData_allRun']]) == "family")]
# bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']] <- bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']][-which(names(bm.opt@options$MARS.binary.earth.earth@args.values[['_allData_allRun']]) == "glm")]
# bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']] <- bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']][-which(names(bm.opt@options$RF.binary.randomForest.randomForest@args.values[['_allData_allRun']]) == "type")]
# bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']] <- bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']][-which(names(bm.opt@options$XGBOOST.binary.xgboost.xgboost@args.values[['_allData_allRun']]) == "objective")]

# bm.opt@models <- sub(".binary", "", bm.opt@models)
# names(bm.opt@options) <- bm.opt@models

# OptionsBigboss <- bm.opt

# usethis::use_data(OptionsBigboss, overwrite = TRUE)
# usethis::use_data(OptionsBigboss, overwrite = TRUE, internal = TRUE)


###################################################################################################
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
#'   \item{PteropusGiganteus}{Presence (1) or Absence (0) for Indian flying fox}
#'   \item{TenrecEcaudatus}{Presence (1) or Absence (0) for tailless tenrec}
#'   \item{VulpesVulpes}{Presence (1) or Absence (0) for red fox}
#' }
#' 
###################################################################################################

"DataSpecies"

# DataSpecies  <-
#   read.csv(
#     "../biomod2_old_inst_folder/external/species/mammals_table.csv",
#     row.names = 1)
# 
# usethis::use_data(DataSpecies, overwrite = TRUE)


###################################################################################################
#' Abundance to build test SDM
#'
#' A dataset covering France with abundance data for 20 bird species. Data comes from the French 
#' Breeding Bird Survey (STOC program) of Vigie-Nature (Fontaine et al. 2020). The original 
#' dataset contains 2,948 sites in which bird species have been observed from 2001 to 2024. 
#' Original raw data is not provided here, but might be made available on request by contacting 
#' \href{benoit.fontaine@mnhn.fr}{Benoît Fontaine}.
#'
#' @format A \code{data.frame} object with 2006 rows and 30 columns (10 variables and 20 species):
#' \describe{
#'   \item{site}{Observation site ID}
#'   \item{period}{yearly average headcounts have been computed over 2 times periods: 
#'   2006-2011 and 2012-2017}
#'   \item{X_WGS84}{Longitude}
#'   \item{Y_WGS84}{Latitude}
#'   \item{temp}{annual mean temperature (CHELSA)}
#'   \item{precip}{annual mean precipitation (CHELSA)}
#'   \item{cover_agri}{percentage of agricultural habitat cover around each point 
#'   (European Union’s Copernicus Land Monitoring Service information)}
#'   \item{cover_water}{percentage of water habitat cover around each point 
#'   (European Union’s Copernicus Land Monitoring Service information)}
#'   \item{cover_wet}{percentage of wetland habitat cover around each point 
#'   (European Union’s Copernicus Land Monitoring Service information)}
#'   \item{sdiv_hab}{Shannon habitat diversity (CORINE)}
#'   \item{Alauda.arvensis}{abundance data for Eurasian skylark}
#'   \item{Cettia.cetti}{abundance data for Cetti's warbler}
#'   \item{Coccothraustes.coccothraustes}{abundance data for hawfinch}
#'   \item{Cuculus.canorus}{abundance data for common cuckoo}
#'   \item{Emberiza.calandra}{abundance data for corn bunting}
#'   \item{Emberiza.citrinella}{abundance data for yellowhammer}
#'   \item{Erithacus.rubecula}{abundance data for European robin}
#'   \item{Fulica.atra}{abundance data for Eurasian coot}
#'   \item{Luscinia.megarhynchos}{abundance data for common nightingale}
#'   \item{Passer.domesticus}{abundance data for house sparrow}
#'   \item{Perdix.perdix}{abundance data for grey partridge}
#'   \item{Periparus.ater}{abundance data for coal tit}
#'   \item{Phylloscopus.bonelli}{abundance data for western Bonelli's warbler}
#'   \item{Phylloscopus.trochilus}{abundance data for willow warbler}
#'   \item{Regulus.regulus}{abundance data for goldcrest}
#'   \item{Serinus.serinus}{abundance data for European serin}
#'   \item{Sitta.europaea}{abundance data for Eurasian nuthatch}
#'   \item{Streptopelia.decaocto}{abundance data for Eurasian collared dove}
#'   \item{Sylvia.melanocephala}{abundance data for Sardinian warbler}
#'   \item{Turdus.philomelos}{abundance data for song thrush}
#' }
#' 
#' @references
#' 
#' \itemize{
#'   \item \href{https://www.vigie-plume.fr/}{STOC EPS} - Vigie-Nature (2025). 
#'   \emph{French Breeding Bird Monitoring Scheme.} Muséum National d’Histoire Naturelle 
#'   (MNHN) - Office Français pour la Biodiversité (OFB) - Ligue pour la Protection des 
#'   Oiseaux (LPO).
#'   \item Fontaine B, Moussy C, Chiffard Carricaburu J, Dupuis J, Corolleur E, Schmaltz L, 
#'   Lorrillière R, Loïs G, Gaudard C. (2020) Suivi des oiseaux communs en France 
#'   1989-2019 : 30 ans de suivis participatifs. MNHN - Centre d'Ecologie et des Sciences 
#'   de la Conservation, LPO BirdLife France - Service Connaissance, Ministère de la 
#'   Transition écologique et solidaire. 46 pp.
#'   \item  Brun P, Zimmermann NE, Hari C, Pellissier L, Karger DN. (2022) 
#'   CHELSA-BIOCLIM+ A novel set of global climate-related predictors at 
#'   kilometre-resolution. EnviDat. 
#'   DOI: \href{https://www.doi.org/10.16904/envidat.332}{10.16904/envidat.332}
#'   \item CORINE Land Cover 2018 (raster 100 m), Europe, 6-yearly - version 2020_20u1, 
#'   May 2020. European Environment Agency. 
#'   DOI: \href{https://doi.org/10.2909/960998c1-1870-4e82-8051-6485205ebbac}{10.2909/960998c1-1870-4e82-8051-6485205ebbac}
#' }
#' 
###################################################################################################

"DataSTOC"

# DataSTOC <- read.csv("inst/external/DATA_biomod2_STOC.csv", header = TRUE, sep = ",")
# DataSTOC <- cbind(DataSTOC[, 1:10], DataSTOC[, sort(colnames(DataSTOC)[11:ncol(DataSTOC)])])
# colnames(DataSTOC)[3:4] = toupper(colnames(DataSTOC)[3:4])
# 
# usethis::use_data(DataSTOC, overwrite = TRUE)
# usethis::use_data(DataSTOC, overwrite = TRUE, internal = TRUE)


###################################################################################################
#' Bioclimatic variables for SDM based on current condition
#'
#' A \code{\link[terra:rast]{SpatRaster}} with 5 bioclimatic variables commonly
#' used for SDM and describing current climate. Additional information available
#' at \href{https://www.worldclim.org/data/bioclim.html}{worldclim}.
#'
#' @format A \code{\link[terra:rast]{SpatRaster}} with 5 layers:
#' \describe{
#'   \item{bio3}{Isothermality}
#'   \item{bio4}{Temperature Seasonality}
#'   \item{bio7}{Temperature Annual Range}
#'   \item{bio11}{Mean Temperature of Coldest Quarter}
#'   \item{bio12}{Annual Precipitation}
#' }
#' 
###################################################################################################

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


###################################################################################################
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
#' 
###################################################################################################

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


###################################################################################################
#' ODMAP empty table
#'
#' A \code{data.frame} containing ODMAP (Overview, Data, Model, Assessment, Prediction) 
#' protocol components.
#'
#' @format A \code{data.frame} object with 84 rows and 4 variables:
#' \describe{
#'   \item{section}{Overview, Data, Model, Assessment or Prediction step}
#'   \item{subsection}{corresponding field}
#'   \item{element}{corresponding field}
#'   \item{value}{to be filled with \code{\link{BIOMOD_Report}} function}
#' }
#' 
###################################################################################################

"ODMAP"

# ODMAP <- data.frame(section = c(rep("Overview", 25)
#                                 , rep("Data", 33)
#                                 , rep("Model", 14)
#                                 , rep("Assessment", 5)
#                                 , rep("Prediction", 7))
#                     , subsection = c(rep("Authorship", 4)
#                                      , rep("Model objective", 2)
#                                      , "Focal Taxon", "Location"
#                                      , rep("Scale of Analysis", 5)
#                                      , rep("Biodiversity data", 2)
#                                      , "Predictors", "Hypotheses", "Assumptions"
#                                      , rep("Algorithms", 3)
#                                      , "Workflow"
#                                      , rep("Software", 3)
#                                      , rep("Biodiversity data", 12)
#                                      , rep("Data partitioning", 3)
#                                      , rep("Predictor variables", 10)
#                                      , rep("Transfer data", 8)
#                                      , "Variable pre-selection", "Multicollinearity"
#                                      , rep("Model settings", 2)
#                                      , rep("Model estimates", 3)
#                                      , rep("Model selection - model averaging - ensembles", 3)
#                                      , rep("Analysis and Correction of non-independence", 3)
#                                      , "Threshold selection"
#                                      , rep("Performance statistics", 3)
#                                      , rep("Plausibility check", 2)
#                                      , rep("Prediction output", 2)
#                                      , rep("Uncertainty quantification", 5))
#                     , element = c("Study title", "Author names", "Contact", "Study link"
#                                   , "Model objective", "Target output", "Focal Taxon"
#                                   , "Location", "Spatial extent", "Spatial resolution"
#                                   , "Temporal extent", "Temporal resolution", "Boundary"
#                                   , "Observation type", "Response data type", "Predictor types"
#                                   , "Hypotheses", "Model assumptions", "Modelling techniques"
#                                   , "Model complexity", "Model averaging", "Model workflow"
#                                   , "Software", "Code availability", "Data availability" ## OVERVIEW
#                                   , "Taxon names", "Taxonomic reference system", "Ecological level"
#                                   , "Data sources", "Sampling design", "Sample size"
#                                   , "Clipping", "Scaling", "Cleaning"
#                                   , "Absence data", "Background data", "Errors and biases"
#                                   , "Training data", "Validation data", "Test data"
#                                   , "Predictor variables"
#                                   , "Data sources", "Spatial extent", "Spatial resolution"
#                                   , "Coordinate reference system"
#                                   , "Temporal extent", "Temporal resolution"
#                                   , "Data processing", "Errors and biases", "Dimension reduction"
#                                   , "Data sources", "Spatial extent", "Spatial resolution"
#                                   , "Temporal extent", "Temporal resolution"
#                                   , "Models and scenarios", "Data processing", "Quantification of Novelty" ## DATA
#                                   , "Variable pre-selection", "Multicollinearity"
#                                   , "Model settings (fitting)", "Model settings (extrapolation)"
#                                   , "Coefficients", "Parameter uncertainty", "Variable importance"
#                                   , "Model selection", "Model averaging", "Model ensembles"
#                                   , "Spatial autocorrelation", "Temporal autocorrelation"
#                                   , "Nested data", "Threshold selection" ## MODEL
#                                   , "Performance on training data", "Performance on validation data"
#                                   , "Performance on test data", "Response shapes", "Expert judgement" ## ASSESSMENT
#                                   , "Prediction unit", "Post-processing"
#                                   , "Algorithmic uncertainty", "Input data uncertainty"
#                                   , "Parameter uncertainty", "Scenario uncertainty"
#                                   , "Novel environments") ## PREDICTION
#                     , value = NA)

# usethis::use_data(ODMAP, overwrite = TRUE)
# usethis::use_data(ODMAP, overwrite = TRUE, internal = TRUE)
