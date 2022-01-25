
## BIOMOD SPECIFIC CAT ----------------------------------------------------------------------------

.bm_cat <- function(x = NULL, ...)
{
  if (is.null(x))
  {
    cat("\n")
    cat(paste(rep("-=", round(.Options$width / 2)), collapse = ""))
    cat("\n")
  } else
  {
    x.length = nchar(x) + 2
    y.length = (.Options$width - x.length) / 2
    cat("\n")
    cat(paste(rep("-=", round(y.length / 2)), collapse = "")
        , x
        , paste(rep("-=", round(y.length / 2)), collapse = ""))
    cat("\n")
  }
}


## TEST PARAMETERS (BIOMOD_ModelingOptions) -------------------------------------------------------
## used in biomod2_classes_1.R file

.fun_testIfInherits <- function(test, objName, objValue, values)
{
  if (!inherits(objValue, values)) {
    stop(paste0("\n", paste0(objName, " must be a '", paste0(values[1:(length(values) -1)], collapse = "', '")
                             , ifelse(length(values) > 1, paste0("' or '", values[length(values)]), "")
                             , "' object")))
    test <- FALSE
  }
  return(test)
}

.fun_testIfIn <- function(test, objName, objValue, values)
{
  if (sum(objValue %in% values) < length(objValue)) {
    stop(paste0("\n", paste0(objName, " must be '", paste0(values[1:(length(values) -1)], collapse = "', '")
                             , ifelse(length(values) > 1, paste0("' or '", values[length(values)]))
                             , "'")))
    test <- FALSE
  }
  return(test)
}

.fun_testIfPosNum <- function(test, objName, objValue)
{
  if (!is.numeric(objValue)) {
    stop(paste0("\n", objName, "must be a numeric"))
    test <- FALSE
  } else if (objValue < 0) {
    stop(paste0("\n", objName, "must be a positive numeric"))
    test <- FALSE
  }
  return(test)
}

.fun_testIf01 <- function(test, objName, objValue)
{
  test <- .fun_testIfPosNum(test, objName, objValue)
  if (test && objValue > 1) {
    stop(paste0("\n", objName, "must be a 0 to 1 numeric"))
    test <- FALSE
  }
  return(test)
}

.fun_testIfPosInt <- function(test, objName, objValue)
{
  if (!is.numeric(objValue)) {
    cat(paste0("\n", objName, "must be a integer"))
    test <- FALSE
  } else if (objValue < 0 | objValue %% 1 != 0) {
    cat(paste0("\n", objName, "must be a positive integer"))
    test <- FALSE
  }
  return(test)
}



## TEMPLATES TO PREDICT MODELS (BIOMOD_Modeling) --------------------------------------------------
## used in biomod2_classes_4.R file

# Fuction to get variables ranges
get_var_type <- function(data) { return(sapply(data, class)) }

get_var_range <- function(data)
{
  get_range <- function(x) {
    if (is.numeric(x)) {
      return(c(min = min(x, na.rm = T), max = max(x, na.rm = T)))
    } else if (is.factor(x)) {
      return(levels(x))
    }
  }
  xx <- lapply(data, get_range)
  names(xx) <- names(data)
  return(xx)
}

# Function to check new data range compatibility with calibrating data #
check_data_range <- function(model, new_data)
{
  ## TODO : remettre en marche cette fonction
  
  #   # get calibration data caracteristics
  #   expl_var_names <- model@expl_var_names
  #   expl_var_type <- model@expl_var_type
  #   expl_var_range <- model@expl_var_range
  #
  #   if(inherits(new_data, "Raster")){ ## raster data case =-=-=-=-=-=-=- #
  #     # check var names compatibility
  #     nd_expl_var_names <- names(new_data)
  #     if(sum(!(expl_var_names %in% nd_expl_var_names) ) > 0 ){
  #       stop("calibration and projections variables names mismatch")
  #     }
  #     # reorder the stack
  #     new_data <- raster::subset(new_data,expl_var_names)
  #     # check var types compatibility (factors)
  #     expl_var_fact <- (expl_var_type=='factor')
  #     nd_expl_var_fact <- is.factor(new_data)
  #     if(sum(! (expl_var_fact==nd_expl_var_fact))>0){
  #       stop("calibration and projections variables class mismatch")
  #     }
  #     # check var range compatibility
  #     ### remove all new factors
  #     if(sum(expl_var_fact)>0){ ## there are factorial variables
  #       for(fact_var_id in which(expl_var_fact)){
  #         ## check if new factors occurs
  #         nd_levels <- levels(raster::subset(new_data,fact_var_id))[[1]]
  #         nd_levels <- as.character(nd_levels[,ncol(nd_levels)])
  #         names(nd_levels) <- levels(raster::subset(new_data,fact_var_id))[[1]]$ID
  #         cd_levels <- as.character(unlist(expl_var_range[[fact_var_id]]))
  #
  #         ## detect new levels
  #         new_levels <- nd_levels[!(nd_levels %in% cd_levels)]
  #
  #         if(length(new_levels)){
  #           for(n_l in new_levels){
  #             # remove points where out of range factors have been detected
  #             new_data[subset(new_data,fact_var_id)[]==as.numeric(names(nd_levels)[which(nd_levels==n_l)])] <- NA
  #           }
  #           warning(paste(nd_expl_var_names[fact_var_id]," new levels have been removed from dataset (",toString(new_levels),")",sep=""))
  #         }
  #       }
  #     }
  #     ## convert data to be sure to get RasterStack output
  # #     new_data <- stack(new_data)
  #   } else{ ## table data case -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  #     # check var names compatibility
  #     nd_expl_var_names <- colnames(new_data)
  #     if(sum(!(expl_var_names %in% nd_expl_var_names) ) > 0 ){
  #       stop("calibration and projections variables names mismatch")
  #     }
  #     # reorder the stack
  #     new_data <- new_data[,expl_var_names, drop=F]
  #     # check var types compatibility (factors)
  #     expl_var_fact <- (expl_var_type=='factor')
  #     nd_expl_var_fact <- sapply(new_data,is.factor)
  #
  #     if(sum(! (expl_var_fact==nd_expl_var_fact))>0){
  #       stop("calibration and projections variables class mismatch")
  #     }
  #     # check var range compatibility
  #     ### remove all new factors
  #     if(sum(expl_var_fact)>0){ ## there are factorial variables
  #       for(fact_var_id in which(expl_var_fact)){
  #         ## check if new factors occurs
  #         nd_levels <- levels(new_data[,fact_var_id])
  #         cd_levels <- as.character(unlist(expl_var_range[[fact_var_id]]))
  #
  #         ## detect new levels
  #         new_levels <- nd_levels[!(nd_levels %in% cd_levels)]
  #
  #         if(length(new_levels)){
  #           # remove points where out of range factors have been detected
  # #           new_data <- new_data[- which(new_data[,fact_var_id] %in% new_levels),]
  #           new_data[which(new_data[,fact_var_id] %in% new_levels),] <- NA
  #           warning(paste(nd_expl_var_names[fact_var_id]," new levels have been removed from dataset (",toString(new_levels),")",sep=""))
  #         }
  #       }
  #     }
  #   }
  return(new_data)
}

.run_pred <- function(object, Prev = 0.5 , dat)
{
  if (is.finite(object$deviance) && 
      is.finite(object$null.deviance) && 
      object$deviance != object$null.deviance)
  {
    if (inherits(dat, 'Raster')) {
      pred <- predict(object = dat, model = object, type = "response")
    } else {
      pred <- predict(object, dat, type = "response")
    }
  }
  
  if (!exists('pred')) {
    if (inherits(dat, 'Raster')) {
      pred <- subset(dat, 1, drop = TRUE)
      if (Prev < 0.5) {
        pred <- reclassify(x = pred, rcl = c(-Inf, Inf, 0))
      } else {
        pred <- reclassify(x = pred, rcl = c(-Inf, Inf, 1))
      }
    } else {
      if (Prev < 0.5) {
        pred <- rep(0, nrow(dat))
      } else {
        pred <- rep(1, nrow(dat))
      }
    }
  }
  
  return(pred)
}

## FOR SINGLE MODELS ----------------------------------------------------------
.template_predict = function(mod, object, newdata, ...)
{
  args <- list(...)
  do_check <- args$do_check
  if (is.null(do_check)) { do_check <- TRUE }
  if (do_check) { newdata <- check_data_range(model = object, new_data = newdata) }
  
  ## Special case for MAXENT.Phillips
  do_raster = FALSE
  newproj = NA
  if (mod == "MAXENT.Phillips" && inherits(newdata, 'Raster')) {
    newraster = newdata[[1]]
    newraster[] = NA
    newdata = na.exclude(rasterToPoints(newdata))
    newraster[cellFromXY(newraster, newdata[, 1:2])] = 1
    do_raster = TRUE
  }
  
  if (inherits(newdata, 'Raster')) {
    eval(parse(text = paste0("res = .predict.", mod, "_biomod2_model.RasterStack(object, newdata, ...)")))
    return(res)
  } else if (inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')) {
    eval(parse(text = paste0("res = .predict.", mod, "_biomod2_model.data.frame(object, newdata, do_raster = do_raster, newraster = newraster, ...)")))
    return(res)
  } else {
    stop("invalid newdata input")
  }
}

.template_predict.RasterStack = function(seedval = NULL, predcommand, object, newdata, ...)
{
  args <- list(...)
  namefile <- args$namefile
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) { overwrite <- TRUE }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  set.seed(seedval)
  eval(parse(text = paste0("proj <- ", predcommand)))
  
  if (length(get_scaling_model(object))) {
    names(proj) <- "pred"
    proj <- .run_pred(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }
  if (on_0_1000) { proj <- round(proj * 1000) }
  
  # save raster on hard drive ?
  if (!is.null(namefile)) {
    cat("\n\t\tWriting projection on hard drive...")
    if (on_0_1000) { ## projections are stored as positive integer
      writeRaster(proj, filename = namefile, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename = namefile, overwrite = overwrite)
    }
    proj <- raster(namefile, RAT = FALSE)
  }
  return(proj)
}

.template_predict.data.frame <- function(seedval = NULL, predcommand, object, newdata, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na
  
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  if (is.null(omit.na)) { omit.na <- FALSE }
  
  ## check if na occurs in newdata cause they are not well supported
  if (omit.na) {
    not_na_rows <- apply(newdata, 1, function(x) { sum(is.na(x)) == 0 })
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }
  
  set.seed(seedval)
  eval(parse(text = paste0("proj <- ", predcommand)))
  
  ## add original NA from formal dataset
  if (sum(!not_na_rows) > 0) {
    tmp <- rep(NA, length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }
  
  if (length(get_scaling_model(object))) {
    proj <- data.frame(pred = proj)
    proj <- .run_pred(object = get_scaling_model(object), Prev = 0.5, dat = proj)
  }
  if (on_0_1000) { proj <- round(proj * 1000) }
  
  return(proj)
}

## FOR ENSEMBLE MODELS --------------------------------------------------------
.template_predictEM = function(mod, object, newdata, formal_predictions, ...)
{
  args <- list(...)
  do_check <- args$do_check
  if (is.null(do_check)) { do_check <- TRUE }
  if (do_check) { newdata <- check_data_range(model = object, new_data = newdata) }
  
  if (inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')) {
    eval(parse(text = paste0("res = .predict.", mod, "_biomod2_model.RasterStack(object, newdata, formal_predictions, ...)")))
    return(res)
  } else if (inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | 
             inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')) {
    eval(parse(text = paste0("res = .predict.", mod, "_biomod2_model.data.frame(object, newdata, formal_predictions, ...)")))
    return(res)
  } else {
    stop("invalid newdata input")
  }
}

.template_predictEM.formal_predictions = function(object, newdata, formal_predictions, ...)
{
  args <- list(...)
  namefile <- args$namefile
  on_0_1000 <- args$on_0_1000
  
  if (is.null(namefile)) { namefile <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (is.null(formal_predictions)) {
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                 function(mod.name, resp_name, modeling.id)
                                 {
                                   ## check if model is loaded on memory
                                   if (is.character(mod.name)) {
                                     mod <- get(load(file.path(resp_name, "models", modeling.id, mod.name)))
                                   }
                                   temp_workdir = NULL
                                   if (length(grep("MAXENT.Phillips$", mod.name)) == 1) {
                                     temp_workdir = mod@model_output_dir
                                   }
                                   return(predict(mod, newdata = newdata, on_0_1000 = on_0_1000
                                                  , temp_workdir = temp_workdir))
                                 }, resp_name = object@resp_name, modeling.id = object@modeling.id)
  }
  
  return(formal_predictions)
}


## PREPARE and DELETE workdir for MAXENT ----------------------------------------------------------
## used in biomod2_classes_4, bm_RunModelsLoop files
.maxent.prepare.workdir <- function(Data, xy, calibLines = NULL, RunName = NULL,
                                    evalData = NULL, evalxy =  NULL,
                                    species.name = NULL, modeling.id = '',
                                    background_data_dir = 'default')
{
  cat('\n\t\tCreating Maxent Temp Proj Data...')
  
  ## initialise output
  MWD <- list()
  class(MWD) <- "maxent_workdir_info"
  
  ## default parameters setting
  if (is.null(RunName)) { RunName <- colnames(Data)[1] }
  if (is.null(species.name)) { species.name <- colnames(Data)[1] }
  if (is.null(calibLines)) { calibLines <- rep(TRUE, nrow(Data)) }
  
  ## define all paths to files needed by MAXENT.Phillips
  nameFolder = file.path(species.name, 'models', modeling.id)
  m_outdir <- file.path(nameFolder, paste0(RunName, '_MAXENT.Phillips_outputs'))
  m_predictDir <- file.path(m_outdir, "Predictions")
  MWD$m_outdir <- m_outdir
  MWD$m_outputFile <- file.path(m_outdir, paste0(RunName, '_Pred_swd.csv'))
  MWD$m_predictDir <- m_predictDir
  
  ## directories creation
  dir.create(m_outdir, showWarnings = FALSE, recursive = TRUE, mode = '777')
  dir.create(m_predictDir, showWarnings = FALSE, recursive = TRUE, mode = '777')
  
  ## Presence Data --------------------------------------------------------------------------------
  presLines <- which((Data[, 1] == 1) & calibLines)
  absLines <- which((Data[, 1] == 0) & calibLines)
  Sp_swd <- cbind(rep(RunName, length(presLines))
                  , xy[presLines, ]
                  , Data[presLines, 2:ncol(Data), drop = FALSE])
  colnames(Sp_swd) <- c('species', 'X', 'Y', colnames(Data)[2:ncol(Data)])
  
  m_speciesFile <- file.path(m_outdir, "Sp_swd.csv")
  write.table(Sp_swd, file = m_speciesFile, quote = FALSE, row.names = FALSE, sep = ",")
  MWD$m_speciesFile <- m_speciesFile
  
  ## Background Data (create background file only if needed) --------------------------------------
  if (background_data_dir == 'default')
  {
    # keep only 0 of calib lines
    Back_swd <- cbind(rep("background", length(absLines))
                      , xy[absLines, ]
                      , Data[absLines, 2:ncol(Data), drop = FALSE])
    colnames(Back_swd) <- c("background", colnames(Back_swd)[-1])
    
    m_backgroundFile <- file.path(m_outdir, "Back_swd.csv")
    write.table(Back_swd, file = m_backgroundFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
    MWD$m_backgroundFile <- m_backgroundFile
  } else { ## use background directory given as an option
    MWD$m_backgroundFile <- background_data_dir
  }
  
  ## Prediction Data ------------------------------------------------------------------------------
  Pred_swd <- cbind(rep("predict", nrow(xy))
                    , xy
                    , Data[, 2:ncol(Data), drop = FALSE])
  colnames(Pred_swd)  <- c("predict", colnames(xy), colnames(Data)[-1])
  
  m_predictFile <- file.path(m_predictDir, "Pred_swd.csv")
  write.table(Pred_swd, file = m_predictFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
  MWD$m_predictFile <- m_predictFile
  
  ## dealing with independent evaluation data -----------------------------------------------------
  if (!is.null(evalData)) {
    Pred_eval_swd <- cbind(rep("predictEval", nrow(evalxy))
                           , evalxy
                           , evalData[, 2:ncol(evalData), drop = FALSE])
    colnames(Pred_eval_swd) <- c("predict", colnames(Back_swd)[-1])
    
    m_predictEvalFile <- file.path(m_predictDir, "PredEval_swd.csv")
    write.table(Pred_eval_swd, file = m_predictEvalFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
    MWD$m_predictEvalFile <- m_predictEvalFile
  }
  
  return(MWD)
}


## CREATE MODEL FORMULA ---------------------------------------------------------------------------
## used in bm_CVnnet, bm_RunModelsLoop files

.scope <- function(enviroTrain, Smoother, degree)
{
  XXX <- enviroTrain
  deg <- degree
  vnames <- names(XXX[])
  step.list <- as.list(vnames)
  names(step.list) <- vnames
  NbVar <- dim(enviroTrain)[2]
  i <- 1
  while (i <= NbVar)
  {
    vname <- names(XXX)[i]
    # loops through independent variable names
    junk <- paste0("1 + ", vname)
    # minimum scope
    if (is.numeric(XXX[, i])) {
      junk <- c(junk, paste0(Smoother, "(", vname, ",", deg, ")"))
      junk <- eval(parse(text = paste("~", paste(junk, collapse = "+"))))
    } else if (is.factor(XXX[, i])) {
      junk <- c(junk, vname)
      junk <- eval(parse(text = paste("~", paste(junk, collapse = "+"))))
    }
    step.list[[vname]] <- junk
    i <- i + 1
  }
  
  return(step.list)
}


.scope_expSyst <- function(enviroTrain, mod)
{
  i <- 1
  junk2 <- c()
  while (i <= dim(enviroTrain)[2])
  {
    vname <- names(enviroTrain)[i]
    
    if (mod %in% c("NNET", "FDA", "GLMs", "CTA", "GBM")) {
      junk <- vname
    } else if (mod == "GLMq") {
      if (is.numeric(enviroTrain[, i])) {
        junk <- paste0(vname, "+I(", vname, "^2)+I(", vname, "^3)")
      } else if (is.factor(enviroTrain[, i])) {
        junk <- vname
      }
    } else if (mod == "GLMp") {
      if (is.numeric(enviroTrain[, i])) {
        junk <- paste0(vname, "+I(", vname, "^2)+I(", vname, "^3)+", "poly(", vname, ",2) + poly(", vname, ",3)")
      } else if(is.factor(enviroTrain[, i])) {
        junk <- vname
      }
    }
    junk2 <- c(junk2, junk)
    i <- i + 1
  }
  
  junk2 <- eval(parse(text = paste("~", paste(junk2, collapse = "+"))))
  return(junk2)
}


## EXTRACT model names according to specific infos ------------------------------------------------
## used in bm_CVnnet, bm_RunModelsLoop files

.extract_modelNamesInfo <- function(model.names, info = 'species')
{
  if (!is.character(model.names)) { stop("model.names must be a character vector") }
  if (!is.character(info) | length(info) != 1 | !(info %in% c('species', 'data.set', 'models', 'run.eval'))) {
    stop("info must be 'species', 'data.set', 'models' or 'run.eval'")
  }
  
  info.tmp <- as.data.frame(strsplit(model.names, "_"))
  
  return(switch(info,
                species = paste(unique(unlist(info.tmp[-c(nrow(info.tmp), nrow(info.tmp) - 1,  nrow(info.tmp) - 2), ])), collapse = "_"), 
                data.set = paste(unique(unlist(info.tmp[(nrow(info.tmp) - 2), ]))), 
                run.eval = paste(unique(unlist(info.tmp[(nrow(info.tmp) - 1), ]))), 
                models = paste(unique(unlist(info.tmp[(nrow(info.tmp)), ])))
  ))
}

