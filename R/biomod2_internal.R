
.getOS = function()
{
  sysinf = Sys.info()
  if (!is.null(sysinf))
  {
    os = sysinf['sysname']
    if (os == 'Darwin') os = "mac"
  } else
  {
    os = .Platform$OS.type
    if (grepl("^darwin", R.version$os)) os = "mac"
    if (grepl("linux-gnu", R.version$os)) os = "linux"
  }
  return(tolower(os))
}


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


## TEST PARAMETERS --------------------------------------------------------------------------------

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
  if (any(! objValue %in% values)) {
    stop(paste0("\n", objName, " must be '", 
                ifelse(length(values) > 1, 
                       paste0(paste0(values[1:(length(values) -1)], collapse = "', '"),
                              "' or '", values[length(values)])
                       , paste0(values,"'"))))
    test <- FALSE
  }
  return(test)
}

.fun_testIfPosNum <- function(test, objName, objValue)
{
  if (!is.numeric(objValue)) {
    stop(paste0("\n", objName, " must be a numeric"))
    test <- FALSE
  } else if (objValue < 0) {
    stop(paste0("\n", objName, " must be a positive numeric"))
    test <- FALSE
  }
  return(test)
}

.fun_testIf01 <- function(test, objName, objValue)
{
  test <- .fun_testIfPosNum(test, objName, objValue)
  if (test && objValue > 1) {
    stop(paste0("\n", objName, " must be a 0 to 1 numeric"))
    test <- FALSE
  }
  return(test)
}

.fun_testIfPosInt <- function(test, objName, objValue)
{
  if (!is.numeric(objValue)) {
    cat(paste0("\n", objName, " must be a integer"))
    test <- FALSE
  } else if (any(objValue < 0) || any(objValue %% 1 != 0)) {
    cat(paste0("\n", objName, " must be a positive integer"))
    test <- FALSE
  }
  return(test)
}


.fun_testMetric <- function(test, objName, objValue, values){
  if (!is.null(objValue)) {
    test <- .fun_testIfInherits(test, objName, objValue, "character")
    if(length(objValue) != 1){
      stop(paste0("`",objName,"` must only have one element"))
    }
    test <- .fun_testIfIn(test, objName, objValue, values)
  }
  return(test)
}

.check_calib.lines_names <- function(calib.lines, expected_PA.names){
  full.names <- colnames(calib.lines)
  if (missing(expected_PA.names)) { # CV only
    expected_CV.names <- c(paste0("_allData_RUN", seq_len(ncol(calib.lines))), "_allData_allRun")
    .fun_testIfIn(TRUE, "colnames(calib.lines)", full.names, expected_CV.names)
  } else {
    err.msg <- "colnames(calib.lines) must follow the following format: '_PAx_RUNy' with x and y integer"
    # check for beginning '_'
    if (!all( substr(full.names, 1, 1) == "_")) {
      stop(err.msg)
    }
    PA.names <- sapply(strsplit(full.names, split = "_"), function(x) x[2])
    CV.names <- sapply(strsplit(full.names, split = "_"), function(x) x[3])
    .fun_testIfIn(TRUE, "Pseudo-absence dataset in colnames(calib.lines)", PA.names, expected_PA.names)
    if (!all( substr(CV.names, 1, 3) == "RUN")) {
      stop(err.msg)
    }
    CV.num <- sapply(strsplit(CV.names, split = "RUN"), function(x) x[2])
    if (any(is.na(as.numeric(CV.num)))) {
      stop(err.msg)
    }
  }
}

# Functions to get variables ranges ---------------------------------
# used bm_RunModelsLoop, BIOMOD_EnsembleModeling
get_var_type <- function(data) { return(sapply(data, function(x) class(x)[1])) }

get_var_range <- function(data)
{
  get_range <- function(x) {
    if (is.numeric(x)) {
      return(c(min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE)))
    } else if (is.factor(x)) {
      return(levels(x))
    }
  }
  xx <- lapply(data, get_range)
  names(xx) <- names(data)
  return(xx)
}

#' @importFrom terra rast classify subset
.run_pred <- function(object, Prev = 0.5, dat, mod.name = NULL)
{
  if (is.finite(object$deviance) && 
      is.finite(object$null.deviance) && 
      object$deviance != object$null.deviance)
  {
    if (inherits(dat, 'SpatRaster')) {
      pred <- predict(object = dat, model = object, type = "response", wopt = list(names = mod.name))
    } else {
      pred <- predict(object, dat, type = "response")
    }
  }
  
  if (!exists('pred')) {
    if (inherits(dat, 'SpatRaster')) {
      pred <- subset(dat, 1)
      if (Prev < 0.5) {
        pred <- classify(x = pred, rcl = c(-Inf, Inf, 0))
      } else {
        pred <- classify(x = pred, rcl = c(-Inf, Inf, 1))
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

## get formal prediction for EM models ------------------------------------

.get_formal_predictions <- function(object, newdata, on_0_1000, seedval)
{
  # make prediction of all models required
  formal_predictions <- sapply(object@model,
                               function(mod.name, dir_name, resp_name, modeling.id)
                               {
                                 ## check if model is loaded on memory
                                 if (is.character(mod.name)) {
                                   mod <- get(load(file.path(dir_name, resp_name, "models", modeling.id, mod.name)))
                                 }
                                 temp_workdir = NULL
                                 if (length(grep("MAXENT$", mod.name)) == 1) {
                                   temp_workdir = mod@model_output_dir
                                 }
                                 return(predict(mod, newdata = newdata, on_0_1000 = on_0_1000
                                                , temp_workdir = temp_workdir, seedval = seedval))
                               }, dir_name = object@dir_name, resp_name = object@resp_name, modeling.id = object@modeling.id)
  
  
  return(formal_predictions)
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


## RESCALE MODEL WITH BINOMIAL GLM ----------------------------------------------------------------
## used in biomod2_classes_4.R, bm_RunModelsLoop files

.scaling_model <- function(dataToRescale, ref = NULL, ...)
{
  args <- list(...)
  prevalence <- args$prevalence
  weights <- args$weights
  
  ## if no weights given, some are created to rise the define prevalence ------
  if (is.null(weights) && ! is.null(prevalence)) {
    nbPres <- sum(ref, na.rm = TRUE)
    nbAbs <- length(ref) - nbPres
    weights <- rep(1, length(ref))
    
    if (nbAbs > nbPres) { # code absences as 1
      weights[which(ref > 0)] <- (prevalence * nbAbs) / (nbPres * (1 - prevalence))
    } else { # code presences as 1
      weights[which(ref == 0 | is.na(ref))] <- (nbPres * (1 - prevalence)) / (prevalence * nbAbs)
    }
    weights = round(weights[]) # to remove glm & gam warnings
    
  } else if (is.null(weights)) { ## only 1 weights vector ---------------------
    weights <- rep(1, length(ref))
  }
  
  ## define a glm to scale predictions from 0 to 1 ----------------------------
  scaling_model <- glm(ref ~ pred, data = data.frame(ref = as.numeric(ref), pred = as.numeric(dataToRescale))
                       , family = binomial(link = probit), x = TRUE, weights = weights)
  
  return(scaling_model)
}


## EXTRACT models outputs according to specific infos ---------------------------------------------
## used in biomod2_classes_3.R

.filter_outputs.df <- function(out, subset.list)
{
  keep_subset <- foreach(sub.i = names(subset.list)) %do%
    {
      keep_lines <- 1:nrow(out)
      if (!is.null(subset.list[[sub.i]])) {
        .fun_testIfIn(TRUE, sub.i, subset.list[[sub.i]], unique(out[, sub.i]))
        keep_lines <- which(out[, sub.i] %in% subset.list[[sub.i]])
      }
      return(keep_lines)
    }
  keep_lines <- Reduce(intersect, keep_subset)
  if (length(keep_lines) == 0) {
    warning(paste0("No information corresponding to the given filter(s) ("
                   , paste0(names(subset.list), collapse = ', '), ")"))
  }
  return(keep_lines)
}

.filter_outputs.vec <- function(out_names, obj.type, subset.list)
{
  keep_layers <- out_names
  if ("full.name" %in% names(subset.list) && !is.null(subset.list[["full.name"]])) {
    .fun_testIfIn(TRUE, "full.name", subset.list[["full.name"]], out_names)
    keep_layers <- which(out_names %in% subset.list[["full.name"]])
  } else {
    keep_subset <- foreach(sub.i = names(subset.list)) %do%
      {
        keep_tmp <- 1:length(out_names)
        if (!is.null(subset.list[[sub.i]])) {
          .fun_testIfIn(TRUE, sub.i, subset.list[[sub.i]], .extract_modelNamesInfo(out_names, obj.type = obj.type, info = sub.i, as.unique = TRUE))
          keep_tmp <- grep(paste(subset.list[[sub.i]], collapse = "|"), out_names)
        }
        return(keep_tmp)
      }
    keep_layers <- Reduce(intersect, keep_subset)
  }
  if (length(keep_layers) == 0) {
    warning(paste0("No information corresponding to the given filter(s) ("
                   , paste0(names(subset.list), collapse = ', '), ")"))
  }
  return(keep_layers)
}

# Format projection output as data.frame with appropriate columns --------------

.format_proj.df <- function(proj, obj.type = "mod"){
  # argument check
  .fun_testIfInherits(TRUE, "proj", proj, "data.frame")
  .fun_testIfIn(TRUE, "obj.type", obj.type, c("mod","em"))
  # function content
  proj$points <- 1:nrow(proj)
  tmp <- melt(proj, id.vars =  "points")
  colnames(tmp) <- c("points", "full.name", "pred")
  tmp$full.name <- as.character(tmp$full.name)
  if(obj.type == "mod"){
    tmp$PA <- .extract_modelNamesInfo(tmp$full.name, obj.type = "mod", info = "PA", as.unique = FALSE)
    tmp$run <- .extract_modelNamesInfo(tmp$full.name, obj.type = "mod", info = "run", as.unique = FALSE)
    tmp$algo <- .extract_modelNamesInfo(tmp$full.name, obj.type = "mod", info = "algo", as.unique = FALSE)
    proj <- tmp[, c("full.name", "PA", "run", "algo", "points", "pred")]
  } else {
    tmp$merged.by.PA <- .extract_modelNamesInfo(tmp$full.name, obj.type = "em", info = "merged.by.PA", as.unique = FALSE)
    tmp$merged.by.run <- .extract_modelNamesInfo(tmp$full.name, obj.type = "em", info = "merged.by.run", as.unique = FALSE)
    tmp$merged.by.algo <- .extract_modelNamesInfo(tmp$full.name, obj.type = "em", info = "merged.by.algo", as.unique = FALSE)
    tmp$filtered.by <- .extract_modelNamesInfo(tmp$full.name, obj.type = "em", info = "filtered.by", as.unique = FALSE)
    tmp$algo <- .extract_modelNamesInfo(tmp$full.name, obj.type = "em", info = "algo", as.unique = FALSE)
    proj <- tmp[, c("full.name", "merged.by.PA", "merged.by.run", "merged.by.algo"
                    , "filtered.by", "algo", "points", "pred")]
    
  }
  proj
}
## TRANSFORM models outputs getting specific slot --------------------
## used in BIOMOD_Modeling.R, BIOMOD_EnsembleModeling.R

.transform_outputs_list = function(obj.type, obj.out, out = 'evaluation')
{
  ## 0. CHECK object type ---------------------------------------------------------------
  .fun_testIfIn(TRUE, "obj.type", obj.type, c("mod", "em"))
  .fun_testIfIn(TRUE, "out", out, c("model", "calib.failure", "models.kept", "pred", "pred.eval", "evaluation", "var.import"))
  
  if (obj.type == "mod") {
    dim_names <- c("PA", "run", "algo")
  } else if (obj.type == "em") {
    dim_names <- c("merged.by", "filtered.by", "algo")
    dim_names.bis <- c("merged.by.PA", "merged.by.run", "merged.by.algo")
  }
  
  if (obj.type == "mod") {
    output <- foreach(i.dim1 = 1:length(obj.out), .combine = "rbind") %do%
      {
        res <- obj.out[[i.dim1]][[out]]
        if (!is.null(res) && length(res) > 0) {
          res <- as.data.frame(res)
          if (out %in% c("model", "calib.failure", "models.kept", "pred", "pred.eval")) {
            colnames(res) <- out
            res[["points"]] <- 1:nrow(res)
            res <- res[, c("points", out)]
          }
          col_names <- colnames(res)
          tmp.full.name <- obj.out[[i.dim1]][["model"]]
          if(out == "calib.failure" | is.null(tmp.full.name)){
            res[["full.name"]] <- NA
            return(res[, c("full.name", col_names)])
          } else {
            res[["full.name"]] <- tmp.full.name
            res[[dim_names[1]]] <- strsplit(tmp.full.name, "_")[[1]][2]
            res[[dim_names[2]]] <- strsplit(tmp.full.name, "_")[[1]][3]
            res[[dim_names[3]]] <- strsplit(tmp.full.name, "_")[[1]][4]
            return(res[, c("full.name", dim_names, col_names)])
          }
        }
      }
  } else if (obj.type == "em") {
    
    ## 1. GET dimension names -------------------------------------------------------------
    ## BIOMOD.models.out            -> dataset / run / algo
    ## BIOMOD.ensemble.models.out   -> mergedBy / filteredBy / algo
    dim1 <- length(obj.out)
    dim2 <- length(obj.out[[1]])
    dim3 <- length(obj.out[[1]][[1]])
    
    ## 2. GET outputs ---------------------------------------------------------------------
    output <- foreach(i.dim1 = 1:dim1, .combine = "rbind") %:%
      foreach(i.dim2 = 1:dim2, .combine = "rbind") %:% 
      foreach(i.dim3 = 1:dim3, .combine = "rbind") %do%
      {
        if (length(obj.out) >= i.dim1 &&
            length(obj.out[[i.dim1]]) >= i.dim2 &&
            length(obj.out[[i.dim1]][[i.dim2]]) >= i.dim3) {
          res <- obj.out[[i.dim1]][[i.dim2]][[i.dim3]][[out]]
          if (!is.null(res) && length(res) > 0) {
            res <- as.data.frame(res)
            if (out %in% c("model", "calib.failure", "models.kept", "pred", "pred.eval")) {
              colnames(res) <- out
              res[["points"]] <- 1:nrow(res)
              res <- res[, c("points", out)]
            }
            col_names <- colnames(res)
            res[[dim_names[1]]] <- names(obj.out)[i.dim1]
            res[[dim_names[2]]] <- names(obj.out[[i.dim1]])[i.dim2]
            res[[dim_names[3]]] <- names(obj.out[[i.dim1]][[i.dim2]])[i.dim3]
            tmp.full.name <- obj.out[[i.dim1]][[i.dim2]][[i.dim3]][["model"]]
            if(out == "calib.failure" | is.null(tmp.full.name)){
              res[["full.name"]] <- NA
            } else {
              res[["full.name"]] <- tmp.full.name
            }
            tmp <- names(obj.out)[i.dim1]
            res[[dim_names.bis[1]]] <- strsplit(tmp, "_")[[1]][1]
            res[[dim_names.bis[2]]] <- strsplit(tmp, "_")[[1]][2]
            res[[dim_names.bis[3]]] <- strsplit(tmp, "_")[[1]][3]
            res[[dim_names[1]]] <- NULL
            return(res[, c("full.name", dim_names.bis, dim_names[-1], col_names)])
          }
        }
      }
  }
  
  if (out %in% c("model", "calib.failure", "models.kept")) {
    if (is.null(output)) { 
      output <- 'none'
    } else { 
      output <- unique(as.character(output[[out]]))
    }
  }
  return(output)
}

## transform predictions data.frame from long to wide with models as columns ----------------
## used in BIOMOD_RangeSize.R
##' @name .transform_model.as.col
##' @author Remi Patin
##' 
##' @title Transform predictions data.frame from long to wide with models as columns
##' 
##' @description This function is used internally in \code{\link{get_predictions}}
##' to ensure back-compatibility with former output of \code{\link{get_predictions}}
##' (i.e. for \pkg{biomod2} version < 4.2-2). It transform a long \code{data.frame}
##' into a wide \code{data.frame} with each column corresponding to a single model.
##' 
##' \emph{Note that the function is intended for internal use but have been 
##' made available for compatibility with \pkg{ecospat}} 
##'
##' @param df a long \code{data.frame}, generally a standard output of 
##' \code{\link{get_predictions}}
##' 
##' @return a wide \code{data.frame}
##' @importFrom stats reshape
##' @export
##' @keywords internal
.transform_model.as.col <- function(df){
  df <- reshape(df[,c("full.name","points","pred")], idvar = "points", timevar = "full.name", direction = "wide")[,-1, drop = FALSE] 
  colnames(df) <- substring(colnames(df), 6)
  df
}

## EXTRACT model names according to specific infos ------------------------------------------------
## used in biomod2_classes_3.R, BIOMOD_LoadModels, BIOMOD_Projection, BIOMOD_EnsembleModeling, bm_PlotRangeSize

.extract_modelNamesInfo <- function(model.names, obj.type, info, as.unique = TRUE)
{
  ## 0. CHECK parameters ----------------------------------------------------------------
  if (!is.character(model.names)) { stop("model.names must be a character vector") }
  .fun_testIfIn(TRUE, "obj.type", obj.type, c("mod", "em"))
  if (obj.type == "mod") {
    .fun_testIfIn(TRUE, "info", info, c("species", "PA", "run", "algo"))
  } else if (obj.type == "em") {
    .fun_testIfIn(TRUE, "info", info, c("species", "merged.by.PA", "merged.by.run", "merged.by.algo", "filtered.by", "algo"))
  }
  
  ## 1. SPLIT model.names ---------------------------------------------------------------
  info.tmp <- strsplit(model.names, "_")
  
  if (obj.type == "mod") {
    res <- switch(info
                  , species = sapply(info.tmp, function(x) x[1])
                  , PA = sapply(info.tmp, function(x) x[2])
                  , run = sapply(info.tmp, function(x) x[3])
                  , algo = sapply(info.tmp, function(x) x[4]))
  } else if (obj.type == "em") {
    res <- switch(info
                  , species = sapply(info.tmp, function(x) x[1])
                  , merged.by.PA = sapply(info.tmp, function(x) x[3])
                  , merged.by.run = sapply(info.tmp, function(x) x[4])
                  , merged.by.algo = sapply(info.tmp, function(x) x[5])
                  , filtered.by = sapply(info.tmp, function(x) sub(".*By", "", x[2]))
                  , algo = sapply(info.tmp, function(x) sub("By.*", "", x[2])))
  }
  
  ## 2. RETURN either the full vector, or only the possible values ----------------------
  if (as.unique == TRUE) {
    return(unique(res))
  } else {
    return(res)
  }
}


## Extract info (out/binary/filtered + metric) from BIOMOD.projection.out  -------------------------------------

.extract_projlinkInfo <- function(obj){
  stopifnot(inherits(obj, 'BIOMOD.projection.out'))
  out.names <- obj@proj.out@link
  sp.name <- obj@sp.name
  sub(pattern     = paste0(".*?_",sp.name), 
      replacement = "",
      x           = out.names)
  out.binary <- grepl(pattern = "bin",
                      x = out.names)
  out.filt <- grepl(pattern = "filt",
                    x = out.names)
  out.type <- ifelse(out.binary,
                     "bin",
                     ifelse(out.filt,
                            "filt", 
                            "out"))
  out.metric <- sapply(seq_len(length(out.names)), 
                       function(i){
                         if (out.type[i] == "out") {
                           return("none")
                         } else {
                           begin.cut <-  gsub(".*_", "", out.names[i]) 
                           return(  sub(paste0(out.type[i],".*"), "", begin.cut) )
                         }
                       })
  data.frame("link"   = out.names,
             "type"   = out.type,
             "metric" = out.metric)
}

# Extract Selected Layers --------------------------------------------
# return relevant layers indices from a BIOMOD.projection.out given metric.binary
# or metric.select
.extract_selected.layers <- function(obj, metric.binary, metric.filter){
  
  ## argument  check ----------------------------------------------------------
  stopifnot(inherits(obj, "BIOMOD.projection.out"))
  
  if (!is.null(metric.binary) & !is.null(metric.filter)) {
    stop("cannot return both binary and filtered projection, please provide either `metric.binary` or `metric.filter` but not both.")
  }
  
  df.info <- .extract_projlinkInfo(obj)
  
  available.metrics.binary <- unique(df.info$metric[which(df.info$type == "bin")])
  available.metrics.filter <- unique(df.info$metric[which(df.info$type == "filt")])
  .fun_testMetric(TRUE, "metric.binary", metric.binary, available.metrics.binary)
  .fun_testMetric(TRUE, "metric.filter", metric.filter, available.metrics.filter)
  
  
  ## select layers  -----------------------------------------------------------
  if (!is.null(metric.binary)) {
    selected.layers <- which(df.info$type == "bin" & df.info$metric == metric.binary)  
  } else if ( !is.null(metric.filter)) {
    selected.layers <- which(df.info$type == "filt" & df.info$metric == metric.filter)  
  } else {
    selected.layers <- which(df.info$type == "out")  
  }
  selected.layers
}

## FILL BIOMOD.models.out elements in bm.mod or bm.em objects -------------------------------------
## used in BIOMOD_Modeling, BIOMOD_EnsembleModeling
.fill_BIOMOD.models.out <- function(objName, objValue, mod.out, inMemory = FALSE, nameFolder = ".")
{
  
  # backup case uses existing @val slot
  if(is.null(objValue)){
    eval(parse(text = paste0("objValue <- mod.out@", objName, "@val")))
  }
  
  save(objValue, file = file.path(nameFolder, objName), compress = TRUE)
  if (inMemory) {
    eval(parse(text = paste0("mod.out@", objName, "@val <- objValue")))
  }
  eval(parse(text = paste0("mod.out@", objName, "@inMemory <- ", inMemory)))
  eval(parse(text = paste0("mod.out@", objName, "@link <- file.path(nameFolder, objName)")))
  return(mod.out)
}


## FUNCTIONS for 'plot' and 'show' objects --------------------------------------------------------
## used in biomod2_classes_1.R file

.clever_cut <- function(x)
{
  nb_col = ceiling(sqrt(x))
  nb_row = ceiling(x / nb_col)
  return(c(nb_row, nb_col))
}

.print_formula <- function(formula = NULL)
{
  ifelse(length(formula) < 1, 'NULL', paste(formula[2], formula[1], formula[3]))
}

.print_control <- function(ctrl)
{
  out <-  paste0(names(ctrl)[1], " = ", ctrl[[1]])
  if (length(ctrl) > 1)
  {
    i = 2
    while (i <= length(ctrl)) {
      if (is.list(ctrl[[i]])) {
        out <- c(out, paste0(", ", names(ctrl)[i], " = list(",
                             paste0(names(ctrl[[i]]), "=", unlist(ctrl[[i]]), collapse = ", "),
                             ")"))
        i <- i + 1
      } else {
        out <- c(out, paste0(", ", names(ctrl)[i], " = ", ctrl[[i]]))
        i <- i + 1
      }
    }
  }
  return(out)
}


## GAM library loading --------------------------------------------------------
## used in biomod2_classes_4.R file for GAM predict template
##' @name .load_gam_namespace
##' 
##' @title Load library for GAM models
##' 
##' @description This function loads library for either GAM and BAM from mgcv
##'   package or for GAM from gam package.
##' 
##' @param model_subclass the subclass of GAM model
##' @keywords internal

.load_gam_namespace <- function(model_subclass = c("GAM_mgcv","BAM_mgcv", "GAM_gam")){
  if (model_subclass %in% c("GAM_mgcv", "BAM_mgcv")) {
    # cat("\n*** unloading gam package / loading mgcv package")
    if (isNamespaceLoaded("gam")) { unloadNamespace("gam") }
    if (!isNamespaceLoaded("mgcv")) { 
      if(!requireNamespace('mgcv', quietly = TRUE)) stop("Package 'mgcv' not found")
    }
  }
  
  if (model_subclass == "GAM_gam") {
    # cat("\n*** unloading mgcv package / loading gam package")
    if (isNamespaceLoaded("mgcv")) {
      if (isNamespaceLoaded("caret")) { unloadNamespace("caret")} ## need to unload caret before car
      if (isNamespaceLoaded("car")) { unloadNamespace("car") } ## need to unload car before mgcv
      unloadNamespace("mgcv")
    }
    if (!isNamespaceLoaded("gam")) { 
      if(!requireNamespace('gam', quietly = TRUE)) stop("Package 'gam' not found")
    }
  }
  invisible(NULL)
}

## Get categorical variable names ---------------------------------------
##' @name .get_categorical_names
##' 
##' @title Get categorical variable names
##' 
##' @description Internal function to get categorical variables name from a data.frame.
##' 
##' @param df data.frame to be checked
##' @return a vector with the name of categorical variables
##' @keywords internal

.get_categorical_names <- function(df){
  unlist(sapply(names(df), function(x) {
    if (is.factor(df[, x])) { return(x) } else { return(NULL) }
  }))
}

## Categorical to numeric ---------------------------------------
##' @name .categorical2numeric
##' 
##' @title Transform categorical into numeric variables
##' 
##' @description Internal function transform categorical variables in a 
##' data.frame into numeric variables. Mostly used with maxent which cannot 
##' read character
##' 
##' @param df data.frame to be transformed
##' @param categorical_var the names of categorical variables in df
##' @return a data.frame without categorical variables
##' @keywords internal

.categorical2numeric <- function(df, categorical_var){
  if(length(categorical_var) > 0){
    for(this_var in categorical_var){
      df[,this_var] <- as.numeric(df[,this_var])
    }
  }
  df
}


## .check_bytes_format for MAXENT options----------------------------
##' @name .check_bytes_format
##' 
##' @title Check bytes formatting 
##' 
##' @description Internal function that check a character string to match a byte
##' format for Java. e.g. 1024M, 1024m, 1024k or 1024K
##' 
##' @param x string to be transformed
##' @return a boolean
##' @keywords internal

.check_bytes_format <- function(this_test, x, varname){
  possible_suffix <- c("k","K","m","M","g","G")
  this_suffix <- substr(x, nchar(x), nchar(x))
  this_number <- substr(x, 1, nchar(x)-1)
  if (! this_suffix %in% possible_suffix) {
    this_test <- FALSE
    cat(paste0("\nMAXENT$",varname," last letter must be among ",
               paste0(possible_suffix[-length(possible_suffix)], collapse = ", "),
               "", " and ", possible_suffix[length(possible_suffix)],". Current value = '", this_suffix,"' was not usable."))
  }
  if (suppressWarnings(is.na(as.numeric(this_number)))) {
    this_test <- FALSE
    cat(paste0("\nMAXENT$",varname," must be a number plus a single letter. Begginning of the given argument was '", this_number, "' and was not convertible into a number."))
  }
  return(this_test)
}

## Tools for SpatRaster ----------------------------
##' @name rast.has.values
##' 
##' @title Check whether SpatRaster is an empty \code{rast()} 
##' 
##' @description Check whether SpatRaster is an empty \code{rast()} 
##' 
##' @param x SpatRaster to be checked
##' @return a boolean
##' @keywords internal

rast.has.values <- function(x){
  stopifnot(inherits(x, 'SpatRaster'))
  suppressWarnings(any(!is.na(values(x))))
}

##' @name check_duplicated_cells
##' 
##' @title Check duplicated cells
##' 
##' @description Identify data that are contained in the same raster cells ; 
##' print warnings if need be and filter those data if asked for.
##' 
##' 
##' @param env a \code{SpatRaster} with environmental data
##' @param xy a \code{data.frame} with coordinates
##' @param sp a \code{vector} with species occurrence data
##' @param filter.raster a \code{boolean} 
##' @return a \code{list}
##' @keywords internal
##' 
##' @importFrom terra cellFromXY

check_duplicated_cells <- function(env, xy, sp, filter.raster, 
                                   PA.user.table = NULL){
  sp.cell <- duplicated(cellFromXY(env, xy))
  if (any(sp.cell)) {
    if (filter.raster) {
      sp <- sp[!sp.cell]
      xy <- xy[!sp.cell,]
      if (!is.null(PA.user.table)) {
        PA.user.table <- PA.user.table[!sp.cell, , drop = FALSE]
      }
      cat("\n !!! Some data are located in the same raster cell. 
          Only the first data in each cell will be kept as `filter.raster = TRUE`.")
    } else {
      cat("\n !!! Some data are located in the same raster cell. 
          Please set `filter.raster = TRUE` if you want an automatic filtering.")
    }
  }
  return(list("sp"  = sp, 
              "xy"  = xy,
              "PA.user.table" = PA.user.table))
  
}

## Get new.env class ----------------------------
##' @name .get_env_class
##' 
##' @title Get class of environmental data provided
##' 
##' @description Get class of environmental data provided
##' 
##' @param new.env object to identify
##' @return a character
##' @keywords internal

.get_env_class <- function(new.env){
  .fun_testIfInherits(TRUE, "new.env", new.env, c('data.frame', 'SpatRaster'))
  if (inherits(new.env,"data.frame")) {
    return("data.frame")
  }
  if (inherits(new.env,"SpatRaster")) {
    return("SpatRaster")
  }
  NULL
}


# Common tools ------------------------------------------------------------

.check_formating_spatial <- function(resp.var, expl.var = NULL, resp.xy = NULL, eval.data = FALSE){
  if (!is.null(resp.xy)) {
    cat("\n      ! XY coordinates of response variable will be ignored because spatial response object is given.")
  }
  
  if (inherits(resp.var, 'SpatialPoints')) { 
    resp.xy <- data.matrix(sp::coordinates(resp.var))
    if (inherits(resp.var, 'SpatialPointsDataFrame')) {
      resp.var <- resp.var@data
    } else {
      cat("\n      ! Response variable is considered as only presences... Is it really what you want?")
      resp.var <- rep(1, nrow(resp.xy))
    }
  }
  
  if (inherits(resp.var, 'SpatVector')) { 
    resp.xy <- data.matrix(crds(resp.var))
    resp.var <- as.data.frame(resp.var)
    if (ncol(resp.var) == 0) {
      if(eval.data){
        stop("eval.resp must have both presences and absences in the data associated to the SpatVector") 
      } else {
        cat("\n      ! Response variable is considered as only presences... Is it really what you want?")
        resp.var <- rep(1, nrow(resp.xy))
      }
    }
  }
  
  if(!eval.data){
    if ( all(!is.na(resp.var)) && 
         all(resp.var == 1, na.rm = TRUE) &&
         !inherits(expl.var, c('Raster','SpatRaster'))) {
      stop("For Presence-Only model based on SpatialPoints or SpatVector, expl.var needs to be a RasterStack or SpatRaster to be able to sample pseudo-absences")
    }
  }
  
  return(
    list(resp.var = resp.var,
         resp.xy = resp.xy)
  )
}

.check_formating_resp.var <- function(resp.var, eval.data = FALSE){
  if (length(which(!(resp.var %in% c(0, 1, NA)))) > 0) {
    cat("\n      ! ", ifelse(eval.data, "Evaluation",""), "Response variable have non-binary values that will be converted into 0 (resp <=0) or 1 (resp > 0).")
    resp.var[which(resp.var > 0)] <- 1
    resp.var[which(resp.var <= 0)] <- 0
  }
  
  if (eval.data) {
    if (!any(resp.var == 1, na.rm = TRUE) || !any(resp.var == 0, na.rm = TRUE))
    {
      stop("Evaluation response data must have both presences and absences")
    }
  }
  
  resp.var
}

.check_formating_table <- function(resp.var){
  resp.var = as.data.frame(resp.var)
  if (ncol(resp.var) > 1) {
    stop("You must give a monospecific response variable (1D object)")
  } else {
    resp.var <- as.numeric(resp.var[, 1])
  }
  resp.var
}

.check_formating_xy <- function(resp.xy, resp.length){
  if (ncol(resp.xy) != 2) {
    stop("If given, resp.xy must be a 2 column matrix or data.frame")
  }
  if (nrow(resp.xy) != resp.length) {
    stop("Response variable and its coordinates don't match")
  }
  as.data.frame(resp.xy)
}

.check_formating_expl.var <- function(expl.var, length.resp.var){
  if (is.matrix(expl.var) | is.numeric(expl.var)) {
    expl.var <- as.data.frame(expl.var)
  }
  
  if (inherits(expl.var, 'Raster')) {
    expl.var <- raster::stack(expl.var, RAT = FALSE)
    if (any(is.factor(expl.var))) {
      expl.var <- .categorical_stack_to_terra(expl.var)
    } else {
      # as of 20/10/2022 the line below does not work if categorical variables
      # are present, hence the trick above. 
      expl.var <- rast(expl.var)
    }
  }
  
  if (inherits(expl.var, 'SpatialPoints')) {
    expl.var <- as.data.frame(expl.var@data)
  }
  if (inherits(expl.var, 'SpatVector')) {
    expl.var <- as.data.frame(expl.var)
  }
  
  if (inherits(expl.var, 'data.frame')) {
    if (nrow(expl.var) != length.resp.var) {
      stop("If explanatory variable is not a raster then dimensions of response variable and explanatory variable must match!")
    }
  }
  expl.var
}


# Categorical Variables Management ----------------------------------------

.categorical_stack_to_terra <- function(myraster, expected_levels = NULL) {
  myTerra <- rast(
    sapply(1:raster::nlayers(myraster), 
           function(thislayer){
             rast(myraster[[thislayer]])
           })
  )
  which.factor <- which(raster::is.factor(myraster))
  for (this.layer in which.factor) {
    this.levels <- raster::levels(myraster)[[this.layer]][[1]]
    ind.levels <- ifelse(ncol(this.levels) > 1, 2, 1)
    if (any(duplicated(this.levels[ , ind.levels]))) {
      stop("duplicated levels in environmental raster")
    }
    if (is.null(expected_levels)) {
      # no check to do, just formatting
      this.levels.df <- data.frame(ID = this.levels[,ind.levels],
                                   value = paste0(this.levels[,ind.levels]))
    } else {
      new.levels <- paste0(this.levels[,ind.levels])
      new.index <- this.levels[,1]
      this.layer.name <- names(myraster)[this.layer]
      fit.levels <- levels(expected_levels[,this.layer.name])
      if (!all(new.levels %in% fit.levels)) {
        cat("\n",
            "!! Levels for layer", colnames(expected_levels)[this.layer],
            " do not match.", new.levels[which(!new.levels %in% fit.levels)], "not found in fit data",
            "\n Fit data levels: ",paste0(fit.levels, collapse = " ; "),
            "\n Projection data levels: ",paste0(new.levels, collapse = " ; "))
        stop(paste0("Levels for ", colnames(expected_levels)[this.layer],
                    " do not match."))
      }
      levels.to.add <- which(!fit.levels %in% new.levels)
      if (length(levels.to.add) > 0) {
        max.index <- max(new.index)
        this.levels.df <- data.frame(ID = c(new.index,
                                            max.index +
                                              seq_along(levels.to.add)),
                                     value = c(new.levels, fit.levels[levels.to.add]))
      } else {
        this.levels.df <- data.frame(ID = new.index,
                                     value = new.levels)
      }
      
    }
    colnames(this.levels.df) <- c("ID", names(myTerra)[this.layer])
    try({myTerra <- categories(myTerra, layer = this.layer, this.levels.df)}, silent = TRUE)
  }
  return(myTerra)
}

.check_env_levels <-  function(new.env, expected_levels) {
  which.factor <- which(sapply(new.env, is.factor))
  if (inherits(new.env, 'SpatRaster')) {
    for (this.layer in which.factor) {
      this.layer.name <- names(new.env)[this.layer]
      this.levels <- cats(new.env)[[this.layer]]
      ind.levels <- ifelse(ncol(this.levels) > 1, 2, 1)
      new.levels <- paste0(this.levels[,ind.levels])
      if( any(duplicated(new.levels)) ) {
        stop("duplicated levels in `new.env`")
      }
      new.index <- this.levels[,1]
      fit.levels <- levels(expected_levels[,this.layer.name])
      if (!all(new.levels %in% fit.levels)) {
        cat("\n",
            "!! Levels for layer", colnames(expected_levels)[this.layer],
            " do not match.", new.levels[which(!new.levels %in% fit.levels)], "not found in fit data",
            "\n Fit data levels: ",paste0(fit.levels, collapse = " ; "),
            "\n Projection data levels: ",paste0(new.levels, collapse = " ; "))
        stop(paste0("Levels for ", colnames(expected_levels)[this.layer],
                    " do not match."))
      }
      levels.to.add <- which(!fit.levels %in% new.levels)
      if (length(levels.to.add) > 0) {
        max.index <- max(new.index)
        this.levels.df <- data.frame(ID = c(new.index,
                                            max.index +
                                              seq_along(levels.to.add)),
                                     value = c(new.levels, fit.levels[levels.to.add]))
        
        colnames(this.levels.df) <- c("ID", names(new.env)[this.layer])
        new.env <- categories(new.env, layer = this.layer, 
                              this.levels.df)
      } 
    }
  } else {
    for (this.layer in which.factor) {
      this.layer.name <- names(new.env)[this.layer]
      this.levels <- levels(new.env[,this.layer.name])
      new.levels <- paste0(this.levels)
      fit.levels <- levels(expected_levels[,this.layer.name])
      if (!all(new.levels %in% fit.levels)) {
        cat("\n",
            "!! Levels for layer", colnames(expected_levels)[this.layer],
            " do not match.", new.levels[which(!new.levels %in% fit.levels)], "not found in fit data",
            "\n Fit data levels: ",paste0(fit.levels, collapse = " ; "),
            "\n Projection data levels: ",paste0(new.levels, collapse = " ; "))
        stop(paste0("Levels for ", colnames(expected_levels)[this.layer],
                    " do not match."))
      }
      levels.to.add <- which(!fit.levels %in% new.levels)
      if (length(levels.to.add) > 0) {
        new.env[ , this.layer.name] <- factor(new.env[ , this.layer.name], 
                                              levels = c(new.levels,
                                                         fit.levels[levels.to.add]))
      }
    }
  }
  new.env 
}
