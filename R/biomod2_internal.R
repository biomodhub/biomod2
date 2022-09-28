
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

# Function to check new data range compatibility with calibrating data
# used in biomod2_classes_4.R file, .template_[...] functions
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
  #     new_data <- new_data[,expl_var_names, drop = FALSE]
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

.run_pred <- function(object, Prev = 0.5, dat)
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
                                   if (length(grep("MAXENT.Phillips$", mod.name)) == 1) {
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


## EXTRACT model names according to specific infos ------------------------------------------------
## used in biomod2_classes_3.R, BIOMOD_LoadModels, BIOMOD_Projection, BIOMOD_EnsembleModeling, bm_PlotRangeSize

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


## FILL BIOMOD.models.out elements in bm.mod or bm.em objects -------------------------------------
## used in BIOMOD_Modeling, BIOMOD_EnsembleModeling
.fill_BIOMOD.models.out <- function(objName, objValue, mod.out, inMemory = FALSE, nameFolder = ".")
{
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


## PLOT MULTIPLE and LEVEL plots ------------------------------------------------------------------
## used in biomod2_classes_3.R file

.level.plot <- function(data.in, XY, color.gradient = 'red'
                        , level.range = c(min(data.in, na.rm = TRUE), max(data.in, na.rm = TRUE))
                        , show.scale = TRUE, SRC = FALSE, plot.title = NULL
                        , AddPresAbs = NULL, cex = 1, PresAbsSymbol = c(cex * 0.8, 16, 4), ...)
{
  extra.args <- list(...)
  multiple.plot <- FALSE
  if (!is.null(extra.args$multiple.plot)) {
    multiple.plot <- extra.args$multiple.plot
  }
  
  .fun_testIfIn(TRUE, "color.gradient", color.gradient, c('grey', 'red', 'blue'))
  if (ncol(XY) != 2) { stop("\n wrong coordinates given in 'XY': there should be two columns \n") }
  if (nrow(XY) != length(data.in)) { stop("\n data and coordinates should be of the same length \n") }
  if (SRC == TRUE && length(unique(data.in)) > 4) {
    cat("\n not possible to render SRC plot -> more than four different values in data ")
    SRC <- FALSE
  }
  
  if (SRC == TRUE) {
    color.system <- c("red", "lightgreen", "grey", "darkgreen")
    col_id <- data.in + 3
    plot.title <- ifelse(is.null(plot.title), "SRC plot", plot.title)
  } else
  {
    color.system = switch(color.gradient
                          , 'grey' = gray(seq(0.95, 0, length.out = 102))
                          , 'blue' = c('grey88'
                                       , rainbow(45, start = 0.5, end = 0.65)
                                       , rainbow(10, start = 0.65, end = 0.7)
                                       , rainbow(45, start = 0.7, end = 0.85)
                                       , 'red')
                          , 'red' = c('grey88'
                                      , c(rep(c(colors()[c(417, 417, 515)]), each = 5)
                                          , rev(rainbow(55, start = 0.13, end = 0.23))
                                          , rev(rainbow(50, start = 0.08, end = 0.13)[seq(1, 50, length.out = 15)])
                                          , rev(rainbow(50, end = 0.08)[seq(1, 50, length.out = 15)]))
                                      , 'brown2')
    )
    
    # define a vector to make correspondence between values and colors
    val_seq <- c(seq(level.range[1], level.range[2], length.out = 101), Inf)
    col_id <- sapply(data.in, function(x) { return(which(x <= val_seq)[1]) })
    plot.title <- ifelse(is.null(plot.title), "Level plot", plot.title)
  }
  
  
  ## PLOT points ------------------------------------------------------------------------
  plot(XY[, 2] ~ XY[, 1],
       col = color.system[col_id],
       cex = cex,
       pch = 19,
       xlab = '',
       ylab = '',
       xaxt = 'n',
       yaxt = 'n',
       main = plot.title)
  
  if (!show.scale && !is.null(AddPresAbs)) { 
    ## Add presence/absence if requested by user:
    points(AddPresAbs[AddPresAbs[, 3] == 1, 1:2], col = "black", pch = pchPres, cex = cex2)
    points(AddPresAbs[AddPresAbs[, 3] == 0, 1:2], col = "black", pch = pchAbs, cex = cex2)
    
  } else if (show.scale) {
    if (!multiple.plot) {
      layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1), heights = c(1, 1))
    }
    
    if (!is.null(AddPresAbs)) {
      ## Add presence/absence if requested by user:
      cex2 <- PresAbsSymbol[1]
      pchPres <- PresAbsSymbol[2]
      pchAbs <- PresAbsSymbol[3]
      
      points(AddPresAbs[AddPresAbs[, 3] == 1, 1:2], col = "black", pch = pchPres, cex = cex2)
      points(AddPresAbs[AddPresAbs[, 3] == 0, 1:2], col = "black", pch = pchAbs, cex = cex2)
    }
    
    ## PLOT -----------------------------------------------------------------------------
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(x = c(-1, 1),
         y = c(0, 1),
         xlim = c(0, 1),
         ylim = c(0, 1),
         type = "n",
         axes = FALSE)
    polygon(x = c(-2, -2, 2, 2),
            y = c(-2, 2, 2, -2),
            col = "#f5fcba",
            border = NA)
    
    if (SRC) {
      legend(0, 0.8, legend = list(' (1) new', ' (0) stable', '(-1) kept', '(-2) lost'), cex = 1, fill = rev(color.system), bty = 'n')
    } else {
      lmin = round(level.range[1], digits = 2)
      if (level.range[1] != min(data.in, na.rm = TRUE)) {
        lmin <- paste0(lmin, " or lower")
      }
      lmax = round(level.range[2], digits = 2)
      if (level.range[2] != max(data.in, na.rm = TRUE)) {
        lmax <- paste0(lmax, " or over")
      }
      
      legend.txt = list(lmax, '', '', '', '',
                        round((3 * level.range[2] + level.range[1]) / 4, digits = 2),
                        '', '', '', '',
                        round(sum(level.range) / 2, digits = 2),
                        '', '', '', '',
                        round((level.range[2] + 3 * level.range[1]) / 4, digits = 2),
                        '', '', '', '',
                        lmin)
      yy = ifelse(multiple.plot, 1.05, 0.92)
      cexx = ifelse(multiple.plot, cex, 1)
      legend(0.2, yy, legend = legend.txt, cex = cexx, bty = 'n'
             , fill = rev(color.system[c(1, seq(2, 101, length.out = 19), 102)]))
    }
  }
}

.multiple.plot <- function(Data, coor
                           , color.gradient = 'red', plots.per.window = 9, cex = 1
                           , AddPresAbs = NULL, PresAbsSymbol = c(cex * 0.8, 16, 4))
{
  .fun_testIfIn(TRUE, "color.gradient", color.gradient, c('grey', 'red', 'blue'))
  if (nrow(coor) != nrow(Data)) { stop("\n data and coordinates should be of the same length \n") }
  
  ## Take off NA data
  Data <- Data[, !is.na(Data[1, ])]
  
  ## FUNCTION plotting color boxes
  pbox <- function(co) {
    plot(x = c(-1, 1), y = c(0, 1), xlim = c(0, 1), ylim = c(0, 1), type = "n", axes = FALSE)
    polygon(x = c(-2,-2, 2, 2), y = c(-2, 2, 2,-2), col = co, border = NA)
  }
  
  ## Calculate the number of windows to open
  NbPlots <- ncol(Data)
  NbWindows <- ceiling(NbPlots / plots.per.window)
  if (NbWindows == 1) { plots.per.window <- NbPlots }
  
  for (W in 1:NbWindows)
  {
    dev.new() 
    
    Wstart <- (W - 1) * plots.per.window + 1
    if (W * plots.per.window > NbPlots) { Wfinal <- NbPlots } else { Wfinal <- W * plots.per.window }
    DataW <- as.data.frame(Data[, Wstart:Wfinal])
    colnames(DataW) <- colnames(Data)[Wstart:Wfinal]
    
    #determine the organisation of the plots on the window
    W.width <- ceiling(sqrt(plots.per.window))
    W.height <- ceiling(plots.per.window / W.width)
    
    #create object for scaling the legend
    legendcex <- 0.64 + 1 / exp(W.height)
    
    #matrix of indexes for ordering the layout
    mat <- c(1,2)
    for (i in 1:(W.width - 1)) { mat <- c(mat, mat[1:2] + 4 * i) }
    mat <- rbind(mat, mat + 2)
    for (i in 1:(W.height - 1)) { mat <- rbind(mat, mat[1:2, ] + W.width * 4 * i) }
    
    layout(mat, widths = rep(c(1, 0.3), W.width), heights = rep(c(0.2, 1), W.height))
    
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    for (i in 1:(Wfinal - Wstart + 1)) {
      pbox("grey98")
      text(x = 0.5, y = 0.8, pos = 1, cex = 1.6, labels = colnames(DataW)[i], col = "#4c57eb")
      pbox("grey98")
      .level.plot(DataW[, i], XY = coor, color.gradient = color.gradient, cex = cex, title = ""
                  , AddPresAbs = AddPresAbs, PresAbsSymbol = PresAbsSymbol, multiple.plot = TRUE)
    }
    
    #fill gaps by grey boxes
    if (W.width * W.height - plots.per.window != 0) {
      for (i in 1:((W.width * W.height - plots.per.window) * 4)) { pbox("grey98") }
    }
  } #W loop
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
    if (!isNamespaceLoaded("mgcv")) { requireNamespace("mgcv", quietly = TRUE) }
  }
  
  if (model_subclass == "GAM_gam") {
    # cat("\n*** unloading mgcv package / loading gam package")
    if (isNamespaceLoaded("mgcv")) {
      if (isNamespaceLoaded("caret")) { unloadNamespace("caret")} ## need to unload caret before car
      if (isNamespaceLoaded("car")) { unloadNamespace("car") } ## need to unload car before mgcv
      unloadNamespace("mgcv")
    }
    if (!isNamespaceLoaded("gam")) { requireNamespace("gam", quietly = TRUE) }
  }
  invisible(NULL)
}