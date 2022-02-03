###################################################################################################
##' @name bm_PlotResponseCurves
##' @author Damien Georges, Maya Gueguen
##' 
##' @title Plot response curves
##' 
##' @description
##' 
##' This function represents response curves of species distribution models, from 
##' \code{BIOMOD.models.out} or \code{BIOMOD.ensemble.models.out} objects that can be obtained 
##' from \code{\link{BIOMOD_Modeling}} or \code{\link{BIOMOD_EnsembleModeling}} functions.
##' Response curves can be represented in either 2 or 3 dimensions (meaning 1 or 2 explanatory 
##' variables at a time, see \href{bm_PlotResponseCurves.html#details}{Details}).
##' 
##'   
##' @param modeling.output a \code{\link{BIOMOD.models.out}} or \code{BIOMOD.ensemble.models.out} 
##' object that can be obtained from \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions
##' @param chosen.models a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' @param new.env a \code{matrix}, \code{data.frame} or \code{\link[raster:stack]{RasterStack}} 
##' object containing the new explanatory variables (in columns or layers, with names matching the 
##' variables names given to the \code{\link{BIOMOD_FormatingData}} function to build 
##' \code{modeling.output}) that will be used to project the species distribution model(s)
##' @param show.variables a \code{vector} containing the names of the explanatory variables 
##' present into \code{new.env} parameter and to be plotted
##' @param do.bivariate (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether the response curves are to be represented in 3 
##' dimensions (meaning 2 explanatory variables at a time) or not (meaning only 1)
##' @param fixed.var.metric a \code{character} corresponding to the statistic to be used to fix as 
##' constant the remaining variables other than the one used to predict response, must be either 
##' \code{mean}, \code{median}, \code{min}, \code{max}
##' @param plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the plot is to be rendered or not
##' @param \ldots some additional arguments (see 
##' \href{bm_PlotResponseCurves.html#details}{Details})
##' 
##' 
##' @return  
##' 
##' A \code{ggplot} object representing response curves.
##' 
##' 
##' @details
##' 
##' This function is an adaptation of the Evaluation Strip method proposed by Elith et al. (2005). 
##' To build the predicted response curves :
##' \itemize{
##'   \item \code{n-1} variables are set constant to a fixed value determined by the 
##'   \code{fixed.var.metric} parameter (in the case of categorical variable, the most represented 
##'   class is taken)
##'   \item the remaining variable is made to vary throughout its range given by the \code{new.env} 
##'   parameter 
##'   \item predicted values are computed with these \code{n-1} fixed variables, and this 
##'   studied variable varying 
##' }
##' If \code{do.bivariate = TRUE}, 2 variables are varying at the same time. \cr \cr
##' 
##' The response curves obtained show the sensibility of the model to the studied variable. \cr 
##' Note that this method does not account for interactions between variables. \cr \cr
##' 
##' \code{...} can take the following values :
##' 
##' \itemize{
##'   \item{\code{main}}{ : a \code{character} corresponding to the graphic title}
##' }
##' 
##' 
##' @references
##' 
##' \itemize{
##'   \item Elith, J., Ferrier, S., Huettmann, FALSE. and Leathwick, J. R. 2005. The evaluation 
##'   strip: A new and robust method for plotting predicted responses from species distribution 
##'   models. \emph{Ecological Modelling}, \bold{186}, 280-289.
##' }
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @keywords ggplot, response curve
##' 
##' 
##' @examples
##' 
##' # species occurrences
##' myFile <- system.file("external/species/mammals_table.csv", package="biomod2")
##' DataSpecies <- read.csv(myFile, row.names = 1)
##' head(DataSpecies)
##' 
##' # the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[, c("X_WGS84", "Y_WGS84")]
##' 
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles = paste0("external/bioclim/current/bio", c(3, 4, 7, 11, 12), ".grd")
##' myExpl = raster::stack(system.file(myFiles[1], package = "biomod2"),
##'                        system.file(myFiles[2], package = "biomod2"),
##'                        system.file(myFiles[3], package = "biomod2"),
##'                        system.file(myFiles[4], package = "biomod2"),
##'                        system.file(myFiles[5], package = "biomod2"))
##' 
##' # 1. Formating Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' # 3. Doing Modelisation
##' myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
##'                                     models = c('SRE','RF'),
##'                                     models.options = myBiomodOption,
##'                                     NbRunEval = 2,
##'                                     DataSplit = 80,
##'                                     VarImport = 3,
##'                                     models.eval.meth = c('TSS','ROC'),
##'                                     do.full.models = FALSE,
##'                                     modeling.id = 'test')
##' 
##' # 4.1 Plotting 2D response curves
##' myRespPlot2D <- bm_PlotResponseCurves(myBiomodModelOut,
##'                                       chosen.models = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
##'                                       do.bivariate = FALSE,
##'                                       fixed.var.metric = 'median')
##'                                       
##' # 4.2 Plotting 3D response curves
##' myRespPlot3D <- bm_PlotResponseCurves(myBiomodModelOut,
##'                                       chosen.models = get_built_models(myBiomodModelOut)[3],
##'                                       do.bivariate = TRUE,
##'                                       fixed.var.metric = 'median')
##'                                       
##'                                       
##' 
##' @importFrom foreach foreach %:%
##' @importFrom reshape2 melt
##' @importFrom ggplot2 ggplot aes_string geom_line geom_rug geom_raster facet_wrap xlab ylab
##' theme element_blank element_rect element_text labs scale_fill_viridis_c
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PlotResponseCurves <- function(modeling.output
                                  , chosen.models = 'all'
                                  , new.env = get_formal_data(modeling.output, 'expl.var')
                                  , show.variables = get_formal_data(modeling.output, 'expl.var.names')
                                  , fixed.var.metric = 'mean'
                                  , do.bivariate = FALSE
                                  , plot = TRUE
                                  , ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  formal_names <- chosen.models
  
  args <- .bm_PlotResponseCurves.check.args(modeling.output, chosen.models, new.env, show.variables, do.bivariate, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## 1. Create output object ----------------------------------------------------------------------
  ref_table <- new.env[1, , drop = FALSE]
  rownames(ref_table) <- NULL
  for (i in 1:ncol(new.env)) {
    if (is.numeric(new.env[, i])) {
      ref_table[, i] <- switch(fixed.var.metric,
                               mean = mean(new.env[data_species == 1, i]),
                               median = median(new.env[data_species == 1, i]),
                               min = min(new.env[data_species == 1, i]),
                               max = max(new.env[data_species == 1, i]))
    } else { # return everytime the majoritary class
      sum_level <- summary(new.env[data_species == 1, i], na.rm = TRUE)
      ref_table[, i] <- names(sum_level)[which.max(sum_level)]
    }
  }
  
  factor_id <- which(sapply(new.env, is.factor))
  
  ## 3. Fill the PLOT object ----------------------------------------------------------------------
  if (!do.bivariate)
  {
    ## NON-BIVARIATE CASE -------------------------------------------------------------------------
   
    PROGRESS = txtProgressBar(min = 0, max = length(show.variables) * length(models), style = 3)
    i.iter = 0
    
    list.out = foreach(vari = show.variables) %do%
      {
        if (is.factor(new.env[, vari])) {
          pts.tmp <- as.factor(levels(new.env[, vari]))
        } else {
          pts.tmp <- seq(min(new.env[, vari], na.rm = TRUE), max(new.env[, vari], na.rm = TRUE), length.out = nb.pts)
        }
        
        ## Creating tmp data --------------------------------------------------------------
        tmp = ref_table[, -which(colnames(ref_table) == vari), drop = FALSE]
        new.env.r.tmp <- eval(parse(text = paste0("cbind(", vari, "= pts.tmp, tmp)")))
        new.env.r.tmp <- new.env.r.tmp[, colnames(ref_table), drop = FALSE]
        if (length(factor_id)) {
          for (f in factor_id) {
            new.env.r.tmp[, f] <- factor(as.character(new.env.r.tmp[, f]), levels = levels(new.env[, f]))
          }
        }
        
        ## Load models --------------------------------------------------------------------
        BIOMOD_LoadModels(modeling.output, full.name = models)
        
        ## Getting predictions for each model ---------------------------------------------
        tab.out = foreach(model = models, .combine = "cbind") %do% 
          {
            mod <- get(model)
            mod.name <- ifelse(use.formal.names, formal_names[which(is.element(models, model))], model)
            
            temp_workdir = NULL
            if (length(grep("MAXENT.Phillips$", mod.name)) == 1) {
              temp_workdir = mod@model_output_dir
            }
            proj.tmp <- predict(mod, new.env.r.tmp, on_0_1000 = on_0_1000, do_check = FALSE, temp_workdir = temp_workdir)
            
            res = data.frame(pts.tmp, proj.tmp)
            colnames(res) = c(vari, mod.name)
            
            i.iter <- i.iter + 1
            setTxtProgressBar(pb = PROGRESS, value = i.iter)
            return(res)
          }
        tab.out = tab.out[, c(vari, models)]
        return(tab.out)
      }
    close(PROGRESS)
    
  } else {
    ## BIVARIATE CASE -----------------------------------------------------------------------------
    
    PROGRESS = txtProgressBar(min = 0, max = ncol(combn(show.variables, 2)) * length(models), style = 3)
    i.iter = 0
    
    list.out = foreach(vari1 = show.variables[-length(show.variables)], .combine = "c") %:%
      foreach(vari2 = show.variables[-(1:which(show.variables == vari1))]) %do%
      {
        pts.tmp1 <- rep(seq(min(new.env[, vari1], na.rm = TRUE), max(new.env[, vari1], na.rm = TRUE)
                            , length.out = sqrt(nb.pts)), each = sqrt(nb.pts))
        pts.tmp2 <- rep(seq(min(new.env[, vari2], na.rm = TRUE), max(new.env[, vari2], na.rm = TRUE)
                            , length.out = sqrt(nb.pts)), sqrt(nb.pts))
        
        ## Creating tmp data --------------------------------------------------------------
        tmp = ref_table[, -which(colnames(ref_table) %in% c(vari1, vari2)), drop = FALSE]
        new.env.r.tmp <- eval(parse(text = paste0("cbind(", vari1," = pts.tmp1, ", vari2, " = pts.tmp2, tmp)")))
        new.env.r.tmp <- new.env.r.tmp[, colnames(ref_table), drop = FALSE]
        if (length(factor_id)) {
          for (f in factor_id) {
            new.env.r.tmp[, f] <- factor(as.character(new.env.r.tmp[, f]), levels = levels(new.env[, f]))
          }
        }
        
        ## Load models --------------------------------------------------------------------
        BIOMOD_LoadModels(modeling.output, full.name = models)
        
        ## Getting predictions for each model ---------------------------------------------
        tab.out = foreach(model = models, .combine = "cbind") %do% 
          {
            mod <- get(model)
            mod.name <- ifelse(use.formal.names, formal_names[which(is.element(models, model))], model)
            
            temp_workdir = NULL
            if (length(grep("MAXENT.Phillips$", mod.name)) == 1) {
              temp_workdir = mod@model_output_dir
            }
            
            proj.tmp <- predict(mod, new.env.r.tmp, on_0_1000 = on_0_1000, do_check = FALSE, temp_workdir = temp_workdir)
            
            res = data.frame(pts.tmp1, pts.tmp2, proj.tmp)
            colnames(res) = c(vari1, vari2, mod.name)
            
            i.iter <- i.iter + 1
            setTxtProgressBar(pb = PROGRESS, value = i.iter)
            return(res)
          }
        tab.out = tab.out[, c(vari1, vari2, models)]
        return(tab.out)
      }
    close(PROGRESS)
  }
  
  # transform list.out into ggplot friendly shape
  ggdat <- .as_ggdat(list.out, do.bivariate)
  
  ## 2. PLOT graphic ------------------------------------------------------------------------------
  if (!do.bivariate) {
    new.env_m <- melt(new.env, variable.name = "expl.name", value.name = "expl.val")
    
    gg <- ggplot(ggdat, aes_string(x = "expl.val", y = "pred.val", color = "pred.name")) +
      geom_line() +
      geom_rug(data = new.env_m, sides = 'b', inherit.aes = FALSE, aes_string(x = "expl.val")) +
      facet_wrap("expl.name", scales = "free_x") +
      xlab("") +
      ylab("") +
      theme(legend.title = element_blank()
            , legend.key = element_rect(fill = "white")
            , legend.position = "bottom"
            , axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    comb.names = sort(unique(ggdat$comb))
    list.ggdat = foreach(combi = comb.names) %do%
      {
        vari1 = strsplit(combi, "_")[[1]][1]
        vari2 = strsplit(combi, "_")[[1]][2]
        tmp = ggdat[which(ggdat$comb == combi), ]
        tmp1 = tmp[which(tmp$expl.name == vari1), c("id", "expl.val", "pred.name", "pred.val", "comb")]
        colnames(tmp1)[which(colnames(tmp1) == "expl.val")] = vari1
        tmp2 = tmp[which(tmp$expl.name == vari2), c("id", "expl.val", "pred.name", "pred.val", "comb")]
        colnames(tmp2)[which(colnames(tmp2) == "expl.val")] = vari2
        tmp = merge(tmp1, tmp2, by = c("id", "pred.name", "pred.val", "comb"))
        return(tmp)
      }
    names(list.ggdat) = comb.names
    
    gg <- ggplot(ggdat, aes_string(fill = "pred.val"))
    for (ii in 1:length(list.ggdat)) {
      combi = names(list.ggdat)[ii]
      vari1 = strsplit(combi, "_")[[1]][1]
      vari2 = strsplit(combi, "_")[[1]][2]
      gg <- gg +
        geom_raster(data = list.ggdat[[ii]], aes_string(x = vari1, y = vari2))
    }
    gg <- gg +
      facet_wrap("pred.name ~ sub('_', ' [ x - y ] ', comb)", scales = "free") +
      xlab("") +
      ylab("") +
      scale_fill_viridis_c("Probability") +
      theme(legend.key = element_rect(fill = "white")
            , legend.position = "bottom"
            , axis.text.x = element_text(angle = 45, hjust = 1))
  }

  if (length(main)) { ## add title
    gg <- gg + labs(title = main)
  }

  if (plot){ print(gg) }
  invisible(gg)
}


###################################################################################################

.bm_PlotResponseCurves.check.args <- function(modeling.output, chosen.models, new.env, show.variables, do.bivariate, ...)
{
  args <- list(...)
  
  ## 1. Check modeling.output argument ----------------------------------------
  .fun_testIfInherits(TRUE, "modeling.output", modeling.output, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 2. Check chosen.models -------------------------------------------------
  if (chosen.models[1] == 'all') {
    chosen.models <- modeling.output@models.computed
  } else {
    chosen.models <- intersect(chosen.models, modeling.output@models.computed)
  }
  if (length(chosen.models) < 1) {
    stop('No models selected')
  }
  models <- chosen.models
  
  ## check that given models exist
  files.check <- paste0(modeling.output@sp.name, '/models/', modeling.output@modeling.id, "/", chosen.models)
  not.checked.files <- grep('MAXENT.Phillips|SRE', files.check)
  if (length(not.checked.files) > 0) {
    files.check <- files.check[-not.checked.files]
  }
  missing.files <- files.check[!file.exists(files.check)]
  if (length(missing.files) > 0) {
    stop(paste0("Projection files missing : ", toString(missing.files)))
    if (length(missing.files) == length(files.check)) {
      stop("Impossible to find any models, might be a problem of working directory")
    }
  }
  
  ## 2. Check args$nb.pts argument ----------------------------------------
  ## defining the number split in each variables range
  nb.pts <- args$nb.pts
  if (!is.null(nb.pts) && do.bivariate == TRUE) {
    nb.pts <- nb.pts ^ 2
  } else if (is.null(nb.pts) && do.bivariate == FALSE) {
    nb.pts <- 100
  } else if (is.null(nb.pts) && do.bivariate == TRUE) {
    nb.pts <- 25 ^ 2
  }
  
  ## 3. Check new.env argument ---------------------------------------------------
  if (inherits(new.env, "Raster")) {
    cat("\n   > Extracting raster infos..")
    DataTmp <- matrix(0, ncol = nlayers(new.env), nrow = nb.pts)
    colnames(DataTmp) <- names(new.env)
    maxVal <- maxValue(new.env)
    minVal <- minValue(new.env)
    for (i in 1:ncol(DataTmp)) {
      DataTmp[, i] <- seq(minVal[i], maxVal[i], length.out = nb.pts)
    }
    new.env <- DataTmp
  }
  
  ## 4. Check show.variables argument -----------------------------------------
  if (length(show.variables) > ncol(new.env) || sum(!(show.variables %in% colnames(new.env)))) {
    stop("columns wanted in show.variables do not match the data \n")
  }
  
  # remove factorial var in do.bivariate case
  if (do.bivariate == TRUE) {
    fact_var <- sapply(new.env[, show.variables, drop = FALSE], is.factor)
    if (sum(fact_var) > 0) {
      cat("\n\tFactorial variables have been automatically removed!")
      show.variables <- show.variables[!fact_var]
    }
  }
  
  ## 5. Check main argument ---------------------------------------------------
  main <- args$main
  if (is.null(main)) { main <- try(paste0("Response curves for ", modeling.output@sp.name, "'s models")) }
  
  ## 6. Check data_species argument -------------------------------------------
  data_species <- args$data_species
  if (is.null(data_species)) {
    data_species <- rep(1, nrow(new.env))
  } else {
    data_species[data_species != 1 | is.na(data_species)] <- 0
  }
  
  return(list(models = models,
              new.env = new.env,
              show.variables = show.variables,
              do.bivariate = do.bivariate,
              nb.pts = nb.pts,
              on_0_1000 = ifelse(is.null(args$on_0_1000), FALSE, args$on_0_1000),
              use.formal.names = ifelse(is.null(args$use.formal.names), FALSE, args$use.formal.names),
              main = main,
              data_species = data_species))
}


###################################################################################################

# .get_biomod2_model_object <- function(mod)
# {
#   
#   tmp.time <- paste0("_AllData_", as.character(format(Sys.time(), "%OS6")))
#   if (inherits(mod, "nnet")) {
#     tmp.class = 'ANN'
#   } else if (inherits(mod, "rpart")) {
#     tmp.class = 'CTA'
#   } else if (inherits(mod, "fda")) {
#     tmp.class = 'FDA'
#   } else if (inherits(mod, "gam")) {
#     tmp.class = 'GAM'
#   } else if (inherits(mod, "gbm")) {
#     tmp.class = 'GBM'
#   } else if (inherits(mod, c("glm", "lm"))) {
#     tmp.class = 'GLM'
#   } else if (inherits(mod, "mars")) {
#     tmp.class = 'MARS'
#   } else if (inherits(mod, "randomForest")) {
#     tmp.class = 'RF'
#   }
#   
#   tmp.subclass = ""
#   if (tmp.class %in% c('ANN', 'RF')) {
#     tmp.resp = ifelse(is.null(mod$terms[[2]]), "species", as.character(mod$terms[[2]]))
#     tmp.expl = ifelse(is.character(attr(mod$terms, "term.labels")), attr(mod$terms, "term.labels"), "")
#   } else if (tmp.class %in% c('CTA', 'FDA', 'GAM', 'GBM', 'GLM')) {
#     tmp.resp = as.character(mod$terms[[2]])
#     tmp.expl = attr(mod$terms, "term.labels")
#     if (tmp.class == 'GAM') {
#       tmp.subclass = ifelse(mod$method == "glm.fit", "GAM_gam", "GAM_mgcv")
#     }
#   } else if (tmp.class == 'MARS') {
#     tmp.resp = "species"
#     tmp.expl = as.character(colnames(mod$factor))
#   } else {
#     stop("Unknown model class")
#   }
#   tmp.name = paste0(tmp.resp, tmp.time, "_", tmp.class)
#   
#   return(new(paste0(tmp.class, "_biomod2_model"),
#              model = mod,
#              model_name = tmp.name,
#              model_class = tmp.class,
#              model_subclass = tmp.subclass,
#              resp_name = tmp.resp,
#              expl_var_names = tmp.expl
#   ))
# }


###################################################################################################

.as_ggdat <- function(list.dat, do.bivariate)
{
  col.expl = c(1, 2)
  if (!do.bivariate) { col.expl = c(1) }
  out_ <- foreach(dat_ = list.dat, .combine = "rbind") %do%
    {
      dat_$id <- as.numeric(rownames(dat_))
      if (do.bivariate) {
        dat_$comb <- paste0(colnames(dat_)[col.expl], collapse = "_")
      }
      id.col = which(colnames(dat_) == "id")
      id.vec = unlist(ifelse(do.bivariate, list(c("id", "comb")), "id"))
      expl.dat_ = melt(dat_[, c(col.expl, id.col)], id.vars = "id", variable.name = "expl.name", value.name = "expl.val")
      pred.dat_ = melt(dat_[, -col.expl], id.vars = id.vec, variable.name = "pred.name", value.name = "pred.val")
      out.dat_ = merge(expl.dat_, pred.dat_, by = "id")
      return(out.dat_)
    }
  return(out_)
}

