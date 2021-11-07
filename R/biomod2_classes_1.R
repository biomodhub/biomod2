# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# BIOMOD objects definition
# Damien Georges, Maya Guéguen
# 09/02/2012, update 18/10/2021
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

requireNamespace("raster", quietly=TRUE)
requireNamespace(".0
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# This file defines the BIOMOD objects and all their methods
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# We choose here to create monospecific objects to make all procedures and parallelising easier
requireNamespacrasterVis", quietly=TRUE)

##' @include biomod2_classes_1.R

###################################################################################################
## 1. BIOMOD.formated.data
###################################################################################################

# 1.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.formated.data",
         representation(sp.name = 'character',
                        coord = "data.frame",
                        data.species = "numeric",
                        data.env.var = "data.frame",
                        data.mask = "RasterStack",
                        has.data.eval = "logical",
                        eval.coord = "data.frame",
                        eval.data.species = "numeric",
                        eval.data.env.var = "data.frame"),
         validity = function(object){ return(TRUE) })


# 1.2 Constructors --------------------------------------------------------------------------------
setGeneric("BIOMOD.formated.data", def = function(sp, env, ...) { standardGeneric("BIOMOD.formated.data") })

setMethod('BIOMOD.formated.data', signature(sp = 'numeric', env = 'data.frame'),
          function(sp, env, xy = NULL, sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE, data.mask = NULL)
          {
            if (is.null(data.mask)) { data.mask <- raster::stack() }
            
            if (is.null(eval.sp)) { ## NO EVALUATION DATA -----------------------------------------
              BFD <- new(
                'BIOMOD.formated.data',
                coord = xy,
                data.species = sp,
                data.env.var = env,
                sp.name = sp.name,
                data.mask = data.mask,
                has.data.eval = FALSE
              )
            } else { ## EVALUATION DATA -----------------------------------------------------------
              BFDeval <- BIOMOD.formated.data(
                sp = eval.sp,
                env = eval.env,
                xy = eval.xy,
                sp.name = sp.name
              )
              
              if (raster::nlayers(BFDeval@data.mask) > 0) {
                data.mask.tmp <- try(raster::addLayer(data.mask, BFDeval@data.mask))
                if (!inherits(data.mask.tmp, "try-error")) {
                  data.mask <- data.mask.tmp
                  names(data.mask) <- c("calibration", "validation")
                }
              }
              
              BFD <- new(
                'BIOMOD.formated.data',
                coord = xy,
                data.species = sp,
                data.env.var = env,
                sp.name = sp.name,
                data.mask = data.mask,
                has.data.eval = TRUE,
                eval.coord = BFDeval@coord,
                eval.data.species = BFDeval@data.species,
                eval.data.env.var = BFDeval@data.env.var
              )
              
              rm('BFDeval')
            }
            
            ## REMOVE NA IF ANY -------------------------------------------------------------------
            if (na.rm) {
              rowToRm <- unique(unlist(lapply(BFD@data.env.var, function(x) { return(which(is.na(x))) })))
              if (length(rowToRm)) {
                cat("\n\t\t\t! Some NAs have been automatically removed from your data")
                BFD@coord <- BFD@coord[-rowToRm, , drop = FALSE]
                BFD@data.species <- BFD@data.species[-rowToRm]
                BFD@data.env.var <- BFD@data.env.var[-rowToRm, , drop = FALSE]
              }
              if (BFD@has.data.eval) {
                rowToRm <- unique(unlist(lapply(BFD@eval.data.env.var, function(x) { return(which(is.na(x))) })))
                if (length(rowToRm)) {
                  cat("\n\t\t\t! Some NAs have been automatically removed from your evaluation data")
                  BFD@eval.coord <- BFD@eval.coord[-rowToRm, , drop = FALSE]
                  BFD@eval.data.species <- BFD@eval.data.species[-rowToRm]
                  BFD@eval.data.env.var <- BFD@eval.data.env.var[-rowToRm, , drop = FALSE]
                }
              }
            }
            return(BFD)
          }
)

setMethod('BIOMOD.formated.data', signature(sp = 'data.frame'),
          function(sp, env, xy = NULL, sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE)
          {
            if (ncol(sp) > 1) { stop("Invalid response variable") }
            sp <- as.numeric(unlist(sp))
            BFD <- BIOMOD.formated.data(sp, env, xy, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm)
            return(BFD)
          }
)

setMethod('BIOMOD.formated.data', signature(sp = 'numeric', env = 'matrix'),
          function(sp, env, xy = NULL, sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE)
          {
            env <- as.data.frame(env)
            BFD <- BIOMOD.formated.data(sp, env, xy, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm)
            return(BFD)
          }
)

setMethod('BIOMOD.formated.data', signature(sp = 'numeric', env = 'RasterStack'),
          function(sp, env, xy = NULL, sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE)
          {
            categorial_var <- names(env)[raster::is.factor(env)]
            
            ## Keep same env variable for eval than calib (+ check for factor)
            if (!is.null(eval.sp) && is.null(eval.env)) {
              eval.env <- as.data.frame(extract(env, eval.xy))
              if (length(categorial_var)) {
                for (cat_var in categorial_var) {
                  eval.env[, cat_var] <- as.factor(eval.env[, cat_var])
                }
              }
            }
            
            if (is.null(xy)) { xy <- as.data.frame(coordinates(env)) }
            
            ## Prepare mask of studied area
            data.mask = reclassify(raster::subset(env, 1, drop = T), c(-Inf, Inf, -1))
            data.mask[cellFromXY(data.mask, xy[which(sp == 1), ])] <- 1
            data.mask[cellFromXY(data.mask, xy[which(sp == 0), ])] <- 0
            data.mask <- raster::stack(data.mask)
            names(data.mask) <- sp.name
            
            ## Keep same env variable for eval than calib (+ check for factor)
            env <- as.data.frame(extract(env, xy, factors = T))
            if (length(categorial_var)) {
              for (cat_var in categorial_var) {
                env[, cat_var] <- as.factor(env[, cat_var])
              }
            }
            
            BFD <- BIOMOD.formated.data(sp, env, xy, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm, data.mask = data.mask)
            return(BFD)
          }
)


# 1.3 Other Functions -----------------------------------------------------------------------------
setMethod('plot', signature(x = 'BIOMOD.formated.data', y = "missing"),
          function(x, coord = NULL, col = NULL)
          {
            if (raster::nlayers(x@data.mask) > 0)
            {
              requireNamespace("rasterVis")
              
              ## check if there is some undefined areas to prevent from strange plotting issues
              if (min(cellStats(x@data.mask, min)) == -1) { # there is undefined area
                my.at <- seq(-1.5, 1.5, by = 1) ## breaks of color key
                my.labs.at <- seq(-1, 1, by = 1) ## labels placed vertically
                my.lab <- c("undifined", "absences", "presences") ## labels
                my.col.regions = c("lightgrey", "red4", "green4") ## colors
                my.cuts <- 2 ## cuts
              } else { # no undefined area.. remove it from plot
                my.at <- seq(-0.5, 1.5, by = 1) ## breaks of color key
                my.labs.at <- seq(0, 1, by = 1) ## labels placed vertically
                my.lab <- c("absences", "presences") ## labels
                my.col.regions = c("red4", "green4") ## colors
                my.cuts <- 1 ## cuts
              }
              
              ## PLOT -----------------------------------------------------------------------------
              levelplot(
                x@data.mask,
                at = my.at,
                cuts = my.cuts,
                margin = T,
                col.regions = my.col.regions,
                main = paste(x@sp.name, "datasets"),
                colorkey = list(labels = list(labels = my.lab, at = my.labs.at))
              )
              
            } else
            {
              # coordinates checking
              if (is.null(coord)) {
                if (sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2]) {
                  stop("coordinates are required to plot your data")
                } else {
                  coord <- x@coord
                }
              }
              
              # colors checking
              if (is.null(col) | length(col) < 3) { col = c('green', 'red', 'grey') }
              
              
              ## PLOT -----------------------------------------------------------------------------
              ## all points (~ mask)
              plot(
                x = x@coord[, 1],
                y = x@coord[, 2],
                col = col[3],
                xlab = 'X',
                ylab = 'Y',
                main = x@sp.name,
                pch = 20
              )
              ## presences
              points(
                x = x@coord[which(x@data.species == 1), 1],
                y = x@coord[which(x@data.species == 1), 2],
                col = col[1],
                pch = 18
              )
              ## true absences
              points(
                x = x@coord[which(x@data.species == 0), 1],
                y = x@coord[which(x@data.species == 0), 2],
                col = col[2],
                pch = 18
              )
            }
          }
)

##' @rdname BIOMOD.formated.data-objects
##' @docType method
##' @aliases show, BIOMOD.formated.data-method
setMethod('show', signature('BIOMOD.formated.data'),
          function(object)
          {
            .bmCat("'BIOMOD.formated.data'")
            cat("\nsp.name = ", object@sp.name, fill = .Options$width)
            cat(
              "\n\t",
              sum(object@data.species, na.rm = TRUE),
              'presences, ',
              sum(object@data.species == 0, na.rm = TRUE),
              'true absences and ',
              sum(is.na(object@data.species), na.rm = TRUE),
              'undifined points in dataset',
              fill = .Options$width
            )
            cat("\n\n\t",
                ncol(object@data.env.var),
                'explanatory variables\n',
                fill = .Options$width)
            print(summary(object@data.env.var))
            
            if (object@has.data.eval) {
              cat("\n\nEvaluation data :", fill = .Options$width)
              cat(
                "\n\t",
                sum(object@eval.data.species, na.rm = TRUE),
                'presences, ',
                sum(object@eval.data.species == 0, na.rm = TRUE),
                'true absences and ',
                sum(is.na(object@eval.data.species), na.rm = TRUE),
                'undifined points in dataset',
                fill = .Options$width
              )
              cat("\n\n", fill = .Options$width)
              print(summary(object@eval.data.env.var))
            }
            
            .bmCat()
          }
)



###################################################################################################
## 2. BIOMOD.formated.data.PA 
## this class inherits from BIOMOD.formated.data and have one more slot 'PA' giving PA selected
###################################################################################################

# 2.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.formated.data.PA",
         contains = "BIOMOD.formated.data",
         representation(PA.strategy = 'character', PA = 'data.frame'),
         validity = function(object) { return(TRUE) })


# 2.2 Constructors --------------------------------------------------------------------------------
# setGeneric("BIOMOD.formated.data.PA", def = function(sp, env, ...) { standardGeneric("BIOMOD.formated.data.PA") })
# 
# setMethod('BIOMOD.formated.data.PA', signature(sp = 'numeric', env = 'data.frame'),
#           function(sp, env, xy = NULL, sp.name = NULL
#                    , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
#                    , na.rm = TRUE, data.mask = NULL)
#           {

BIOMOD.formated.data.PA <-  function(sp, env, xy, sp.name
                                     , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                                     , PA.NbRep = 1, PA.strategy = 'random', PA.nb.absences = NULL
                                     , PA.dist.min = 0, PA.dist.max = NULL
                                     , PA.sre.quant = 0.025, PA.table = NULL
                                     , na.rm = TRUE)
{
  
  categorial_var <- NULL
  if (inherits(env, 'Raster')) { categorial_var <- names(env)[raster::is.factor(env)] }
  
  ## Keep same env variable for eval than calib (+ check for factor)
  if (!is.null(eval.sp) && is.null(eval.env)) {
    if (inherits(env, 'Raster')) {
      eval.env <- as.data.frame(extract(env, eval.xy))
      if (length(categorial_var)) {
        for (cat_var in categorial_var) {
          eval.env[, cat_var] <- as.factor(eval.env[, cat_var])
        }
      }
    } else { stop("No evaluation explanatory variable given") }
  }
  
  # Convert sp in SpatialPointsDataFrame
  if (is.numeric(sp)) {
    if (is.null(xy)) {
      sp <- SpatialPointsDataFrame(matrix(0, ncol = 2, nrow = length(sp)), data.frame(sp), match.ID = FALSE)
    } else {
      sp <- SpatialPointsDataFrame(data.matrix(xy), data.frame(sp), match.ID = FALSE)
    }
  }
  
  pa.data.tmp <- .pseudo.absences.sampling(sp = sp,
                                           env = env,
                                           nb.repet = PA.NbRep,
                                           strategy = PA.strategy,
                                           nb.points = PA.nb.absences,
                                           distMin = PA.dist.min,
                                           distMax = PA.dist.max,
                                           quant.SRE = PA.sre.quant,
                                           PA.table = PA.table)
  
  if (!is.null(pa.data.tmp))
  {
    ## Keep same env variable for eval than calib (+ check for factor)
    if (length(categorial_var)) {
      for (cat_var in categorial_var) {
        pa.data.tmp$env[, cat_var] <- as.factor(pa.data.tmp$env[, cat_var])
      }
    }
    
    ## REMOVE NA IF ANY ---------------------------------------------------------------------------
    if (na.rm) {
      rowToRm <- unique(unlist(lapply(pa.data.tmp$env, function(x) { return(which(is.na(x))) })))
      if (length(rowToRm)) {
        cat("\n\t\t\t! Some NAs have been automatically removed from your data")
        pa.data.tmp$xy <- pa.data.tmp$xy[-rowToRm, , drop = FALSE]
        pa.data.tmp$sp <- pa.data.tmp$sp[-rowToRm, drop = FALSE]
        pa.data.tmp$env <- pa.data.tmp$env[-rowToRm, , drop = FALSE]
        pa.data.tmp$pa.tab <- pa.data.tmp$pa.tab[-rowToRm, , drop = FALSE]
      }
    }
    
    BFD <- BIOMOD.formated.data(sp = pa.data.tmp$sp,
                                env = pa.data.tmp$env,
                                xy = as.data.frame(pa.data.tmp$xy),
                                sp.name = sp.name, 
                                eval.sp = eval.sp,
                                eval.env = eval.env,
                                eval.xy = eval.xy,
                                na.rm = na.rm)
    
    if(inherits(env,'Raster'))
    {
      ## create data.mask for ploting
      data.mask.tmp <- reclassify(raster::subset(env, 1), c(-Inf, Inf, -1))
      data.mask <- stack(data.mask.tmp)
      xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp == 1), , drop = FALSE]
      xy_abs <- pa.data.tmp$xy[which(pa.data.tmp$sp == 0), , drop = FALSE]
      if (nrow(xy_pres)) { data.mask[cellFromXY(data.mask.tmp, xy_pres)] <- 1 }
      if (nrow(xy_abs)) { data.mask[cellFromXY(data.mask.tmp, xy_abs)] <- 0 }
      names(data.mask) <- "input_data"
      
      ## add eval data
      if(BFD@has.data.eval){ } ### TO DO
      
      ## add pa data
      for(pa in 1:ncol(as.data.frame(pa.data.tmp$pa.tab)))
      {
        data.mask.tmp2 <- data.mask.tmp
        
        ind.pa <- as.data.frame(pa.data.tmp$pa.tab)[, pa]
        xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp == 1 & ind.pa), , drop = FALSE]
        xy_abs <- pa.data.tmp$xy[which((pa.data.tmp$sp != 1 | is.na(pa.data.tmp$sp)) & ind.pa), , drop = FALSE]
        if (nrow(xy_pres)) {
          id_pres <- cellFromXY(data.mask.tmp, xy_pres)
          data.mask.tmp2[id_pres] <- 1
        }
        if (nrow(xy_abs)) {
          id_abs <- cellFromXY(data.mask.tmp, xy_abs)
          data.mask.tmp2[id_abs] <- 0
        }
        
        data.mask <- addLayer(data.mask, data.mask.tmp2)
      }
      names(data.mask) <- c("input_data", colnames(as.data.frame(pa.data.tmp$pa.tab)))
      
    } else {  data.mask <- stack() }
    
    BFDP <- new('BIOMOD.formated.data.PA',
                sp.name = BFD@sp.name,
                coord = BFD@coord,
                data.env.var = BFD@data.env.var,
                data.species = BFD@data.species,
                data.mask = data.mask,
                has.data.eval = BFD@has.data.eval,
                eval.coord = BFD@eval.coord,
                eval.data.species = BFD@eval.data.species,
                eval.data.env.var = BFD@eval.data.env.var,
                PA = as.data.frame(pa.data.tmp$pa.tab),
                PA.strategy = PA.strategy)
    
    rm(list='BFD')
  } else
  {
    cat("\n   ! PA selection not done", fill = .Options$width)
    
    BFDP <- BIOMOD.formated.data(sp = sp,
                                 env = env,
                                 xy = xy,
                                 sp.name = sp.name, 
                                 eval.sp = eval.sp,
                                 eval.env = eval.env,
                                 eval.xy = eval.xy)
  }
  rm(list = "pa.data.tmp" )
  return(BFDP)
}


# 2.3 other functions -----------------------------------------------------------------------------
setMethod('plot', signature(x = 'BIOMOD.formated.data.PA', y = "missing"),
          function(x, coord = NULL, col = NULL)
          {
            if (raster::nlayers(x@data.mask) > 0)
            {
              requireNamespace("rasterVis")
              
              ## check if there is some undefined areas to prevent from strange plotting issues
              if (min(cellStats(x@data.mask, min)) == -1) { # there is undefined area
                my.at <- seq(-1.5, 1.5, by = 1) ## breaks of color key
                my.labs.at <- seq(-1, 1, by = 1) ## labels placed vertically
                my.lab <- c("undifined", "absences", "presences") ## labels
                my.col.regions = c("lightgrey", "red4", "green4") ## colors
                my.cuts <- 2 ## cuts
              } else { # no undefined area.. remove it from plot
                my.at <- seq(-0.5, 1.5, by = 1) ## breaks of color key
                my.labs.at <- seq(0, 1, by = 1) ## labels placed vertically
                my.lab <- c("absences", "presences") ## labels
                my.col.regions = c("red4", "green4") ## colors
                my.cuts <- 1 ## cuts
              }
              
              ## PLOT -----------------------------------------------------------------------------
              levelplot(
                x@data.mask,
                at = my.at,
                cuts = my.cuts,
                margin = T,
                col.regions = my.col.regions,
                main = paste(x@sp.name, "datasets"),
                colorkey = list(labels = list(labels = my.lab, at = my.labs.at))
              )
            } else
            {
              # coordinates checking
              if (is.null(coord)) {
                if (sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2]) {
                  stop("coordinates are required to plot your data")
                } else {
                  coord <- x@coord
                }
              }
              
              # colors checking
              if (is.null(col) | length(col) < 3) { col = c('green', 'red', 'orange', 'grey') }
              
              ## PLOT -----------------------------------------------------------------------------
              par(mfrow = c(.clever_cut(ncol(x@PA) + 1)))
              
              # all points (~mask)
              plot(
                x = x@coord[, 1],
                y = x@coord[, 2],
                col = col[4],
                xlab = 'X',
                ylab = 'Y',
                main = paste(x@sp.name, " original data", sep = ""),
                pch = 20
              )
              # presences
              points(
                x = x@coord[which(x@data.species == 1), 1],
                y = x@coord[which(x@data.species == 1), 2],
                col = col[1],
                pch = 18
              )
              # true absences
              points(
                x = x@coord[which(x@data.species == 0), 1],
                y = x@coord[which(x@data.species == 0), 2],
                col = col[2],
                pch = 18
              )
              # PA data
              for(i in 1:ncol(x@PA))
              {
                # all points (~mask)
                plot(
                  x = x@coord[, 1],
                  y = x@coord[, 2],
                  col = col[4],
                  xlab = 'X',
                  ylab = 'Y',
                  main = paste0(x@sp.name, " Pseudo Absences ", i),
                  pch = 20
                )
                # presences
                points(
                  x = x@coord[(x@data.species == 1) & x@PA[, i], 1],
                  y = x@coord[(x@data.species == 1) & x@PA[, i], 2],
                  col = col[1],
                  pch = 18
                )
                # true absences
                points(
                  x = x@coord[(x@data.species == 0) & x@PA[, i], 1],
                  y = x@coord[(x@data.species == 0) & x@PA[, i], 2],
                  col = col[2],
                  pch = 18
                )
                # PA
                points(
                  x = x@coord[is.na(x@data.species) & x@PA[, i], 1],
                  y = x@coord[is.na(x@data.species) & x@PA[, i], 2],
                  col = col[3],
                  pch = 18
                )
              }
            }
          }
)

##' @rdname BIOMOD.formated.data.PA-objects
##' @docType method
##' @aliases show, BIOMOD.formated.data.PA-method
setMethod('show', signature('BIOMOD.formated.data.PA'),
          function(object)
          {
            .bmCat("'BIOMOD.formated.data.PA'")
            cat("\nsp.name = ", object@sp.name, fill = .Options$width)
            cat(
              "\n\t",
              sum(object@data.species, na.rm = TRUE),
              'presences, ',
              sum(object@data.species == 0, na.rm = TRUE),
              'true absences and ',
              sum(is.na(object@data.species), na.rm = TRUE),
              'undifined points in dataset',
              fill = .Options$width
            )
            cat("\n\n\t",
                ncol(object@data.env.var),
                'explanatory variables\n',
                fill = .Options$width)
            print(summary(object@data.env.var))
            
            if (object@has.data.eval) {
              cat("\n\nEvaluation data :", fill = .Options$width)
              cat(
                "\n\t",
                sum(object@eval.data.species, na.rm = TRUE),
                'presences, ',
                sum(object@eval.data.species == 0, na.rm = TRUE),
                'true absences and ',
                sum(is.na(object@eval.data.species), na.rm = TRUE),
                'undifined points in dataset',
                fill = .Options$width
              )
              cat("\n\n", fill = .Options$width)
              print(summary(object@eval.data.env.var))
            }
            
            cat(
              "\n\n",
              ncol(object@PA),
              'Pseudo Absences dataset available (',
              colnames(object@PA),
              ") with ",
              sum(object@PA[, 1], na.rm = T) - sum(object@data.species, na.rm = TRUE),
              'absences in each (true abs + pseudo abs)',
              fill = .Options$width
            )
            .bmCat()
          }
)



###################################################################################################
## 3. BIOMOD.Model.Options
###################################################################################################


##' @rdname BIOMOD.Model.Options-objects
##' @name BIOMOD.Model.Options-class
##' @docType class
##' @aliases   BIOMOD.Model.Options-class
##' 
##' @title BIOMOD_ModelingOptions outputs objects class
##' 
##' @description
##' BIOMOD.Model.Options objects are created, used and returned
##' by BIOMOD functions. These objects will contains for each
##' model support within \pkg{biomod2}, a set of options that
##' users can change. Please refer to 
##' \code{\link[biomod2]{BIOMOD_ModelingOptions}} for further
##'  details. 
##'   
##' - output of: \code{\link[biomod2]{BIOMOD_ModelingOptions}}
##' - input of:  \code{\link[biomod2]{BIOMOD_Modeling}}
##' 
##' @param object init list of options
##' 
##' @details  
##' Please refer to \code{\link[biomod2]{BIOMOD_ModelingOptions}}
##' for each model arguments supported.
##' 
##' @slot GLM "list", list of GLM supported options
##' @slot GBM "list", list of GBM supported options 
##' @slot GAM "list", list of GAM supported options
##' @slot CTA "list", list of CTA supported options
##' @slot ANN "list", list of ANN supported options
##' @slot SRE "list", list of SRE supported options
##' @slot FDA "list", list of FDA supported options
##' @slot MARS "list", list of MARS supported options
##' @slot RF "list", list of RF supported options
##' @slot MAXENT.Phillips "list", list of MAXENT.Phillips
##'   supported options
##' @slot MAXENT.Phillips.2 "list", list of maxnet
##'   supported options
##'   
##' @author Damien Georges
##' @seealso \code{\link[biomod2]{BIOMOD_ModelingOptions}}
##' @keywords models
##' @keywords options
##' 
##' @examples
##' showClass("BIOMOD.Model.Options")
##' 


setClass("BIOMOD.Model.Options",
         representation(GLM = "list",
                        GBM = "list",
                        GAM = "list",
                        CTA = "list",
                        ANN = "list",
                        SRE = "list",
                        FDA = "list",
                        MARS = "list",
                        RF = "list",
                        MAXENT.Phillips = "list",
                        MAXENT.Phillips.2 = "list"),
         prototype(GLM = list(type = 'quadratic',
                              interaction.level = 0,
                              myFormula = NULL,
                              test = 'AIC',
                              family = binomial(link = 'logit'),
                              mustart = 0.5,
                              control = glm.control(maxit = 50)),
                   GBM = list(distribution = 'bernoulli',
                              n.trees = 2500,
                              interaction.depth = 7,
                              n.minobsinnode = 5,
                              shrinkage = 0.001,
                              bag.fraction = 0.5,
                              train.fraction = 1,
                              cv.folds = 3,
                              keep.data = FALSE,
                              verbose = FALSE,
                              # class.stratify.cv = 'bernoulli',
                              perf.method = 'cv',
                              n.cores = 1),
                   GAM = list(algo = "GAM_mgcv",
                              type = "s_smoother",
                              k = NULL,
                              interaction.level = 0,
                              myFormula = NULL,
                              family = binomial(link = 'logit'),
                              control = list(epsilon = 1e-06, trace = FALSE, maxit = 100),
                              method = "GCV.Cp",
                              optimizer = c("outer", "newton"),
                              select = FALSE,
                              knots = NULL,
                              paraPen = NULL),
                   CTA = list(method = 'class',
                              parms = 'default',
                              # control = rpart.control(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25),
                              control = list(xval = 5, 
                                             minbucket = 5, 
                                             minsplit = 5,
                                             cp = 0.001, 
                                             maxdepth = 25),
                              cost = NULL),
                   ANN = list(NbCV = 5,
                              size = NULL,
                              decay = NULL,
                              rang = 0.1,
                              maxit = 200),
                   SRE = list(quant = 0.025),
                   FDA = list(method = 'mars', add_args = NULL),
                   MARS = list(type = 'simple',
                               interaction.level = 0,
                               myFormula = NULL,
                               # degree = 1,
                               nk = NULL,
                               penalty = 2,
                               thresh = 0.001,
                               nprune = NULL,
                               pmethod = 'backward'),
                   RF = list(do.classif = TRUE,
                             ntree = 500,
                             mtry = 'default',
                             nodesize = 5,
                             maxnodes= NULL),
                   MAXENT.Phillips = list(path_to_maxent.jar = getwd(),
                                          memory_allocated = 512,
                                          background_data_dir = 'default',
                                          maximumbackground = 'default',
                                          maximumiterations = 200,
                                          visible = FALSE,
                                          linear = TRUE,
                                          quadratic = TRUE,
                                          product = TRUE,
                                          threshold = TRUE,
                                          hinge = TRUE,
                                          lq2lqptthreshold = 80,
                                          l2lqthreshold = 10,
                                          hingethreshold = 15,
                                          beta_threshold = -1.0,
                                          beta_categorical = -1.0,
                                          beta_lqp = -1.0,
                                          beta_hinge = -1.0,
                                          betamultiplier = 1,
                                          defaultprevalence = 0.5),
                   MAXENT.Phillips.2 = list(myFormula = NULL,
                                            regmult = 1,
                                            regfun = maxnet::maxnet.default.regularization)
         ),
         validity = function(object)
         {
           test <- TRUE
           
           ## GLM ##
           test <- .fun_testIfIn(test, "GLM$type", object@GLM$type, c('simple', 'quadratic', 'polynomial', 'user.defined'))
           test <- .fun_testIfPosInt(test, "GLM$interaction.level", object@GLM$interaction.level)
           if(!is.null(object@GLM$myFormula)) if(!inherits(object@GLM$myFormula, "formula")){ cat("\nGLM$myFormula must be NULL or a formula object"); test <- FALSE }
           test <- .fun_testIfIn(test, "GLM$test", object@GLM$test, c('AIC', 'BIC', 'none'))
           fam <- 'none'
           if(!inherits(object@GLM$family, "family")){ cat("\nGLM$family must be a valid family object"); test <- FALSE }
           if(!is.list(object@GLM$control)){cat("\nGLM$control must be a list like that returned by glm.control"); test <- FALSE}
           
           ## GBM ##
           test <- .fun_testIfIn(test, "GBM$distribution", object@GBM$distribution, c("bernoulli", "huberized", "multinomial", "adaboost"))
           # test <- .fun_testIfPosInt(test, "GBM$n.trees", object@GBM$n.trees)
           if(!is.numeric(object@GBM$n.trees)){ cat("\nGBM$n.trees must be a integer"); test <- FALSE } else{
             if(object@GBM$n.trees < 0 | floor(object@GBM$n.trees) != object@GBM$n.trees){ cat("\nGBM$n.trees must be a positive integer"); test <- FALSE }
           }
           test <- .fun_testIfPosInt(test, "GBM$interaction.depth", object@GBM$interaction.depth)
           test <- .fun_testIfPosInt(test, "GBM$n.minobsinnode", object@GBM$n.minobsinnode)
           test <- .fun_testIfPosNum(test, "GBM$shrinkage", object@GBM$shrinkage)
           test <- .fun_testIf01(test, "GBM$bag.fraction", object@GBM$bag.fraction)
           test <- .fun_testIf01(test, "GBM$train.fraction", object@GBM$train.fraction)
           test <- .fun_testIfPosInt(test, "GBM$cv.folds", object@GBM$cv.folds)
           if(!is.logical(object@GBM$keep.data)){ cat("\nGBM$keep.data must be a logical"); test <- FALSE }
           if(!is.logical(object@GBM$verbose)){ cat("\nGBM$verbose must be a logical"); test <- FALSE }
           # test <- .fun_testIfIn(test, "GBM$class.stratify.cv", object@GBM$class.stratify.cv, c('bernoulli', 'multinomial'))
           test <- .fun_testIfIn(test, "GBM$perf.method", object@GBM$perf.method, c('OOB', 'test', 'cv'))
           
           ## GAM ##
           test <- .fun_testIfIn(test, "GAM$algo", object@GAM$algo, c('GAM_mgcv', 'GAM_gam', 'BAM_mgcv'))
           test <- .fun_testIfIn(test, "GAM$type", object@GAM$type, c('s_smoother', 's', 'lo', 'te'))
           if(! is.null(object@GAM$k)){
             if(! is.numeric(object@GAM$k)  ){ cat("\nGAM$k must be a integer"); test <- FALSE } else{
               if(object@GAM$k < -1 | object@GAM$k%%1!=0){ cat("\nGAM$k must be > -1"); test <- FALSE }
             }
           }
           test <- .fun_testIfPosInt(test, "GAM$interaction.level", object@GAM$interaction.level)
           if(!is.null(object@GAM$myFormula)) if(!inherits(object@GAM$myFormula, "formula")){ cat("\nGAM$myFormula must be NULL or a formula object"); test <- FALSE }
           if(!inherits(object@GAM$family, "family")){ cat("\nGAM$family must be a valid family object"); test <- FALSE }
           if(!is.list(object@GAM$control)){cat("\nGAM$control must be a list like that returned by gam.control"); test <- FALSE}
           test <- .fun_testIfIn(test, "GAM$method", object@GAM$method, c('GCV.Cp', 'GACV.Cp', 'REML', 'P-REML', 'ML', 'P-ML'))
           if(sum(! object@GAM$optimizer %in% c('perf','outer', 'newton', 'bfgs', 'optim', 'nlm', 'nlm.fd')) > 0 ){cat("\nGAM$optimizer bad definition (see ?mgcv::gam)") ; test <- FALSE}
           if(!is.logical(object@GAM$select)){ cat("\nGAM$select must be a logical"); test <- FALSE }
           #            knots=NULL,
           #            paraPen=NULL
           
           ## CTA ##
           test <- .fun_testIfIn(test, "CTA$method", object@CTA$method, c('anova', 'poisson', 'class', 'exp'))
           #parms = 'default',
           if(!is.list(object@CTA$control)){cat("\nCTA$control must be a list like that returned by rpart.control"); test <- FALSE}
           if(length(object@CTA$cost)){
             if(!is.numeric(object@CTA$cost)){cat("\nCTA$cost must be a non negative cost vector"); test <- FALSE}
           }
           
           ## ANN ##
           test <- .fun_testIfPosInt(test, "ANN$NbCV", object@ANN$NbCV)
           if( ( is.null(object@ANN$size) | length(object@ANN$size)>1 ) & object@ANN$NbCV <= 0){ cat("\nANN$size has to be defined as a single integer if ANN$NbCV=0"); test <- FALSE } else{
             if(!is.null(object@ANN$size)) if( !is.numeric(object@ANN$size) | !all( object@ANN$size > 0 ) | !all( object@ANN$size %% 1 == 0 ) ){ cat("\nANN$size must be NULL or a positive (vector of) integer"); test <- FALSE }
           }
           if( ( is.null(object@ANN$decay) | length(object@ANN$decay)>1 ) & object@ANN$NbCV <= 0){ cat("\nANN$decay has to be defined as a single number if ANN$NbCV=0"); test <- FALSE } else{
             if(!is.null(object@ANN$decay)) if( !is.numeric(object@ANN$decay) | !all( object@ANN$decay > 0 ) ){ cat("\nANN$decay must be NULL or a positive (vector of) number"); test <- FALSE }
           }
           test <- .fun_testIfPosNum(test, "ANN$rang", object@ANN$rang)
           test <- .fun_testIfPosInt(test, "ANN$maxit", object@ANN$maxit)
           
           ## FDA ##
           test <- .fun_testIfIn(test, "FDA$method", object@FDA$method, c('polyreg', 'mars', 'bruto'))
           if(!is.null(object@FDA$add_args)){ if(!is.list(object@FDA$add_args)) {cat("\nFDA$add_args must be a list or NULL"); test <- FALSE } }
           
           ## SRE ##
           if(!is.numeric(object@SRE$quant)){ cat("\nSRE$quant must be a numeric"); test <- FALSE } else{
             if(object@SRE$quant >= 0.5 | object@SRE$quant < 0){ cat("\nSRE$quant must between 0 and 0.5"); test <- FALSE }
           }
           
           ## MARS ##
           test <- .fun_testIfIn(test, "MARS$type", object@MARS$type,  c('simple', 'quadratic', 'polynomial', 'user.defined'))
           test <- .fun_testIfPosInt(test, "MARS$interaction.level", object@MARS$interaction.level)
           if(!is.null(object@MARS$myFormula)) if(!inherits(object@MARS$myFormula, "formula")){ cat("\nMARS$myFormula must be NULL or a formula object"); test <- FALSE }
           # test <- .fun_testIfPosInt(test, "MARS$degree", object@MARS$degree)
           if(!is.null(object@MARS$nk)){
             if(object@MARS$nk < 0 | object@MARS$nk%%1!=0){ cat("\nMARS$nk must be a positive integer or NULL if you want to use default parameter"); test <- FALSE }
           }
           test <- .fun_testIfPosInt(test, "MARS$penalty", object@MARS$penalty)
           test <- .fun_testIfPosNum(test, "MARS$thresh", object@MARS$thresh)
           if(!is.null(object@MARS$nprune)){ if(!is.numeric(object@MARS$nprune)){ cat("\nMARS$nprune must be a numeric or NULL"); test <- FALSE }}
           supported.pmethod <- c('backward', 'none', 'exhaustive', 'forward', 'seqrep', 'cv')
           if(!is.element(object@MARS$pmethod, supported.pmethod)){cat("\nMARS$pmethod must be a one of", supported.pmethod); test <- FALSE }
           
           ## RF ##
           if(!is.logical(object@RF$do.classif)){ cat("\nRF$do.classif must be a logical"); test <- FALSE }
           test <- .fun_testIfPosInt(test, "RF$ntree", object@RF$ntree)
           if (object@RF$mtry != 'default') { test <- .fun_testIfPosInt(test, "RF$mtry", object@RF$mtry) }
           test <- .fun_testIfPosInt(test, "RF$nodesize", object@RF$nodesize)
           if(length(object@RF$maxnodes)) { test <- .fun_testIfPosInt(test, "RF$maxnodes", object@RF$maxnodes) }
           
           ## MAXENT.Phillips ##
           if(!is.character(object@MAXENT.Phillips$path_to_maxent.jar)){ cat("\nMAXENT.Phillips$path_to_maxent.jar must be a character"); test <- FALSE }
           if(!is.null(object@MAXENT.Phillips$memory_allocated)){
             if(!is.numeric(object@MAXENT.Phillips$memory_allocated)){
               cat("\nMAXENT.Phillips$memory_allocated must be a positive integer or NULL for unlimited memory allocation"); test <- FALSE }
           }
           if(!is.character(object@MAXENT.Phillips$background_data_dir)){ cat("\nMAXENT.Phillips$background_data_dir must be 'default' (=> use the same pseudo absences than other models as background) or a path to the directory where your environmental layer are stored"); test <- FALSE }
           tt <- is.character(object@MAXENT.Phillips$maximumbackground) | is.numeric(object@MAXENT.Phillips$maximumbackground)
           if(is.character(object@MAXENT.Phillips$maximumbackground)) if(object@MAXENT.Phillips$maximumbackground != 'default') tt <- FALSE
           if(!tt){ cat("\nMAXENT.Phillips$maximumbackground must be 'default' or numeric"); test <- FALSE }
           test <- .fun_testIfPosInt(test, "MAXENT.Phillips$maximumiterations", object@MAXENT.Phillips$maximumiterations)
           if(!is.logical(object@MAXENT.Phillips$visible)){ cat("\nMAXENT.Phillips$visible must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$linear)){ cat("\nMAXENT.Phillips$linear must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$quadratic)){ cat("\nMAXENT.Phillips$quadratic must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$product)){ cat("\nMAXENT.Phillips$product must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$threshold)){ cat("\nMAXENT.Phillips$threshold must be a logical"); test <- FALSE }
           if(!is.logical(object@MAXENT.Phillips$hinge)){ cat("\nMAXENT.Phillips$hinge must be a logical"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$lq2lqptthreshold)){ cat("\nMAXENT.Phillips$lq2lqptthreshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$l2lqthreshold)){ cat("\nMAXENT.Phillips$l2lqthreshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$lq2lqptthreshold)){ cat("\nMAXENT.Phillips$lq2lqptthreshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$hingethreshold)){ cat("\nMAXENT.Phillips$hingethreshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$beta_threshold)){ cat("\nMAXENT.Phillips$beta_threshold must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$beta_categorical)){ cat("\nMAXENT.Phillips$beta_categorical must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$beta_lqp)){ cat("\nMAXENT.Phillips$beta_lqp must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$beta_hinge)){ cat("\nMAXENT.Phillips$beta_hinge must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$betamultiplier)){ cat("\nMAXENT.Phillips$betamultiplier must be a numeric"); test <- FALSE }
           if(!is.numeric(object@MAXENT.Phillips$defaultprevalence)){ cat("\nMAXENT.Phillips$defaultprevalence must be a numeric"); test <- FALSE }
           
           ## MAXENT.Phillips.2 (MAXENT.Tsuruoka)
           ### TO BE DONE ===
           # 		       if(!is.numeric(object@MAXENT.Tsuruoka$l1_regularizer)){ cat("\nMAXENT.Tsuruoka$l1_regularizer must be a numeric"); test <- FALSE }
           # 		       if(!is.numeric(object@MAXENT.Tsuruoka$l2_regularizer)){ cat("\nMAXENT.Tsuruoka$l2_regularizer must be a numeric"); test <- FALSE }
           # 		       if(!is.logical(object@MAXENT.Tsuruoka$use_sgd)){ cat("\nMAXENT.Tsuruoka$use_sgd must be a logical"); test <- FALSE }
           # 		       if(!is.numeric(object@MAXENT.Tsuruoka$set_heldout)){ cat("\nMAXENT.Tsuruoka$set_heldout must be a numeric"); test <- FALSE }
           # 		       if(!is.logical(object@MAXENT.Tsuruoka$verbose)){ cat("\nMAXENT.Tsuruoka$verbose must be a logical"); test <- FALSE }
           
           return(test)
         }
)



##' @rdname BIOMOD.Model.Options-objects
setMethod('show', signature('BIOMOD.Model.Options'),
          function(object)
          {
            .bmCat(" 'BIOMOD.Model.Options' ")
            cat("\n")
            
            ## GLM options
            cat("\nGLM = list( type = '", object@GLM$type, "',", sep = "")
            cat("\n            interaction.level = ", object@GLM$interaction.level, ",", sep = "")
            cat("\n            myFormula = ",  ifelse(length(object@GLM$myFormula) < 1, 'NULL', paste(object@GLM$myFormula[2],
                                                                                                      object@GLM$myFormula[1],
                                                                                                      object@GLM$myFormula[3])), ",", sep = "")
            cat("\n            test = '", object@GLM$test, "',", sep = "")
            cat("\n            family = ", object@GLM$family$family, "(link = '", object@GLM$family$link, "'),", sep = "")
            cat("\n            mustart = ", object@GLM$mustart, ",", sep = "")
            cat("\n            control = glm.control(", .print.control(object@GLM$control), ") ),", sep = "", fill = .Options$width)
            
            ## GBM options
            cat("\n")
            cat("\nGBM = list( distribution = '", object@GBM$distribution, "',", sep = "")
            cat("\n            n.trees = ", object@GBM$n.trees, ",", sep = "")
            cat("\n            interaction.depth = ", object@GBM$interaction.depth, ",", sep = "")
            cat("\n            n.minobsinnode = ", object@GBM$n.minobsinnode, ",", sep = "")
            cat("\n            shrinkage = ", object@GBM$shrinkage, ",", sep = "")
            cat("\n            bag.fraction = ", object@GBM$bag.fraction, ",", sep = "")
            cat("\n            train.fraction = ", object@GBM$train.fraction, ",", sep = "")
            cat("\n            cv.folds = ", object@GBM$cv.folds, ",", sep = "")
            cat("\n            keep.data = ", object@GBM$keep.data, ",", sep = "")
            cat("\n            verbose = ", object@GBM$verbose, ",", sep = "")
            #             cat("\n            class.stratify.cv = '", object@GBM$class.stratify.cv, "',", sep="")
            cat("\n            perf.method = '", object@GBM$perf.method, "',", sep = "")
            cat("\n            n.cores = ", ifelse(length(object@GBM$n.cores), object@GBM$n.cores, 'NULL'), "),", sep = "")
            
            ## GAM options
            cat("\n")
            cat("\nGAM = list( algo = '", object@GAM$algo, "',", sep = "")
            cat("\n            type = '", object@GAM$type, "',", sep = "")
            cat("\n            k = ", ifelse(length(object@GAM$k) < 1, 'NULL', object@GAM$k), ",", sep = "")
            cat("\n            interaction.level = ", object@GAM$interaction.level, ",", sep = "")
            cat("\n            myFormula = ", ifelse(length(object@GAM$myFormula) < 1, 'NULL', paste(object@GAM$myFormula[2],
                                                                                                     object@GAM$myFormula[1],
                                                                                                     object@GAM$myFormula[3])), ",", sep = "")
            cat("\n            family = ", object@GAM$family$family, "(link = '", object@GAM$family$link, "'),", sep = "")
            if (object@GAM$algo == 'GAM_mgcv') {
              cat("\n            method = '", object@GAM$method, "', ", sep="")
              cat("\n            optimizer = c('", paste(object@GAM$optimizer, collapse = "','"), "'),", sep = "")
              cat("\n            select = ", object@GAM$select, ",", sep = "")
              cat("\n            knots = ",  ifelse(length(object@GLM$knots) < 1, 'NULL', "'user.defined'"), ",", sep = "")
              cat("\n            paraPen = ",  ifelse(length(object@GLM$paraPen) < 1, 'NULL', "'user.defined'"), ",", sep = "")
            }
            cat("\n            control = list(", .print.control(object@GAM$control), ") ),", sep = "", fill = .Options$width)
            
            ## CTA options
            cat("\n")
            cat("\nCTA = list( method = '", object@CTA$method, "',", sep = "")
            cat("\n            parms = '", object@CTA$parms, "',", sep = "")
            cat("\n            cost = ", ifelse(length(object@CTA$cost) < 1, 'NULL', object@CTA$cost), ",", sep = "")
            cat("\n            control = list(", .print.control(object@CTA$control), ") ),", sep = "", fill = .Options$width)
            
            ## ANN options
            cat("\n")
            cat("\nANN = list( NbCV = ", object@ANN$NbCV, ",", sep = "")
            cat("\n            size = ", ifelse(length(object@ANN$size) < 1, 'NULL', object@ANN$size), ",", sep = "")
            cat("\n            decay = ", ifelse(length(object@ANN$decay) < 1, 'NULL', object@ANN$decay), ",", sep = "")
            cat("\n            rang = ", object@ANN$rang, ",", sep = "")
            cat("\n            maxit = ", object@ANN$maxit, "),", sep = "")
            
            ## SRE options
            cat("\n")
            cat("\nSRE = list( quant = ", object@SRE$quant, "),", sep = "")
            
            ## FDA options
            cat("\n")
            cat("\nFDA = list( method = '", object@FDA$method, "',", sep = "")
            cat("\n            add_args = ", ifelse(length(object@FDA$add_args) < 1, 'NULL'
                                                    , paste("list(", paste(.print.control(object@FDA$add_args), collapse = "")
                                                            , ")", sep = "")), "),", sep = "")
            
            ## MARS options
            cat("\n")
            cat("\nMARS = list( type = '", object@MARS$type, "',", sep = "")
            cat("\n             interaction.level = ", object@MARS$interaction.level, ",", sep = "")
            cat("\n             myFormula = ",  ifelse(length(object@MARS$myFormula) < 1, 'NULL', paste(object@GLM$myFormula[2],
                                                                                                        object@GLM$myFormula[1],
                                                                                                        object@GLM$myFormula[3])), ",", sep = "")
            #             cat("\n             degree = ", object@MARS$degree, ",", sep="")
            cat("\n             nk = ", ifelse(length(object@MARS$nk) < 1, 'NULL', object@MARS$nk), ",", sep = "")
            cat("\n             penalty = ", object@MARS$penalty, ",", sep = "")
            cat("\n             thresh = ", object@MARS$thresh, ",", sep = "")
            cat("\n             nprune = ", ifelse(length(object@MARS$nprune) < 1, 'NULL', object@MARS$nprune), ",", sep = "")
            cat("\n             pmethod = '", object@MARS$pmethod, "'),", sep = "")
            
            ## RF options
            cat("\n")
            cat("\nRF = list( do.classif = ", object@RF$do.classif, ",", sep = "")
            cat("\n           ntree = ", object@RF$ntree, ",", sep = "")
            cat("\n           mtry = '", object@RF$mtry, "',", sep = "")
            cat("\n           nodesize = ", object@RF$nodesize, ",", sep = "")
            cat("\n           maxnodes = ", ifelse(length(object@RF$maxnodes) < 1, 'NULL', object@RF$maxnodes),  "),", sep = "")
            
            ## MAXENT.Phillips options
            cat("\n")
            cat("\nMAXENT.Phillips = list( path_to_maxent.jar = '", object@MAXENT.Phillips$path_to_maxent.jar, "', ", sep="")
            cat("\n               memory_allocated = ", ifelse(length(object@MAXENT.Phillips$memory_allocated) < 1, 'NULL'
                                                               , object@MAXENT.Phillips$memory_allocated), ",", sep = "")
            cat("\n               background_data_dir = ", ifelse(is.character(object@MAXENT.Phillips$background_data_dir), "'", "")
                , object@MAXENT.Phillips$background_data_dir, ifelse(is.character(object@MAXENT.Phillips$background_data_dir), "'", ""), ",", sep = "")
            cat("\n               maximumbackground = ", ifelse(is.character(object@MAXENT.Phillips$maximumbackground), "'", "")
                , object@MAXENT.Phillips$maximumbackground, ifelse(is.character(object@MAXENT.Phillips$maximumbackground), "'", ""), ",", sep = "")
            cat("\n               maximumiterations = ", object@MAXENT.Phillips$maximumiterations, ",", sep = "")
            cat("\n               visible = ", object@MAXENT.Phillips$visible, ",", sep = "")
            cat("\n               linear = ", object@MAXENT.Phillips$linear, ",", sep = "")
            cat("\n               quadratic = ", object@MAXENT.Phillips$quadratic, ",", sep = "")
            cat("\n               product = ", object@MAXENT.Phillips$product, ",", sep = "")
            cat("\n               threshold = ", object@MAXENT.Phillips$threshold, ",", sep = "")
            cat("\n               hinge = ", object@MAXENT.Phillips$hinge, ",", sep = "")
            cat("\n               lq2lqptthreshold = ", object@MAXENT.Phillips$lq2lqptthreshold, ",", sep = "")
            cat("\n               l2lqthreshold = ", object@MAXENT.Phillips$l2lqthreshold, ",", sep = "")
            cat("\n               hingethreshold = ", object@MAXENT.Phillips$hingethreshold, ",", sep = "")
            cat("\n               beta_threshold = ", object@MAXENT.Phillips$beta_threshold, ",", sep = "")
            cat("\n               beta_categorical = ", object@MAXENT.Phillips$beta_categorical, ",", sep = "")
            cat("\n               beta_lqp = ", object@MAXENT.Phillips$beta_lqp, ",", sep = "")
            cat("\n               beta_hinge = ", object@MAXENT.Phillips$beta_hinge, ",", sep = "")
            cat("\n               betamultiplier = ", object@MAXENT.Phillips$betamultiplier, ",", sep = "")
            cat("\n               defaultprevalence = ", object@MAXENT.Phillips$defaultprevalence, "),", sep = "")
            
            ## MAXENT.Phillips.2 options
            cat("\n")
            cat("\n MAXENT.Phillips.2 = list( myFormula = ", .print.formula(object@MAXENT.Phillips.2$myFormula), ",", sep = "")
            cat("\n     regmult = ", object@MAXENT.Phillips.2$regmult, ",", sep = "")
            cat("\n     regfun = <function> )")
            cat("\n)")
            
            # ## MAXENT.Tsuruoka
            # cat("\n")
            # cat("\nMAXENT.Tsuruoka = list( l1_regularizer = ", object@MAXENT.Tsuruoka$l1_regularizer, ",", sep="")
            # cat("\n                        l2_regularizer = ", object@MAXENT.Tsuruoka$l2_regularizer, ",", sep="")
            # cat("\n                        use_sgd = ", object@MAXENT.Tsuruoka$use_sgd, ",", sep="")
            # cat("\n                        set_heldout = ", object@MAXENT.Tsuruoka$set_heldout, ",", sep="")
            # cat("\n                        verbose = ", object@MAXENT.Tsuruoka$verbose, ")", sep="")
            
            .bmCat()
          }
)


###################################################################################################
## USED IN BIOMOD_Modeling FUNCTION
###################################################################################################

setGeneric(".Models.prepare.data", def = function(data, ...) { standardGeneric(".Models.prepare.data") })

setMethod('.Models.prepare.data', signature('BIOMOD.formated.data'),
          function(data, NbRunEval, DataSplit, Yweights = NULL, Prevalence = NULL
                   , do.full.models = TRUE, DataSplitTable = NULL)
          {
            list.out <- list()
            name <- paste0(data@sp.name, '_AllData')
            xy <- data@coord
            dataBM <- bind_cols(tibble(!!data@sp.name := data@data.species), data@data.env.var)
            
            ## dealing with evaluation data
            if (data@has.data.eval) {
              evalDataBM <- data.frame(cbind(data@eval.data.species, data@eval.data.env.var[, , drop = FALSE]))
              colnames(evalDataBM)[1] <- data@sp.name
              eval.xy <- data@eval.coord
            } else {
              evalDataBM <- eval.xy <- NULL
            }
            
            ## Calib/Valid lines
            if (!is.null(DataSplitTable)) {
              calibLines <- DataSplitTable
              colnames(calibLines) <- paste('_RUN', 1:ncol(calibLines), sep = '')
            } else {
              if (NbRunEval == 0) { # take all available data
                calibLines <- matrix(rep(TRUE, length(data@data.species)), ncol = 1)
                colnames(calibLines) <- '_Full'
              } else {
                calibLines <- .sample_mat(data.sp = data@data.species,
                                          dataSplit = DataSplit,
                                          nbRun = NbRunEval,
                                          data.env = data@data.env.var)
                if (do.full.models) {
                  calibLines <- cbind(calibLines, rep(TRUE, length(data@data.species)))
                  colnames(calibLines)[NbRunEval + 1] <- '_Full'
                }
              }
            }
            ## force calib.lines object to be 3D array
            if (length(dim(calibLines)) < 3) {
              dn_tmp <- dimnames(calibLines) ## keep track of dimnames
              dim(calibLines) <- c(dim(calibLines), 1)
              dimnames(calibLines) <- list(dn_tmp[[1]], dn_tmp[[2]], "_AllData")
            }
            
            if (is.null(Yweights)) { # 1 for all points
              if (!is.null(Prevalence)) {
                cat("\n\t> Automatic weights creation to rise a", Prevalence, "prevalence")
                Yweights <- .automatic_weights_creation(data@data.species , prev = Prevalence)
              } else {
                cat("\n\t> No weights : all observations will have the same weight")
                Yweights <- rep(1, length(data@data.species))
              }
            }
            
            list.out[[name]] <- list(name = name,
                                     xy = xy,
                                     dataBM = dataBM, 
                                     calibLines = calibLines,
                                     Yweights = Yweights,
                                     evalDataBM = evalDataBM,
                                     eval.xy = eval.xy)
            return(list.out)
          }
)

setMethod('.Models.prepare.data', signature(data='BIOMOD.formated.data.PA'),
          function(data, NbRunEval, DataSplit, Yweights = NULL, Prevalence = NULL
                   , do.full.models = TRUE, DataSplitTable = NULL)
          {
            list.out <- list()
            formal_weights <- Yweights
            for (pa in 1:ncol(data@PA))
            {
              Yweights <- formal_weights
              name <- paste0(data@sp.name, "_", colnames(data@PA)[pa])
              xy <- data@coord[data@PA[, pa], ]
              resp <- data@data.species[data@PA[, pa]] # response variable (with pseudo absences selected)
              resp[is.na(resp)] <- 0
              dataBM <- data.frame(cbind(resp, data@data.env.var[data@PA[, pa], , drop = FALSE]))
              colnames(dataBM)[1] <- data@sp.name
              
              ## Calib/Valid lines
              if (!is.null(DataSplitTable))
              {
                if (length(dim(DataSplitTable)) == 2) {
                  calibLines <- DataSplitTable
                } else {
                  calibLines <- asub(DataSplitTable, pa, 3, drop = TRUE)
                }
                colnames(calibLines) <- paste0('_RUN', 1:ncol(calibLines))
                calibLines[which(!data@PA[, pa]), ] <- NA
              } else {
                if (NbRunEval == 0) { # take all available data
                  calibLines <- matrix(NA, nrow = length(data@data.species), ncol = 1)
                  calibLines[data@PA[, pa], 1] <- TRUE
                  colnames(calibLines) <- '_Full'
                } else {
                  calibLines <- matrix(NA, nrow = length(data@data.species), ncol = NbRunEval)
                  sampled.mat <- .sample_mat(
                    data.sp = data@data.species[data@PA[, pa]],
                    dataSplit = DataSplit,
                    nbRun = NbRunEval,
                    data.env = data@data.env.var[data@PA[, pa], , drop = FALSE]
                  )
                  calibLines[data@PA[, pa], ] <- sampled.mat
                  colnames(calibLines) <- colnames(sampled.mat)
                  if (do.full.models) {
                    calibLines <- cbind(calibLines, rep(NA, length(data@data.species)))
                    calibLines[data@PA[, pa], NbRunEval + 1] <- TRUE
                    colnames(calibLines)[NbRunEval + 1] <- '_Full'
                  }
                }
              }
              
              ## force calib.lines object to be 3D array
              if (length(dim(calibLines)) < 3) {
                dn_tmp <- dimnames(calibLines) ## keep track of dimnames
                dim(calibLines) <- c(dim(calibLines), 1)
                dimnames(calibLines) <- list(dn_tmp[[1]], dn_tmp[[2]], paste0("_PA", pa))
              }
              
              # dealing with evaluation data
              if (data@has.data.eval) {
                evalDataBM <- data.frame(cbind(data@eval.data.species, data@eval.data.env.var))
                colnames(evalDataBM)[1] <- data@sp.name
                eval.xy <- data@eval.coord
              } else {
                evalDataBM <- eval.xy <- NULL
              }
              
              if (is.null(Yweights)) { # prevalence of 0.5... may be parametrize
                if (is.null(Prevalence)) { Prevalence <- 0.5 }
                cat("\n\t\t\t! Weights where automatically defined for", name, "to rise a", Prevalence, "prevalence !")
                Yweights <- rep(NA, length(data@data.species))
                Yweights[data@PA[, pa]] <- .automatic_weights_creation(as.numeric(dataBM[, 1]) , prev = Prevalence)
              } else { # remove useless weights
                Yweights[!data@PA[, pa]] <- NA
              }
              
              list.out[[name]] <- list(name = name,
                                       xy = xy, 
                                       dataBM = dataBM,
                                       calibLines = calibLines,
                                       Yweights = Yweights,
                                       evalDataBM = evalDataBM,
                                       eval.xy = eval.xy)
            }
            return(list.out)
          }
)

