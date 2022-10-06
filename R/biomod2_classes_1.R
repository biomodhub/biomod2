
## --------------------------------------------------------------------------- #
## 1. BIOMOD.formated.data ---------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.formated.data
##' @aliases BIOMOD.formated.data-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_FormatingData()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_FormatingData}}, and used by 
##' \code{\link{BIOMOD_Tuning}}, \code{\link{BIOMOD_CrossValidation}} and 
##' \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @param dir.name a \code{character} corresponding to the modeling folder
##' @param sp.name a \code{character} corresponding to the species name
##' 
##' @param sp a \code{vector}, \code{\link[sp]{SpatialPoints}} (\emph{if presence-only}) or 
##' \code{\link[sp]{SpatialPointsDataFrame}} object containing binary data (\code{0} : absence, 
##' \code{1} : presence, \code{NA} : indeterminate) for a single species that will be used to 
##' build the species distribution model(s)
##' @param env a \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} 
##' or \code{\link[raster:stack]{RasterStack}} object containing the explanatory variables 
##' (in columns or layers) that will be used to build the species distribution model(s)
##' @param xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be used to build the 
##' species distribution model(s)
##' @param eval.sp (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector}, \code{\link[sp]{SpatialPoints}} (\emph{if presence-only}) or 
##' \code{\link[sp]{SpatialPointsDataFrame}} object containing binary data (\code{0} : absence, 
##' \code{1} : presence, \code{NA} : indeterminate) for a single species that will be used to 
##' evaluate the species distribution model(s) with independent data
##' @param eval.env (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} or 
##' \code{\link[raster:stack]{RasterStack}} object containing the explanatory variables (in 
##' columns or layers) that will be used to evaluate the species distribution model(s) with 
##' independent data
##' @param eval.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be used to evaluate 
##' the species distribution model(s) with independent data
##' 
##' @param na.rm (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether points having one or several missing values for 
##' explanatory variables should be removed from the analysis or not
##' 
##' @param data.mask a \code{\link[raster:stack]{RasterStack}} object containing the mask of the 
##' studied area
##' 
##' @param coord a 2-columns \code{data.frame} containing \code{X} and \code{Y} coordinates for plot
##' @param col a \code{vector} containing colors for plot (default : \code{c('green', 'red', 
##' 'orange', 'grey')})
##' @param x a \code{\link{BIOMOD.formated.data.PA}} object
##' @param object a \code{\link{BIOMOD.formated.data.PA}} object
##' 
##' 
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot coord a 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates
##' @slot data.species a \code{vector} containing the species observations (\code{0}, \code{1} or 
##' \code{NA})
##' @slot data.env.var a \code{data.frame} containing explanatory variables
##' @slot data.mask a \code{\link[raster:stack]{RasterStack}} object containing the mask of the 
##' studied area
##' @slot has.data.eval a \code{logical} value defining whether evaluation data is given
##' @slot eval.coord (\emph{optional, default} \code{NULL}) \cr 
##' A 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates for evaluation data
##' @slot eval.data.species (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing the species observations (\code{0}, \code{1} or \code{NA}) for 
##' evaluation data
##' @slot eval.data.env.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{data.frame} containing explanatory variables for evaluation data
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_Tuning}}, 
##' \code{\link{BIOMOD_CrossValidation}}, \code{\link{BIOMOD_Modeling}}, 
##' \code{\link{bm_RunModelsLoop}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.formated.data")
##' 
##' ## ----------------------------------------------------------------------- #
##' 
##' # Load species occurrences (6 species available)
##' myFile <- system.file('external/species/mammals_table.csv', package = 'biomod2')
##' DataSpecies <- read.csv(myFile, row.names = 1)
##' head(DataSpecies)
##' 
##' # Select the name of the studied species
##' myRespName <- 'GuloGulo'
##' 
##' # Get corresponding presence/absence data
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' # Get corresponding XY coordinates
##' myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
##' 
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles <- paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
##' myExpl <- raster::stack(system.file(myFiles, package = 'biomod2'))
##' 
##' \dontshow{
##' myExtent <- raster::extent(0,30,45,70)
##' myExpl <- raster::stack(raster::crop(myExpl, myExtent))
##' }
##' 
##' ## ----------------------------------------------------------------------- #
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' myBiomodData
##' plot(myBiomodData)
##' 
##' 
##' 
NULL

##' @name BIOMOD.formated.data-class
##' @rdname BIOMOD.formated.data
##' @importFrom raster stack nlayers addLayer is.factor subset extract cellStats cellFromXY
##' @export
##' 

# 1.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.formated.data",
         representation(dir.name = 'character',
                        sp.name = 'character',
                        coord = "data.frame",
                        data.species = "numeric",
                        data.env.var = "data.frame",
                        data.mask = "RasterStack",
                        has.data.eval = "logical",
                        eval.coord = "data.frame",
                        eval.data.species = "numeric",
                        eval.data.env.var = "data.frame"),
         validity = function(object){ return(TRUE) })


# 1.2 Constructors -------------------------------------------------------------
setGeneric("BIOMOD.formated.data", def = function(sp, env, ...) { standardGeneric("BIOMOD.formated.data") })

## BIOMOD.formated.data(sp = numeric, env = data.frame) -----------------------------------
##' 
##' @rdname BIOMOD.formated.data
##' @export
##' 

setMethod('BIOMOD.formated.data', signature(sp = 'numeric', env = 'data.frame'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE, data.mask = NULL)
          {
            if (is.null(data.mask)) { data.mask <- stack() }
            
            if (is.null(eval.sp)) { ## NO EVALUATION DATA
              BFD <- new(
                'BIOMOD.formated.data',
                coord = xy,
                data.species = sp,
                data.env.var = env,
                dir.name = dir.name,
                sp.name = sp.name,
                data.mask = data.mask,
                has.data.eval = FALSE
              )
            } else { ## EVALUATION DATA
              BFDeval <- BIOMOD.formated.data(
                sp = eval.sp,
                env = eval.env,
                xy = eval.xy,
                dir.name = dir.name,
                sp.name = sp.name
              )
              
              if (nlayers(BFDeval@data.mask) == 1) {
                if (nlayers(data.mask) == 1) {
                  data.mask.tmp <- try(addLayer(data.mask, BFDeval@data.mask))
                  if (!inherits(data.mask.tmp, "try-error")) {
                    data.mask <- data.mask.tmp
                    names(data.mask) <- c("calibration", "validation")
                  }
                } else if (nlayers(data.mask) == 0) {
                  # in this case the data.mask from calibration is added later
                  data.mask <- BFDeval@data.mask
                  names(data.mask) <- c("validation")
                }
              }
              
              BFD <- new(
                'BIOMOD.formated.data',
                coord = xy,
                data.species = sp,
                data.env.var = env,
                dir.name = dir.name,
                sp.name = sp.name,
                data.mask = data.mask,
                has.data.eval = TRUE,
                eval.coord = BFDeval@coord,
                eval.data.species = BFDeval@data.species,
                eval.data.env.var = BFDeval@data.env.var
              )
              
              rm('BFDeval')
            }
            
            ## REMOVE NA IF ANY
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

## BIOMOD.formated.data(sp = data.frame) -----------------------------------
##' 
##' @rdname BIOMOD.formated.data
##' @export
##' 

setMethod('BIOMOD.formated.data', signature(sp = 'data.frame'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE)
          {
            if (ncol(sp) > 1) { stop("Invalid response variable") }
            sp <- as.numeric(unlist(sp))
            BFD <- BIOMOD.formated.data(sp, env, xy, dir.name, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm)
            return(BFD)
          }
)

## BIOMOD.formated.data(sp = numeric, env = matrix) ----------------------------
##' 
##' @rdname BIOMOD.formated.data
##' @export
##' 

setMethod('BIOMOD.formated.data', signature(sp = 'numeric', env = 'matrix'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE)
          {
            env <- as.data.frame(env)
            BFD <- BIOMOD.formated.data(sp, env, xy, dir.name, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm)
            return(BFD)
          }
)

## BIOMOD.formated.data(sp = numeric, env = RasterStack) -----------------------
##' 
##' @rdname BIOMOD.formated.data
##' @export
##' 

setMethod('BIOMOD.formated.data', signature(sp = 'numeric', env = 'RasterStack'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE)
          {
            categorical_var <- names(env)[is.factor(env)]
            
            ## Keep same env variable for eval than calib (+ check for factor)
            if (!is.null(eval.sp) && is.null(eval.env)) {
              eval.env <- as.data.frame(extract(env, eval.xy))
              if (length(categorical_var)) {
                for (cat_var in categorical_var) {
                  eval.env[, cat_var] <- as.factor(eval.env[, cat_var])
                }
              }
            }
            
            if (is.null(xy)) { xy <- as.data.frame(coordinates(env)) }
            
            ## Prepare mask of studied area
            data.mask = reclassify(subset(env, 1, drop = TRUE), c(-Inf, Inf, -1))
            data.mask[cellFromXY(data.mask, xy[which(sp == 1), ])] <- 1
            data.mask[cellFromXY(data.mask, xy[which(sp == 0), ])] <- 0
            data.mask <- stack(data.mask)
            names(data.mask) <- sp.name
            
            ## Keep same env variable for eval than calib (+ check for factor)
            env <- as.data.frame(extract(env, xy, factors = TRUE))
            if (length(categorical_var)) {
              for (cat_var in categorical_var) {
                env[, cat_var] <- as.factor(env[, cat_var])
              }
            }
            
            BFD <- BIOMOD.formated.data(sp, env, xy, dir.name, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm, data.mask = data.mask)
            return(BFD)
          }
)


# 1.3 Other Functions -----------------------------------------------------------------------------

### plot.BIOMOD.formated.data  --------------------------------------------------
##' 
##' @rdname BIOMOD.formated.data
##' @export
##' 

setMethod('plot', signature(x = 'BIOMOD.formated.data', y = "missing"),
          function(x, coord = NULL, col = NULL)
          {
            if (nlayers(x@data.mask) > 0)
            {
              requireNamespace("rasterVis")
              
              ## check if there is some undefined areas to prevent from strange plotting issues
              if (min(cellStats(x@data.mask, min)) == -1) { # there is undefined area
                my.at <- seq(-1.5, 1.5, by = 1) ## breaks of color key
                my.labs.at <- seq(-1, 1, by = 1) ## labels placed vertically
                my.lab <- c("undefined", "absences", "presences") ## labels
                my.col.regions = c("lightgrey", "red4", "green4") ## colors
                my.cuts <- 2 ## cuts
              } else { # no undefined area.. remove it from plot
                my.at <- seq(-0.5, 1.5, by = 1) ## breaks of color key
                my.labs.at <- seq(0, 1, by = 1) ## labels placed vertically
                my.lab <- c("absences", "presences") ## labels
                my.col.regions = c("red4", "green4") ## colors
                my.cuts <- 1 ## cuts
              }
              
              ## PLOT
              rasterVis::levelplot(
                x@data.mask,
                at = my.at,
                cuts = my.cuts,
                margin = TRUE,
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
              
              
              ## PLOT
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

### show.BIOMOD.formated.data  --------------------------------------------------
##' 
##' @rdname BIOMOD.formated.data
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('BIOMOD.formated.data'),
          function(object)
          {
            .bm_cat("BIOMOD.formated.data")
            cat("\ndir.name = ", object@dir.name, fill = .Options$width)
            cat("\nsp.name = ", object@sp.name, fill = .Options$width)
            cat("\n\t",
                sum(object@data.species, na.rm = TRUE),
                'presences, ',
                sum(object@data.species == 0, na.rm = TRUE),
                'true absences and ',
                sum(is.na(object@data.species), na.rm = TRUE),
                'undefined points in dataset',
                fill = .Options$width)
            cat("\n\n\t",
                ncol(object@data.env.var),
                'explanatory variables\n',
                fill = .Options$width)
            print(summary(object@data.env.var))
            
            if (object@has.data.eval) {
              cat("\n\nEvaluation data :", fill = .Options$width)
              cat("\n\t",
                  sum(object@eval.data.species, na.rm = TRUE),
                  'presences, ',
                  sum(object@eval.data.species == 0, na.rm = TRUE),
                  'true absences and ',
                  sum(is.na(object@eval.data.species), na.rm = TRUE),
                  'undefined points in dataset',
                  fill = .Options$width)
              cat("\n\n", fill = .Options$width)
              print(summary(object@eval.data.env.var))
            }
            .bm_cat()
          }
)



## --------------------------------------------------------------------------- #
## 2. BIOMOD.formated.data.PA ------------------------------------------------
## this class inherits from BIOMOD.formated.data and have one more slot 'PA', 
## giving PA selected
## --------------------------------------------------------------------------- #

##' @name BIOMOD.formated.data.PA
##' @aliases BIOMOD.formated.data.PA-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_FormatingData()} output object class (with pseudo-absences)
##' 
##' @description Class returned by \code{\link{BIOMOD_FormatingData}}, and used by 
##' \code{\link{BIOMOD_Tuning}}, \code{\link{BIOMOD_CrossValidation}} and 
##' \code{\link{BIOMOD_Modeling}}
##' 
##' @param dir.name a \code{character} corresponding to the modeling folder
##' @param sp.name a \code{character} corresponding to the species name
##' 
##' @param sp a \code{vector}, \code{\link[sp]{SpatialPoints}} (\emph{if presence-only}) or 
##' \code{\link[sp]{SpatialPointsDataFrame}} object containing binary data (\code{0} : absence, 
##' \code{1} : presence, \code{NA} : indeterminate) for a single species that will be used to 
##' build the species distribution model(s)
##' @param env a \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} 
##' or \code{\link[raster:stack]{RasterStack}} object containing the explanatory variables 
##' (in columns or layers) that will be used to build the species distribution model(s)
##' @param xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be used to build the 
##' species distribution model(s)
##' @param eval.sp (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector}, \code{\link[sp]{SpatialPoints}} (\emph{if presence-only}) or 
##' \code{\link[sp]{SpatialPointsDataFrame}} object containing binary data (\code{0} : absence, 
##' \code{1} : presence, \code{NA} : indeterminate) for a single species that will be used to 
##' evaluate the species distribution model(s) with independent data
##' @param eval.env (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} or 
##' \code{\link[raster:stack]{RasterStack}} object containing the explanatory variables (in 
##' columns or layers) that will be used to evaluate the species distribution model(s) with 
##' independent data
##' @param eval.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be used to evaluate 
##' the species distribution model(s) with independent data
##' 
##' @param na.rm (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether points having one or several missing values for 
##' explanatory variables should be removed from the analysis or not
##' @param na.rm (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether points having one or several missing values for 
##' explanatory variables should be removed from the analysis or not
##' 
##' @param PA.nb.rep (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection, an \code{integer} corresponding to the number of sets 
##' (repetitions) of pseudo-absence points that will be drawn
##' @param PA.strategy (\emph{optional, default} \code{NULL}) \cr 
##' If pseudo-absence selection, a \code{character} defining the strategy that will be used to 
##' select the pseudo-absence points. Must be \code{random}, \code{sre}, \code{disk} or 
##' \code{user.defined} (see Details)
##' @param PA.nb.absences (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection, and \code{PA.strategy = 'random'} or \code{PA.strategy = 'sre'} 
##' or \code{PA.strategy = 'disk'}, an \code{integer} corresponding to the number of pseudo-absence 
##' points that will be selected for each pseudo-absence repetition (true absences included)
##' @param PA.sre.quant (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'sre'}, a \code{numeric} between \code{0} 
##' and \code{0.5} defining the half-quantile used to make the \code{sre} pseudo-absence selection 
##' (see Details)
##' @param PA.dist.min (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'disk'}, a \code{numeric} defining the 
##' minimal distance to presence points used to make the \code{disk} pseudo-absence selection 
##' (in meters, see Details)
##' @param PA.dist.max (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'disk'}, a \code{numeric} defining the 
##' maximal distance to presence points used to make the \code{disk} pseudo-absence selection 
##' (in meters, see Details)
##' @param PA.user.table (\emph{optional, default} \code{NULL}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'user.defined'}, a \code{matrix} or 
##' \code{data.frame} with as many rows as \code{resp.var} values, as many columns as 
##' \code{PA.nb.rep}, and containing \code{TRUE} or \code{FALSE} values defining which points 
##' will be used to build the species distribution model(s) for each repetition (see Details)
##' 
##' @param coord a 2-columns \code{data.frame} containing \code{X} and \code{Y} coordinates for plot
##' @param col a \code{vector} containing colors for plot (default : \code{c('green', 'red', 
##' 'orange', 'grey')})
##' @param x a \code{\link{BIOMOD.formated.data.PA}} object
##' @param object a \code{\link{BIOMOD.formated.data.PA}} object
##' 
##' 
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot coord a 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates
##' @slot data.species a \code{vector} containing the species observations (\code{0}, \code{1} or 
##' \code{NA})
##' @slot data.env.var a \code{data.frame} containing explanatory variables
##' @slot data.mask a \code{\link[raster:stack]{RasterStack}} object containing the mask of the 
##' studied area
##' @slot has.data.eval a \code{logical} value defining whether evaluation data is given
##' @slot eval.coord (\emph{optional, default} \code{NULL}) \cr 
##' A 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates for evaluation data
##' @slot eval.data.species (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing the species observations (\code{0}, \code{1} or \code{NA}) for 
##' evaluation data
##' @slot eval.data.env.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{data.frame} containing explanatory variables for evaluation data
##' @slot PA.strategy a \code{character} corresponding to the pseudo-absence selection strategy
##' @slot PA.table a \code{data.frame} containing the corresponding table of selected 
##' pseudo-absences (indicated by \code{TRUE} or \code{FALSE}) from the \code{pa.tab} list 
##' element returned by the \code{\link{bm_PseudoAbsences}} function
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{bm_PseudoAbsences}}, 
##' \code{\link{BIOMOD_Tuning}}, \code{\link{BIOMOD_CrossValidation}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{bm_RunModelsLoop}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.formated.data.PA")
##' 
##' ## ----------------------------------------------------------------------- #
##' 
##' # Load species occurrences (6 species available)
##' myFile <- system.file('external/species/mammals_table.csv', package = 'biomod2')
##' DataSpecies <- read.csv(myFile, row.names = 1)
##' head(DataSpecies)
##' 
##' # Select the name of the studied species
##' myRespName <- 'GuloGulo'
##' 
##' # Get corresponding presence/absence data
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' # Get corresponding XY coordinates
##' myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
##' 
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles <- paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
##' myExpl <- raster::stack(system.file(myFiles, package = 'biomod2'))
##' 
##' \dontshow{
##' myExtent <- raster::extent(0,30,45,70)
##' myExpl <- raster::stack(raster::crop(myExpl, myExtent))
##' }
##' 
##' ## ----------------------------------------------------------------------- #
##' # Format Data with pseudo-absences : random method
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName,
##'                                      PA.nb.rep = 0,
##'                                      PA.strategy = 'random',
##'                                      PA.nb.absences = 1000)
##' myBiomodData
##' plot(myBiomodData)
##' 
##' 
NULL

##' @name BIOMOD.formated.data.PA-class
##' @rdname BIOMOD.formated.data.PA
##' @importFrom raster stack nlayers addLayer is.factor subset cellFromXY cellStats
## @importFrom rasterVis levelplot
##' 
##' @export
##' 

# 2.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.formated.data.PA",
         contains = "BIOMOD.formated.data",
         representation(PA.strategy = 'character', PA.table = 'data.frame'),
         validity = function(object) { return(TRUE) })


# 2.2 Constructors -------------------------------------------------------------
setGeneric("BIOMOD.formated.data.PA", 
           def = function(sp, env, ...) { 
             standardGeneric("BIOMOD.formated.data.PA") 
           })

### BIOMOD.formated.data.PA(sp = numeric, env = data.frame) --------------------
##' 
##' @rdname BIOMOD.formated.data.PA
##' @export
##' 

setMethod('BIOMOD.formated.data.PA', signature(sp = 'numeric', env = 'data.frame'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , PA.nb.rep = 1, PA.strategy = 'random', PA.nb.absences = NULL
                   , PA.dist.min = 0, PA.dist.max = NULL
                   , PA.sre.quant = 0.025, PA.user.table = NULL
                   , na.rm = TRUE) {
            .BIOMOD.formated.data.PA(sp, env, xy, dir.name, sp.name
                                     , eval.sp, eval.env, eval.xy
                                     , PA.nb.rep, PA.strategy, PA.nb.absences
                                     , PA.dist.min, PA.dist.max
                                     , PA.sre.quant, PA.user.table
                                     , na.rm)
          })

### BIOMOD.formated.data.PA(sp = numeric, env = RasterStack) -------------------
##' 
##' @rdname BIOMOD.formated.data.PA
##' @export
##' 

setMethod('BIOMOD.formated.data.PA', signature(sp = 'numeric', env = 'RasterStack'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , PA.nb.rep = 1, PA.strategy = 'random', PA.nb.absences = NULL
                   , PA.dist.min = 0, PA.dist.max = NULL
                   , PA.sre.quant = 0.025, PA.user.table = NULL
                   , na.rm = TRUE) {
            .BIOMOD.formated.data.PA(sp, env, xy, dir.name, sp.name
                                     , eval.sp, eval.env, eval.xy
                                     , PA.nb.rep, PA.strategy, PA.nb.absences
                                     , PA.dist.min, PA.dist.max
                                     , PA.sre.quant, PA.user.table
                                     , na.rm)
          })

### .BIOMOD.formated.data.PA ---------------------------------------------------
.BIOMOD.formated.data.PA <-  function(sp, env, xy, dir.name, sp.name
                                      , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                                      , PA.nb.rep = 1, PA.strategy = 'random', PA.nb.absences = NULL
                                      , PA.dist.min = 0, PA.dist.max = NULL
                                      , PA.sre.quant = 0.025, PA.user.table = NULL
                                      , na.rm = TRUE)
{
  
  categorical_var <- NULL
  if (inherits(env, 'Raster')) { categorical_var <- names(env)[is.factor(env)] }
  
  ## Keep same env variable for eval than calib (+ check for factor)
  if (!is.null(eval.sp) && is.null(eval.env)) {
    if (inherits(env, 'Raster')) {
      eval.env <- as.data.frame(extract(env, eval.xy))
      if (length(categorical_var)) {
        for (cat_var in categorical_var) {
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
  
  pa.data.tmp <- bm_PseudoAbsences(resp.var = sp,
                                   expl.var = env,
                                   nb.rep = PA.nb.rep,
                                   strategy = PA.strategy,
                                   nb.absences = PA.nb.absences,
                                   sre.quant = PA.sre.quant,
                                   dist.min = PA.dist.min,
                                   dist.max = PA.dist.max,
                                   user.table = PA.user.table)
  
  if (!is.null(pa.data.tmp)) {
    ## Keep same env variable for eval than calib (+ check for factor)
    if (length(categorical_var)) {
      for (cat_var in categorical_var) {
        pa.data.tmp$env[, cat_var] <- as.factor(pa.data.tmp$env[, cat_var])
      }
    }
    
    ## REMOVE NA IF ANY
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
                                dir.name = dir.name,
                                sp.name = sp.name, 
                                eval.sp = eval.sp,
                                eval.env = eval.env,
                                eval.xy = eval.xy,
                                na.rm = na.rm)
    
    if (inherits(env,'Raster')) {
      ## create data.mask for ploting
      data.mask.tmp <- reclassify(subset(env, 1), c(-Inf, Inf, -1))
      data.mask <- stack(data.mask.tmp)
      xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp == 1), , drop = FALSE]
      xy_abs <- pa.data.tmp$xy[which(pa.data.tmp$sp == 0), , drop = FALSE]
      if (nrow(xy_pres)) { data.mask[cellFromXY(data.mask.tmp, xy_pres)] <- 1 }
      if (nrow(xy_abs)) { data.mask[cellFromXY(data.mask.tmp, xy_abs)] <- 0 }
      names(data.mask) <- "input_data"
      
      ## add eval data
      if (BFD@has.data.eval) {
        if (nlayers(BFD@data.mask) == 1 && names(BFD@data.mask) == "validation") {
          data.mask.eval.tmp <- try(addLayer(data.mask, BFD@data.mask))
          if (!inherits(data.mask.eval.tmp, "try-error")) {
            data.mask <- data.mask.eval.tmp
            names(data.mask) <- c("input_data", "validation")
          }
        }
      } 
      
      ## add pa data
      data.mask.names.tmp <- names(data.mask)
      for (pa in 1:ncol(as.data.frame(pa.data.tmp$pa.tab))) {
        data.mask.tmp2 <- data.mask.tmp
        
        ind.pa <- as.data.frame(pa.data.tmp$pa.tab)[, pa]
        xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp == 1 & ind.pa == TRUE), , drop = FALSE]
        xy_abs <- pa.data.tmp$xy[which((pa.data.tmp$sp != 1 | is.na(pa.data.tmp$sp)) & ind.pa == TRUE), , drop = FALSE]
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
      
      names(data.mask) <- c(data.mask.names.tmp,
                            colnames(as.data.frame(pa.data.tmp$pa.tab)))
      
    } else {  data.mask <- stack() }
    
    BFDP <- new('BIOMOD.formated.data.PA',
                dir.name = BFD@dir.name,
                sp.name = BFD@sp.name,
                coord = BFD@coord,
                data.env.var = BFD@data.env.var,
                data.species = BFD@data.species,
                data.mask = data.mask,
                has.data.eval = BFD@has.data.eval,
                eval.coord = BFD@eval.coord,
                eval.data.species = BFD@eval.data.species,
                eval.data.env.var = BFD@eval.data.env.var,
                PA.strategy = PA.strategy,
                PA.table = as.data.frame(pa.data.tmp$pa.tab))
    
    rm(list = 'BFD')
  } else {
    cat("\n   ! PA selection not done", fill = .Options$width)
    
    BFDP <- BIOMOD.formated.data(sp = as.vector(sp@data),
                                 env = env,
                                 xy = xy,
                                 dir.name = dir.name,
                                 sp.name = sp.name, 
                                 eval.sp = eval.sp,
                                 eval.env = eval.env,
                                 eval.xy = eval.xy)
  }
  rm(list = "pa.data.tmp" )
  return(BFDP)
}


# 2.3 other functions ----------------------------------------------------------

### plot.BIOMOD.formated.data.PA -----------------------------------------------
##' 
##' @rdname BIOMOD.formated.data.PA
##' @export
##' 

setMethod('plot', signature(x = 'BIOMOD.formated.data.PA', y = "missing"),
          function(x, coord = NULL, col = NULL)
          {
            if (nlayers(x@data.mask) > 0)
            {
              requireNamespace("rasterVis")
              
              ## check if there is some undefined areas to prevent from strange plotting issues
              if (min(cellStats(x@data.mask, min)) == -1) { # there is undefined area
                my.at <- seq(-1.5, 1.5, by = 1) ## breaks of color key
                my.labs.at <- seq(-1, 1, by = 1) ## labels placed vertically
                my.lab <- c("undefined", "absences", "presences") ## labels
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
              rasterVis::levelplot(
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
              par(mfrow = c(.clever_cut(ncol(x@PA.table) + 1)))
              
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
              for (i in 1:ncol(x@PA.table))
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
                  x = x@coord[(x@data.species == 1) & x@PA.table[, i], 1],
                  y = x@coord[(x@data.species == 1) & x@PA.table[, i], 2],
                  col = col[1],
                  pch = 18
                )
                # true absences
                points(
                  x = x@coord[(x@data.species == 0) & x@PA.table[, i], 1],
                  y = x@coord[(x@data.species == 0) & x@PA.table[, i], 2],
                  col = col[2],
                  pch = 18
                )
                # PA
                points(
                  x = x@coord[is.na(x@data.species) & x@PA.table[, i], 1],
                  y = x@coord[is.na(x@data.species) & x@PA.table[, i], 2],
                  col = col[3],
                  pch = 18
                )
              }
            }
          }
)

### show.BIOMOD.formated.data.PA -----------------------------------------------
##' 
##' @rdname BIOMOD.formated.data.PA
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('BIOMOD.formated.data.PA'),
          function(object)
          {
            .bm_cat("BIOMOD.formated.data.PA")
            cat("\ndir.name = ", object@dir.name, fill = .Options$width)
            cat("\nsp.name = ", object@sp.name, fill = .Options$width)
            cat(
              "\n\t",
              sum(object@data.species, na.rm = TRUE),
              'presences, ',
              sum(object@data.species == 0, na.rm = TRUE),
              'true absences and ',
              sum(is.na(object@data.species), na.rm = TRUE),
              'undefined points in dataset',
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
                'undefined points in dataset',
                fill = .Options$width
              )
              cat("\n\n", fill = .Options$width)
              print(summary(object@eval.data.env.var))
            }
            
            cat(
              "\n\n",
              ncol(object@PA.table),
              'Pseudo Absences dataset available (',
              colnames(object@PA.table),
              ") with ",
              sum(object@PA.table[, 1], na.rm = TRUE) - sum(object@data.species, na.rm = TRUE),
              'absences in each (true abs + pseudo abs)',
              fill = .Options$width
            )
            .bm_cat()
          }
)



## --------------------------------------------------------------------------- #
## 3. BIOMOD.models.options --------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.models.options
##' @aliases BIOMOD.models.options-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_ModelingOptions()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_ModelingOptions}}, and used by 
##' \code{\link{BIOMOD_Tuning}} and \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @slot GLM a \code{list} containing GLM options
##' @slot GBM a \code{list} containing GBM options
##' @slot GAM a \code{list} containing GAM options
##' @slot CTA a \code{list} containing CTA options
##' @slot ANN a \code{list} containing ANN options
##' @slot SRE a \code{list} containing SRE options
##' @slot FDA a \code{list} containing FDA options
##' @slot MARS a \code{list} containing MARS options
##' @slot RF a \code{list} containing RF options
##' @slot MAXENT.Phillips a \code{list} containing MAXENT.Phillips options
##' @slot MAXENT.Phillips.2 a \code{list} containing MAXENT.Phillips options
##' 
##' @param object a \code{\link{BIOMOD.models.options}} object
##' 
##' 
##' @seealso \code{\link{BIOMOD_ModelingOptions}}, \code{\link{BIOMOD_Tuning}}, 
##' \code{\link{BIOMOD_Modeling}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.models.options")
##' 
##' ## ----------------------------------------------------------------------- #
##' ## default BIOMOD.models.options object
##' myBiomodOptions <- BIOMOD_ModelingOptions()
##'
##' ## print the object
##' myBiomodOptions
##'
##' 
NULL

##' @name BIOMOD.models.options-class
##' @rdname BIOMOD.models.options
##' @export
##' 

setClass("BIOMOD.models.options",
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
                             sampsize = NULL,
                             nodesize = 5,
                             maxnodes = NULL),
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
         validity = function(object) {
           test <- TRUE
           
           ## GLM ##
           test <- .fun_testIfIn(test, "GLM$type", object@GLM$type, c('simple', 'quadratic', 'polynomial', 'user.defined'))
           test <- .fun_testIfPosInt(test, "GLM$interaction.level", object@GLM$interaction.level)
           if (!is.null(object@GLM$myFormula)) if (!inherits(object@GLM$myFormula, "formula")) { cat("\nGLM$myFormula must be NULL or a formula object"); test <- FALSE }
           test <- .fun_testIfIn(test, "GLM$test", object@GLM$test, c('AIC', 'BIC', 'none'))
           fam <- 'none'
           if (!inherits(object@GLM$family, "family")) { cat("\nGLM$family must be a valid family object"); test <- FALSE }
           if (!is.list(object@GLM$control)) {cat("\nGLM$control must be a list like that returned by glm.control"); test <- FALSE}
           
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
           if(!is.null(object@RF$sampsize)) { test <- .fun_testIfPosInt(test, "RF$sampsize", object@RF$sampsize) }
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

### show.BIOMOD.models.options -------------------------------------------------
##'
##' @rdname BIOMOD.models.options
##' @importMethodsFrom methods show
##' @export
##'

setMethod('show', signature('BIOMOD.models.options'),
          function(object)
          {
            .bm_cat("BIOMOD.models.options")
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
            cat("\n            control = glm.control(", .print_control(object@GLM$control), ") ),", sep = "", fill = .Options$width)
            
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
              cat("\n            method = '", object@GAM$method, "', ", sep = "")
              cat("\n            optimizer = c('", paste(object@GAM$optimizer, collapse = "','"), "'),", sep = "")
              cat("\n            select = ", object@GAM$select, ",", sep = "")
              cat("\n            knots = ",  ifelse(length(object@GLM$knots) < 1, 'NULL', "'user.defined'"), ",", sep = "")
              cat("\n            paraPen = ",  ifelse(length(object@GLM$paraPen) < 1, 'NULL', "'user.defined'"), ",", sep = "")
            }
            cat("\n            control = list(", .print_control(object@GAM$control), ") ),", sep = "", fill = .Options$width)
            
            ## CTA options
            cat("\n")
            cat("\nCTA = list( method = '", object@CTA$method, "',", sep = "")
            cat("\n            parms = '", object@CTA$parms, "',", sep = "")
            cat("\n            cost = ", ifelse(length(object@CTA$cost) < 1, 'NULL', object@CTA$cost), ",", sep = "")
            cat("\n            control = list(", .print_control(object@CTA$control), ") ),", sep = "", fill = .Options$width)
            
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
                                                    , paste("list(", paste(.print_control(object@FDA$add_args), collapse = "")
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
            cat("\n           sampsize = ", ifelse(length(object@RF$sampsize) < 1, 'NULL', object@RF$sampsize), ",", sep = "")
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
            cat("\n MAXENT.Phillips.2 = list( myFormula = ", .print_formula(object@MAXENT.Phillips.2$myFormula), ",", sep = "")
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
            
            .bm_cat()
          }
)

