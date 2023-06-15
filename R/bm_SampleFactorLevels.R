###################################################################################################
##' @name bm_SampleFactorLevels
##' @aliases bm_SampleFactorLevels
##' @aliases bm_SampleFactorLevels.raster
##' @aliases bm_SampleFactorLevels.data.frame
##' @author Damien Georges
##' 
##' @title Tool to ensure the sampling of all levels of a factorial variable 
##' 
##' @description This internal \pkg{biomod2} function samples randomly an element of each level of 
##' all the factorial variables contained in a \code{raster*} or \code{data.frame} object.
##' 
##' @param expl.var a \code{data.frame} or \code{\link[terra:rast]{SpatRaster}}
##' object containing the explanatory variables (in columns or layers)
##' @param mask.out a \code{data.frame} or \code{\link[terra:rast]{SpatRaster}}
##' object containing the area that has already been sampled (\emph{factor 
##' levels within this mask will not be sampled})
##' @param mask.in a \code{data.frame} or \code{\link[terra:rast]{SpatRaster}}
##' object containing areas where factor levels are to be sampled in priority. 
##' \emph{Note that if after having explored these masks, some factor levels 
##' remain unsampled, they will be sampled in the reference input object \code{expl.var}.}
##' 
##' 
##' @return  
##' 
##' A \code{numeric vector} containing point IDs (either cell number for \code{raster*} objects, 
##' or row number for \code{data.frame}), each refering to a single level of a single factorial 
##' variable.
##' 
##' In case any factorial variable is found in the input object, \code{NULL} is returned.
##' 
##'
##' @details 
##' 
##' The \code{expl.var}, \code{mask.out} and \code{mask.in} parameters must be coherent in terms of 
##' dimensions :
##' \itemize{
##'   \item same number of rows for \code{data.frame} objects
##'   \item same resolution, projection system and number of cells for \code{raster*} objects 
##'   \cr \cr
##' }
##' 
##' If \code{mask.in} contains several masks (either it is a 
##' \code{\link[raster:stack]{RasterStack}} object or a multi-columns \code{data.frame}), then 
##' the order of masks / columns matters : they will be considered successively to sample missing 
##' factor levels. \cr \cr
##' 
##' \itemize{
##'   \item \code{raster*} masks will be understood as :
##'   \itemize{
##'     \item \code{NA} : out of mask
##'     \item \code{not NA} : in mask
##'   }
##'   \item \code{data.frame} masks will be understood as :
##'   \itemize{
##'     \item \code{FALSE} : out of mask
##'     \item \code{TRUE} : in mask
##'   }
##' }
##' 
##' 
##' @seealso \code{\link{bm_PseudoAbsences}}, \code{\link{bm_RunModelsLoop}}, 
##' \code{\link{BIOMOD_Modeling}}
##' @family Secundary functions
##' 
##'   
##' @examples
##' 
##' library(terra)
##' 
##' ## Create raster data
##' ras.1 <- ras.2 <- mask.out <- rast(nrows = 10, ncols = 10)
##' ras.1[] <- as.factor(rep(c(1, 2, 3, 4, 5), each = 20))
##' ras.1 <- as.factor(ras.1)
##' ras.2[] <- rnorm(100)
##' stk <- c(ras.1, ras.2)
##' names(stk) <- c("varFact", "varNorm")
##' 
##' ## define a mask for already sampled points
##' mask.out[1:40] <- 1
##' 
##' ## define a list of masks where we want to sample in priority
##' mask.in <- list(ras.1, ras.1)
##' mask.in[[1]][1:80] <- NA ## only level 5 should be sampled in this mask
##' mask.in[[1]][21:80] <- NA ## only levels 1 and 5 should be sampled in this mask
##' 
##' ## Sample all factor levels
##' samp1 <- bm_SampleFactorLevels(expl.var = stk, mask.out = mask.out)
##' samp2 <- bm_SampleFactorLevels(expl.var = stk, mask.in = mask.in)
##' samp3 <- bm_SampleFactorLevels(expl.var = stk, mask.out = mask.out, mask.in = mask.in)
##' 
##' 
##' @importFrom terra rast cats mask subset is.factor values
##' @export
##' 
##' 
###################################################################################################

bm_SampleFactorLevels <- function(expl.var, mask.out = NULL, mask.in = NULL)
{
  if (inherits(expl.var, 'SpatRaster')) {
    fact.level.cells <- bm_SampleFactorLevels.SpatRaster(expl.var, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else if (inherits(expl.var, 'data.frame')) {
    fact.level.cells <- bm_SampleFactorLevels.data.frame(expl.var, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else {
    warning(paste0("\nunsupported input data.",
                   "\nexpl.var should be a Raster* object or a data.frame.",
                   "\n NULL returned"))
    return(NULL)
  }
}

bm_SampleFactorLevels.SpatRaster <- function(expl.var, mask.out = NULL, mask.in = NULL)
{
  ## check if some factorial variables are in the input data
  fact.var <- which(is.factor(expl.var))
  if(any(fact.var))
  { ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f)
    {
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      
      ## get the factor levels on the full dataset
      fact.level.names <- cats(subset(expl.var, f))[[1]][,2]
      fact.level <- fact.level.original <- cats(subset(expl.var, f))[[1]][,1]
      cat("\n\t> fact.level for",  names(expl.var)[f], ":\t", paste(fact.level, fact.level.names, sep = ":", collapse = "\t"))
      ## mask containing points that have already been sampled ------------------------------------
      if (!is.null(mask.out))
      {
        ## check the factor levels that have already been sampled
        fact.levels.sampled <- 
          unique(na.omit(
            values((mask(subset(expl.var, f), mask.out, maskvalues = c(1), inverse = TRUE)))
          ))
        ## update levels names (lost during mask conversion)
        cat("\n\t - according to mask.out levels", fact.levels.sampled, "have already been sampled")
        ## update the list of factor levels to sample
        fact.level <- setdiff(fact.level, fact.levels.sampled)
      }
      ## if there still is some levels to sample --------------------------------------------------
      ## take a random value of them in the full dataset 
      if (length(fact.level) > 0) {
        cat("\n\t - levels", fact.level, "will be sampled in the original raster")
        selected.cells <- c(selected.cells, sapply(fact.level, function(fl) {
          ind = which(subset(expl.var, f)[] == fl)
          if (length(ind) > 0) {
            return(sample(ind, 1))
          }
        }))
      }
      
      return(unlist(selected.cells))
    })))
    return(unique(fact.level.cells))
  } else { ## no factorial variable
    return(NULL)
  }
}

bm_SampleFactorLevels.data.frame <- function(expl.var, mask.out = NULL, mask.in = NULL)
{
  ## check if some factorial variables are in the input data
  fact.var <- which(sapply(expl.var, is.factor))
  if(any(fact.var))
  { ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f)
    {
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      
      ## get the factor levels on the full dataset
      fact.level <- fact.level.original <- levels(expl.var[, f])
      cat("\n\t> fact.level for",  colnames(expl.var)[f], ":\t", paste(1:length(fact.level), fact.level, sep = ":", collapse = "\t"))
      
      ## mask containing points that have already been sampled ------------------------------------
      if (!is.null(mask.out))
      {
        ## check the factor levels that have already been sampled
        fact.levels.sampled <- unique(na.omit(as.character(expl.var[mask.out[, 1], f])))
        ## remove already sampled points from candidates
        expl.var[mask.out[, 1], ] <- NA
        cat("\n\t - according to mask.out levels", fact.levels.sampled, "have already been sampled")
        ## update the list of factor levels to sample
        fact.level <- setdiff(fact.level, fact.levels.sampled)
      }
      
      ## if there still is some levels to sample --------------------------------------------------
      if(length(fact.level))
      {
        ## a. try first to sample factors in the given masks -------------------
        if (!is.null(mask.in))
        { ## list of mask we want to sample in (order matter!)
          for (mask.in.id in 1:ncol(mask.in))
          {
            if (length(fact.level) > 0) ## if there still is some levels to sample
            {
              ## update the masked version of the factorial raster
              x.f.masked <- as.character(expl.var[, f])
              x.f.masked[!mask.in[, mask.in.id]] <- NA
              x.f.levels <- unique(na.omit(x.f.masked))
              ## get the list of levels that could be sampled in this mask
              fact.levels.in.m.in <- intersect(fact.level, x.f.levels)
              if (length(fact.levels.in.m.in) > 0) {
                cat("\n\t - levels", fact.levels.in.m.in, "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, function(fl){
                  candidate.cells <- na.omit(which(x.f.masked[] == fl))
                  selected.cell <- NULL
                  if (length(candidate.cells) == 1) { ## single candidate cell
                    selected.cell <- candidate.cells
                  } else if (length(candidate.cells) > 1) { ## multi candidate cells
                    selected.cell <- sample(candidate.cells, 1)
                  }
                  return(selected.cell)
                }))
                ## update the list of factor levels to sample
                fact.level <- setdiff(fact.level, fact.levels.in.m.in)
              }
            } 
          } ## end loop over mask.in
        }
        
        ## if there still is some levels to sample ------------------------------------------------
        ## b. take a random value of them in the full dataset 
        ## !! this should be tricky if mask.in arg is given because the value will be picked out of 
        ## mask.in but is necessary to ensure that models will run smoothly
        if (length(fact.level) > 0){
          cat("\n\t - levels", fact.level, "will be sampled in the original data.frame")
          selected.cells <- c(selected.cells, sapply(fact.level, function(fl) {
            candidate.cells <- na.omit(which(expl.var[, f] == fl))
            selected.cell <- NULL
            if (length(candidate.cells) <= 1) { ## single candidate cell
              selected.cell <- candidate.cells
            } else if (length(candidate.cells) > 1) { ## multi candidate cells
              selected.cell <- sample(candidate.cells, 1)
            }
            return(selected.cell)
          }))
        }
      }
      return(unlist(selected.cells))
    })))
    return(unique(fact.level.cells))
  } else { ## no factorial variable
    return(NULL)
  }
}

