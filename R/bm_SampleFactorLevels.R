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
##' @param x a \code{data.frame} or \code{\link[raster:stack]{RasterStack}} object containing the 
##' explanatory variables (in columns or layers)
##' @param mask.out a \code{data.frame} or \code{\link[raster:raster]{raster}} object containing 
##' the area that has already been sampled (\emph{factor levels within this mask will not be 
##' sampled})
##' @param mask.in a \code{data.frame}, \code{\link[raster:raster]{raster}} or 
##' \code{\link[raster:stack]{RasterStack}} object containing areas where factor levels are to be 
##' sampled in priority. \emph{Note that if after having explored these masks, some factor levels 
##' remain unsampled, they will be sampled in the reference input object \code{x}.}
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
##' The \code{x}, \code{mask.out} and \code{mask.in} parameters must be coherent in terms of 
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
##' @examples
##' 
##' ## Example with raster* object ----------------------------------
##' 
##' library(raster)
##' 
##' ## create a factorial raster
##' r1 <- raster()
##' r1[] <- 1; r1[1] <- 2; r1[2:3] <- 3
##' r1 <- as.factor(r1)
##' 
##' ## create a continuous raster
##' r2 <- raster()
##' r2[] <- rnorm(ncell(r2))
##' 
##' ## put the raster into a RasterStack
##' stk <- stack(r1, r2)
##' is.factor(stk)
##' 
##' ## define a mask for already sampled points
##' mask.out <- r1
##' mask.out[] <- NA; mask.out[2:3] <- 1
##' 
##' ## define a list of masks where we want to sample in priority
##' mask.in.1 <- mask.in.2 <- r1
##' mask.in.1[1:10] <- NA ## only level 1 should be sampled in this mask
##' mask.in.2[1] <- NA ## only levels 1 and 3 should be sampled in this mask
##' mask.in <- list(mask.in.1 = mask.in.1, 
##'                 mask.in.2 = mask.in.2)
##' 
##' ## test different version of the function
##' bm_SampleFactorLevels(stk, mask.out = mask.out)
##' bm_SampleFactorLevels(stk, mask.in = mask.in)
##' bm_SampleFactorLevels(stk, mask.out = mask.out, mask.in = mask.in)
##' 
##' 
##' @importFrom raster is.factor as.factor levels subset mask
##' 
##' @export
##' 
##' 
###################################################################################################

bm_SampleFactorLevels <- function(x, mask.out = NULL, mask.in = NULL)
{
  ## make some checking of given parameters
  ## TODO(damien)
  
  if (inherits(x, 'Raster')) {
    fact.level.cells <- bm_SampleFactorLevels.raster(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else if (inherits(x, 'data.frame')) {
    fact.level.cells <- bm_SampleFactorLevels.data.frame(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else {
    warning(paste0("\nunsupported input data.",
                   "\nx should be a Raster* object or a data.frame.",
                   "\n NULL returned"))
    return(NULL)
  }
}

bm_SampleFactorLevels.raster <- function(x, mask.out = NULL, mask.in = NULL)
{
  ## check if some factorial variables are in the input data
  fact.var <- which(is.factor(x))
  if(any(fact.var))
  { ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f)
    {
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      
      ## get the factor levels on the full dataset
      fact.level <- fact.level.original <- unlist(levels(subset(x, f)))
      cat("\n\t> fact.level for",  names(x)[f], ":\t", paste(fact.level, names(fact.level), sep = ":", collapse = "\t"))
      
      ## mask containing points that have already been sampled ------------------------------------
      if (!is.null(mask.out))
      {
        ## check the factor levels that have already been sampled
        fact.levels.sampled <- unlist(levels(as.factor(mask(subset(x, f), mask.out))))
        ## update levels names (lost during mask conversion)
        attr(fact.levels.sampled, "names") <- attr(fact.level.original, "names")[fact.levels.sampled]
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
          for (mask.in.id in 1:length(mask.in))
          {
            if (length(fact.level)) ## if there still is some levels to sample
            {
              ## update the masked version of the factorial raster
              x.f.masked <- as.factor(mask(subset(x, f), mask.in[[mask.in.id]]))
              x.f.levels <- unlist(levels(x.f.masked))
              ## update levels names (lost during mask conversion)
              attr(x.f.levels, "names") <- attr(fact.level.original, "names")[x.f.levels]              
              ## get the list of levels that could be sampled in this mask
              fact.levels.in.m.in <- intersect(fact.level, x.f.levels)
              if (length(fact.levels.in.m.in)) {
                cat("\n\t - levels", fact.levels.in.m.in, "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, function(fl){
                  sample(which(x.f.masked[] == fl), 1)
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
        if (length(fact.level)){
          cat("\n\t - levels", fact.level, "will be sampled in the original raster")
          selected.cells <- c(selected.cells, sapply(fact.level, function(fl) {
            sample(which(subset(x, f)[] == fl), 1)
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

bm_SampleFactorLevels.data.frame <- function(x, mask.out = NULL, mask.in = NULL)
{
  ## check if some factorial variables are in the input data
  fact.var <- which(sapply(x, is.factor))
  if(any(fact.var))
  { ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f)
    {
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      
      ## get the factor levels on the full dataset
      fact.level <- fact.level.original <- levels(x[, f])
      cat("\n\t> fact.level for",  colnames(x)[f], ":\t", paste(1:length(fact.level), fact.level, sep = ":", collapse = "\t"))
      
      ## mask containing points that have already been sampled ------------------------------------
      if (!is.null(mask.out))
      {
        ## check the factor levels that have already been sampled
        fact.levels.sampled <- unique(na.omit(as.character(x[mask.out[, 1], f])))
        ## remove already sampled points from candidates
        x[mask.out[, 1], ] <- NA
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
            if (length(fact.level)) ## if there still is some levels to sample
            {
              ## update the masked version of the factorial raster
              x.f.masked <- as.character(x[, f])
              x.f.masked[!mask.in[, mask.in.id]] <- NA
              x.f.levels <- unique(na.omit(x.f.masked))
              ## get the list of levels that could be sampled in this mask
              fact.levels.in.m.in <- intersect(fact.level, x.f.levels)
              if (length(fact.levels.in.m.in)) {
                cat("\n\t - levels", fact.levels.in.m.in, "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, function(fl){
                  candidate.cells <- na.omit(which(x.f.masked[] == fl))
                  selected.cell <- NULL
                  if (length(candidate.cells) == 1) { ## single candidate cell
                    selected.cell <- candidate.cells
                  } else { ## multi candidate cells
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
        if (length(fact.level)){
          cat("\n\t - levels", fact.level, "will be sampled in the original data.frame")
          selected.cells <- c(selected.cells, sapply(fact.level, function(fl) {
            candidate.cells <- na.omit(which(x[, f] == fl))
            selected.cell <- NULL
            if (length(candidate.cells) <= 1) { ## single candidate cell
              selected.cell <- candidate.cells
            } else { ## multi candidate cells
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

