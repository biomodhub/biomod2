##' 
##' @name sample.factor.levels
##' @aliases sample.factor.levels
##' 
##' @title Tool to ensure the sampling of all levels of a factorial variable 
##' @description This function will sample randomly an element of each level
##'   of all the factorial variables contains in a Raster* object or a data.frame
##' @author damien g.
##' 
##' @param x         a Raster* object or a data.frame
##' @param mask.out  a Raster/data.frame mask containing area that have already 
##'   been sampled. The factor levels within this mask will not be sampled.
##' @param mask.in   a Raster/list of Raster/data.frame mask (potentially a stack of 
##'   masks) indicating areas were we want to sample our factor level in priority.
##'   Note that if after having explore this masks some levels of the considered
##'   factorial varialble remains unsampled, this levels will be sampled in the 
##'   reference input object (here 'x')
##'      
##' @note 
##'   - The x/mask.out/mask.in should be coherent in term of dimention (same number of 
##'     rows for data.frame and same number of rows, column, identic resolution 
##'     and projection coordinates system for Raster* objects)
##'   - If mask.in contains several masks (RasterStack or multi-column data.frame)
##'     then the order of the mask matter. The mask will be considered successively.
##'     The first will be use prioritarly to sample our variable factor levels and 
##'     so on.
##'   - Raster* masks will be understood as: 
##'       - NA: out of of mask
##'       - not NA: in mask
##'   - data.frame masks will be understood as:
##'       - FALSE: out of mask
##'       - TRUE: in mask
##' 
##' @details In case any factorial variable is found in the input object then 
##'   NULL is returned.
##'   
##' @return a numeric vector the number (cell number for Raster* objects or row 
##'   number for data.frame) where each will refer to a single level of a single
##'   factorial variable.
##'   
##' @examples
##' ## example with raster* object ---------- 
##' library(raster)
##' ## create a factorial raster
##' r1 <- raster()
##' r1[] <- 1; r1[1] <- 2; r1[2:3] <- 3
##' r1 <- as.factor(r1)
##' ## create a continuous raster
##' r2 <- raster()
##' r2[] <- rnorm(ncell(r2))
##' ## pull the raster into a RasterStack
##' stk <- stack(r1, r2)
##' is.factor(stk)
##' 
##' ## define a mask for already sampled points
##' mask.out <- r1
##' mask.out[] <- NA; mask.out[2:3] <- 1
##' 
##' ## define a list of mask where we want to sample in priority
##' mask.in.1 <- mask.in.2 <- r1
##' mask.in.1[1:10] <- NA ## only level 1 should be sampled in this mask
##' mask.in.2[1] <- NA ## only levels 1 and 3 should be sampled in this mask
##' mask.in <- list(mask.in.1 = mask.in.1, 
##'                 mask.in.2 = mask.in.2)
##' 
##' ## test different version of the function
##' sample.factor.levels(stk, mask.out = mask.out)
##' sample.factor.levels(stk, mask.in = mask.in)
##' sample.factor.levels(stk, mask.out = mask.out, mask.in = mask.in)
##' 
sample.factor.levels <- function(x, mask.out = NULL, mask.in = NULL){
  ## make some checking of given parameters
  ## TODO(damien)
  if(inherits(x, 'Raster')){
    fact.level.cells <- .sample.factor.levels.raster(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else if(inherits(x, 'data.frame')){
    fact.level.cells <- .sample.factor.levels.data.frame(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else {
    warning(paste0("\nunsupported input data.",
                   "\nx should be a Raster* object or a data.frame.",
                   "\n NULL returned"))
    return(NULL)
  }
}

.sample.factor.levels.raster <- function(x, mask.out = NULL, mask.in = NULL){
  ## identificate the factorial variables
  fact.var <- which(is.factor(x))
  ## check if some factorial variables are in the input data
  if(any(fact.var)){ ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f){
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      ## get the levels of the factor on the full dataset
      fact.level.original <- unlist(raster::levels(subset(x, f)))
      fact.level <- fact.level.original
      cat("\n\t> fact.level for",  names(x)[f], ":\t", paste(fact.level, names(fact.level), sep = ":", collapse = "\t"))
      if(!is.null(mask.out)){ ## mask containing points that have already been sampled
        ## check the levels of the fector that have been already sampled
        fact.levels.sampled <- unlist(levels(as.factor(mask(subset(x, f), mask.out))))
        ## update levels names (lost during mask conversion)
        attr(fact.levels.sampled, "names") <- attr(fact.level.original, "names")[fact.levels.sampled]
        cat("\n\t - according to mask.out levels", fact.levels.sampled, "have already been sampled")
        ## update the list of factor levels to sample
        fact.level <- fact.level[!is.element(fact.level, fact.levels.sampled)]
      }
      if(length(fact.level)){
        ## try first to sample factors in the given masks
        if(!is.null(mask.in)){ ## list of mask we want to sample in (order matter!)
          for(mask.in.id in 1:length(mask.in)){
            ## check that some levels remains to be sampled
            if(length(fact.level)){
              ## update the masked version of the factorial raster
              x.f.masked <- as.factor(mask(subset(x, f), mask.in[[mask.in.id]]))
              x.f.levels <- unlist(levels(x.f.masked))
              ## update levels names (lost during mask conversion)
              attr(x.f.levels, "names") <- attr(fact.level.original, "names")[x.f.levels]              
              ## get the list of levels that coulb be sampled in this mask
              fact.levels.in.m.in <- fact.level[is.element(fact.level, x.f.levels)]
              if(length(fact.levels.in.m.in)){
                cat("\n\t - levels", fact.levels.in.m.in, "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, function(fl){
                  sample(which(x.f.masked[] == fl), 1)
                }))
                ## update the list of factor levels to sample
                fact.level <- fact.level[!is.element(fact.level, fact.levels.in.m.in)]
              }
            } 
          } ## end loop over mask.in
        }
        ## @note if some levels remains unsampled then we will take a random value of
        ## them in the full dataset => !! this should be tricky if mask.in arg is given
        ## because the value will be picked out of mask.in but is necessary to 
        ## ensure that models will run smoothly
        if(length(fact.level)){
          cat("\n\t - levels", fact.level, "will be sampled in the original raster")
          selected.cells <- c(selected.cells, sapply(fact.level, function(fl){
            sample(which(subset(x, f)[] == fl), 1)}))
        }
      }
      return(unlist(selected.cells))
    })))
    return(fact.level.cells)
  } else { ## no factorial variable
    return(NULL)
  }
}

.sample.factor.levels.data.frame.old <- function(x){
  ## identificate the factorial variables
  fact.var <- which(sapply(x, is.factor))
  ## check if some factorial variables are in the input data
  if(any(fact.var)){ ## some factorial variables present
    fact.level.cells <- as.numeric(sapply(fact.var, function(f){
      fact.level <- levels(x[, f])
      cat("\n fact.level for",  f, ":", fact.level)
      sapply(fact.level, function(fl){
        sample(which(x[,f] == fl), 1)
      })
    }))
    return(unique(fact.level.cells))
  } else { ## no factorial variable
    return(NULL)
  }
}

.sample.factor.levels.data.frame <- function(x, mask.out = NULL, mask.in = NULL){
  ## identificate the factorial variables
  fact.var <- which(sapply(x, is.factor))
  ## check if some factorial variables are in the input data
  if(any(fact.var)){ ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f){
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      ## get the levels of the factor on the full dataset
      fact.level.original <- levels(x[, f])
      fact.level <- fact.level.original
      cat("\n> fact.level for",  colnames(x)[f], ":\t", paste(1:length(fact.level), fact.level, sep = ":", collapse = "\t"))
      if(!is.null(mask.out)){ ## mask containing points that have already been sampled
        ## check the levels of the fector that have been already sampled
        fact.levels.sampled <- unique(na.omit(as.character(x[mask.out[, 1], f])))
        ## remove already sampled points from candidates
        x[mask.out[, 1], ] <- NA
#         ## update levels names (lost during mask conversion)
#         attr(fact.levels.sampled, "names") <- attr(fact.level.original, "names")[fact.levels.sampled]
        cat("\n - according to mask.out levels", fact.levels.sampled, "have already been sampled")
        ## update the list of factor levels to sample
        fact.level <- setdiff(fact.level, fact.levels.sampled)
      }
      if(length(fact.level)){
        ## try first to sample factors in the given masks
        if(!is.null(mask.in)){ ## list of mask we want to sample in (order matter!)
          for(mask.in.id in 1:ncol(mask.in)){
            ## check that some levels remains to be sampled
            if(length(fact.level)){
              ## update the masked version of the factorial raster
              x.f.masked <- as.character(x[, f])
              x.f.masked[!mask.in[, mask.in.id]] <- NA
              x.f.levels <- unique(na.omit(x.f.masked))
              ## get the list of levels that coulb be sampled in this mask
              fact.levels.in.m.in <- intersect(fact.level, x.f.levels)
              if(length(fact.levels.in.m.in)){
                cat("\n - levels", fact.levels.in.m.in, "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, function(fl){
                  candidate.cells <- na.omit(which(x.f.masked[] == fl))
                  selected.cell <- NULL
                  if(length(candidate.cells) == 1){ ## single candiate cell
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
        ## @note if some levels remains unsampled then we will take a random value of
        ## them in the full dataset => !! this should be tricky if mask.in arg is given
        ## because the value will be picked out of mask.in but is necessary to 
        ## ensure that models will run smoothly
        if(length(fact.level)){
          cat("\n - levels", fact.level, "will be sampled in the original data.frame")
          selected.cells <- c(selected.cells, sapply(fact.level, function(fl){
            candidate.cells <- na.omit(which(x[, f] == fl))
            selected.cell <- NULL
            if(length(candidate.cells) <= 1){ ## single candidate cell
              selected.cell <- candidate.cells
            } else { ## multi candidate cells
              selected.cell <- sample(candidate.cells, 1)
            }
            return(selected.cell)
          }))
        }
      }
      return(selected.cells)
    })))
    return(unique(fact.level.cells))
  } else { ## no factorial variable
    return(NULL)
  }
}
