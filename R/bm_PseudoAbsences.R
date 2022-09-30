###################################################################################################
##' @name bm_PseudoAbsences
##' @aliases bm_PseudoAbsences
##' @aliases bm_PseudoAbsences_random
##' @aliases bm_PseudoAbsences_sre
##' @aliases bm_PseudoAbsences_disk
##' @aliases bm_PseudoAbsences_user.defined
##' @author Wilfried Thuiller, Damien Georges
##' 
##' @title Select pseudo-absences
##' 
##' @description This internal \pkg{biomod2} function allows to select pseudo-absences according 
##' to 4 different methods : \code{random}, \code{sre}, \code{disk} or \code{user.defined} (see Details).
##' 
##' 
##' @param resp.var a \code{vector}, \code{\link[sp]{SpatialPoints}} or 
##' \code{\link[sp]{SpatialPointsDataFrame}} object containing binary data (\code{0} : absence, 
##' \code{1} : presence, \code{NA} : indeterminate) for a single species that will be used to 
##' find the pseudo-absences
##' @param expl.var a \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} 
##' or \code{\link[raster:stack]{RasterStack}} object containing the explanatory variables (in 
##' columns or layers) that will be used to find the pseudo-absences
##' @param \ldots (\emph{optional, one or several of the following arguments depending on the selected 
##' method)}) 
##' 
##' @param nb.rep an \code{integer} corresponding to the number of sets (repetitions) of 
##' pseudo-absence points that will be drawn
##' @param strategy a \code{character} corresponding to the pseudo-absence selection strategy, 
##' must be among \code{random}, \code{sre}, \code{disk} or \code{user.defined}
##' @param nb.absences (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'random'} or \code{strategy = 'sre'} or \code{strategy = 'disk'}, an 
##' \code{integer} corresponding to the number of pseudo-absence points that will be selected for 
##' each pseudo-absence repetition (true absences included)
##' @param sre.quant (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'sre'}, a \code{numeric} between \code{0} and \code{0.5} defining the 
##' half-quantile used to make the \code{sre} pseudo-absence selection (see \code{\link{bm_SRE}})
##' @param dist.min (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'disk'}, a \code{numeric} defining the minimal distance to presence points 
##' used to make the \code{disk} pseudo-absence selection (in meters)
##' @param dist.max (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'disk'}, a \code{numeric} defining the maximal distance to presence points 
##' used to make the \code{disk} pseudo-absence selection (in meters)
##' @param user.table (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'user.defined'}, a \code{matrix} or \code{data.frame} with as many rows as 
##' \code{resp.var} values, as many columns as \code{nb.rep}, and containing \code{TRUE} or 
##' \code{FALSE} values defining which points will be used to build the species distribution 
##' model(s) for each repetition
##' 
##' 
##' @return 
##' 
##' A \code{list} containing the following elements :
##' \itemize{
##'   \item{\code{xy} : }{the coordinates of the species observations}
##'   \item{\code{sp} : }{the values of the species observations (\code{0}, \code{1} or \code{NA})}
##'   \item{\code{env} : }{the explanatory variables}
##'   \item{\code{pa.tab} : }{the corresponding table of selected pseudo-absences (indicated by 
##'   \code{TRUE} or \code{FALSE})}
##' }
##' 
##'
##' @details
##' 
##' \bold{Concerning random selection :}
##' 
##' The idea is to select pseudo-absences randomly in spatial locations where the species has not 
##' been sampled. This method is the simplest one and the most appropriate if lacking information 
##' about the presence sampling (non-exhaustive, biased sampling, etc). \cr \cr
##' 
##' \bold{Concerning SRE selection :}
##' 
##' The idea is to select pseudo-absences in spatial locations whose environmental conditions are 
##' different from those of the presence points. This method is appropriate when most of the 
##' environmental space of the species has been sampled. \cr \cr
##' 
##' \bold{Concerning disk selection :}
##' 
##' The idea is to select pseudo-absences, not too close from presence points, but not too far 
##' away either. This method is appropriate when most of the spatial range of the species has 
##' been sampled. \cr \cr
##' 
##' \bold{Concerning user defined selection :}
##' 
##' The user can provide pseudo-absences locations through a table containing spatial locations 
##' in rows, pseudo-absences repetitions in columns, and \code{TRUE/FALSE} values indicating 
##' whether each point is to be considered as pseudo-absence or not for each dataset.
##' 
##'
##' @keywords pseudo-absence random SRE disk
##' 
##' 
##' @seealso \code{\link{BIOMOD.formated.data.PA}}, \code{\link{BIOMOD_FormatingData}}
##' @family Secundary functions
##'
##' 
##' @importFrom raster extract coordinates subset reclassify mask cellFromXY xyFromCell distance 
##' sampleRandom
##' @importFrom sp SpatialPointsDataFrame coordinates
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PseudoAbsences <- function(resp.var, expl.var, nb.rep = 1, strategy = 'random', dist.min = 0, dist.max = NULL
                              , nb.absences = NULL, sre.quant = 0, user.table = NULL)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PseudoAbsences.check.args(resp.var, expl.var, nb.rep, strategy, dist.min, dist.max, nb.absences, sre.quant)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  if ((nb.rep == 0 | nb.absences <= 0) & strategy != 'user.defined') {
    out <- NULL
  } else {
    out <- switch(strategy,
                  user.defined = bm_PseudoAbsences_user.defined(resp.var, expl.var, user.table),
                  random = bm_PseudoAbsences_random(resp.var, expl.var, nb.absences, nb.rep),
                  sre = bm_PseudoAbsences_sre(resp.var, expl.var, sre.quant, nb.absences, nb.rep),
                  disk = bm_PseudoAbsences_disk(resp.var, expl.var, dist.min, dist.max, nb.absences, nb.rep))
  }
  
  return(out)
}


###################################################################################################

.bm_PseudoAbsences.check.args <- function(resp.var, expl.var, nb.rep, strategy, dist.min, dist.max, nb.absences, sre.quant)
{
  cat('\n\nChecking Pseudo-absence selection arguments...\n')
  
  ## 1. Check resp.var argument -----------------------------------------------
  if (is.vector(resp.var)) {
    resp.var <- SpatialPointsDataFrame(matrix(0, ncol = 2, nrow = length(resp.var)), data.frame(resp.var))
  }
  .fun_testIfInherits(TRUE, "resp.var", resp.var, "SpatialPoints")
  
  ## 2. Check expl.var argument -----------------------------------------------
  if (is.matrix(expl.var) | is.data.frame(expl.var)) {
    if (nrow(expl.var) != nrow(resp.var)) {
      stop("Species and Explanatory must have same dimensions")
    }
    expl.var <- SpatialPointsDataFrame(coordinates(resp.var), as.data.frame(expl.var))
  }
  .fun_testIfInherits(TRUE, "expl.var", expl.var, c("SpatialPoints", "Raster"))
  
  ## 3. Check strategy argument -----------------------------------------------
  availableStrategies <- c("random", "sre", "disk", "user.defined")
  if (!(strategy %in% availableStrategies) || sum(abs(coordinates(resp.var))) == 0) {
    # no coordinates or unknown strategy
    strategy <- "random"
    cat("\n   ! Random strategy was automatically selected (that can be due to points coordinates lack or unavailable strategy choosen)")
  }
  
  ## 4. Check nb.absences argument --------------------------------------------
  if (strategy != "user.defined") {
    if (is.null(nb.absences)) {
      stop("You must give the number of pseudo absences you want")
    } else{
      nbTrueAbs <- .get_nb_true_abs(resp.var)
      if (nbTrueAbs >= nb.absences) {
        cat("\n    ! There is more 'true absences' than desired pseudo absences. No pseudo absences selection done.")
        nb.absences = 0
      } else {
        nb.absences = nb.absences - nbTrueAbs
      }
    }
  }
  
  ## 5. Check dist.min and dist.max arguments ---------------------------------
  if (!is.null(dist.min) && dist.min < 0) {
    dist.min <- 0
  }
  if (!is.null(dist.max) && dist.max < 0) {
    dist.max <- NULL
  }
  if (!is.null(dist.max) && !is.null(dist.min) && dist.min >= dist.max) {
    stop("dist.min >= dist.max")
  }
  
  ## 6. Check sre.quant argument ----------------------------------------------
  if (strategy == 'SRE' && (sre.quant >= 0.5 || sre.quant < 0)) {
    stop("\n    ! SRE Quant should be a value between 0 and 0.5 ")
  }
  
  return(list(resp.var = resp.var,
              expl.var = expl.var,
              nb.rep = nb.rep,
              strategy = strategy,
              dist.min = dist.min,
              dist.max = dist.max,
              nb.absences = nb.absences,
              sre.quant = sre.quant))
}

###################################################################################################

.get_nb_true_abs <- function(sp)
{
  if (is.vector(sp)) {
    return(sum(sp == 0, na.rm = TRUE))
  } else if (inherits(sp, 'SpatialPoints')) {
    return(sum(sp@data == 0, na.rm = TRUE))
  } else if (inherits(sp, 'Raster')) {
    return(sum(sp[] == 0, na.rm = TRUE))
  }
}

.get_nb_available_pa_cells <- function(data, PA.flag = NA)
{
  if (is.vector(data) || is.data.frame(data) || is.matrix(data)) {
    return(ifelse(is.na(PA.flag), sum(is.na(data)), sum(data == PA.flag, na.rm = TRUE)))
  } else if(inherits(data, 'SpatialPoints')){
    return(ifelse(is.na(PA.flag), sum(is.na(data@data)), sum(data@data == PA.flag, na.rm = TRUE)))
  } else if(inherits(data, 'Raster')){
    return(ifelse(is.na(PA.flag), sum(is.na(data[])), sum(data[] == PA.flag, na.rm = TRUE)))
  }
}

.add_pa_rownames <- function(xy)
{
  rn <- row.names(xy)
  # if (length(rn) == 0) {
  #   rn = rep("", nrow(xy))
  # }
  missing_rn <- which(rn == "")
  if (length(missing_rn)) {
    rn[missing_rn] <- paste0("pa", 1:length(missing_rn))
  }
  rownames(xy) <- rn
  return(xy)
}


###################################################################################################

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_user.defined",
           def = function(resp.var, expl.var, ...) {
             standardGeneric( "bm_PseudoAbsences_user.defined")
           })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_user.defined', signature(expl.var = "SpatialPointsDataFrame"),
          function(resp.var, expl.var, user.table) {
            cat("\n   > User defined pseudo absences selection")
            return(list(xy = coordinates(resp.var),
                        sp = as.numeric(unlist(resp.var@data, use.names = FALSE)),
                        env = as.data.frame(expl.var@data),
                        pa.tab = user.table))
          })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_user.defined', signature(expl.var = "RasterStack"), 
          function(resp.var, expl.var, user.table) {
            cat("\n   > User defined pseudo absences selection")
            expl.var <- as.data.frame(extract(expl.var, coordinates(resp.var)))
            return(list(xy = coordinates(resp.var),
                        sp = as.numeric(unlist(resp.var@data, use.names = FALSE)), 
                        env = as.data.frame(expl.var),
                        pa.tab = user.table))
          })


###################################################################################################

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_random",
           def = function(resp.var, expl.var, ...) {
             standardGeneric( "bm_PseudoAbsences_random")
           })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_random', signature(expl.var = "SpatialPointsDataFrame"),
          function(resp.var, expl.var, nb.absences, nb.rep)
          {
            cat("\n   > random pseudo absences selection")
            
            # 1. Check if NA are present in resp.var observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(resp.var)
            if (nb.cells > 0) { # PA will be taken into response variable
              
              # 2. If nb NA < nb.absences, select all NA cells
              if (nb.cells <= nb.absences) {
                nb.rep <- 1
                nb.absences <- nb.cells
                cat("\n   > All available cells have been selected (", nb.absences, "cells )")
              }
              
              # 3. Select always the presences and the true absences
              pa.tab <- matrix(FALSE, ncol = nb.rep, nrow = nrow(resp.var))
              colnames(pa.tab) <- paste0("PA", 1:nb.rep)
              pa.tab[c(which(resp.var@data[, 1] == 1), which(resp.var@data[, 1] == 0)),] <- TRUE
              
              # 4. For each repetition, select among NA cells
              cand.cells <- which(is.na(resp.var@data[, 1]))
              for (j in 1:ncol(pa.tab)) {
                ## force to get at least one value of each factorial variable
                fact.level.cells <- bm_SampleFactorLevels(expl.var = as.data.frame(expl.var),
                                                          mask.out = pa.tab[, j, drop = FALSE])
                if (length(fact.level.cells)) {
                  pa.tab[fact.level.cells, j] <- TRUE
                  cand.cells <- setdiff(cand.cells, fact.level.cells)
                }
                pa.tab[sample(x = cand.cells,
                              size = nb.absences - length(fact.level.cells),
                              replace = FALSE), j] <- TRUE
              }
              
              return(list(xy = coordinates(resp.var),
                          sp = as.numeric(unlist(resp.var@data, use.names = FALSE)),
                          env = as.data.frame(expl.var@data),
                          pa.tab = pa.tab))
            } else {
              cat("\nUnsupported case yet!")
              return(NULL)
            }
          })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_random', signature(expl.var = "RasterStack"),
          function(resp.var, expl.var, nb.absences, nb.rep)
          {
            cat("\n   > random pseudo absences selection")
            
            # 1. Check if NA are present in resp.var observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(resp.var)
            if (nb.cells > 0) { # PA will be taken into response variable
              
              # 2. If nb NA < nb.absences, select all NA cells
              if (nb.cells <= nb.absences) {
                nb.rep <- 1
                nb.absences <- nb.cells
                cat("\n   > All available cells have been selected (", nb.absences, "cells )")
              }
              
              # 3. Select always the presences and the true absences
              pa.tab <- matrix(FALSE, ncol = nb.rep, nrow = nrow(resp.var))
              colnames(pa.tab) <- paste0("PA", 1:nb.rep)
              pa.tab[c(which(resp.var@data[, 1] == 1), which(resp.var@data[, 1] == 0)),] <- TRUE
              
              # 4. For each repetition, select among NA cells
              cand.cells <- which(is.na(resp.var@data[, 1]))
              for (j in 1:ncol(pa.tab)) {
                pa.tab[sample(x = cand.cells, size = nb.absences, replace = FALSE), j] <- TRUE
              }
              
              expl.var <- as.data.frame(extract(expl.var, coordinates(resp.var)))
              return(list(xy = coordinates(resp.var),
                          sp = as.numeric(unlist(resp.var@data, use.names = FALSE)), 
                          env = as.data.frame(expl.var),
                          pa.tab = as.data.frame(pa.tab)))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              
              # create a mask containing all not already sampled points (presences and absences)
              mask.env <- mask.out <- subset(expl.var, 1, drop = TRUE)
              mask.env <- reclassify(mask.env, c(-Inf, Inf, -1)) ## the area we want to sample
              mask.out[] <- NA
              
              # add presences and true absences in our mask
              in.cells <- cellFromXY(mask.env, coordinates(resp.var))
              mask.env[in.cells] <- NA
              mask.out[in.cells] <- 1
              
              # checking of nb candidates
              nb.cells <- .get_nb_available_pa_cells(mask.env, PA.flag = -1)
              
              # 2. If nb NA < nb.absences, select all NA cells
              if (nb.cells <= nb.absences) {
                nb.rep <- 1
                nb.absences <- nb.cells
                cat("\n   > All available cells have been selected (", nb.absences, "cells )")
              }
              
              if (nb.absences == 0) {
                cat("\n   > No cells are available (0 pseudo absences selected )")
                return(NULL)
              } else {
                # 4. For each repetition, select among raster cells
                pa.tab.tmp <- matrix(NA, ncol = nb.rep, nrow = nb.absences)
                for (j in 1:ncol(pa.tab.tmp)) {
                  SR <- NULL ## initialise the vector of sample cells
                  mask.env.tmp <- mask.env ## define a copy of the sampling mask
                  
                  ## force to get at least one value of each factorial variable
                  fact.level.cells <- bm_SampleFactorLevels(expl.var = expl.var, mask.out = mask.out)
                  if (length(fact.level.cells)) {
                    SR <- c(SR, fact.level.cells)
                    mask.env.tmp[SR] <- NA ## update the mask by removing already selected cells
                  }
                  
                  SR <- c(SR, sampleRandom(x = mask.env.tmp,
                                           size = nb.absences - length(SR),
                                           cells = TRUE,
                                           na.rm = TRUE)[, "cell", drop = TRUE])
                  ## repeat sampling until having the right number of points
                  while(length(SR)<nb.absences){
                    mask.env.tmp[SR] <- NA ## update the mask by removing already selected cells
                    SR <- c(SR, sampleRandom(x = mask.env.tmp,
                                             size = nb.absences - length(SR),
                                             cells = TRUE,
                                             na.rm = TRUE)[, "cell", drop = TRUE])
                  }
                  pa.tab.tmp[, j] <- SR
                }
                
                # putting cells in good format
                selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
                pa.tab <- matrix(FALSE, ncol = nb.rep, nrow = length(selected.cells))
                colnames(pa.tab) <- paste0("PA", 1:nb.rep)
                for (j in 1:ncol(pa.tab)) {
                  pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
                }
                
                # putting presences, true absences and pseudo absences together
                xy = coordinates(resp.var)
                rownames(xy) = paste0("pres", 1:nrow(xy))
                xy <- rbind(xy, xyFromCell(mask.env, selected.cells))
                xy <- .add_pa_rownames(xy)
                resp.var <- as.numeric(unlist(c(resp.var@data[, 1], rep(NA, length(selected.cells))), use.names = FALSE))
                expl.var <- extract(expl.var, xy)
                pa.tab <- rbind(matrix(TRUE, nrow = (nrow(xy) - length(selected.cells)),
                                       ncol = ncol(pa.tab)), pa.tab)
                
                return(list(xy = xy,
                            sp = resp.var,
                            env = as.data.frame(expl.var),
                            pa.tab = as.data.frame(pa.tab)))
              }
            }
          })


###################################################################################################

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_sre",
           def = function(resp.var, expl.var, ...) {
             standardGeneric("bm_PseudoAbsences_sre")
           })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_sre', signature(expl.var = "SpatialPointsDataFrame"), 
          function(resp.var, expl.var, sre.quant, nb.absences, nb.rep)
          {
            cat("\n   > SRE pseudo absences selection")
            
            # 0. calculate SRE to determine available pixels
            mask.in <- bm_SRE(resp.var = resp.var, expl.var = expl.var, new.env = expl.var@data, quant = sre.quant)
            mask.in <- data.frame(mask.in = !as.logical(mask.in)) ## revert the mask to sample PA out of SRE
            
            # 1. Check if NA are present in resp.var observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(mask.in$mask.in, PA.flag = TRUE)
            
            # 2. If nb NA < nb.absences, select all NA cells
            if (nb.cells <= nb.absences) {
              nb.rep <- 1
              nb.absences <- nb.cells
              cat("\n   > All available cells have been selected (", nb.absences, "cells )")
            }
            
            # 3. Select always the presences and the true absences
            pa.tab <- matrix(FALSE, ncol = nb.rep, nrow = nrow(resp.var))
            colnames(pa.tab) <- paste0("PA", 1:nb.rep)
            pa.tab[c(which(resp.var@data[, 1] == 1), which(resp.var@data[, 1] == 0)),] <- TRUE
            
            # 4. For each repetition, select among NA cells
            cand.cells <- which(!mask.in$mask.in)
            for (j in 1:ncol(pa.tab)) {
              ## force to get at least one value of each factorial variable
              fact.level.cells <- bm_SampleFactorLevels(expl.var = as.data.frame(expl.var),
                                                        mask.out = pa.tab[, j, drop = FALSE],
                                                        mask.in = mask.in)
              pa.tab[c(fact.level.cells,
                       sample(x = setdiff(cand.cells, fact.level.cells),
                              size = nb.absences - length(fact.level.cells),
                              replace = FALSE)), j] <- TRUE
            }
            
            return(list(xy = coordinates(resp.var),
                        sp = as.numeric(unlist(resp.var@data, use.names = FALSE)),
                        env = as.data.frame(expl.var@data),
                        pa.tab = pa.tab))
          })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_sre', signature(expl.var = "RasterStack"), 
          function(resp.var, expl.var, sre.quant, nb.absences, nb.rep)
          {
            cat("\n   > SRE pseudo absences selection")
            
            # 0. calculate SRE to determine available pixels
            mask.in <- bm_SRE(resp.var = resp.var, expl.var = expl.var, new.env = expl.var, quant = sre.quant)
            mask.in[mask.in[] > 0] <- NA ## remove points that are in SRE
            
            ## mask of already sampled points (presences/absences)
            mask.out <- subset(expl.var, 1)
            mask.out[] <- NA
            
            # 1. Check if NA are present in resp.var observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(mask.in, PA.flag = 0)
            
            # 2. If nb NA < nb.absences, select all NA cells
            if (nb.cells <= nb.absences) {
              nb.rep <- 1
              nb.absences <- nb.cells
              cat("\n   > All available cells have been selected (", nb.absences, "cells )")
            }
            
            # 4. For each repetition, select among raster cells
            pa.tab.tmp <- matrix(NA, ncol = nb.rep, nrow = nb.absences)
            for (j in 1:ncol(pa.tab.tmp)) {
              SR <- NULL ## initialise the vector of sample cells
              mask.in.tmp <- mask.in ## define a copy of the sampling mask
              
              ## force to get at least one value of each factorial variable
              fact.level.cells <- bm_SampleFactorLevels(expl.var = expl.var, mask.out = mask.out, mask.in = mask.in)
              if (length(fact.level.cells)) {
                SR <- c(SR, fact.level.cells)
                mask.in.tmp[SR] <- NA ## update the mask by removing already selected cells
              }
              
              SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                       size = nb.absences - length(SR),
                                       cells = TRUE,
                                       na.rm = TRUE)[, "cell", drop = TRUE])
              ## repeat sampling until having the right number of points
              while (length(SR) < nb.absences) {
                mask.in.tmp[SR] <- NA ## update the mask by removing already selected cells
                SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                         size = nb.absences - length(SR),
                                         cells = TRUE,
                                         na.rm = TRUE)[, "cell", drop = TRUE])
              }
              pa.tab.tmp[, j] <- SR
            }
            
            # putting cells in good format
            selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
            pa.tab <- matrix(FALSE, ncol = nb.rep, nrow = length(selected.cells))
            colnames(pa.tab) <- paste0("PA", 1:nb.rep)
            for (j in 1:ncol(pa.tab)) {
              pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
            }
            
            # putting presences, true absences and pseudo absences together
            xy <- rbind(coordinates(resp.var)[which(!is.na(resp.var@data[, 1])), ],
                        xyFromCell(mask.in, selected.cells))
            xy <- .add_pa_rownames(xy)
            resp.var <- as.numeric(unlist(c(na.omit(resp.var@data[, 1]), rep(NA, length(selected.cells))), use.names = FALSE))
            expl.var <- extract(expl.var, xy)
            pa.tab <- rbind(matrix(TRUE, nrow = (nrow(xy) - length(selected.cells)),
                                   ncol = ncol(pa.tab)), pa.tab)
            
            return(list(xy = xy,
                        sp = resp.var,
                        env = as.data.frame(expl.var),
                        pa.tab = as.data.frame(pa.tab)))
          })


###################################################################################################

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_disk",
           def = function(resp.var, expl.var, ...) {
             standardGeneric("bm_PseudoAbsences_disk")
           })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_disk', signature(expl.var = "SpatialPointsDataFrame"),
          function(resp.var, expl.var, dist.min, dist.max, nb.absences, nb.rep)
          {
            cat("\n   > Disk pseudo absences selection")
            
            # 1. determining area which can be selected
            coor <- coordinates(resp.var)
            pres <- which(resp.var@data[, 1] == 1)
            true.abs <- which(resp.var@data[, 1] == 0)
            tmp.abs <- which(is.na(resp.var@data[, 1]))
            outside <- rep(0, length(abs))
            inside <- rep(0, length(abs))
            
            for (i in 1:length(pres)) {
              # removing points too close from presences
              inside <- inside + (sqrt((coor[tmp.abs, 1] - coor[pres[i], 1]) ^ 2 + 
                                         (coor[tmp.abs, 2] - coor[pres[i], 2])^2) > dist.min)
              # keeping points not to far from presences
              if (!is.null(dist.max)) {
                outside <- outside + (sqrt((coor[tmp.abs, 1] - coor[pres[i], 1]) ^ 2 + 
                                             (coor[tmp.abs, 2] - coor[pres[i], 2]) ^ 2) < dist.max )
              }
            }
            if (is.null(dist.max)) { # no cells are too far
              outside <- outside + 1
            }
            selected.abs <- tmp.abs[(inside == length(pres)) & (outside > 0)]
            
            # 2. adding presences and true absences and selecting randomly pseudo absences
            return(bm_PseudoAbsences_random(resp.var[c(pres, true.abs, selected.abs), ],
                                            expl.var[c(pres, true.abs, selected.abs), ],
                                            nb.absences, nb.rep))
          })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_disk', signature(expl.var = "RasterStack"),
          function(resp.var, expl.var, dist.min, dist.max, nb.absences, nb.rep)
          {
            cat("\n   > Disk pseudo absences selection")
            
            # 1. Check if NA are present in resp.var observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(resp.var)
            if (nb.cells > 0) { # PA will be taken into response variable
              expl.var.tmp <- SpatialPointsDataFrame(coords = coordinates(resp.var),
                                                     data = as.data.frame(extract(expl.var, coordinates(resp.var))))
              return(bm_PseudoAbsences_disk(resp.var, expl.var.tmp, dist.min, dist.max, nb.absences, nb.rep))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              
              # create a mask
              dist.mask <- subset(expl.var, 1, drop = TRUE)
              dist.mask[] <- NA
              
              pres.xy <- coordinates(resp.var[which(resp.var@data[, 1] == 1), ])
              dist.mask[cellFromXY(dist.mask, pres.xy)] <- 1
              
              dist.mask <- distance(dist.mask)
              dist.mask <- mask(dist.mask, subset(expl.var, 1, drop = TRUE))
              
              if (is.null(dist.max)) { dist.max <- Inf }
              mask.in <- reclassify(dist.mask, c(-Inf, dist.min, NA, dist.min, dist.max, 1, dist.max, Inf, NA))
              mask.in[cellFromXY(mask.in, pres.xy)] <- 1
              mask.in = stack(expl.var * mask.in)
              names(mask.in) = names(expl.var)
              
              # 2. selecting randomly pseudo absences
              return(bm_PseudoAbsences_random(resp.var, expl.var = mask.in, nb.absences, nb.rep))
            }
          })

