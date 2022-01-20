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
##' @description 
##' 
##' This internal \pkg{biomod2} function allows to select pseudo-absences according to 4 different 
##' methods : \code{random}, \code{sre}, \code{disk} or \code{user.defined} (see 
##' \href{bm_PseudoAbsences/html#details}{Details}).
##' 
##' @param sp the species observations
##' @param env the explanatory variables
##' @param nb.repet the number of repetitions
##' @param strategy a \code{character} corresponding to the pseudo-absence selection strategy, 
##' must be among \code{random}, \code{sre}, \code{disk} or \code{user.defined}
##' @param nb.points (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy} is \code{random}, \code{sre} or \code{disk}, the number of pseudo-absences 
##' to be selected
##' @param quant.SRE (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'sre'}, a \code{numeric} between \code{0} and \code{0.5} corresponding to 
##' the most extreme value for each variable not to be taken into account for determining the 
##' tolerance boundaries of the considered species (see \code{\link{bm_SRE}})
##' @param distMin (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'disk'}, the minimum distance away from presence points to select 
##' pseudo-absences
##' @param distMax (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'disk'}, the maximum distance away from presence points to select 
##' pseudo-absences
##' @param PA.table (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'user.defined'}, a \code{matrix} or \code{data.frame} with as many rows as 
##' \code{resp.var} values, as many columns as \code{PA.nb.rep}, and containing \code{TRUE} or 
##' \code{FALSE} values defining which points will be used to build the species distribution 
##' model(s) for each repetition (see \code{\link{BIOMOD_FormatingData}})
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
##' \bold{Concerning SRE selection :}
##' 
##' \bold{Concerning disk selection :}
##' 
##' \bold{Concerning user defined selection :}
##' 
##'
##' @keywords pseudo-absence, random, SRE, disk
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}
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


bm_PseudoAbsences <- function(sp, env, nb.repet = 1, strategy = 'random', distMin = 0, distMax = NULL
                              , nb.points = NULL, quant.SRE = 0, PA.table = NULL)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PseudoAbsences.check.args(sp, env, nb.repet, strategy, distMin, distMax, nb.points, quant.SRE)
  sp <- args$sp
  env <- args$env
  nb.repet <- args$nb.repet
  strategy <- args$strategy
  distMin <- args$distMin
  distMax <- args$distMax
  nb.points <- args$nb.points
  quant.SRE <- args$quant.SRE
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  if ((nb.repet == 0 | nb.points <= 0) & strategy != 'user.defined') {
    out <- NULL
  } else {
    out <- switch(strategy,
                  user.defined = bm_PseudoAbsences_user.defined(sp, env, PA.table),
                  random = bm_PseudoAbsences_random(sp, env, nb.points, nb.repet),
                  sre = bm_PseudoAbsences_sre(sp, env, quant.SRE, nb.points, nb.repet),
                  disk = bm_PseudoAbsences_disk(sp, env, distMin, distMax, nb.points, nb.repet))
  }
  
  return(out)
}


###################################################################################################

.bm_PseudoAbsences.check.args <- function(sp, env, nb.repet, strategy, distMin, distMax, nb.points, quant.SRE)
{
  cat('\n\nChecking Pseudo-absence selection arguments...\n')

  ## 1. Check sp argument -----------------------------------------------------
  if (is.vector(sp)) {
    sp <- SpatialPointsDataFrame(matrix(0, ncol = 2, nrow = length(sp)), data.frame(sp))
  }
  .fun_testIfInherits(TRUE, "sp", sp, "SpatialPoints")
  
  ## 2. Check env argument ----------------------------------------------------
  if (is.matrix(env) | is.data.frame(env)) {
    if (nrow(env) != nrow(sp)) {
      stop("Species and Explanatory must have same dimensions")
    }
    env <- SpatialPointsDataFrame(coordinates(sp), as.data.frame(env))
  }
  .fun_testIfInherits(TRUE, "env", env, c("SpatialPoints", "Raster"))
  
  ## 3. Check strategy argument -----------------------------------------------
  availableStrategies <- c("random", "sre", "disk", "user.defined")
  if (!(strategy %in% availableStrategies) || sum(abs(coordinates(sp))) == 0) {
    # no coordinates or unknown strategy
    strategy <- "random"
    cat("\n   ! Random strategy was automatically selected (that can be due to points coordinates lack or unavailable strategy choosen)")
  }
  
  ## 4. Check nb.points argument ----------------------------------------------
  if (strategy != "user.defined") {
    if (is.null(nb.points)) {
      stop("You must give the number of pseudo absences you want")
    } else{
      nbTrueAbs <- .get_nb_true_abs(sp)
      if (nbTrueAbs >= nb.points) {
        cat("\n    ! There is more 'true absences' than desired pseudo absences. No pseudo absences selection done.")
        nb.points = 0
      } else {
        nb.points = nb.points - nbTrueAbs
      }
    }
  }
  
  ## 5. Check distMin and distMax arguments -----------------------------------
  if (!is.null(distMin) && distMin < 0) {
    distMin <- 0
  }
  if (!is.null(distMax) && distMax < 0) {
    distMax <- NULL
  }
  if (!is.null(distMax) && !is.null(distMin) && distMin >= distMax) {
    stop("distMin >= distMax")
  }
  
  ## 6. Check quant.SRE argument ----------------------------------------------
  if (strategy == 'SRE' && (quant.SRE >= 0.5 || quant.SRE < 0)) {
    stop("\n    ! SRE Quant should be a value between 0 and 0.5 ")
  }
  
  return(list(sp = sp,
              env = env,
              nb.repet = nb.repet,
              strategy = strategy,
              distMin = distMin,
              distMax = distMax,
              nb.points = nb.points,
              quant.SRE = quant.SRE))
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
           def = function(sp, env, ...) {
             standardGeneric( "bm_PseudoAbsences_user.defined")
           })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_user.defined', signature(env = "SpatialPointsDataFrame"),
          function(sp, env, pa.table) {
            cat("\n   > User defined pseudo absences selection")
            return(list(xy = coordinates(sp),
                        sp = as.vector(sp@data),
                        env = as.data.frame(env@data),
                        pa.tab = pa.table))
          })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_user.defined', signature(env = "RasterStack"), 
          function(sp, env, pa.table) {
            cat("\n   > User defined pseudo absences selection")
            env <- as.data.frame(extract(env, coordinates(sp)))
            return(list(xy = coordinates(sp),
                        sp = as.numeric(unlist(sp@data, use.names = FALSE)), 
                        env = as.data.frame(env),
                        pa.tab = pa.table))
          })


###################################################################################################

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_random",
           def = function(sp, env, ...) {
             standardGeneric( "bm_PseudoAbsences_random")
           })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_random', signature(env = "SpatialPointsDataFrame"),
          function(sp, env, nb.points, nb.repet)
          {
            cat("\n   > random pseudo absences selection")
            
            # 1. Check if NA are present in sp observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(sp)
            if (nb.cells > 0) { # PA will be taken into response variable
              
              # 2. If nb NA < nb.points, select all NA cells
              if (nb.cells <= nb.points) {
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "cells )")
              }
              
              # 3. Select always the presences and the true absences
              pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = nrow(sp))
              colnames(pa.tab) <- paste0("PA", 1:nb.repet)
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              
              # 4. For each repetition, select among NA cells
              cand.cells <- which(is.na(sp@data))
              for (j in 1:ncol(pa.tab)) {
                ## force to get at least one value of each factorial variable
                fact.level.cells <- bm_SampleFactorLevels(x = as.data.frame(env),
                                                          mask.out = pa.tab[, j, drop = FALSE])
                if (length(fact.level.cells)) {
                  pa.tab[fact.level.cells, j] <- TRUE
                  cand.cells <- setdiff(cand.cells, fact.level.cells)
                }
                pa.tab[sample(x = cand.cells,
                              size = nb.points - length(fact.level.cells),
                              replace = FALSE), j] <- TRUE
              }
              
              return(list(xy = coordinates(sp),
                          sp = as.vector(sp@data),
                          env = as.data.frame(env@data),
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

setMethod('bm_PseudoAbsences_random', signature(env = "RasterStack"),
          function(sp, env, nb.points, nb.repet)
            {
            cat("\n   > random pseudo absences selection")
            
            # 1. Check if NA are present in sp observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(sp)
            if (nb.cells > 0) { # PA will be taken into response variable
              
              # 2. If nb NA < nb.points, select all NA cells
              if (nb.cells <= nb.points) {
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "cells )")
              }
              
              # 3. Select always the presences and the true absences
              pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = nrow(sp))
              colnames(pa.tab) <- paste0("PA", 1:nb.repet)
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              
              # 4. For each repetition, select among NA cells
              cand.cells <- which(is.na(sp@data))
              for (j in 1:ncol(pa.tab)) {
                pa.tab[sample(x = cand.cells, size = nb.points, replace = FALSE), j] <- TRUE
              }
              
              env <- as.data.frame(extract(env, coordinates(sp)))
              return(list(xy = coordinates(sp),
                          sp = as.numeric(unlist(sp@data, use.names = FALSE)), 
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              
              # create a mask containing all not already sampled points (presences and absences)
              mask.env <- mask.out <- subset(env, 1, drop = TRUE)
              mask.env <- reclassify(mask.env, c(-Inf, Inf, -1)) ## the area we want to sample
              mask.out[] <- NA
              
              # add presences and true absences in our mask
              in.cells <- cellFromXY(mask.env, coordinates(sp))
              mask.env[in.cells] <- NA
              mask.out[in.cells] <- 1
              
              # checking of nb candidates
              nb.cells <- .get_nb_available_pa_cells(mask.env, PA.flag = -1)
              
              # 2. If nb NA < nb.points, select all NA cells
              if (nb.cells <= nb.points) {
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "cells )")
              }
              
              if (nb.points == 0) {
                cat("\n   > No cells are available (0 pseudo absences selected )")
                return(NULL)
              } else {
                # 4. For each repetition, select among raster cells
                pa.tab.tmp <- matrix(NA, ncol = nb.repet, nrow = nb.points)
                for (j in 1:ncol(pa.tab.tmp)) {
                  SR <- NULL ## initialise the vector of sample cells
                  mask.env.tmp <- mask.env ## define a copy of the sampling mask
                  
                  ## force to get at least one value of each factorial variable
                  fact.level.cells <- bm_SampleFactorLevels(env, mask.out = mask.out)
                  if (length(fact.level.cells)) {
                    SR <- c(SR, fact.level.cells)
                    mask.env.tmp[SR] <- NA ## update the mask by removing already selected cells
                  }
                  
                  SR <- c(SR, sampleRandom(x = mask.env.tmp,
                                           size = nb.points - length(SR),
                                           cells = TRUE,
                                           na.rm = TRUE)[, "cell", drop = TRUE])
                  ## repeat sampling until having the right number of points
                  while(length(SR)<nb.points){
                    mask.env.tmp[SR] <- NA ## update the mask by removing already selected cells
                    SR <- c(SR, sampleRandom(x = mask.env.tmp,
                                             size = nb.points - length(SR),
                                             cells = TRUE,
                                             na.rm = TRUE)[, "cell", drop = TRUE])
                  }
                  pa.tab.tmp[, j] <- SR
                }
                
                # putting cells in good format
                selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
                pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
                colnames(pa.tab) <- paste0("PA", 1:nb.repet)
                for (j in 1:ncol(pa.tab)) {
                  pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
                }
                
                # putting presences, true absences and pseudo absences together
                xy <- rbind(coordinates(sp), xyFromCell(mask.env, selected.cells))
                xy <- .add_pa_rownames(xy)
                sp <- as.numeric(unlist(c(as.vector(sp@data), rep(NA, length(selected.cells))), use.names = FALSE))
                env <- extract(env, xy)
                pa.tab <- rbind(matrix(TRUE, nrow = (nrow(xy) - length(selected.cells)),
                                       ncol = ncol(pa.tab)), pa.tab)
                
                return(list(xy = xy,
                            sp = sp,
                            env = as.data.frame(env),
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
           def = function(sp, env, ...) {
             standardGeneric("bm_PseudoAbsences_sre")
           })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_sre', signature(env = "SpatialPointsDataFrame"), 
          function(sp, env, quant.SRE, nb.points, nb.repet)
          {
            cat("\n   > SRE pseudo absences selection")
            
            # 0. calculate SRE to determine available pixels
            mask.in <- bm_SRE(Response = sp, Explanatory = env, NewData = env@data, Quant = quant.SRE)
            mask.in <- data.frame(mask.in = !as.logical(mask.in)) ## revert the mask to sample PA out of SRE
            
            # 1. Check if NA are present in sp observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(mask.in$mask.in, PA.flag = TRUE)
            
            # 2. If nb NA < nb.points, select all NA cells
            if (nb.cells <= nb.points) {
              nb.repet <- 1
              nb.points <- nb.cells
              cat("\n   > All available cells have been selected (", nb.points, "cells )")
            }
            
            # 3. Select always the presences and the true absences
            pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = nrow(sp))
            colnames(pa.tab) <- paste0("PA", 1:nb.repet)
            pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
            
            # 4. For each repetition, select among NA cells
            cand.cells <- which(!mask.in$mask.in)
            for (j in 1:ncol(pa.tab)) {
              ## force to get at least one value of each factorial variable
              fact.level.cells <- bm_SampleFactorLevels(as.data.frame(env),
                                                        mask.out = pa.tab[, j, drop = FALSE],
                                                        mask.in = mask.in)
              pa.tab[c(fact.level.cells,
                       sample(x = setdiff(cand.cells, fact.level.cells),
                              size = nb.points - length(fact.level.cells),
                              replace = FALSE)), j] <- TRUE
            }
            
            return(list(xy = coordinates(sp),
                        sp = as.vector(sp@data),
                        env = as.data.frame(env@data),
                        pa.tab = pa.tab))
          })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_sre', signature(env = "RasterStack"), 
          function(sp, env, quant.SRE, nb.points, nb.repet)
          {
            cat("\n   > SRE pseudo absences selection")
            
            # 0. calculate SRE to determine available pixels
            mask.in <- bm_SRE(Response = sp, Explanatory = env, NewData = env, Quant = quant.SRE)
            mask.in[mask.in[] > 0] <- NA ## remove points that are in SRE
            
            ## mask of already sampled points (presences/absences)
            mask.out <- subset(env, 1)
            mask.out[] <- NA
            
            # 1. Check if NA are present in sp observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(mask.in, PA.flag = 0)
            
            # 2. If nb NA < nb.points, select all NA cells
            if (nb.cells <= nb.points) {
              nb.repet <- 1
              nb.points <- nb.cells
              cat("\n   > All available cells have been selected (", nb.points, "cells )")
            }
            
            # 4. For each repetition, select among raster cells
            pa.tab.tmp <- matrix(NA, ncol = nb.repet, nrow = nb.points)
            for (j in 1:ncol(pa.tab.tmp)) {
              SR <- NULL ## initialise the vector of sample cells
              mask.in.tmp <- mask.in ## define a copy of the sampling mask
              
              ## force to get at least one value of each factorial variable
              fact.level.cells <- bm_SampleFactorLevels(env, mask.out = mask.out, mask.in = mask.in)
              if (length(fact.level.cells)) {
                SR <- c(SR, fact.level.cells)
                mask.in.tmp[SR] <- NA ## update the mask by removing already selected cells
              }

              SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                       size = nb.points - length(SR),
                                       cells = TRUE,
                                       na.rm = TRUE)[, "cell", drop = TRUE])
              ## repeat sampling until having the right number of points
              while (length(SR) < nb.points) {
                mask.in.tmp[SR] <- NA ## update the mask by removing already selected cells
                SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                         size = nb.points - length(SR),
                                         cells = TRUE,
                                         na.rm = TRUE)[, "cell", drop = TRUE])
              }
              pa.tab.tmp[, j] <- SR
            }
              
            # putting cells in good format
            selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
            pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
            colnames(pa.tab) <- paste0("PA", 1:nb.repet)
            for (j in 1:ncol(pa.tab)) {
              pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
            }
            
            # putting presences, true absences and pseudo absences together
            xy <- rbind(coordinates(sp)[which(!is.na(as.vector(sp@data))), ],
                        xyFromCell(mask.in, selected.cells))
            xy <- .add_pa_rownames(xy)
            sp <- as.numeric(unlist(c(na.omit(as.vector(sp@data)), rep(NA, length(selected.cells))), use.names = FALSE))
            env <- extract(env, xy)
            pa.tab <- rbind(matrix(TRUE, nrow = (nrow(xy) - length(selected.cells)),
                                   ncol = ncol(pa.tab)), pa.tab)
            
            return(list(xy = xy,
                        sp = sp,
                        env = as.data.frame(env),
                        pa.tab = as.data.frame(pa.tab)))
          })


###################################################################################################

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_disk",
           def = function(sp, env, ...) {
             standardGeneric("bm_PseudoAbsences_disk")
           })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_disk', signature(env = "SpatialPointsDataFrame"),
          function(sp, env, distMin, distMax, nb.points, nb.repet)
            {
            cat("\n   > Disk pseudo absences selection")
            
            # 1. determining area which can be selected
            coor <- coordinates(sp)
            pres <- which(sp@data[, 1] == 1)
            true.abs <- which(sp@data[, 1] == 0)
            tmp.abs <- which(is.na(sp@data[, 1]))
            outside <- rep(0, length(abs))
            inside <- rep(0, length(abs))
            
            for (i in 1:length(pres)) {
              # removing points too close from presences
              inside <- inside + (sqrt((coor[tmp.abs, 1] - coor[pres[i], 1]) ^ 2 + 
                                         (coor[tmp.abs, 2] - coor[pres[i], 2])^2) > distMin)
              # keeping points not to far from presences
              if (!is.null(distMax)) {
                outside <- outside + (sqrt((coor[tmp.abs, 1] - coor[pres[i], 1]) ^ 2 + 
                                             (coor[tmp.abs, 2] - coor[pres[i], 2]) ^ 2) < distMax )
              }
            }
            if (is.null(distMax)) { # no cells are too far
              outside <- outside + 1
            }
            selected.abs <- tmp.abs[(inside == length(pres)) & (outside > 0)]
            
            # 2. adding presences and true absences and selecting randomly pseudo absences
            return(bm_PseudoAbsences_random(sp[c(pres, true.abs, selected.abs), ],
                                               env[c(pres, true.abs, selected.abs), ],
                                               nb.points, nb.repet))
          })

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_disk', signature(env = "RasterStack"),
          function(sp, env, distMin, distMax, nb.points, nb.repet)
          {
            cat("\n   > Disk pseudo absences selection")
            
            # 1. Check if NA are present in sp observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(sp)
            if (nb.cells > 0) { # PA will be taken into response variable
              env.tmp <- SpatialPointsDataFrame(coords = coordinates(sp),
                                                data = as.data.frame(extract(env, coordinates(sp))))
              return(bm_PseudoAbsences_disk(sp, env.tmp, distMin, distMax, nb.points, nb.repet))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              
              # create a mask
              dist.mask <- subset(env, 1, drop = TRUE)
              dist.mask[] <- NA
              
              pres.xy <- coordinates(sp[which(sp@data[, 1] == 1), ])
              dist.mask[cellFromXY(dist.mask, pres.xy)] <- 1
              
              dist.mask <- distance(dist.mask)
              dist.mask <- mask(dist.mask, subset(env, 1, drop = TRUE))

              if (is.null(distMax)) { distMax <- Inf }
              mask.in <- reclassify(dist.mask, c(-Inf, distMin, NA, distMin, distMax, -1, distMax, Inf, NA))

              # 2. selecting randomly pseudo absences
              return(bm_PseudoAbsences_random(sp, env = mask.in, nb.points, nb.repet))
            }
          })

