# bm_PseudoAbsences doc ---------------------------------------------------------
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
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables (in 
##' columns or layers) that will be used to find the pseudo-absences
##' @param \ldots (\emph{optional, one or several of the following arguments depending on the selected 
##' method}) 
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
##' @examples 
##' 
##' library(terra)
##' # Load species occurrences (6 species available)
##' data(DataSpecies)
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
##' data(bioclim_current)
##' myExpl <- terra::rast(bioclim_current)
##' 
##' \dontshow{
##' myExtent <- terra::ext(0,30,45,70)
##' myExpl <- terra::crop(myExpl, myExtent)
##' }
##' 
##' 
##' # --------------------------------------------------------------- #
##' # Create the different pseudo-absence datasets
##' 
##' # Transform true absences into potential pseudo-absences
##' myResp.PA <- ifelse(myResp == 1, 1, NA)
##' myResp.PA.vect <- vect(cbind(myRespXY, myResp.PA), geom = c("X_WGS84","Y_WGS84"))
##' 
##' # random method
##' PA.r <- bm_PseudoAbsences(resp.var = myResp.PA,
##'                           expl.var = myExpl,
##'                           nb.rep = 4,
##'                           nb.absences = 1000,
##'                           strategy = 'random')
##' 
##' # disk method
##' PA.d <- bm_PseudoAbsences(resp.var = myResp.PA,
##'                           expl.var = myExpl,
##'                           nb.rep = 4,
##'                           nb.absences = 500,
##'                           strategy = 'disk',
##'                           dist.min = 5,
##'                           dist.max = 35)
##' 
##' # SRE method
##' PA.s <- bm_PseudoAbsences(resp.var = myResp.PA,
##'                           expl.var = myExpl,
##'                           nb.rep = 4,
##'                           nb.absences = 1000,
##'                           strategy = 'sre',
##'                           sre.quant = 0.025)
##' 
##' # user.defined method
##' myPAtable <- data.frame(PA1 = ifelse(myResp == 1, TRUE, FALSE),
##'                         PA2 = ifelse(myResp == 1, TRUE, FALSE))
##' for (i in 1:ncol(myPAtable)) myPAtable[sample(which(myPAtable[, i] == FALSE), 500), i] = TRUE
##' PA.u <- bm_PseudoAbsences(resp.var = myResp.PA,
##'                           expl.var = myExpl,
##'                           strategy = 'user.defined',
##'                           user.table = myPAtable)
##' 
##' str(PA.r)
##' head(PA.r$pa.tab)
##' apply(PA.r$pa.tab, 2, table)
##' 
##' head(PA.d$pa.tab)
##' apply(PA.d$pa.tab, 2, table)
##' 
##' head(PA.s$pa.tab)
##' apply(PA.s$pa.tab, 2, table)
##' 
##' tail(PA.u$pa.tab)
##' apply(PA.u$pa.tab, 2, table)
##' 
##' 
##' # random method : different number of PA
##' PA.r_mult <- bm_PseudoAbsences(resp.var = myResp.PA.vect,
##'                                expl.var = myExpl,
##'                                nb.rep = 4,
##'                                nb.absences = c(1000, 500, 500, 200),
##'                                strategy = 'random')
##' 
##' str(PA.r_mult)
##' head(PA.r_mult$pa.tab)
##' apply(PA.r_mult$pa.tab, 2, table)
##'
##' 
##' @importFrom terra rast vect freq spatSample values extract
##' @importFrom utils packageVersion
##'
##' @export
##' 
##' 
## --------------------------------------------------------------------------- #


bm_PseudoAbsences <- function(resp.var, expl.var, nb.rep = 1, strategy = 'random', nb.absences = NULL
                              , sre.quant = 0, dist.min = 0, dist.max = NULL, user.table = NULL)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PseudoAbsences.check.args(resp.var, expl.var, nb.rep, strategy, nb.absences
                                        , sre.quant, dist.min, dist.max, user.table)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  if ((nb.rep == 0 || any(nb.absences <= 0)) & strategy != 'user.defined') {
    out <- NULL
  } else {
    if (length(nb.absences) == 1) {
      out <- switch(strategy,
                    user.defined = bm_PseudoAbsences_user.defined(resp.var, expl.var, user.table),
                    random = bm_PseudoAbsences_random(resp.var, expl.var, nb.absences, nb.rep),
                    sre = bm_PseudoAbsences_sre(resp.var, expl.var, sre.quant, nb.absences, nb.rep),
                    disk = bm_PseudoAbsences_disk(resp.var, expl.var, dist.min, dist.max, nb.absences, nb.rep))
    } else if (length(nb.absences) == nb.rep) {
      out.list = foreach(i.abs = unique(nb.absences)) %do% 
        {
          i.rep = which(nb.absences == i.abs)
          cat("\n > Set ", paste0(i.rep, collapse = ", "), " (", i.abs, " pseudo absences wanted)", sep = "")
          
          out <- switch(strategy,
                        user.defined = bm_PseudoAbsences_user.defined(resp.var, expl.var, user.table),
                        random = bm_PseudoAbsences_random(resp.var, expl.var, i.abs, length(i.rep)),
                        sre = bm_PseudoAbsences_sre(resp.var, expl.var, sre.quant, i.abs, length(i.rep)),
                        disk = bm_PseudoAbsences_disk(resp.var, expl.var, dist.min, dist.max, i.abs, length(i.rep)))
          
          ## CASE where all available cells have been selected :
          ## give back only one dataset, even if several were asked
          if (ncol(out$pa.tab) == 1 && length(i.rep) > 1) {
            col1 <- out$pa.tab
            for (j in 2:length(i.rep)) {
              out$pa.tab <- cbind(out$pa.tab, col1)
            }
            # i.rep = i.rep[1] ## NOT working with summary and plot functions
          }
          
          colnames(out$pa.tab) = paste0("PA", i.rep)
          return(out)
        }
      
      ## GET XY -----------------------------------------------------
      
      ## Get coordinates of presences
      ind.pres = which(out.list[[1]][["sp"]] == 1)
      xy.pres = out.list[[1]][["xy"]][ind.pres, ]
      nb.pres = length(ind.pres)
      
      ## Get coordinates of pseudo-absences
      out.xy = foreach(i = 1:length(out.list)) %do%
        {
          ind.keep = 1:nrow(out.list[[i]][["xy"]])
          ind.keep = ind.keep[-which(ind.keep %in% ind.pres)]
          res = out.list[[i]][["xy"]][ind.keep, ]
          res = cbind(res, ind.keep)
          return(res)
        }
      
      ## Merge all coordinates of pseudo-absences (may be duplicates)
      out.xy = Reduce(function(x, y) merge(x, y, by = c("x", "y"), all = TRUE), out.xy)
      
      ## Get indexes of merged PA coordinates for each set
      ## To be used to rebuild env and pa.tab 
      out.index = out.xy[, -which(colnames(out.xy) %in% c("x", "y"))]
      out.order = !is.na(out.index)
      out.order = t(apply(out.order, 1, cumsum))
      
      ## Keep presences + pseudo-absences coordinates
      out.xy = out.xy[, c("x", "y")]
      out.xy = rbind(xy.pres, out.xy)
      out.sp = c(rep(1, nb.pres), rep(NA, nrow(out.xy) - nb.pres))
      
      ## GET ENV & PA.TAB -------------------------------------------
      
      ## Initialize env matrix
      out.env = matrix(NA, nrow = nrow(out.xy), ncol = ncol(out.list[[1]][["env"]]))
      out.env = as.data.frame(out.env)
      colnames(out.env) = colnames(out.list[[1]][["env"]])
      out.env[1:nb.pres, ] = out.list[[1]][["env"]][1:nb.pres, ]
      
      ## Initialize pa.tab matrix
      out.pa.tab = matrix(NA, nrow = nrow(out.xy), ncol = nb.rep)
      out.pa.tab = as.data.frame(out.pa.tab)
      colnames(out.pa.tab) = paste0("PA", 1:nb.rep)
      out.pa.tab[1:nb.pres, ] = TRUE
      
      ind.start = 1
      ind.end = ncol(out.list[[1]][["pa.tab"]])
      
      ## Fill first column
      ind = which(out.order[, 1] == 1)
      out.env[ind + nb.pres, ] = out.list[[1]][["env"]][out.index[ind, 1], ]
      out.pa.tab[ind + nb.pres, ind.start:ind.end] = out.list[[1]][["pa.tab"]][out.index[ind, 1], ]
      
      ## Fill all other columns
      for (j in 2:ncol(out.order)) {
        ind = which(out.order[, j] != out.order[, j-1])
        if (length(ind) > 0) {
          
          ## For env
          out.env[ind + nb.pres, ] = out.list[[j]][["env"]][out.index[ind, j], ]
          
          ## For pa.tab
          ind.start = ind.end + 1
          ind.end = ind.start + ncol(out.list[[j]][["pa.tab"]]) - 1
          out.pa.tab[ind + nb.pres, ind.start:ind.end] = out.list[[j]][["pa.tab"]][out.index[ind, j], ]
        }
      }
      
      ## GET everything ---------------------------------------------
      out = list(xy = out.xy,
                 sp = out.sp,
                 env = out.env,
                 pa.tab = out.pa.tab)
    }
  }
  cat("\n")
  return(out)
}


# Argument Check --------------------------------------------------------------

.bm_PseudoAbsences.check.args <- function(resp.var, expl.var, nb.rep, strategy, nb.absences, sre.quant, dist.min, dist.max, user.table)
{
  cat('\n\nChecking Pseudo-absence selection arguments...\n')
  ## 1. Check resp.var argument -----------------------------------------------
  if (is.vector(resp.var)) {
    resp.var <- vect(data.frame(x = 0,
                                y = 0,
                                resp = resp.var),
                     geom = c("x","y"))
    if (!is.null(nb.absences) && length(nb.absences) > 1) {
      stop("Selection of multiple number of pseudo-absences depends on coordinates. Please provide some.")
    }
  }
  .fun_testIfInherits(TRUE, "resp.var", resp.var, "SpatVector")
  
  ## 2. Check expl.var argument -----------------------------------------------
  if (inherits(expl.var, c("matrix","data.frame"))) {
    if (nrow(expl.var) != nrow(resp.var)) {
      stop("Species and Explanatory must have same dimensions")
    }
    # transform expl.var into SpatVector
    rownames(expl.var) <- NULL
    expl.var <- vect(
      cbind(data.frame(x = crds(resp.var)[,1],
                       y = crds(resp.var)[,2]),
            as.data.frame(expl.var)
      ),
      geom = c("x","y")
    )
  }
  .fun_testIfInherits(TRUE, "expl.var", expl.var, c("SpatVector", "SpatRaster"))
  
  args <- .BIOMOD.formated.data.check.args(sp = resp.var, env = expl.var)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 3. Check strategy argument -----------------------------------------------
  availableStrategies <- c("random", "sre", "disk", "user.defined")
  if (is.null(strategy) || !(strategy %in% availableStrategies) || 
      (strategy != 'user.defined' && all(crds(resp.var) == 0))) {
    # no coordinates or unknown strategy
    strategy <- "random"
    cat("\n   ! Random strategy was automatically selected (that can be due to points coordinates lack or unavailable strategy choosen)")
  }
  
  ## 4. Check sre.quant argument ----------------------------------------------
  if (strategy == 'SRE' && (sre.quant >= 0.5 || sre.quant < 0)) {
    stop("\n    ! SRE Quant should be a value between 0 and 0.5 ")
  }
  
  ## 5. Check dist.min and dist.max arguments ---------------------------------
  if (strategy == 'disk') {
    if (!is.null(dist.min) && dist.min < 0) {
      dist.min <- 0
    }
    if (!is.null(dist.max) && dist.max < 0) {
      dist.max <- NULL
    }
    if (!is.null(dist.max) && !is.null(dist.min) && dist.min >= dist.max) {
      stop("dist.min >= dist.max")
    }
  }
  
  ## 6. Check nb.absences argument --------------------------------------------
  if (strategy != "user.defined") {
    if (is.null(nb.absences)) {
      stop("You must give the number of pseudo absences you want")
    } else {
      if (length(nb.absences) > 1) {
        if (length(nb.absences) != nb.rep) {
          stop("You must give one value for pseudo absences, or as many as the number of repetitions")
        } else if (length(unique(nb.absences)) == 1) {
          nb.absences = unique(nb.absences)
        }
      }
      nbTrueAbs <- .get_nb_true_abs(resp.var)
      if (nbTrueAbs) {
        stop("Your dataset contains true absences. This should not be mixed with pseudo absences selection")
      }
    }
  }
  
  ## 7. Check user.table argument --------------------------------------------
  if (strategy == "user.defined") {
    if (is.null(user.table)) {
      stop("You must give a table defining the pseudo absences you want")
    } else {
      if (!(is.matrix(user.table) | is.data.frame(user.table))) {
        stop("\n PA.user.table must be a matrix or a data.frame")
      }
      if (nrow(user.table) != length(resp.var)) {
        stop("\n PA.user.table must have as many row than the number of observation of your response variable")
      }
      colnames(user.table) <- paste0("PA", 1:ncol(user.table))
      nb.absences <- nrow(user.table)
    }
  }
  
  return(list(resp.var = resp.var,
              expl.var = expl.var,
              nb.rep = nb.rep,
              strategy = strategy,
              nb.absences = nb.absences,
              sre.quant = sre.quant,
              dist.min = dist.min,
              dist.max = dist.max,
              user.table = user.table))
}

# Additionnal tools ------------------------------------------------------------

.get_nb_true_abs <- function(sp)
{
  if (is.vector(sp)) {
    return(length(which(sp == 0)))
  } else if (inherits(sp, 'SpatVector')) {
    return(length(which(values(sp)[,1] == 0)))
  } else if (inherits(sp, 'SpatRaster')) {
    return(freq(sp, value = 0)$count)
  }
}

.get_nb_available_pa_cells <- function(data, PA.flag = NA)
{
  if (is.vector(data) || inherits(data, c("matrix","data.frame"))) {
    return(
      ifelse(is.na(PA.flag),
             length(which(is.na(data))), 
             length(which(data == PA.flag)))
    )
  } else if (inherits(data, 'SpatVector')) {
    return(
      ifelse(is.na(PA.flag), 
             length(which(is.na(as.data.frame(data)[,1]))), 
             length(which(values(data)[,1] == PA.flag)))
    )
  } else if (inherits(data, 'SpatRaster')) {
    return(
      freq(data, value = PA.flag)$count
    )
  }
}

.add_pa_rownames <- function(xy)
{
  rn <- row.names(xy)
  # if (length(rn) == 0) {
  #   rn = rep("", nrow(xy))
  # }
  missing_rn <- which(rn == "")
  if (length(missing_rn) > 0) {
    rn[missing_rn] <- paste0("pa", 1:length(missing_rn))
  }
  rownames(xy) <- rn
  return(xy)
}


# bm_PseudoAbsences user-defined methods --------------------------------------

##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_user.defined",
           def = function(resp.var, expl.var, ...) {
             standardGeneric( "bm_PseudoAbsences_user.defined")
           })

## bm_PseudoAbsences user-defined SpatVector methods --------------------------------------
##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_user.defined', signature(expl.var = "SpatVector"),
          function(resp.var, expl.var, user.table) {
            cat("\n   > User defined pseudo absences selection")
            return(list(xy = crds(resp.var),
                        sp = as.numeric(unlist(values(resp.var), use.names = FALSE)),
                        env = as.data.frame(expl.var),
                        pa.tab = user.table))
          })

## bm_PseudoAbsences user-defined SpatRaster methods --------------------------------------
##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_user.defined', signature(expl.var = "SpatRaster"), 
          function(resp.var, expl.var, user.table) {
            cat("\n   > User defined pseudo absences selection")
            expl.var <- extract(expl.var, resp.var, ID = FALSE)
            return(list(xy = crds(resp.var),
                        sp = as.numeric(unlist(values(resp.var), use.names = FALSE)), 
                        env = expl.var,
                        pa.tab = user.table))
          })


# bm_PseudoAbsences random methods --------------------------------------


##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_random",
           def = function(resp.var, expl.var, ...) {
             standardGeneric( "bm_PseudoAbsences_random")
           })

## bm_PseudoAbsences random SPDF methods --------------------------------------
##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_random', signature(expl.var = "SpatVector"),
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
              pa.tab[c(which(values(resp.var)[, 1] == 1),
                       which(values(resp.var)[, 1] == 0)),] <- TRUE
              
              # 4. For each repetition, select among NA cells
              cand.cells <- which(is.na(values(resp.var)[, 1]))
              for (j in 1:ncol(pa.tab)) {
                ## force to get at least one value of each factorial variable
                fact.level.cells <- bm_SampleFactorLevels(expl.var = as.data.frame(expl.var),
                                                          mask.out = pa.tab[, j, drop = FALSE])
                if (length(fact.level.cells) > 0) {
                  pa.tab[fact.level.cells, j] <- TRUE
                  cand.cells <- setdiff(cand.cells, fact.level.cells)
                }
                pa.tab[sample(x = cand.cells,
                              size = nb.absences - length(fact.level.cells),
                              replace = FALSE), j] <- TRUE
              }
              
              return(list(xy = crds(resp.var),
                          sp = as.numeric(unlist(values(resp.var), use.names = FALSE)),
                          env = as.data.frame(expl.var),
                          pa.tab = pa.tab))
            } else {
              cat("\nUnsupported case yet!")
              return(NULL)
            }
          })

## bm_PseudoAbsences random SpatRaster methods --------------------------------------
##'
##' @rdname bm_PseudoAbsences
##' @export
##' @importFrom terra rast xyFromCell cellFromXY classify crds extract subset
##'

setMethod('bm_PseudoAbsences_random', signature(expl.var = "SpatRaster"),
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
              pa.tab[c(which(values(resp.var)[, 1] == 1),
                       which(values(resp.var)[, 1] == 0)),] <- TRUE
              
              # 4. For each repetition, select among NA cells
              cand.cells <- which(is.na(values(resp.var)[, 1]))
              for (j in 1:ncol(pa.tab)) {
                pa.tab[sample(x = cand.cells, size = nb.absences, replace = FALSE), j] <- TRUE
              }
              
              expl.var <- extract(expl.var, resp.var, ID = FALSE)
              return(list(xy = crds(resp.var),
                          sp = as.numeric(unlist(values(resp.var), use.names = FALSE)), 
                          env = as.data.frame(expl.var),
                          pa.tab = as.data.frame(pa.tab)))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              # create a mask containing all not already sampled points (presences and absences)
              mask.env <- mask.out <- subset(expl.var, 1)
              ## the area we want to sample
              mask.env <- classify(mask.env, matrix(c(-Inf, Inf, -1), 
                                                    ncol = 3, byrow = TRUE)) 
              mask.out[] <- NA
              
              # add presences and true absences in our mask
              in.cells <- cellFromXY(mask.env, crds(resp.var))
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
                  if (length(fact.level.cells) > 0) {
                    SR <- c(SR, fact.level.cells)
                    mask.env.tmp[SR] <- NA ## update the mask by removing already selected cells
                  }
                  ## repeat sampling until having the right number of points
                  ## spatSample with na.rm = TRUE may return less points than asked
                  while(length(SR) < nb.absences){
                    SR <- c(SR, spatSample(x = mask.env.tmp,
                                           method = "random",
                                           size = nb.absences - length(SR),
                                           cells = TRUE,
                                           na.rm = TRUE)[, "cell", drop = TRUE])
                    mask.env.tmp[SR] <- NA ## update the mask by removing already selected cells
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
                xy = crds(resp.var)
                rownames(xy) = paste0("pres", 1:nrow(xy))
                xy <- rbind(xy, xyFromCell(mask.env, selected.cells))
                xy <- .add_pa_rownames(xy)
                resp.var <- as.numeric(unlist(c(values(resp.var)[, 1], 
                                                rep(NA, length(selected.cells))),
                                              use.names = FALSE))
                expl.var <- extract(expl.var, xy)
                pa.tab <- rbind(matrix(TRUE, 
                                       nrow = (nrow(xy) - length(selected.cells)),
                                       ncol = ncol(pa.tab)), 
                                pa.tab)
                
                return(list(xy = xy,
                            sp = resp.var,
                            env = as.data.frame(expl.var),
                            pa.tab = as.data.frame(pa.tab)))
              }
            }
          })


# bm_PseudoAbsences SRE methods --------------------------------------


##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_sre",
           def = function(resp.var, expl.var, ...) {
             standardGeneric("bm_PseudoAbsences_sre")
           })

## bm_PseudoAbsences SRE SpatVector methods ----------------------------------
##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_sre', signature(expl.var = "SpatVector"), 
          function(resp.var, expl.var, sre.quant, nb.absences, nb.rep)
          {
            cat("\n   > SRE pseudo absences selection")
            
            # 0. calculate SRE to determine available pixels
            mask.in <- bm_SRE(resp.var = resp.var, expl.var = expl.var, new.env = values(expl.var), quant = sre.quant)
            mask.in <- data.frame(mask.in = !as.logical(mask.in)) ## revert the mask to sample PA out of SRE
            mask.in$mask.in[c(which(values(resp.var)[, 1] == 1)
                              , which(values(resp.var)[, 1] == 0))] <- FALSE
            
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
            pa.tab[c(which(values(resp.var)[, 1] == 1),
                     which(values(resp.var)[, 1] == 0)),] <- TRUE
            
            # 4. For each repetition, select among NA cells
            cand.cells <- which(mask.in$mask.in == TRUE)
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
            
            return(list(xy = crds(resp.var),
                        sp = as.numeric(unlist(values(resp.var), use.names = FALSE)),
                        env = values(expl.var),
                        pa.tab = pa.tab))
          })

## bm_PseudoAbsences SRE SpatRaster methods -----------------------------------
##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_sre', signature(expl.var = "SpatRaster"), 
          function(resp.var, expl.var, sre.quant, nb.absences, nb.rep)
          {
            cat("\n   > SRE pseudo absences selection")
            
            # 0. calculate SRE to determine available pixels
            mask.in <- bm_SRE(resp.var = resp.var, expl.var = expl.var, new.env = expl.var, quant = sre.quant)
            mask.in[mask.in[] > 0] <- NA ## remove points that are in SRE
            mask.in[cellFromXY(mask.in, crds(resp.var)[which(!is.na(values(resp.var)[, 1])), ])] <- NA
            
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
              if (length(fact.level.cells) > 0) {
                SR <- c(SR, fact.level.cells)
                mask.in.tmp[SR] <- NA ## update the mask by removing already selected cells
              }
              
              ## repeat sampling until having the right number of points
              while (length(SR) < nb.absences) {
                SR <- c(SR, spatSample(x = mask.in.tmp,
                                       method = "random",
                                       size = nb.absences - length(SR),
                                       cells = TRUE,
                                       na.rm = TRUE)[, "cell", drop = TRUE])
                mask.in.tmp[SR] <- NA ## update the mask by removing already selected cells
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
            xy <- rbind(crds(resp.var)[which(!is.na(values(resp.var)[, 1])), ],
                        xyFromCell(mask.in, selected.cells))
            xy <- .add_pa_rownames(xy)
            resp.var <- as.numeric(unlist(c(na.omit(values(resp.var)[, 1]),
                                            rep(NA, length(selected.cells))),
                                          use.names = FALSE))
            expl.var <- extract(expl.var, xy)
            pa.tab <- rbind(matrix(TRUE, nrow = (nrow(xy) - length(selected.cells)),
                                   ncol = ncol(pa.tab)), pa.tab)
            
            return(list(xy = xy,
                        sp = resp.var,
                        env = as.data.frame(expl.var),
                        pa.tab = as.data.frame(pa.tab)))
          })


# bm_PseudoAbsences disk methods --------------------------------------


##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setGeneric("bm_PseudoAbsences_disk",
           def = function(resp.var, expl.var, ...) {
             standardGeneric("bm_PseudoAbsences_disk")
           })

## bm_PseudoAbsences disk SpatVector methods --------------------------------------
##'
##' @rdname bm_PseudoAbsences
##' @export
##'

setMethod('bm_PseudoAbsences_disk', signature(expl.var = "SpatVector"),
          function(resp.var, expl.var, dist.min, dist.max, nb.absences, nb.rep) {
            cat("\n   > Disk pseudo absences selection")
            # 1. determining area which can be selected
            coor <- crds(resp.var)
            pres <- which(values(resp.var)[, 1] == 1)
            true.abs <- which(values(resp.var)[, 1] == 0)
            tmp.abs <- which(is.na(values(resp.var)[, 1]))
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

## bm_PseudoAbsences disk SpatRaster methods --------------------------------------
##'
##' @rdname bm_PseudoAbsences
##' @export
##' @importFrom terra rast distance subset mask crds cellFromXY extract
##'

setMethod('bm_PseudoAbsences_disk', signature(expl.var = "SpatRaster"),
          function(resp.var, expl.var, dist.min, dist.max, nb.absences, nb.rep)
          {
            cat("\n   > Disk pseudo absences selection")
            
            # 1. Check if NA are present in resp.var observations or not to determine which dataset to use
            nb.cells <- .get_nb_available_pa_cells(resp.var)
            if (nb.cells > 0) { # PA will be taken into response variable
              expl.var.tmp <- extract(expl.var, resp.var[,-1], bind = TRUE, ID = FALSE)
              return(
                bm_PseudoAbsences_disk(resp.var, expl.var.tmp, 
                                       dist.min, dist.max,
                                       nb.absences, nb.rep)
              )
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              
              # create a mask
              dist.mask <- subset(expl.var, 1)
              dist.mask[] <- NA
              
              pres.xy <- crds(resp.var[which(values(resp.var)[, 1] == 1), ])
              dist.mask[cellFromXY(dist.mask, pres.xy)] <- 1
              
              dist.mask <- distance(dist.mask)
              dist.mask <- mask(dist.mask, subset(expl.var, 1))
              
              if (is.null(dist.max)) { 
                dist.max <- Inf 
                }
              mask.in <- classify(dist.mask, 
                                  matrix(c(-Inf, dist.min, NA,
                                               dist.min, dist.max, 1, 
                                               dist.max, Inf, NA),
                                         ncol = 3, byrow = TRUE))
              mask.in[cellFromXY(mask.in, pres.xy)] <- 1
              mask.in = expl.var * mask.in
              names(mask.in) = names(expl.var)

              # 2. selecting randomly pseudo absences
              return(bm_PseudoAbsences_random(resp.var, expl.var = mask.in, nb.absences, nb.rep))
            }
          })

