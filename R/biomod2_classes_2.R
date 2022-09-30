
##' @name BIOMOD.stored.data
##' @aliases BIOMOD.stored.data-class
##' @aliases BIOMOD.stored.array-class
##' @aliases BIOMOD.stored.data.frame-class
##' @aliases BIOMOD.stored.raster.stack-class
##' @aliases BIOMOD.stored.files-class
##' @aliases BIOMOD.stored.formated.data-class
##' @aliases BIOMOD.stored.models.options-class
##' @aliases BIOMOD.stored.models.out-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_EnsembleModeling()} output object class
##' 
##' @description Classes used by \code{\link{BIOMOD_Modeling}} and 
##' \code{\link{BIOMOD_EnsembleModeling}} to build their output object (see 
##' \code{\link{BIOMOD.models.out}} objects)
##' 
##' 
##' @slot inMemory a \code{logical} defining whether the \code{val} slot has been loaded in 
##' memory or not
##' @slot link a \code{character} containing the file name of the saved \code{val} slot
##' @slot val an object of type depending on the \code{BIOMOD.stored.[...]} class (see Details)
##' 
##' @details 
##' 
##' \code{BIOMOD.stored.data} is the basic object containing the slots \code{inMemory} and 
##' \code{link}. \cr All listed classes below are derived from \code{BIOMOD.stored.data}, and 
##' contain a \code{val} slot of specific type :
##' 
##' \itemize{
##'   \item{\code{BIOMOD.stored.array} : }{\code{val} is an \code{array}}
##'   \item{\code{BIOMOD.stored.data.frame} : }{\code{val} is a \code{data.frame}}
##'   \item{\code{BIOMOD.stored.raster.stack} : }{\code{val} is a 
##'   \code{\link[raster:stack]{RasterStack}}}
##'   \item{\code{BIOMOD.stored.files} : }{\code{val} is a \code{character}}
##'   \item{\code{BIOMOD.stored.formated.data} : }{\code{val} is a 
##'   \code{\link{BIOMOD.formated.data}} object}
##'   \item{\code{BIOMOD.stored.models.options} : }{\code{val} is a 
##'   \code{\link{BIOMOD.models.options}} object}
##'   \item{\code{BIOMOD.stored.models.out} : }{\code{val} is a 
##'   \code{\link{BIOMOD.models.out}} object}
##' }
##' 
##' 
##' @seealso \code{\link{BIOMOD.formated.data}}, \code{\link{BIOMOD.models.options}},
##' \code{\link{BIOMOD.models.out}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleForecasting}}
##' @family Toolbox objects
##' 
##' 
##' @examples 
##' 
##' showClass("BIOMOD.stored.data")
##' showClass("BIOMOD.stored.array") 
##' showClass("BIOMOD.stored.data.frame") 
##' showClass("BIOMOD.stored.raster.stack") 
##' showClass("BIOMOD.stored.files") 
##' showClass("BIOMOD.stored.formated.data") 
##' showClass("BIOMOD.stored.models.options") 
##' showClass("BIOMOD.stored.models.out") 
##' 
##' 
NULL

##' @name BIOMOD.stored.data-class
##' @rdname BIOMOD.stored.data
##' @importFrom raster stack
##' 
##' @export
##' 

## --------------------------------------------------------------------------- #
## BIOMOD.stored[...] objects ------------------------------------------------
## --------------------------------------------------------------------------- #

setClass("BIOMOD.stored.data",
         representation(inMemory = 'logical', link = 'character'),
         prototype(inMemory = FALSE, link = ''),
         validity = function(object){ return(TRUE) })

##' @name BIOMOD.stored.array-class
##' @rdname BIOMOD.stored.data
##' 
setClass("BIOMOD.stored.array",
         contains = "BIOMOD.stored.data",
         representation(val = 'array'),
         prototype(val = array()),
         validity = function(object){ return(TRUE) })

##' @name BIOMOD.stored.data.frame-class
##' @rdname BIOMOD.stored.data
##' 
setClass("BIOMOD.stored.data.frame",
         contains = "BIOMOD.stored.data",
         representation(val = 'data.frame'),
         prototype(val = data.frame()),
         validity = function(object){ return(TRUE) })

##' @name BIOMOD.stored.raster.stack-class
##' @rdname BIOMOD.stored.data
##' 
setClass("BIOMOD.stored.raster.stack",
         contains = "BIOMOD.stored.data",
         representation(val = 'RasterStack'),
         prototype(val = stack()),
         validity = function(object){ return(TRUE) })

##' @name BIOMOD.stored.files-class
##' @rdname BIOMOD.stored.data
##' 
setClass("BIOMOD.stored.files",
         contains = "BIOMOD.stored.data",
         representation(val = 'character'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) })

##' @name BIOMOD.stored.formated.data-class
##' @rdname BIOMOD.stored.data
##' 
setClass("BIOMOD.stored.formated.data",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.formated.data'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) })

##' @name BIOMOD.stored.models.options-class
##' @rdname BIOMOD.stored.data
##' 
setClass("BIOMOD.stored.models.options",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.options'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) })

# setClass("BIOMOD.stored.models.out",
#          contains = "BIOMOD.stored.data",
#          representation(val = 'BIOMOD.models.out'),
#          prototype(val = NULL),
#          validity = function(object){ return(TRUE) } )


###################################################################################################
## LOAD BIOMOD.stored[...] objects
###################################################################################################

##' @name load_stored_object
##' @author Damien Georges
##' 
##' @title Functions to load \code{\link{BIOMOD.stored.data}} objects
##' 
##' @description This functions allow the user to load \code{\link{BIOMOD.stored.data}} objects 
##' into memory.
##' 
##' @param obj a \code{\link{BIOMOD.stored.data}} object
##' 
##' @seealso \code{\link{BIOMOD.stored.data}}
##' @family Toolbox functions
##' 
##' @export
##' 

setGeneric("load_stored_object", function(obj) { standardGeneric("load_stored_object") })

##' 
##' @rdname load_stored_object
##' @export
##' 

setMethod("load_stored_object", "BIOMOD.stored.data",
          function(obj)
          {
            if(obj@inMemory){ return(obj@val) }
            
            # different comportement with raster
            if (inherits(obj, "BIOMOD.stored.raster.stack")) {
              if (length(obj@link) == 1 & all(grepl(".RData", obj@link))) {
                return(get(load(obj@link)))
              } else if (all(grepl(".grd", obj@link) | grepl(".img", obj@link))) {
                out <- stack(x = obj@link, RAT = FALSE)
                ## rename layer in case of individual projections
                if (all(grepl("individual_projections", obj@link))) {
                  # remove directories arch and extention
                  xx <- sub("[:.:].+$", "", sub("^.+individual_projections/", "", obj@link))
                  # remove projection name
                  to_rm <- unique(sub("[^_]+[:_:][^_]+[:_:][^_]+[:_:][^_]+$", "", xx))
                  xx <- sub(to_rm, "", xx)
                  names(out) <- xx
                }
                return(out)
              }
            } else { # for all other stored objects
              if (obj@link != '') {
                return(get(load(obj@link)))
              } else {
                warning("No link provided for this object")
                return(NA)
              }
            }
          }
)

