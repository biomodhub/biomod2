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



###################################################################################################
## BIOMOD.stored[...] objects
###################################################################################################

setClass("BIOMOD.stored.data",
         representation(inMemory = 'logical', link = 'character'),
         prototype(inMemory = FALSE, link = ''),
         validity = function(object){ return(TRUE) })

setClass("BIOMOD.stored.array",
         contains = "BIOMOD.stored.data",
         representation(val = 'array'),
         prototype(val = array()),
         validity = function(object){ return(TRUE) })

setClass("BIOMOD.stored.data.frame",
         contains = "BIOMOD.stored.data",
         representation(val = 'data.frame'),
         prototype(val = data.frame()),
         validity = function(object){ return(TRUE) })

setClass("BIOMOD.stored.raster.stack",
         contains = "BIOMOD.stored.data",
         representation(val = 'RasterStack'),
         prototype(val = stack()),
         validity = function(object){ return(TRUE) })

setClass("BIOMOD.stored.files",
         contains = "BIOMOD.stored.data",
         representation(val = 'character'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) })

setClass("BIOMOD.stored.formated.data",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.formated.data'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) })

setClass("BIOMOD.stored.models.options",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.Model.Options'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) })

setClass("BIOMOD.stored.models.out",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.out'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) } )


###################################################################################################
## LOAD BIOMOD.stored[...] objects
###################################################################################################

setGeneric("load_stored_object", function(obj, ...) { standardGeneric("load_stored_object") })

setMethod("load_stored_object", "BIOMOD.stored.data",
          function(obj)
          {
            if(obj@inMemory){ return(obj@val) }
            
            # different comportement with raster
            if (inherits(obj, "BIOMOD.stored.raster.stack")) {
              if (length(obj@link) == 1 & all(grepl(".RData", obj@link))) {
                return(get(load(obj@link)))
              } else if (all(grepl(".grd", obj@link) | grepl(".img", obj@link))) {
                out <- raster::stack(x = obj@link, RAT = FALSE)
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
              return(get(load(obj@link)))
            }
          }
)

