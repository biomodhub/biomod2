# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# BIOMOD objects definition
# Damien Georges
# 09/02/2012
# v2.0
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# This file defines the BIOMOD objects and all their methods 
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# We choose here to create monospecific objects to make all procedures and parallelising easier
requireNamespace("sp", quietly=TRUE)
requireNamespace("raster", quietly=TRUE)
requireNamespace("rasterVis", quietly=TRUE)

# 0. Generic Functions definition -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
setGeneric("get_predictions",
           function(obj, ...){
             standardGeneric("get_predictions")
           })

setGeneric("get_projected_models",
           function(obj, ...){
             standardGeneric("get_projected_models")
           })

setGeneric("get_evaluations",
           function(obj, ...){
             standardGeneric("get_evaluations")
           })

setGeneric("get_calib_lines",
           function(obj, ...){
             standardGeneric("get_calib_lines")
           })

setGeneric("get_variables_importance",
           function(obj, ...){
             standardGeneric("get_variables_importance")
           })

setGeneric("get_options",
           function(obj, ...){
             standardGeneric("get_options")
           })

setGeneric("get_formal_data",
           function(obj, ...){
             standardGeneric("get_formal_data")
           })

setGeneric("get_built_models",
           function(obj, ...){
             standardGeneric("get_built_models")
           })

setGeneric("get_needed_models",
           function(obj, ...){
             standardGeneric("get_needed_models")
           })

setGeneric("get_kept_models",
           function(obj, ...){
             standardGeneric("get_kept_models")
           })

setGeneric("load_stored_object",
           function(obj, ...){
             standardGeneric("load_stored_object")
           })

setGeneric("RemoveProperly",
           function(obj, obj.name=deparse(substitute(obj)), ...){
             standardGeneric("RemoveProperly")
           })

setGeneric("free",
           function(obj, ...){
             standardGeneric("free")
           })

setGeneric( ".Models.prepare.data", 
            def = function(data, ...){
              standardGeneric( ".Models.prepare.data" )
            } )


# 1. The BIOMOD.formated.data -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
# this object is the basic one

# 1.1 Class Definition
setClass("BIOMOD.formated.data",
         representation(sp.name = 'character',
                        coord = "data.frame",
                        data.species = "numeric",
                        #                         data.counting = "matrix",
                        data.env.var = "data.frame",
                        data.mask = "RasterStack",
                        has.data.eval = "logical",
                        eval.coord = "data.frame",
                        eval.data.species = "numeric",
                        eval.data.env.var = "data.frame"),
         validity = function(object){ return(TRUE) } )

# 1.2 Constructors
# if( !isGeneric( "BIOMOD.formated.data" ) ) {
setGeneric( "BIOMOD.formated.data", 
            def = function(sp, env, ...){
              standardGeneric( "BIOMOD.formated.data" )
            } )
# }

setMethod('BIOMOD.formated.data', signature(sp='numeric', env='data.frame' ), 
          function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE, data.mask=NULL ){
            if(is.null(data.mask)) data.mask <- raster::stack()
            
            if(is.null(eval.sp)){
              BFD <- new('BIOMOD.formated.data', 
                         coord=xy, 
                         data.species=sp, 
                         data.env.var=env, 
                         sp.name=sp.name,
                         data.mask=data.mask,
                         has.data.eval=FALSE)
            } else{
              BFDeval <- BIOMOD.formated.data(sp=eval.sp,
                                              env=eval.env,
                                              xy=eval.xy,
                                              sp.name=sp.name)
              
              if(raster::nlayers(BFDeval@data.mask)>0){
                data.mask.tmp <- try(raster::addLayer(data.mask,BFDeval@data.mask))
                if( !inherits(data.mask.tmp,"try-error")){
                  data.mask <- data.mask.tmp
                  names(data.mask) <- c("calibration","validation")
                }
              }
              
              BFD <- new('BIOMOD.formated.data', 
                         coord=xy, 
                         data.species=sp, 
                         data.env.var=env, 
                         sp.name=sp.name,
                         data.mask=data.mask,
                         has.data.eval=TRUE,
                         eval.coord = BFDeval@coord,
                         eval.data.species = BFDeval@data.species,
                         eval.data.env.var = BFDeval@data.env.var )
              
              
              rm('BFDeval')
            }
            if(na.rm){
              rowToRm <- unique(unlist(lapply(BFD@data.env.var,function(x){return(which(is.na(x)))})))
              if(length(rowToRm)){
                cat("\n\t\t\t! Some NAs have been automaticly removed from your data")
                BFD@coord <- BFD@coord[-rowToRm,,drop=FALSE]
                BFD@data.species <- BFD@data.species[-rowToRm]
                BFD@data.env.var <- BFD@data.env.var[-rowToRm,,drop=FALSE]
              }
              if(BFD@has.data.eval){
                rowToRm <- unique(unlist(lapply(BFD@eval.data.env.var,function(x){return(which(is.na(x)))})))
                if(length(rowToRm)){
                  cat("\n\t\t\t! Some NAs have been automaticly removed from your evaluation data")
                  BFD@eval.coord <- BFD@eval.coord[-rowToRm,,drop=FALSE]
                  BFD@eval.data.species <- BFD@eval.data.species[-rowToRm]
                  BFD@eval.data.env.var <- BFD@eval.data.env.var[-rowToRm,,drop=FALSE]
                }      
              }
            }
            
            
            
            # count data occutances
            #     BFD@data.counting <- matrix(c(sum(BFD@data.species, na.rm=TRUE),sum(BFD@data.species==0, na.rm=TRUE)),
            #                             ncol=1,nrow=2, dimnames=list(c("nb.pres","nb.abs"),c("data.species") ) )
            #     
            #     if(BFD@has.data.eval){
            #       BFD@data.counting <- cbind(BFD@data.counting,c(sum(BFD@eval.data.species, na.rm=TRUE),sum(BFD@eval.data.species==0, na.rm=TRUE)))
            #       colnames(BFD@data.counting)[ncol(BFD@data.counting)] <- "eval.data.species"
            #     }
            
            return(BFD)
          }
)

setMethod('BIOMOD.formated.data', signature(sp='data.frame'), 
          function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE){
            if(ncol(sp) > 1 ){
              stop("Invalid response variable")
            }
            sp <- as.numeric(unlist(sp))
            BFD <- BIOMOD.formated.data(sp,env,xy,sp.name, eval.sp, eval.env, eval.xy, na.rm=na.rm)
            return(BFD)
          }
)

setMethod('BIOMOD.formated.data', signature(sp='numeric', env='matrix' ), 
          function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE){
            env <- as.data.frame(env)
            BFD <- BIOMOD.formated.data(sp,env,xy,sp.name, eval.sp, eval.env, eval.xy, na.rm=na.rm)
            return(BFD)
          }
)

setMethod('BIOMOD.formated.data', signature(sp='numeric', env='RasterStack' ), 
          function(sp,env,xy=NULL,sp.name=NULL, eval.sp=NULL, eval.env=NULL, eval.xy=NULL, na.rm=TRUE){
            categorial_var <- names(env)[raster::is.factor(env)]  
            
            # take the same eval environemental variables than calibrating ones 
            if(!is.null(eval.sp)){
              if(is.null(eval.env)){
                #         eval.env_levels <- levels(eval.env)
                eval.env <- as.data.frame(extract(env,eval.xy))
                if(length(categorial_var)){
                  for(cat_var in categorial_var){
                    eval.env[,cat_var] <- as.factor(eval.env[,cat_var])
                  }
                }
              }
            }

            if(is.null(xy)) xy <- as.data.frame(coordinates(env))
                        
            data.mask = reclassify(raster::subset(env,1,drop=T), c(-Inf,Inf,-1))
            data.mask[cellFromXY(data.mask,xy[which(sp==1),])] <- 1
            data.mask[cellFromXY(data.mask,xy[which(sp==0),])] <- 0
            data.mask <- raster::stack(data.mask)
            names(data.mask) <- sp.name
                        
            #     env_levels <- levels(env)
            env <- as.data.frame(extract(env,xy, factors=T))
                        
            if(length(categorial_var)){
              for(cat_var in categorial_var){
                env[,cat_var] <- as.factor(env[,cat_var])
              }
            }
            
            BFD <- BIOMOD.formated.data(sp,env,xy,sp.name,eval.sp, eval.env, eval.xy, na.rm=na.rm, data.mask=data.mask)
            return(BFD)
          }
)

# 1.3 Other Functions
# if( !isGeneric( "plot" ) ) {
#   setGeneric( "plot", 
#               def = function(x, ...){
#   	                  standardGeneric( "plot" )
#                       } )
# }

setMethod('plot', signature(x='BIOMOD.formated.data', y="missing"),
          function(x,coord=NULL,col=NULL){
            if(raster::nlayers(x@data.mask)>0){
              requireNamespace("rasterVis")
              
              ## check that there is some undefined areas to prevent from strange plotting issues
              if(min(cellStats(x@data.mask,min)) == -1){ # there is undifined area
                ## define the breaks of the color key
                my.at <- seq(-1.5,1.5,by=1)
                ## the labels will be placed vertically centered
                my.labs.at <- seq(-1,1,by=1)
                ## define the labels
                my.lab <- c("undifined","absences","presences")
                ## define colors
                my.col.regions = c("lightgrey","red4","green4")
                ## defined cuts
                my.cuts <- 2
              } else{ # no undefined area.. remove it from plot 
                ## define the breaks of the color key
                my.at <- seq(-0.5,1.5,by=1)
                ## the labels will be placed vertically centered
                my.labs.at <- seq(0,1,by=1)
                ## define the labels
                my.lab <- c("absences","presences")
                ## define colors
                my.col.regions = c("red4","green4")
                ## defined cuts
                my.cuts <- 1
              }
              
              
              levelplot(x@data.mask, at=my.at, cuts=my.cuts, margin=T, col.regions=my.col.regions,
                        main=paste(x@sp.name,"datasets"),
                        colorkey=list(labels=list(
                          labels=my.lab,
                          at=my.labs.at)))
              
            } else{
              # coordinates checking
              if(is.null(coord)){
                if( sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2] ){
                  stop("coordinates are required to plot your data")
                } else {
                  coord <- x@coord
                }
              }
              
              # colors checking
              if(is.null(col) | length(col) < 3){
                col = c('green', 'red', 'grey')
              }
              
              # plot data
              # all points (~mask) 
              
              plot(x=x@coord[,1], y=x@coord[,2], col=col[3], xlab = 'X', ylab = 'Y',
                   main = paste(x@sp.name, sep=""), pch=20 )
              # presences 
              points(x=x@coord[which(x@data.species == 1),1],
                     y=x@coord[which(x@data.species == 1),2],
                     col=col[1],pch=18)
              # true absences
              points(x=x@coord[which(x@data.species == 0),1],
                     y=x@coord[which(x@data.species == 0),2],
                     col=col[2],pch=18)
              
            }
            
            
            
          })

setMethod('show', signature('BIOMOD.formated.data'),
          function(object){
            .bmCat("'BIOMOD.formated.data'")
            cat("\nsp.name = ", object@sp.name, fill=.Options$width)
            cat("\n\t", sum(object@data.species, na.rm=TRUE), 'presences, ',
                sum(object@data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
            cat("\n\n\t", ncol(object@data.env.var), 'explanatory variables\n', fill=.Options$width)
            print(summary(object@data.env.var))
            
            if(object@has.data.eval){
              cat("\n\nEvaluation data :", fill=.Options$width)
              cat("\n\t", sum(object@eval.data.species, na.rm=TRUE), 'presences, ',
                  sum(object@eval.data.species==0, na.rm=TRUE), 'true absences and ', 
                  sum(is.na(object@eval.data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
              cat("\n\n", fill=.Options$width)
              print(summary(object@eval.data.env.var))
            }
            
            .bmCat()
          })

# 2. The BIOMOD.formated.data.PA =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=- #
# this class inherits from BIOMOD.formated.data and have one more slot 'PA' giving PA selected

# 2.1 Class Definition
setClass("BIOMOD.formated.data.PA",
         contains = "BIOMOD.formated.data",
         representation(PA.strategy='character', PA = 'data.frame'),
         validity = function(object){
           return(TRUE)
         })

# 2.2 Constructors
# if( !isGeneric( "BIOMOD.formated.data.PA" ) ){
#   setGeneric( "BIOMOD.formated.data.PA", 
#               def = function(sp, env, PA.NbRep, ...){
#                 standardGeneric( "BIOMOD.formated.data.PA" )
#               })
# }

# setMethod('BIOMOD.formated.data.PA',signature(sp='ANY',
#                                               env='ANY',
#                                               PA.NbRep='integer'),

BIOMOD.formated.data.PA <-  function(sp, env, xy, sp.name,
                                     eval.sp=NULL, eval.env=NULL, eval.xy=NULL,
                                     PA.NbRep=1,
                                     PA.strategy='random',
                                     PA.nb.absences = NULL,
                                     PA.dist.min = 0,
                                     PA.dist.max = NULL,
                                     PA.sre.quant = 0.025,
                                     PA.table = NULL,
                                     na.rm=TRUE){
  
  if(inherits(env,'Raster')){
    categorial_var <- names(env)[raster::is.factor(env)]
  }  else categorial_var <- NULL
  
  
  # take the same eval environemental variables than calibrating ones 
  if(!is.null(eval.sp)){
    if(is.null(eval.env)){
      if(inherits(env,'Raster')){
        eval.env <- as.data.frame(extract(env,eval.xy))
        if(length(categorial_var)){
          for(cat_var in categorial_var){
            eval.env[,cat_var] <- as.factor(eval.env[,cat_var])
          }
        }
      } else{
        stop("No evaluation explanatory variable given")
      }
    }
  }
  
  # convert sp in spatial obj
  if(is.numeric(sp)){
    if(is.null(xy)){
      sp <- SpatialPointsDataFrame(matrix(0,ncol=2,nrow=length(sp)), data.frame(sp),match.ID=FALSE)
    } else{
      sp <- SpatialPointsDataFrame(data.matrix(xy), data.frame(sp),match.ID=FALSE)
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
  
  if(!is.null(pa.data.tmp)){
    
    if(length(categorial_var)){
      for(cat_var in categorial_var){
        pa.data.tmp$env[,cat_var] <- as.factor(pa.data.tmp$env[,cat_var])
      }
    }
    
    
    if(na.rm){
      rowToRm <- unique(unlist(lapply(pa.data.tmp$env,function(x){return(which(is.na(x)))})))
      if(length(rowToRm)){
        cat("\n\t\t\t! Some NAs have been automaticly removed from your data")
        pa.data.tmp$xy <- pa.data.tmp$xy[-rowToRm,,drop=FALSE]        
        pa.data.tmp$sp <- pa.data.tmp$sp[-rowToRm, drop=FALSE]
        pa.data.tmp$env <- pa.data.tmp$env[-rowToRm,,drop=FALSE]
        pa.data.tmp$pa.tab <- pa.data.tmp$pa.tab[-rowToRm,,drop=FALSE]
        #         cat("\n***\n")
        #         cat("\ndim(xy) <-", dim(pa.data.tmp$xy))
        #         cat("\nclass(sp) <- ", class(pa.data.tmp$sp))
        #         cat("\ndim(sp) <-", dim(pa.data.tmp$sp))
        #         cat("\ndim(env) <-", dim(pa.data.tmp$env))
        #         cat("\ndim(pa.tab) <-", dim(pa.data.tmp$pa.tab))
        
      }      
    }
    
    
    # data counting
    #     pa.data.tmp$data.counting <- apply(pa.data.tmp$pa.tab,2,function(x){nbPres <- sum(pa.data.tmp$sp[x],na.rm=T) ; return(c(nbPres,sum(x)-nbPres))})
    #     colnames(pa.data.tmp$data.counting) <- colnames(pa.data.tmp$pa.tab)
    
    BFD <- BIOMOD.formated.data(sp=pa.data.tmp$sp,
                                env=pa.data.tmp$env,
                                xy=as.data.frame(pa.data.tmp$xy),
                                sp.name=sp.name,
                                eval.sp=eval.sp,
                                eval.env=eval.env,
                                eval.xy=eval.xy,
                                na.rm=na.rm) # because check is already done
    
    if(inherits(env,'Raster')){
      
      ## create data.mask for ploting
      data.mask.tmp <- reclassify(raster::subset(env,1), c(-Inf,Inf,-1))
      data.mask <- stack(data.mask.tmp)
      xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp==1),]
      xy_abs <- pa.data.tmp$xy[which(pa.data.tmp$sp==0),]
      if(nrow(xy_pres)){
        data.mask[cellFromXY(data.mask.tmp, xy_pres)] <- 1
      }
      if(nrow(xy_abs)){
        data.mask[cellFromXY(data.mask.tmp, xy_abs)] <- 0
      }
      names(data.mask) <- "input_data"
      
      ## add eval data
      if(BFD@has.data.eval){
        ### TO DO
        
      }
      
      for(pa in 1:ncol(as.data.frame(pa.data.tmp$pa.tab))){
        data.mask.tmp2 <- data.mask.tmp
        
        xy_pres <- pa.data.tmp$xy[which(pa.data.tmp$sp==1 & as.data.frame(pa.data.tmp$pa.tab)[,pa]) ,]
        xy_abs <- pa.data.tmp$xy[which( (pa.data.tmp$sp!=1 | is.na(pa.data.tmp$sp)) & as.data.frame(pa.data.tmp$pa.tab)[,pa]) ,]
        
        if(nrow(xy_pres)){
          id_pres <- cellFromXY(data.mask.tmp, xy_pres)
          data.mask.tmp2[id_pres] <- 1
        }
        
        if(nrow(xy_abs)){
          id_abs <- cellFromXY(data.mask.tmp, xy_abs)
          data.mask.tmp2[id_abs] <- 0
        }
        
        data.mask <- addLayer(data.mask, data.mask.tmp2)
      }
      
      names(data.mask) <- c("input_data", colnames(as.data.frame(pa.data.tmp$pa.tab)))
      
    } else{
      data.mask <- stack()
    }
    
    BFDP <- new('BIOMOD.formated.data.PA',
                sp.name = BFD@sp.name,
                coord = BFD@coord,
                #                 data.counting = cbind(BFD@data.counting,pa.data.tmp$data.counting) ,
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
  } else {
    cat("\n   ! PA selection not done", fill=.Options$width)
    
    BFDP <- BIOMOD.formated.data(sp=sp,
                                 env=env,
                                 xy=xy,
                                 sp.name=sp.name,
                                 eval.sp=eval.sp,
                                 eval.env=eval.env,
                                 eval.xy=eval.xy)
    
  }
  
  rm(list = "pa.data.tmp" )
  
  return(BFDP)
  
}


# 2.3 other functions
setMethod('plot', signature(x='BIOMOD.formated.data.PA', y="missing"),
          function(x,coord=NULL,col=NULL){

            if(raster::nlayers(x@data.mask)>0){
              requireNamespace("rasterVis")
              
              ## check that there is some undefined areas to prevent from strange plotting issues
              if(min(cellStats(x@data.mask,min)) == -1){ # there is undifined area
                ## define the breaks of the color key
                my.at <- seq(-1.5,1.5,by=1)
                ## the labels will be placed vertically centered
                my.labs.at <- seq(-1,1,by=1)
                ## define the labels
                my.lab <- c("undifined","absences","presences")
                ## define colors
                my.col.regions = c("lightgrey","red4","green4")
                ## defined cuts
                my.cuts <- 2
              } else{ # no undefined area.. remove it from plot 
                ## define the breaks of the color key
                my.at <- seq(-0.5,1.5,by=1)
                ## the labels will be placed vertically centered
                my.labs.at <- seq(0,1,by=1)
                ## define the labels
                my.lab <- c("absences","presences")
                ## define colors
                my.col.regions = c("red4","green4")
                ## defined cuts
                my.cuts <- 1
              }
              
              
              levelplot(x@data.mask, at=my.at, cuts=my.cuts, margin=T, col.regions=my.col.regions,
                        main=paste(x@sp.name,"datasets"),
                        colorkey=list(labels=list(
                          labels=my.lab,
                          at=my.labs.at)))
              
            } else{
              # coordinates checking
              if(is.null(coord)){
                if( sum(is.na(x@coord)) == dim(x@coord)[1] * dim(x@coord)[2] ){
                  stop("coordinates are required to plot your data")
                } else {
                  coord <- x@coord
                }
              }
              
              # colors checking
              if(is.null(col) | length(col) < 3){
                col = c('green', 'red', 'orange', 'grey')
              }
              
              # plot data
              par(mfrow=c(.CleverCut(ncol(x@PA)+1)))
              # all points (~mask)
              plot(x=x@coord[,1], y=x@coord[,2], col=col[4], xlab = 'X', ylab = 'Y',
                   main = paste(x@sp.name," original data", sep=""), pch=20 )
              # presences 
              points(x=x@coord[which(x@data.species == 1),1],
                     y=x@coord[which(x@data.species == 1),2],
                     col=col[1],pch=18)
              # true absences
              points(x=x@coord[which(x@data.species == 0),1],
                     y=x@coord[which(x@data.species == 0),2],
                     col=col[2],pch=18)
              
              # PA data
              for(i in 1:ncol(x@PA)){
                # all points (~mask)
                plot(x=x@coord[,1], y=x@coord[,2], col=col[4], xlab = 'X', ylab = 'Y',
                     main = paste(x@sp.name," Pseudo Absences ", i, sep=""), pch=20 )
                # presences 
                points(x=x@coord[(x@data.species == 1) & x@PA[,i],1],
                       y=x@coord[(x@data.species == 1) & x@PA[,i],2],
                       col=col[1],pch=18)
                # true absences
                points(x=x@coord[(x@data.species == 0) & x@PA[,i],1],
                       y=x@coord[(x@data.species == 0) & x@PA[,i],2],
                       col=col[2],pch=18)
                # PA
                points(x=x@coord[is.na(x@data.species) & x@PA[,i],1],
                       y=x@coord[is.na(x@data.species) & x@PA[,i],2],
                       col=col[3],pch=18)
              }
            } })

setMethod('show', signature('BIOMOD.formated.data.PA'),
          function(object){
            .bmCat("'BIOMOD.formated.data.PA'")
            cat("\nsp.name = ", object@sp.name,fill=.Options$width)
            cat("\n\t", sum(object@data.species, na.rm=TRUE), 'presences, ',
                sum(object@data.species==0, na.rm=TRUE), 'true absences and ', 
                sum(is.na(object@data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
            cat("\n\n\t", ncol(object@data.env.var), 'explanatory variables\n', fill=.Options$width)
            print(summary(object@data.env.var))
            
            if(object@has.data.eval){
              cat("\n\nEvaluation data :", fill=.Options$width)
              cat("\n\t", sum(object@eval.data.species, na.rm=TRUE), 'presences, ',
                  sum(object@eval.data.species==0, na.rm=TRUE), 'true absences and ', 
                  sum(is.na(object@eval.data.species), na.rm=TRUE),'undifined points in dataset', fill=.Options$width)
              cat("\n\n", fill=.Options$width)
              print(summary(object@eval.data.env.var))
            }
            
            cat("\n\n", ncol(object@PA), 'Pseudo Absences dataset available (', colnames(object@PA),") with ",
                sum(object@PA[,1], na.rm=T) - sum(object@data.species, na.rm=TRUE), 'absences in each (true abs + pseudo abs)', fill=.Options$width)
            .bmCat()
          })


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
                        MAXENT.Tsuruoka = "list"),
         
         prototype(GLM = list( type = 'quadratic',
                               interaction.level = 0,
                               myFormula = NULL,
                               test = 'AIC',
                               family = binomial(link='logit'),
                               mustart = 0.5,
                               control = glm.control(maxit = 50)),
                   
                   GBM = list(  distribution = 'bernoulli',
                                n.trees = 2500,
                                interaction.depth = 7,
                                n.minobsinnode = 5,
                                shrinkage = 0.001,
                                bag.fraction = 0.5,
                                train.fraction = 1,
                                cv.folds = 3,
                                keep.data = FALSE,
                                verbose = FALSE,
                                #                                 class.stratify.cv = 'bernoulli',
                                perf.method = 'cv'),
                   
                   GAM = list( algo = "GAM_mgcv",
                               type = "s_smoother",
                               k = NULL,
                               interaction.level = 0,
                               myFormula = NULL,
                               family = binomial(link='logit'),
                               control = list(epsilon = 1e-06, trace = FALSE ,maxit = 100),
                               method = "GCV.Cp",
                               optimizer=c("outer","newton"),
                               select=FALSE,
                               knots=NULL,
                               paraPen=NULL), 
                   
                   CTA = list(method = 'class',
                              parms = 'default',
                              #                               control = rpart.control(xval = 5, minbucket = 5, minsplit = 5,
                              #                                                       cp = 0.001, maxdepth = 25),
                              control = list(xval = 5, minbucket = 5, minsplit = 5,
                                             cp = 0.001, maxdepth = 25),                              
                              cost = NULL ),
                   
                   ANN = list(NbCV = 5,
                              size = NULL,
                              decay = NULL,
                              rang = 0.1,
                              maxit = 200),
                   
                   SRE = list(quant = 0.025),
                   
                   FDA = list(method = 'mars',
                              add_args = NULL),
                   
                   MARS = list(type = 'simple',
                               interaction.level = 0,
                               myFormula = NULL,
#                                degree = 1,
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
                   
                   MAXENT.Tsuruoka = list(l1_regularizer = 0.0,
                                          l2_regularizer = 0.0,
                                          use_sgd = FALSE,
                                          set_heldout = 0,
                                          verbose = FALSE)
                   
         ),
         validity = function(object){
           test <- TRUE
           ## GLM ##
           if(!(object@GLM$type %in% c('simple','quadratic','polynomial','user.defined'))){ cat("\nGLM$type must be 'simple',  'quadratic', 'polynomial' or 'user.defined'"); test <- FALSE}
           if(!is.numeric(object@GLM$interaction.level)){ cat("\nGLM$interaction.level must be a integer"); test <- FALSE } else{
             if(object@GLM$interaction.level < 0 | object@GLM$interaction.level%%1!=0){ cat("\nGLM$interaction.level must be a positive integer"); test <- FALSE }
           }
           if(!is.null(object@GLM$myFormula)) if(class(object@GLM$myFormula) != "formula"){ cat("\nGLM$myFormula must be NULL or a formula object"); test <- FALSE }
           if(!(object@GLM$test %in% c('AIC','BIC','none'))){ cat("\nGLM$test must be 'AIC','BIC' or 'none'"); test <- FALSE}
           fam <- 'none'
           if(class(object@GLM$family) != "family"){ cat("\nGLM$family must be a valid family object"); test <- FALSE }
           if(!is.list(object@GLM$control)){cat("\nGLM$control must be a list like that returned by glm.control"); test <- FALSE}
           
           
           ## GBM ##
           if(! object@GBM$distribution %in% c("bernoulli","huberized", "multinomial", "adaboost")){cat("\nGBM$distribution must be 'bernoulli', 'huberized', 'multinomial'or  'adaboost'") }
           
           if(!is.numeric(object@GBM$n.trees)){ cat("\nGBM$n.trees must be a integer"); test <- FALSE } else{
             if(object@GBM$n.trees < 0 | floor(object@GBM$n.trees) != object@GBM$n.trees){ cat("\nGBM$n.trees must be a positive integer"); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$interaction.depth)){ cat("\nGBM$interaction.depth must be a integer"); test <- FALSE } else{
             if(object@GBM$interaction.depth < 0 | object@GBM$interaction.depth%%1!=0){ cat("\nGBM$interaction.depth must be a positive integer"); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$n.minobsinnode)){ cat("\nGBM$n.minobsinnode must be a integer"); test <- FALSE } else{
             if(object@GBM$n.minobsinnode < 0 | object@GBM$n.minobsinnode%%1!=0){ cat("\nGBM$n.minobsinnode must positive "); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$shrinkage)){ cat("\nGBM$shrinkage must be a numeric"); test <- FALSE } else{
             if(object@GBM$shrinkage < 0 ){ cat("\nGBM$shrinkage must positive "); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$bag.fraction)){ cat("\nGBM$bag.fraction must be a numeric"); test <- FALSE } else{
             if(object@GBM$bag.fraction < 0 | object@GBM$bag.fraction > 1){ cat("\nGBM$bag.fraction must be a 0 to 1 numeric"); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$train.fraction)){ cat("\nGBM$train.fraction must be a numeric"); test <- FALSE } else{
             if(object@GBM$train.fraction < 0 | object@GBM$train.fraction > 1){ cat("\nGBM$train.fraction must be a 0 to 1 numeric"); test <- FALSE }
           }
           
           if(!is.numeric(object@GBM$cv.folds)){ cat("\nGBM$cv.folds must be a integer"); test <- FALSE } else{
             if(object@GBM$cv.folds < 0 | object@GBM$cv.folds%%1!=0){ cat("\nGBM$cv.folds must be a positive integer"); test <- FALSE }
           }
           
           if(!is.logical(object@GBM$keep.data)){ cat("\nGBM$keep.data must be a logical"); test <- FALSE }
           
           if(!is.logical(object@GBM$verbose)){ cat("\nGBM$verbose must be a logical"); test <- FALSE }
           
           #            if(! object@GBM$class.stratify.cv %in% c("bernoulli", "multinomial")){cat("\nGBM$class.stratify.cv must be 'bernoulli' or 'multinomial'") }
           
           if(! object@GBM$perf.method %in% c('OOB', 'test', 'cv')){cat("\nGBM$perf.method must be 'OOB', 'test' or 'cv'"); test <- FALSE }
           
           
           ## GAM ##
           if(! object@GAM$algo %in% c('GAM_mgcv','GAM_gam', 'BAM_mgcv')){cat("\nGAM$algo must be 'GAM_mgcv','GAM_gam' or  'BAM_mgcv'"); test <- FALSE }
           
           if(! object@GAM$type %in% c('s_smoother','s', 'lo', 'te')){cat("\nGAM$type must be c('s_smoother','s', 'lo' or 'te'"); test <- FALSE }
           
           if(! is.null(object@GAM$k)){
             if(! is.numeric(object@GAM$k)  ){ cat("\nGAM$k must be a integer"); test <- FALSE } else{
               if(object@GAM$k < -1 | object@GAM$k%%1!=0){ cat("\nGAM$k must be > -1"); test <- FALSE }
             }
           }
           
           
           if(!is.numeric(object@GAM$interaction.level)){ cat("\nGAM$interaction.level must be a integer"); test <- FALSE } else{
             if(object@GAM$interaction.level < 0 | object@GAM$interaction.level%%1!=0){ cat("\nGAM$interaction.level must be a positive integer"); test <- FALSE }
           }
           
           if(!is.null(object@GAM$myFormula)) if(class(object@GAM$myFormula) != "formula"){ cat("\nGAM$myFormula must be NULL or a formula object"); test <- FALSE }
           
           if(class(object@GAM$family) != "family"){ cat("\nGAM$family must be a valid family object"); test <- FALSE }
           
           if(!is.list(object@GAM$control)){cat("\nGAM$control must be a list like that returned by gam.control"); test <- FALSE}
           if(! object@GAM$method %in% c('GCV.Cp','GACV.Cp','REML', 'P-REML', 'ML', 'P-ML')){cat("\nGAM$method must be 'GCV.Cp','GACV.Cp','REML', 'P-REML', 'ML'or 'P-ML'"); test <- FALSE}
           
           if(sum(! object@GAM$optimizer %in% c('perf','outer', 'newton', 'bfgs', 'optim', 'nlm', 'nlm.fd')) > 0 ){cat("\nGAM$optimizer bad definition (see ?mgcv::gam)") ; test <- FALSE}
           
           if(!is.logical(object@GAM$select)){ cat("\nGAM$select must be a logical"); test <- FALSE }
           
           #            knots=NULL,
           #            paraPen=NULL
           
           
           ## CTA ##
           if(! object@CTA$method %in% c( 'anova', 'poisson', 'class', 'exp')){cat("\nCTA$method must be 'anova', 'poisson', 'class' or 'exp'"); test <- FALSE }
           
           #parms = 'default',
           
           if(!is.list(object@CTA$control)){cat("\nCTA$control must be a list like that returned by rpart.control"); test <- FALSE}                           
           if(length(object@CTA$cost)){
             if(!is.numeric(object@CTA$cost)){cat("\nCTA$cost must be a non negative cost vector"); test <- FALSE}
           }
           
           
           
           ## ANN ##
           if(!is.numeric(object@ANN$NbCV)){ cat("\nANN$NbCV must be a integer"); test <- FALSE } else{
             if(object@ANN$NbCV < 0 | object@ANN$NbCV%%1!=0){ cat("\nANN$NbCV must be a positive integer"); test <- FALSE }
           }
           
           if( ( is.null(object@ANN$size) | length(object@ANN$size)>1 ) & object@ANN$NbCV <= 0){ cat("\nANN$size has to be defined as a single integer if ANN$NbCV=0"); test <- FALSE } else{
             if(!is.null(object@ANN$size)) if( !is.numeric(object@ANN$size) | !all( object@ANN$size > 0 ) | !all( object@ANN$size %% 1 == 0 ) ){ cat("\nANN$size must be NULL or a positive (vector of) integer"); test <- FALSE }
           }

           if( ( is.null(object@ANN$decay) | length(object@ANN$decay)>1 ) & object@ANN$NbCV <= 0){ cat("\nANN$decay has to be defined as a single number if ANN$NbCV=0"); test <- FALSE } else{
             if(!is.null(object@ANN$decay)) if( !is.numeric(object@ANN$decay) | !all( object@ANN$decay > 0 ) ){ cat("\nANN$decay must be NULL or a positive (vector of) number"); test <- FALSE }
           }
           
           if(!is.numeric(object@ANN$rang)){ cat("\nANN$rang must be a numeric"); test <- FALSE } else{
             if(object@ANN$rang < 0 ){ cat("\nANN$rang must be positive"); test <- FALSE }
           }
           
           if(!is.numeric(object@ANN$maxit)){ cat("\nANN$maxit must be a integer"); test <- FALSE } else{
             if(object@ANN$maxit < 0 | object@ANN$maxit%%1!=0){ cat("\nANN$maxit must be a positive integer"); test <- FALSE }
           }
           
           
           
           ## FDA ##
           if(! object@FDA$method %in% c( 'polyreg', 'mars', 'bruto')){cat("\nFDA$method must be 'polyreg', 'mars' or 'bruto'"); test <- FALSE }
           if(!is.null(object@FDA$add_args)){ if(!is.list(object@FDA$add_args)) {cat("\nFDA$add_args must be a list or NULL"); test <- FALSE } }
           
           
           ## SRE ##
           if(!is.numeric(object@SRE$quant)){ cat("\nSRE$quant must be a numeric"); test <- FALSE } else{
             if(object@SRE$quant >= 0.5 | object@SRE$quant < 0){ cat("\nSRE$quant must between 0 and 0.5"); test <- FALSE }
           }
           
           
           
           ## MARS ##
           if(!(object@MARS$type %in% c('simple','quadratic','polynomial','user.defined'))){ cat("\nMARS$type must be 'simple',  'quadratic', 'polynomial' or 'user.defined'"); test <- FALSE}
           if(!is.numeric(object@MARS$interaction.level)){ cat("\nMARS$interaction.level must be a integer"); test <- FALSE } else{
             if(object@MARS$interaction.level < 0 | object@MARS$interaction.level%%1!=0){ cat("\nMARS$interaction.level must be a positive integer"); test <- FALSE }
           }
           if(!is.null(object@MARS$myFormula)) if(class(object@MARS$myFormula) != "formula"){ cat("\nMARS$myFormula must be NULL or a formula object"); test <- FALSE }
#            if(!is.numeric(object@MARS$degree)){ cat("\nMARS$degree must be a integer"); test <- FALSE } else{
#              if(object@MARS$degree < 0 | object@MARS$degree%%1!=0){ cat("\nMARS$degree must be a positive integer"); test <- FALSE }
#            }
           if(!is.null(object@MARS$nk)){ 
             if(object@MARS$nk < 0 | object@MARS$nk%%1!=0){ cat("\nMARS$nk must be a positive integer or NULL if you want to use default parameter"); test <- FALSE }
           }
           if(!is.numeric(object@MARS$penalty)){ cat("\nMARS$penalty must be a integer"); test <- FALSE } else{
             if(object@MARS$penalty < 0 | object@MARS$penalty%%1!=0){ cat("\nMARS$penalty must be a positive integer"); test <- FALSE }
           }
           if(!is.numeric(object@MARS$thresh)){ cat("\nMARS$thresh must be a numeric"); test <- FALSE } else{
             if(object@MARS$thresh < 0 ){ cat("\nMARS$thresh must be positive"); test <- FALSE }
           }
           if(!is.null(object@MARS$nprune)){ if(!is.numeric(object@MARS$nprune)){ cat("\nMARS$nprune must be a numeric or NULL"); test <- FALSE }}
           supported.pmethod <- c('backward', 'none', 'exhaustive', 'forward', 'seqrep', 'cv')
           if(!is.element(object@MARS$pmethod, supported.pmethod)){cat("\nMARS$pmethod must be a one of", supported.pmethod); test <- FALSE }

           
           ## RF ##
           if(!is.logical(object@RF$do.classif)){ cat("\nRF$do.classif must be a logical"); test <- FALSE }
           
           if(!is.numeric(object@RF$ntree)){ cat("\nRF$ntree must be a integer"); test <- FALSE } else{
             if(object@RF$ntree < 0 | object@RF$ntree%%1!=0){ cat("\nRF$ntree must be a positive integer"); test <- FALSE }
           }
           
           if( object@RF$mtry != 'default'){
             if(!is.numeric(object@RF$mtry)){ cat("\nRF$mtry must be a integer"); test <- FALSE } else{
               if(object@RF$mtry < 0 | object@RF$mtry%%1!=0){ cat("\nRF$mtry must be a positive integer"); test <- FALSE }
             }
           }
           
           if(!is.numeric(object@RF$nodesize)){ cat("\nRF$nodesize must be a integer"); test <- FALSE } else{
             if(object@RF$nodesize < 0 | object@RF$nodesize%%1!=0){ cat("\nRF$nodesize must be a positive integer"); test <- FALSE }
           }
           
           if(length(object@RF$maxnodes)){
             if(!is.numeric(object@RF$maxnodes)){ cat("\nRF$maxnodes must be a integer"); test <- FALSE } else{
               if(object@RF$maxnodes < 0 | object@RF$maxnodes%%1!=0){ cat("\nRF$maxnodes must be a positive integer"); test <- FALSE }
             }             
           }
           
           
           
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
           if(!is.numeric(object@MAXENT.Phillips$maximumiterations)){ cat("\nMAXENT.Phillips$maximumiterations must be a integer"); test <- FALSE } else{
             if(object@MAXENT.Phillips$maximumiterations < 0 | object@MAXENT.Phillips$maximumiterations%%1!=0){ cat("\nMAXENT.Phillips$maximumiterations must be a positive integer"); test <- FALSE }
           }
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
           
           ## MAXENT.Tsuruoka
		       if(!is.numeric(object@MAXENT.Tsuruoka$l1_regularizer)){ cat("\nMAXENT.Tsuruoka$l1_regularizer must be a numeric"); test <- FALSE }		   
		       if(!is.numeric(object@MAXENT.Tsuruoka$l2_regularizer)){ cat("\nMAXENT.Tsuruoka$l2_regularizer must be a numeric"); test <- FALSE }
		       if(!is.logical(object@MAXENT.Tsuruoka$use_sgd)){ cat("\nMAXENT.Tsuruoka$use_sgd must be a logical"); test <- FALSE }
		       if(!is.numeric(object@MAXENT.Tsuruoka$set_heldout)){ cat("\nMAXENT.Tsuruoka$set_heldout must be a numeric"); test <- FALSE }
		       if(!is.logical(object@MAXENT.Tsuruoka$verbose)){ cat("\nMAXENT.Tsuruoka$verbose must be a logical"); test <- FALSE }
		   
           return(test)
         })

setMethod('show', signature('BIOMOD.Model.Options'),
          function(object){
            .bmCat(" 'BIOMOD.Model.Options' ")
            cat("\n")
            
            ## GLM options
            cat("\nGLM = list( type = '", object@GLM$type, "',", sep="")
            cat("\n            interaction.level = ", object@GLM$interaction.level, ",", sep="")
            cat("\n            myFormula = ",  ifelse(length(object@GLM$myFormula) < 1,'NULL',paste(object@GLM$myFormula[2],object@GLM$myFormula[1],object@GLM$myFormula[3])), ",", sep="")
            cat("\n            test = '", object@GLM$test, "',", sep="")
            cat("\n            family = ", object@GLM$family$family,"(link = '",object@GLM$family$link,"'),", sep="")
            cat("\n            mustart = ", object@GLM$mustart, ",", sep="")
            cat("\n            control = glm.control(", .print.control(object@GLM$control), ") ),", sep="", fill=.Options$width)
            
            ## GBM options
            cat("\n")
            cat("\nGBM = list( distribution = '", object@GBM$distribution, "',", sep="")
            cat("\n            n.trees = ", object@GBM$n.trees, ",", sep="")
            cat("\n            interaction.depth = ", object@GBM$interaction.depth, ",", sep="")
            cat("\n            n.minobsinnode = ", object@GBM$n.minobsinnode, ",", sep="")
            cat("\n            shrinkage = ", object@GBM$shrinkage, ",", sep="")
            cat("\n            bag.fraction = ", object@GBM$bag.fraction, ",", sep="")
            cat("\n            train.fraction = ", object@GBM$train.fraction, ",", sep="")
            cat("\n            cv.folds = ", object@GBM$cv.folds, ",", sep="")
            cat("\n            keep.data = ", object@GBM$keep.data, ",", sep="")
            cat("\n            verbose = ", object@GBM$verbose, ",", sep="")
            #             cat("\n            class.stratify.cv = '", object@GBM$class.stratify.cv, "',", sep="")
            cat("\n            perf.method = '", object@GBM$perf.method, "'),", sep="")
            
            ## GAM options
            cat("\n")
            cat("\nGAM = list( algo = '", object@GAM$algo, "',", sep="")
            cat("\n            type = '", object@GAM$type, "',", sep="")
            cat("\n            k = ", ifelse(length(object@GAM$k) < 1,'NULL',object@GAM$k), ",", sep="")
            cat("\n            interaction.level = ", object@GAM$interaction.level, ",", sep="")
            cat("\n            myFormula = ", ifelse(length(object@GAM$myFormula) < 1,'NULL',paste(object@GAM$myFormula[2],object@GAM$myFormula[1],object@GAM$myFormula[3])), ",", sep="")
            cat("\n            family = ", object@GAM$family$family,"(link = '",object@GAM$family$link,"'),", sep="")
            
            if(object@GAM$algo=='GAM_mgcv'){
              cat("\n            method = '", object@GAM$method, "',", sep="")
              cat("\n            optimizer = c('", paste(object@GAM$optimizer,collapse="','"), "'),", sep="")
              cat("\n            select = ", object@GAM$select, ",", sep="")
              cat("\n            knots = ",  ifelse(length(object@GLM$knots) < 1,'NULL',"'user.defined'"), ",", sep="")
              cat("\n            paraPen = ",  ifelse(length(object@GLM$paraPen) < 1,'NULL',"'user.defined'"), ",", sep="")
            }
            
            cat("\n            control = list(", .print.control(object@GAM$control), ") ),", sep="", fill=.Options$width)
            
            
            
            ## CTA options
            cat("\n")
            cat("\nCTA = list( method = '", object@CTA$method, "',", sep="")
            cat("\n            parms = '", object@CTA$parms, "',", sep="")
            cat("\n            cost = ", ifelse(length(object@CTA$cost)<1,'NULL',object@CTA$cost), ",", sep="")
            cat("\n            control = list(", .print.control(object@CTA$control), ") ),", sep="", fill=.Options$width)
            
            ## ANN options
            cat("\n")
            cat("\nANN = list( NbCV = ", object@ANN$NbCV, ",", sep="")
            cat("\n            size = ", ifelse(length(object@ANN$size)<1,'NULL',object@ANN$size), ",", sep="")
            cat("\n            decay = ", ifelse(length(object@ANN$decay)<1,'NULL',object@ANN$decay), ",", sep="")
            cat("\n            rang = ", object@ANN$rang, ",", sep="")
            cat("\n            maxit = ", object@ANN$maxit, "),", sep="")
            
            ## SRE options
            cat("\n")
            cat("\nSRE = list( quant = ", object@SRE$quant, "),", sep="")
            
            ## FDA options
            cat("\n")
            cat("\nFDA = list( method = '", object@FDA$method, "',", sep="")
            cat("\n            add_args = ", ifelse(length(object@FDA$add_args)<1,
                                                    'NULL', 
                                                    paste("list(", paste(.print.control(object@FDA$add_args), collapse=""), ")", sep="")), "),",sep="")
            
            ## MARS options
            cat("\n")
            cat("\nMARS = list( type = '", object@MARS$type, "',", sep="")
            cat("\n             interaction.level = ", object@MARS$interaction.level, ",", sep="")
            cat("\n             myFormula = ",  ifelse(length(object@MARS$myFormula) < 1,'NULL',paste(object@GLM$myFormula[2],object@GLM$myFormula[1],object@GLM$myFormula[3])), ",", sep="") 
#             cat("\n             degree = ", object@MARS$degree, ",", sep="")
            cat("\n             nk = ", ifelse(length(object@MARS$nk) < 1,'NULL',object@MARS$nk), ",", sep="")
            cat("\n             penalty = ", object@MARS$penalty, ",", sep="")
            cat("\n             thresh = ", object@MARS$thresh, ",", sep="")
            cat("\n             nprune = ", ifelse(length(object@MARS$nprune) < 1,'NULL',object@MARS$nprune), ",", sep="")
            cat("\n             pmethod = '", object@MARS$pmethod, "'),", sep="")
            
            ## RF options
            cat("\n")
            cat("\nRF = list( do.classif = ", object@RF$do.classif, ",", sep="")
            cat("\n           ntree = ", object@RF$ntree, ",", sep="")
            cat("\n           mtry = '", object@RF$mtry, "',", sep="")
            cat("\n           nodesize = ", object@RF$nodesize, ",", sep="")
            cat("\n           maxnodes = ", ifelse(length(object@RF$maxnodes) < 1,'NULL',object@RF$maxnodes), "),", sep="")
            
            ## MAXENT.Phillips options
            cat("\n")
            cat("\nMAXENT.Phillips = list( path_to_maxent.jar = '", object@MAXENT.Phillips$path_to_maxent.jar, "',", sep="")
            cat("\n               memory_allocated = ", ifelse(length(object@MAXENT.Phillips$memory_allocated) < 1,'NULL',object@MAXENT.Phillips$memory_allocated), ",", sep="")
            cat("\n               background_data_dir = ", ifelse(is.character(object@MAXENT.Phillips$background_data_dir), "'", ""), object@MAXENT.Phillips$background_data_dir, ifelse(is.character(object@MAXENT.Phillips$background_data_dir), "'", ""), ",", sep="")
            cat("\n               maximumbackground = ", ifelse(is.character(object@MAXENT.Phillips$maximumbackground), "'", ""), object@MAXENT.Phillips$maximumbackground, ifelse(is.character(object@MAXENT.Phillips$maximumbackground), "'", ""), ",", sep="")
            cat("\n               maximumiterations = ", object@MAXENT.Phillips$maximumiterations, ",", sep="")
            cat("\n               visible = ", object@MAXENT.Phillips$visible, ",", sep="")
            cat("\n               linear = ", object@MAXENT.Phillips$linear, ",", sep="")
            cat("\n               quadratic = ", object@MAXENT.Phillips$quadratic, ",", sep="")
            cat("\n               product = ", object@MAXENT.Phillips$product, ",", sep="")
            cat("\n               threshold = ", object@MAXENT.Phillips$threshold, ",", sep="")
            cat("\n               hinge = ", object@MAXENT.Phillips$hinge, ",", sep="")
            cat("\n               lq2lqptthreshold = ", object@MAXENT.Phillips$lq2lqptthreshold, ",", sep="")
            cat("\n               l2lqthreshold = ", object@MAXENT.Phillips$l2lqthreshold, ",", sep="")
            cat("\n               hingethreshold = ", object@MAXENT.Phillips$hingethreshold, ",", sep="")
            cat("\n               beta_threshold = ", object@MAXENT.Phillips$beta_threshold, ",", sep="")
            cat("\n               beta_categorical = ", object@MAXENT.Phillips$beta_categorical, ",", sep="")
            cat("\n               beta_lqp = ", object@MAXENT.Phillips$beta_lqp, ",", sep="")
            cat("\n               beta_hinge = ", object@MAXENT.Phillips$beta_hinge, ",", sep="")
            cat("\n               betamultiplier = ", object@MAXENT.Phillips$betamultiplier, ",", sep="")
            cat("\n               defaultprevalence = ", object@MAXENT.Phillips$defaultprevalence, "),", sep="")

            ## MAXENT.Tsuruoka
            cat("\n")
            cat("\nMAXENT.Tsuruoka = list( l1_regularizer = ", object@MAXENT.Tsuruoka$l1_regularizer, ",", sep="")
            cat("\n                        l2_regularizer = ", object@MAXENT.Tsuruoka$l2_regularizer, ",", sep="")
            cat("\n                        use_sgd = ", object@MAXENT.Tsuruoka$use_sgd, ",", sep="")
            cat("\n                        set_heldout = ", object@MAXENT.Tsuruoka$set_heldout, ",", sep="")
            cat("\n                        verbose = ", object@MAXENT.Tsuruoka$verbose, ")", sep="")
            
            .bmCat()
          })

.print.control <- function(ctrl){
  out <-  paste(names(ctrl)[1], " = ", ctrl[[1]], sep="")
  
  if(length(ctrl) > 1){
    i=2
    while(i <= length(ctrl)){
      if(is.list(ctrl[[i]])){
        out <- c(out, paste(", ", names(ctrl)[i], " = list(",
                            paste(names(ctrl[[i]]), "=",unlist(ctrl[[i]]), sep="", collapse=", "),")", sep=""))
        #         i <- i+length(ctrl[[i]])
        i <- i+1
      } else {
        out <- c(out, paste(", ", names(ctrl)[i], " = ", ctrl[[i]], sep=""))
        i <- i+1
      }
    }    
  }
  #   return(toString(out))
  return(out)
  
}
####################################################################################################
### BIOMOD Storing Results Objects #################################################################
####################################################################################################
setClass("BIOMOD.stored.data",
         representation(inMemory = 'logical',
                        link = 'character'),
         prototype(inMemory=FALSE,
                   link = ''),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.array",
         contains = "BIOMOD.stored.data",
         representation(val = 'array'),
         prototype(val = array()),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.data.frame",
         contains = "BIOMOD.stored.data",
         representation(val = 'data.frame'),
         prototype(val = data.frame()),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.raster.stack",
         contains = "BIOMOD.stored.data",
         representation(val = 'RasterStack'),
         prototype(val = stack()),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.files",
         contains = "BIOMOD.stored.data",
         representation(val = 'character'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.formated.data",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.formated.data'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.models.options",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.Model.Options'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setMethod("load_stored_object", "BIOMOD.stored.data",
          function(obj){
            if(obj@inMemory){
              return(obj@val)
            }
            
            # different comportement with raster
            if(inherits(obj, "BIOMOD.stored.raster.stack")){
              if( length(obj@link) == 1 & all(grepl(".RData", obj@link)) ){
                return(get(load(obj@link)))
              } else if(all(grepl(".grd", obj@link) | grepl(".img", obj@link))){
                out <- raster::stack(x = obj@link, RAT=FALSE)
                ## rename layer in case of individual projections
                if(all(grepl("individual_projections",obj@link))){
                  # remove directories arch and extention
                  xx <- sub("[:.:].+$", "", sub("^.+individual_projections/", "",obj@link))
                  # remove projection name
                  to_rm <- unique(sub("[^_]+[:_:][^_]+[:_:][^_]+[:_:][^_]+$", "", xx))                  
                  xx <- sub(to_rm,"", xx)
                  names(out) <- xx
                }
                return(out)
              } # else {
              #                 filesToLoad <- list.files(path=sub("/individual_projections","", obj@link), full.names=T)
              #                 toMatch <- c('.grd$','.img$')
              #                 filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)  
              #                 if(length(filesToLoad)){
              #                   return(raster::stack(filesToLoad[1]))
              #                 } else {
              #                   filesToLoad <- list.files(path=obj@link, full.names=T)
              #                   toMatch <- c('.grd$','.img$')
              #                   filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
              #                   toMatch <- obj@models.projected
              #                   filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
              #                   proj <- raster::stack(filesToLoad)
              #                   toMatch <- c(obj@link,".img$",'.grd$', .Platform$file.sep)
              #                   names(proj) <- gsub(pattern=paste(toMatch,collapse="|"), "", filesToLoad)
              #                   return(proj)
              #                 }   
              #               } 
            } else { # for all other stored objects
              return(get(load(obj@link)))
            }
          }
          
)



setClass("BIOMOD.models.out",
         representation(modeling.id = 'character', 
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.computed = 'character',
                        models.failed = 'character',
                        has.evaluation.data = 'logical',
                        rescal.all.models = 'logical',
                        models.evaluation = 'BIOMOD.stored.array',
                        variables.importances = 'BIOMOD.stored.array',
                        models.prediction = 'BIOMOD.stored.array',
                        models.prediction.eval = 'BIOMOD.stored.array',
                        formated.input.data = 'BIOMOD.stored.formated.data',
                        calib.lines = 'BIOMOD.stored.array',
                        models.options = 'BIOMOD.stored.models.options',
                        link = 'character'),
         prototype(modeling.id = as.character(format(Sys.time(), "%s")),
                   sp.name='',
                   expl.var.names = '',
                   models.computed='',
                   models.failed='',
                   has.evaluation.data=FALSE,
                   rescal.all.models=TRUE, 
                   models.evaluation = new('BIOMOD.stored.array'),
                   variables.importances = new('BIOMOD.stored.array'),
                   models.prediction = new('BIOMOD.stored.array'),
                   models.prediction.eval = new('BIOMOD.stored.array'),
                   formated.input.data = new('BIOMOD.stored.formated.data'),
                   calib.lines = new('BIOMOD.stored.array'),
                   models.options = new('BIOMOD.stored.models.options'),
                   link=''),
         validity = function(object){
           return(TRUE)
         })

setClass("BIOMOD.stored.models.out",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.out'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

setMethod('show', signature('BIOMOD.models.out'),
          function(object){
            .bmCat("BIOMOD.models.out")
            cat("\nModeling id :", object@modeling.id, fill=.Options$width)
            cat("\nSpecies modeled :", object@sp.name, fill=.Options$width)
            cat("\nConsidered variables :", object@expl.var.names, fill=.Options$width)
            
            cat("\n\nComputed Models : ", object@models.computed, fill=.Options$width)
            cat("\n\nFailed Models : ", object@models.failed, fill=.Options$width)
            .bmCat()
          })


setClass("BIOMOD.stored.models.out",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.out'),
         prototype(val = NULL),
         validity = function(object){
           return(TRUE)
         })

### GETTEURS ###

setMethod("get_predictions", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE, evaluation = FALSE){
            # check evaluation data avialability
            if( evaluation & (! obj@has.evaluation.data) ){
              warning("calibration data returned because no evaluation data available")
              evaluation = FALSE
            }
            
            # select calibration or eval data
            if(evaluation) pred <- obj@models.prediction.eval else pred <- obj@models.prediction
            
            if(!as.data.frame){
              if(pred@inMemory ){
                return(pred@val)
              } else{
                if(pred@link != ''){
                  return(get(load(pred@link)))
                } else{ return(NULL) }
              }              
            } else {
              if(pred@inMemory ){
                mod.pred <- as.data.frame(pred@val)
                names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                                 function(x){
                                                   x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                                   data.set.id <- x.rev[1]
                                                   cross.valid.id <- x.rev[2]
                                                   algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                                   model.id <- paste(obj@sp.name,
                                                                     data.set.id,
                                                                     cross.valid.id,
                                                                     algo.id, sep="_")
                                                   return(model.id)
                                                 }))
                return(mod.pred)
              } else{
                if(pred@link != ''){
                  mod.pred <- as.data.frame(get(load(pred@link)))    
                  names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                                   function(x){
                                                     x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                                     data.set.id <- x.rev[1]
                                                     cross.valid.id <- x.rev[2]
                                                     algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                                     model.id <- paste(obj@sp.name,
                                                                       data.set.id,
                                                                       cross.valid.id,
                                                                       algo.id, sep="_")
                                                     return(model.id)
                                                   }))
                  return(mod.pred)
                } else{ return(NULL) }
              }
              
            }
          }
)

setMethod("get_evaluations", "BIOMOD.models.out",
          function(obj, ...){
            args <- list(...)
            
            ## fill some additional parameters 
            as.data.frame <- ifelse(!is.null(args$as.data.frame), args$as.data.frame, FALSE)
            
            out <- NULL
            if(obj@models.evaluation@inMemory ){
              out <- obj@models.evaluation@val
            } else{
              if(obj@models.evaluation@link != ''){
                out <- get(load(obj@models.evaluation@link))          
              }
            }
            
            ## transform into data.frame object if needed
            if(as.data.frame){
              tmp <- reshape::melt.array(out,varnames=c("eval.metric", "test","m","r","d"))
              model_names <- unique(apply(tmp[,c("m","r","d"), drop=F], 1, paste, collapse="_"))
              out <- data.frame() #NULL
              for(mod in model_names){
                m = unlist(strsplit(mod,"_"))[1]
                r = unlist(strsplit(mod,"_"))[2]
                d = unlist(strsplit(mod,"_"))[3]
                eval.met = as.character(unique(tmp[which( tmp$m == m & tmp$r == r & tmp$d == d), "eval.metric", drop=T]))
                for(em in eval.met){
                  out <- rbind(out, 
                               data.frame( Model.name = mod,
                                           Eval.metric = em,
                                           Testing.data = as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Testing.data"), "value", drop=T]),
                                           Evaluating.data = ifelse("Evaluating.data" %in% tmp$test, as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Evaluating.data"), "value", drop=T]), NA ),
                                           Cutoff = as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Cutoff"), "value", drop=T]),
                                           Sensitivity = as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Sensitivity"), "value", drop=T]),
                                           Specificity = as.numeric( tmp[which( tmp$m == m & tmp$r == r & tmp$d == d & tmp$eval.metric == em & tmp$test == "Specificity"), "value", drop=T]) )
                  ) 
                } # end loop on eval metric
              } # end loop on models names
              
            }
            
            return(out)
          }
)


setMethod("get_calib_lines", "BIOMOD.models.out",
   function(obj, as.data.frame = FALSE, ...){
     calib_lines <- load_stored_object(obj@calib.lines)
     return(calib_lines)
   }
          
)


setMethod("get_variables_importance", "BIOMOD.models.out",
          function(obj, ...){
            if(obj@variables.importances@inMemory ){
              return(obj@variables.importances@val)
            } else{
              if(obj@variables.importances@link != ''){
                return(get(load(obj@variables.importances@link)))
              } else{ return(NA) }
            }
          }
)



setMethod("get_options", "BIOMOD.models.out",
          function(obj){
            if(obj@models.options@inMemory ){
              return(obj@models.options@val)
            } else{
              if(obj@models.options@link != ''){
                return(get(load(obj@models.options@link)))                
              } else{ return(NA) }
            }
          }
)

setMethod("get_formal_data", "BIOMOD.models.out",
          function(obj, subinfo = NULL){
            if(is.null(subinfo)){
              if(obj@formated.input.data@inMemory ){
                return(obj@formated.input.data@val)
              } else{
                if(obj@formated.input.data@link != ''){
                  data <- get(load(obj@formated.input.data@link))
                  return(data)
                } else{ return(NA) }
              }              
            } else if(subinfo == 'MinMax'){
              return(apply(get_formal_data(obj, "expl.var"),2, function(x){
                if(is.numeric(x)){
                  return( list(min = min(x,na.rm=T), max = max(x, na.rm=T) ) )
                } else if(is.factor(x)){
                  return(list(levels = levels(x)))
                }
              }) )
            } else if(subinfo == 'expl.var'){
              return(as.data.frame(get_formal_data(obj)@data.env.var))
            } else if(subinfo == 'expl.var.names'){
              return(obj@expl.var.names)
            } else if(subinfo == 'resp.var'){
              return(as.numeric(get_formal_data(obj)@data.species))
            } else if(subinfo == 'eval.resp.var'){
              return(as.numeric(get_formal_data(obj)@eval.data.species))
            } else if(subinfo == 'eval.expl.var'){
              return(as.data.frame(get_formal_data(obj)@eval.data.env.var))
            } else{
              stop("Unknow subinfo tag")
            }
            
          }
)

setMethod("get_built_models", "BIOMOD.models.out",
          function(obj, ...){
            return(obj@models.computed)
          }
)


setMethod("RemoveProperly", "BIOMOD.models.out",
          function(obj, obj.name=deparse(substitute(obj))){
            cat("\n\t> Removing .BIOMOD_DATA files...")
            unlink(file.path(obj@sp.name, ".BIOMOD_DATA", obj@modeling.id), recursive=T, force=TRUE)
            cat("\n\t> Removing models...")
            unlink(file.path(obj@sp.name, "models", obj@modeling.id), recursive=T, force=TRUE)
            cat("\n\t> Removing object hard drive copy...")
            unlink(obj@link, recursive=T, force=TRUE)
            cat("\n\t> Removing object from memory...")
            rm(list=obj.name,envir=sys.frame(-2))
            cat("\nCompleted!")
          })


####################################################################################################
### BIOMOD Storing Projection Objects ##############################################################
####################################################################################################
# setClass("BIOMOD.projection",
#          representation(proj.names = 'character',
#                         sp.name = 'character',
#                         expl.var.names = 'character',
#                         models.computed = 'character',
#                         models.failed = 'character',
#                         models.thresholds = 'BIOMOD.stored.array',
#                         models.prediction = 'BIOMOD.stored.array',
#                         formated.input.data = 'BIOMOD.stored.formated.data',
#                         calib.lines = 'BIOMOD.stored.array',
#                         models.options = 'BIOMOD.stored.models.options'),
#          prototype(sp.name='',
#                    expl.var.names = '',
#                    models.computed='',
#                    models.failed='',
#                    models.evaluation = new('BIOMOD.stored.array'),
#                    variables.importances = new('BIOMOD.stored.array'),
#                    models.prediction = new('BIOMOD.stored.array'),
#                    formated.input.data = new('BIOMOD.stored.formated.data'),
#                    calib.lines = new('BIOMOD.stored.array'),
#                    models.options = new('BIOMOD.stored.models.options')),
#          validity = function(object){
#            return(TRUE)
#            })

setClass("BIOMOD.projection.out",
         representation(proj.names = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.projected = 'character',
                        scaled.models = 'logical',
                        modeling.object = 'BIOMOD.stored.data',
                        modeling.object.id = 'character',
                        type = 'character',
                        proj = 'BIOMOD.stored.data',
                        xy.coord = 'matrix'),
         prototype(proj.names = '',
                   sp.name='',
                   expl.var.names='',
                   models.projected='',
                   scaled.models=TRUE,
                   modeling.object.id='',
                   type='',
                   xy.coord = matrix()),
         validity = function(object){
           return(TRUE)
         })

setMethod("get_projected_models", "BIOMOD.projection.out",
          function(obj){
            return(obj@models.projected)
          })

setMethod("get_predictions", "BIOMOD.projection.out",
          function(obj, as.data.frame=FALSE, full.name=NULL, model=NULL, run.eval=NULL, data.set=NULL){
            models_selected <- get_projected_models(obj)
            if(length(full.name)){
              models_selected <- intersect(full.name, models_selected)
            } else if(length(model) | length(run.eval) | length(data.set)){
              # models subselection according to model, run.eval and sata.set parameters
              if(length(model)) grep_model <- paste("(",paste(model,collapse="|"),")", sep="") else grep_model = "*"
              if(length(run.eval)) grep_run.eval <- paste("(",paste(run.eval,collapse="|"),")", sep="") else grep_run.eval = "*"
              if(length(data.set)) grep_data.set <- paste("(",paste(data.set,collapse="|"),")", sep="") else grep_data.set = "*"
              grep_full <- paste(grep_data.set,"_",grep_run.eval,"_",grep_model,"$",sep="")
              
              models_selected <- grep(pattern=grep_full, models_selected, value=T)
            }
            
            #             cat("\n*** models_selected = ", models_selected)
            
            if (length(models_selected)){
              proj <- load_stored_object(obj@proj)
              names(proj) <- obj@models.projected
              if(inherits(proj,'Raster')){
                proj <- raster::subset(proj,models_selected,drop=FALSE)
              } else {
                if(length(dim(proj)) == 4){ ## 4D arrays
                  proj <- proj[,.extractModelNamesInfo(model.names=models_selected,info='models'),
                               .extractModelNamesInfo(model.names=models_selected,info='run.eval'),
                               .extractModelNamesInfo(model.names=models_selected,info='data.set'), drop=FALSE]
                } else{ ## matrix (e.g. from ensemble models projections)
                  proj <- proj[,models_selected,drop=FALSE]
                }
                
              }
              
              if(as.data.frame){
                proj <- as.data.frame(proj)
                ## set correct names
                
                if(obj@type == 'array' & sum(!(names(proj) %in% models_selected))>0 ){ ## from array & not valid names
                  names(proj) <- unlist(lapply(strsplit(names(proj),".", fixed=TRUE), 
                                                   function(x){
                                                     x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                                     data.set.id <- x.rev[1]
                                                     cross.valid.id <- x.rev[2]
                                                     algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                                     model.id <- paste(obj@sp.name,
                                                                       data.set.id,
                                                                       cross.valid.id,
                                                                       algo.id, sep="_")
                                                     return(model.id)
                                                   }))
                }
                
                # reorder the data.frame
                proj <- proj[,models_selected]
                
              }  
            } else{
              proj <- NULL
            }
            
            return(proj)          
          })




# 2.3 other functions
setMethod('plot', signature(x='BIOMOD.projection.out', y="missing"),
          function(x,col=NULL, str.grep=NULL){
            models_selected <- x@models.projected 
            if(length(str.grep)){
              models_selected <- grep(paste(str.grep,collapse="|"), models_selected,value=T)
            } 
            
            if(!length(models_selected)) stop("invalid str.grep arg")
            
            
            
            if(class(x@proj) == "BIOMOD.stored.raster.stack"){
              requireNamespace("rasterVis")
              
              ## define the breaks of the color key
              my.at <- seq(0,1000,by=100)
              ## the labels will be placed vertically centered
              my.labs.at <- seq(0,1000,by=250)
              ## define the labels
              my.lab <- seq(0,1000,by=250)
              ## define colors
              #               my.col <- colorRampPalette(c("red4","orange4","yellow4","green4"))(100)
              my.col <- colorRampPalette(c("grey90","yellow4","green4"))(100)
              
              ## try to use levelplot function
              try_plot <- try(
                levelplot(get_predictions(x, full.name=models_selected),
                          at=my.at, margin=T, col.regions=my.col,
                          main=paste(x@sp.name,x@proj.names,"projections"),
                          colorkey=list(labels=list(
                            labels=my.lab,
                            at=my.labs.at)))
              ) 
              if(! inherits(try_plot,"try-error")){ ## produce plot
                print(try_plot)
              } else{## try classical plot
                cat("\nrasterVis' levelplot() function failed. Try to call standard raster plotting function.",
                    "It can lead to unooptimal representations.",
                    "You should try to do it by yourself extracting predicions (see : get_predictions() function)", fill=options()$width)
                try_plot <- try(
                  plot(get_predictions(x, full.name=models_selected))
                )
              } 
              
              if(inherits(try_plot,"try-error")){ # try classical plot
                cat("\n Plotting function failed.. You should try to do it by yourself!")
              }
              
            } else if(class(x@proj) == "BIOMOD.stored.array"){
              if(ncol(x@xy.coord) != 2){
                cat("\n ! Impossible to plot projections because xy coordinates are not available !")
              } else {
                multiple.plot(Data = get_predictions(x, full.name=models_selected, as.data.frame=T), coor = x@xy.coord)
              }
              
            } else {cat("\n !  Biomod Projection plotting issue !", fill=.Options$width)}
            
          })

setMethod('show', signature('BIOMOD.projection.out'),
          function(object){
            .bmCat("'BIOMOD.projection.out'")
            cat("\nProjection directory :", paste(object@sp.name,"/",object@proj.names, sep=""), fill=.Options$width)
            cat("\n")
            cat("\nsp.name :", object@sp.name, fill=.Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill=.Options$width)
            cat("\n")
            cat("\nmodeling id :", object@modeling.object.id ,"(",object@modeling.object@link ,")", fill=.Options$width)
            cat("\nmodels projected :", toString(object@models.projected), fill=.Options$width)
            
            .bmCat()
          })


setMethod(f='free', 
          signature='BIOMOD.projection.out',
          definition = function(obj){
            if(inherits(obj@proj,"BIOMOD.stored.array")){
              obj@proj@val = array()
            } else if(inherits(obj@proj,"BIOMOD.stored.raster.stack")){
              obj@proj@val = stack()
            } else{
              obj@proj@val = NULL
            }
            obj@proj@inMemory = FALSE
            
            return(obj)
          })

####################################################################################################
### BIOMOD Storing Ensemble Modeling Objects #######################################################
####################################################################################################   
setClass("BIOMOD.EnsembleModeling.out",
         representation(sp.name = 'character',
                        expl.var.names = 'character',
                        models.out.obj = 'BIOMOD.stored.models.out',
                        eval.metric = 'character',
                        eval.metric.quality.threshold = 'numeric',
                        em.computed = 'character',
                        em.by = 'character',
                        em.models = 'ANY',
                        modeling.id = 'character',
                        link = 'character'),
         #                         em.models.kept = 'list',
         #                         em.prediction = 'BIOMOD.stored.array',
         #                         em.evaluation = 'BIOMOD.stored.array',
         #                         em.res = 'list',
         #                         em.ci.alpha = 'numeric',
         #                         em.weight = 'list',
         #                         em.bin.tresh = 'list'),
         prototype( sp.name = '',
                    expl.var.names = '',
                    models.out.obj = new('BIOMOD.stored.models.out'),
                    eval.metric = '',
                    eval.metric.quality.threshold = 0,
                    em.models = NULL,
                    em.computed = character(),
                    modeling.id = '.',
                    link = ''),
         #                     em.models.kept = NULL,
         #                     em.prediction = NULL,
         #                     #                     em.evaluation = NULL,
         #                     em.res = list(),
         #                     em.ci.alpha = 0.05,
         #                     em.weight = list(),
         #                     em.bin.tresh = list()),
         validity = function(object){
           return(TRUE)
         })


setMethod('show', signature('BIOMOD.EnsembleModeling.out'),
          function(object){
            .bmCat("'BIOMOD.EnsembleModeling.out'")
            cat("\nsp.name :", object@sp.name, fill=.Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill=.Options$width)
            cat("\n")
            cat("\nmodels computed:", toString(object@em.computed), fill=.Options$width)
            
            .bmCat()
          })

setMethod("get_needed_models", "BIOMOD.EnsembleModeling.out",
          function(obj, subset='all', ...){
            add.args <- list(...)
            needed_models <- lapply(obj@em.models, function(x){
              return(x@model)
            })
            needed_models <- unique(unlist(needed_models))
            return(needed_models)
          }
)

setMethod("get_needed_models", "BIOMOD.EnsembleModeling.out",
          function(obj, selected.models='all', ...){           ### MODIFIED ROBIN : ici j'ai renom? "subset" en "selected.models" pour que ?a soit le m?me nom que dans les autres fonctions
            add.args <- list(...)  
            if(selected.models[[1]]=="all") selected.index <- c(1:length(obj@em.models)) else selected.index <- which(names(obj@em.models) %in% selected.models) ### MODIFIED ROBIN : ici je selectionne uniquement le sous-ensemble des mod?les que l'utilisateur a sp?cifi?.
            needed_models <- lapply(obj@em.models[selected.index], function(x) return(x@model))  ### MODIFIED ROBIN : ici je selectionne uniquement le sous-ensemble des mod?les que l'utilisateur a sp?cifi?.
            needed_models <- unique(unlist(needed_models))
            return(needed_models)
          }
)


setMethod("get_kept_models", "BIOMOD.EnsembleModeling.out",
          function(obj, model, ...){
            if(is.character(model) | is.numeric(model)){
              return(obj@em.models[[model]]@model)
            } else{
              kept_mod <- lapply(obj@em.models, function(x){return(x@model)})
              names(kept_mod) <- names(obj@em.models)
              return(kept_mod)
            }
            
          }
)

setMethod("get_evaluations", "BIOMOD.EnsembleModeling.out",
          function(obj, ...){
            args <- list(...)
            
            ## fill some additional parameters 
            as.data.frame <- ifelse(!is.null(args$as.data.frame), args$as.data.frame, FALSE)
            
            out <- list()
            
            ## list of computed models
            models <- obj@em.computed
            
            ## extract evaluation scores as a list
            for(mod in models){ 
              out[[mod]] <- obj@em.models[[mod]]@model_evaluation[,,drop=F]
            }
            
            ## transform into data.frame object if needed
            if(as.data.frame){
              tmp <- melt(out, varnames=c("eval.metric", "test"))
              tmp$model.name <- sapply(tmp$L1, function(x){paste(unlist(strsplit(x, "_"))[-1], collapse="_")})
              out <- data.frame() #NULL
              for(mod in unique(tmp$model.name)){
                eval.met = as.character(unique(tmp[which( tmp$model.name == mod ), "eval.metric", drop=T]))
                for(em in eval.met){
                  out <- rbind(out, 
                               data.frame( Model.name = mod,
                                           Eval.metric = em,
                                           Testing.data = as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Testing.data"), "value", drop=T]),
                                           Evaluating.data = ifelse("Evaluating.data" %in% tmp$test, as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Evaluating.data"), "value", drop=T]), NA ),
                                           Cutoff = as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Cutoff"), "value", drop=T]),
                                           Sensitivity = as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Sensitivity"), "value", drop=T]),
                                           Specificity = as.numeric( tmp[which( tmp$model.name == mod & tmp$eval.metric == em & tmp$test == "Specificity"), "value", drop=T]) )
                  ) 
                } # end loop on eval metric
              } # end loop on models names
            } # end as.data.frame == TRUE
            
            return(out)
          }
)


setMethod("get_variables_importance", "BIOMOD.EnsembleModeling.out",
          function(obj, ...){
            vi <- NULL
            for (mod in get_built_models(obj) ){
              (vi_tmp <- obj@em.models[[mod]]@model_variables_importance)
              vi <- abind::abind(vi, vi_tmp, along=3)
            }
            dimnames(vi)[[3]] <- get_built_models(obj)
            
            return(vi)
          }
)

setMethod("get_built_models", "BIOMOD.EnsembleModeling.out",
          function(obj, ...){
            return(obj@em.computed)
          })

setMethod("get_predictions", "BIOMOD.EnsembleModeling.out",
  function(obj, ...){
    ## note: ensemble models predicitons are stored within the directory 
    ##  <sp.name>/.BIOMOD_DATA/<modelling.id>/ensemble.models/ensemble.models.projections/
    ##  This function is just a friendly way to load this data
    
    ## get the path to projections files we want to load
    files.to.load <- file.path(obj@sp.name, ".BIOMOD_DATA", obj@modeling.id, "ensemble.models", 
                              "ensemble.models.predictions", paste0(obj@em.computed, ".predictions"))
    ## load and merge projection files within a data.frame
    bm.pred <- do.call(cbind, lapply(files.to.load, function(ftl) get(load(ftl))))
    colnames(bm.pred) <- obj@em.computed
    return(bm.pred)
  }
)

"/home/georgeda/Work/BIOMOD/RForge/tests/workdir/GuloGulo/.BIOMOD_DATA/test//GuloGulo_EMmeanByTSS_mergedAlgo_mergedRun_mergedData.predictions"
####################################################################################################
### BIOMOD Storing Ensemble Forecasting Objects ####################################################
####################################################################################################

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# if( !isGeneric( ".Models.prepare.data" ) ) {
# }

setMethod('.Models.prepare.data', signature(data='BIOMOD.formated.data'),
          function(data, NbRunEval, DataSplit, Yweights=NULL, Prevalence=NULL, do.full.models=TRUE, DataSplitTable=NULL){
            list.out <- list()
            name <- paste(data@sp.name,'_AllData',sep="")
            xy <- data@coord
            dataBM <- data.frame(cbind(data@data.species,data@data.env.var))
            colnames(dataBM)[1] <- data@sp.name
            
            # dealing with evaluation data
            if(data@has.data.eval){
              evalDataBM <- data.frame(cbind(data@eval.data.species,data@eval.data.env.var[,,drop=FALSE]))
              colnames(evalDataBM)[1] <- data@sp.name
              eval.xy <- data@eval.coord
            } else{ evalDataBM <- eval.xy <- NULL }
            
            ### Calib/Valid lines
            if(!is.null(DataSplitTable)){
              calibLines <- DataSplitTable
              colnames(calibLines) <- paste('_RUN',1:ncol(calibLines), sep='')
            } else {
              if(NbRunEval == 0){ # take all available data
                calibLines <- matrix(rep(TRUE,length(data@data.species)),ncol=1)
                colnames(calibLines) <- '_Full'
              } else {
                calibLines <- .SampleMat(data.sp = data@data.species, 
                                         dataSplit = DataSplit, 
                                         nbRun = NbRunEval, 
                                         data.env = data@data.env.var)                    
                if(do.full.models){
                  calibLines <- cbind(calibLines, rep(TRUE,length(data@data.species)))
                  colnames(calibLines)[NbRunEval+1] <- '_Full'
                }
              }
            }
            ## force calib.lines object to be 3D array
            if(length(dim(calibLines)) < 3 ){
              dn_tmp <- dimnames(calibLines) ## keep track of dimnames
              dim(calibLines) <- c(dim(calibLines),1)
              dimnames(calibLines) <- list(dn_tmp[[1]], dn_tmp[[2]], "_AllData")
            }
            
            if(is.null(Yweights)){ # 1 for all points
              if(!is.null(Prevalence)){
                cat("\n\t> Automatic weights creation to rise a", Prevalence,"prevalence")
                Yweights <- .automatic_weights_creation(data@data.species ,prev=Prevalence)
              } else{
                cat("\n\t> No weights : all observations will have the same weight")
                Yweights <- rep(1,length(data@data.species))
              }
              
            }
            list.out[[name]] <- list(name=name,
                                     xy=xy,
                                     dataBM=dataBM,
                                     calibLines=calibLines,
                                     Yweights = Yweights,
                                     evalDataBM = evalDataBM,
                                     eval.xy = eval.xy)
            return(list.out)
          })

setMethod('.Models.prepare.data', signature(data='BIOMOD.formated.data.PA'),
          function(data, NbRunEval, DataSplit, Yweights=NULL, Prevalence=NULL, do.full.models=TRUE, DataSplitTable=NULL){
            list.out <- list()
            formal_weights <- Yweights
            for(pa in 1:ncol(data@PA)){
              Yweights <- formal_weights
              name <- paste(data@sp.name,"_",colnames(data@PA)[pa],sep="")
              xy <- data@coord[data@PA[,pa],]
              resp <- data@data.species[data@PA[,pa]] # response variable (with pseudo absences selected)
              resp[is.na(resp)] <- 0
              dataBM <- data.frame(cbind(resp,
                                         data@data.env.var[data@PA[,pa],,drop=FALSE]))
              colnames(dataBM)[1] <- data@sp.name
              
              ### Calib/Valid lines
              if(!is.null(DataSplitTable)){                
                if(length(dim(DataSplitTable))==2){
                  calibLines <- DataSplitTable
                } else {
                  calibLines <- asub(DataSplitTable,pa,3,drop=TRUE)
                }
                colnames(calibLines) <- paste('_RUN',1:ncol(calibLines), sep='')
                calibLines[which(!data@PA[,pa]),] <- NA
              } else{
                if(NbRunEval == 0){ # take all available data
                  calibLines <- matrix(NA,nrow=length(data@data.species),ncol=1)
                  calibLines[data@PA[,pa],1] <- TRUE
                  colnames(calibLines) <- '_Full'
                } else {
                  calibLines <- matrix(NA,nrow=length(data@data.species),ncol=NbRunEval)
                  sampled.mat <- .SampleMat(data.sp = data@data.species[data@PA[,pa]],
                                            dataSplit = DataSplit,
                                            nbRun = NbRunEval,
                                            data.env = data@data.env.var[data@PA[,pa], , drop = FALSE]) 
                  calibLines[data@PA[,pa],] <- sampled.mat
                  colnames(calibLines) <- colnames(sampled.mat)
                  if(do.full.models){
                    calibLines <- cbind(calibLines, rep(NA,length(data@data.species)))
                    calibLines[data@PA[,pa],NbRunEval+1] <- TRUE
                    colnames(calibLines)[NbRunEval+1] <- '_Full'
                  }                
                }
              }

              ## force calib.lines object to be 3D array
              if(length(dim(calibLines)) < 3 ){
                dn_tmp <- dimnames(calibLines) ## keep track of dimnames
                dim(calibLines) <- c(dim(calibLines),1)
                dimnames(calibLines) <- list(dn_tmp[[1]], dn_tmp[[2]], paste("_PA",pa, sep=""))
              }
                
              # dealing with evaluation data
              if(data@has.data.eval){
                evalDataBM <- data.frame(cbind(data@eval.data.species,data@eval.data.env.var))
                colnames(evalDataBM)[1] <- data@sp.name
                eval.xy <- data@eval.coord
              } else{ evalDataBM <- eval.xy <- NULL }
              
              if(is.null(Yweights)){ # prevalence of 0.5... may be parametrize
                if(is.null(Prevalence)) Prevalence <- 0.5
                
                cat("\n\t\t\t! Weights where automaticly defined for", name, "to rise a", Prevalence, "prevalence !")
                
                
                Yweights <- rep(NA, length(data@data.species))
                Yweights[data@PA[,pa]] <- .automatic_weights_creation(as.numeric(dataBM[,1]) ,prev=Prevalence)#, subset=data@PA[,pa])
              } else{
                # remove useless weights
                Yweights[!data@PA[,pa]] <- NA
              }
              
              list.out[[name]] <- list(name=name,
                                       xy=xy,
                                       dataBM=dataBM,
                                       calibLines=calibLines,
                                       Yweights = Yweights,
                                       evalDataBM = evalDataBM,
                                       eval.xy = eval.xy)
            }
            return(list.out)
          })


.automatic_weights_creation <- function(resp,prev=0.5, subset=NULL){
  if(is.null(subset)) subset<- rep(TRUE, length(resp))
  
  nbPres <- sum(resp[subset], na.rm=TRUE)
  nbAbsKept <- sum(subset, na.rm=T) - sum(resp[subset], na.rm=TRUE) # The number of true absences + pseudo absences to maintain true value of prevalence
  Yweights <- rep(1,length(resp))
  
  if(nbAbsKept > nbPres){ # code absences as 1
    Yweights[which(resp>0)] <- (prev * nbAbsKept) / (nbPres * (1-prev))
  } else{ # code presences as 1
    Yweights[which(resp==0 | is.na(resp))] <- (nbPres * (1-prev)) / (prev * nbAbsKept)
  }
  Yweights = round(Yweights[])
  Yweights[!subset] <- 0
  
  return(Yweights)
}
