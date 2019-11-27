.pseudo.absences.sampling <-
function(sp, env, nb.repet=1, strategy='random', distMin=0, distMax=NULL, nb.points=NULL, quant.SRE = 0, PA.table = NULL){

  # 1. Parameters checking
  args <- .check.params.pseudo.absences.sampling(sp, env, nb.repet, strategy, distMin, distMax, nb.points, quant.SRE)

  sp <- args$sp
  env <- args$env
  nb.repet <- args$nb.repet
  strategy <- args$strategy
  distMin <- args$distMin
  distMax <- args$distMax
  nb.points <- args$nb.points
  quant.SRE <- args$quant.SRE

  rm("args")

  if( (nb.repet == 0 | nb.points <= 0) & strategy != 'user.defined'){
    out <- NULL
  } else {
    out <- switch(strategy,
                   user.defined = user.defined.pseudo.abs.selection(sp, env, PA.table),
                   random = random.pseudo.abs.selection( sp, env, nb.points, nb.repet ),
                   sre = sre.pseudo.abs.selection(sp, env, quant.SRE, nb.points, nb.repet),
                   disk = disk.pseudo.abs.selection(sp, env, distMin, distMax, nb.points, nb.repet))
  }

  return(out)

#   # 2. Check if NA are present in sp or not to determine which dataset to use
#   if(sum(is.na(sp@data)) > 0 ){ # PA will be taken into response variable
#     cat("\n*** PA selection")
#     pa.tab <- switch(strategy,
#                      random = random.pseudo.abs.selection(data=sp, nb.points=nb.points, nb.repet=nb.repet),
#                      sre = sre.pseudo.abs.selection(sp),
#                      disk = disk.pseudo.abs.selection(sp))
#     .arranging.pa.table()
#   } else{ # PA will be taken into explanatory variables
#     if(inherits(env, 'Raster')){ # Raster env var case
#
#     } else if(inherits(env, 'SpatialPoints')){ # spatial data.frame case
#
#     }
#
#   }
}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.check.params.pseudo.absences.sampling <- function(sp, env, nb.repet, strategy, distMin, distMax, nb.points, quant.SRE){
  cat("\n   > Pseudo Absences Selection checkings...")

  # define here the implemented strategies
  availableStrategies <- c("random", "sre", "disk", "user.defined")

  # 1. sp input checking
  if(is.vector(sp)){
    sp <- SpatialPointsDataFrame(matrix(0,ncol=2,nrow=length(sp)), data.frame(sp))
  }

  if(!(inherits(sp, 'SpatialPoints'))){
    stop("species input must be a SpatialPointsDataFrame object")
  }

  # 2. env input checking
  if(is.matrix(env) | is.data.frame(env)){
    if(nrow(env) != nrow(sp)){
      stop("Species and Explanatory must have same dimensions")
    }
    env <- SpatialPointsDataFrame(coordinates(sp), as.data.frame(env))
  }

  if(!inherits(env, 'SpatialPoints') & !inherits(env, 'Raster')){
    stop("Explanatory variables input must be a SpatialPointsDataFrame or a RasterStack object")
  }

  # 3. Strategy checking
  if( ! (strategy %in% c("random", "sre", "user.defined")) ){
    if( ( sum(abs(coordinates(sp))) == 0 ) | !( strategy %in% availableStrategies ) ){ # no coordinates or unknow strategy
      strategy <- "random"
      cat("\n   ! Random strategy was automatically selected (that can be due to points coordinates lack or unavailable strategy choosen)")
    }
  }

  # 4. Nb points checking
  if(strategy != "user.defined"){
    if(is.null(nb.points)){
      stop("You must give the number of pseudo absences you want")
    } else{
      nbTrueAbs <- .get.nb.true.abs(sp)
      if(nbTrueAbs >= nb.points){
        cat("\n    ! There is more 'true absences' than desired pseudo absences. No pseudo absences selection done.")
        nb.points = 0
        #       #### Return a flag that tell to function that no PA selected
        #       return(NULL)
      } else {
        nb.points = nb.points - nbTrueAbs
      }
    }
  }


  # 4. Nb repetition checking

  # 5. Distances checking
  if(!is.null(distMin)){
    if(distMin < 0){
        distMin <- 0
    }
  }

  if(!is.null(distMax)){
    if(distMax < 0){
        distMax <- NULL
    }
  }

  if(!is.null(distMax) & !is.null(distMin)){
    if(distMin >= distMax){
      stop("distMin >= distMax")
    }
  }

  # 6. SRE quantil checking
  if(strategy == 'SRE'){
    if( quant.SRE >= 0.5 | quant.SRE <0 ){
      stop("\n    ! SRE Quant should be a value between 0 and 0.5 ")
    }
  }

  # 7. return checked params
  return(list(sp = sp,
              env = env,
              nb.repet = nb.repet,
              strategy = strategy,
              distMin = distMin,
              distMax = distMax,
              nb.points = nb.points,
              quant.SRE = quant.SRE))

}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.get.nb.true.abs <- function(sp){
  if(is.vector(sp)) return(sum(sp==0, na.rm=TRUE))

  if(inherits(sp, 'SpatialPoints')) return(sum(sp@data==0, na.rm=TRUE))

  if(inherits(sp, 'Raster')) return(sum(sp[]==0, na.rm=TRUE))
}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.nb.available.pa.cells <- function(data, PA.flag = NA){
  if(is.vector(data)){
    return(ifelse(is.na(PA.flag), sum(is.na(data)), sum(data == PA.flag, na.rm = TRUE)))
  }
  if(is.data.frame(data) | is.matrix(data)){
    return(ifelse(is.na(PA.flag), sum(is.na(data)), sum(data == PA.flag, na.rm = TRUE)))
  }
  if(inherits(data, 'SpatialPoints')){
    return(ifelse(is.na(PA.flag), sum(is.na(data@data)), sum(data@data == PA.flag, na.rm = TRUE)))
  }
  if(inherits(data, 'Raster')){
    return(ifelse(is.na(PA.flag), sum(is.na(data[])), sum(data[] == PA.flag, na.rm = TRUE)))
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.rand.pseudo.abs.selection <- function(data, nb.points){
  if(is.vector(data)){ return(sample(1:length(data), nb.points, replace=FALSE)) }

  if(inherits(data, 'SpatialPoints')){ return(sample(1:nrow(data@data), nb.points, replace=FALSE))}

  if(inherits(data, 'Raster')){ return(sort(sampleRandom(x=data, size=nb.points, cells=T)[,"cell"]))}
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# if( !isGeneric( "random.pseudo.abs.selection" ) ) {
  setGeneric( "random.pseudo.abs.selection",
              def = function(sp,env, ...){
                      standardGeneric( "random.pseudo.abs.selection" )
                      } )
# }

setMethod('random.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function( sp, env, nb.points, nb.repet ){
            cat("\n   > random pseudo absences selection")

            # 1. Check if NA are present in sp or not to determine which dataset to use
            if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
              nb.cells <- .nb.available.pa.cells(sp)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }
              pa.tab <- matrix(FALSE, ncol=nb.repet, nrow=nrow(sp))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              # select always the presences and the true absences
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              # and a subset of candidates cells
              cand.cells <- which(is.na(sp@data))
              for(j in 1:ncol(pa.tab)){
                ## force to get at least one value of each factorial variable
                fact.level.cells <- sample.factor.levels(x = as.data.frame(env),
                                                         mask.out = pa.tab[, j, drop = FALSE])
                if(length(fact.level.cells)){
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

setMethod('random.pseudo.abs.selection', signature(env="RasterStack"),
          function( sp, env, nb.points, nb.repet ){
#             require('raster',quietly=T)
            cat("\n   > random pseudo absences selection")

            # 1. Check if NA are present in sp or not to determine which dataset to use
            if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
              nb.cells <- .nb.available.pa.cells(sp)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }
              pa.tab <- matrix(FALSE, ncol=nb.repet, nrow=nrow(sp))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              # select always the presences and the true absences
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              # and a subset of candidates cells
              cand.cells <- which(is.na(sp@data))
              for(j in 1:ncol(pa.tab)){
                pa.tab[sample(x=cand.cells,size=nb.points,replace=FALSE),j] <- TRUE
              }
              env <- as.data.frame(extract(env, coordinates(sp)))

              return(list(xy = coordinates(sp),
                          sp = as.numeric(unlist(sp@data, use.names=FALSE)),
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")
              # create a mask containing all not already sampled points (presences and absences)
              mask.env <- mask.out <- raster::subset(env, 1, drop = TRUE)
              mask.env <- raster::reclassify(mask.env, c(-Inf, Inf, -1)) ## the area we want to sample
              mask.out[] <- NA

              # add presences and true absences in our mask
              in.cells <- cellFromXY(mask.env, coordinates(sp))
              mask.env[in.cells] <- NA
              mask.out[in.cells] <- 1

              # checking of nb candidates
              nb.cells <- .nb.available.pa.cells(mask.env, PA.flag = -1)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }

              # select cells into raster
              pa.tab.tmp <- matrix(NA, ncol = nb.repet, nrow = nb.points)
              for( j in 1:ncol(pa.tab.tmp)){
                ## initialise the vector of sample cells
                SR <- NULL
                ## define a compy of the sampling mask
                mask.env.tmp <- mask.env
                ## force to get at least one value of each factorial variable
                fact.level.cells <- sample.factor.levels(env, mask.out = mask.out)
                if(length(fact.level.cells)){
                  SR <- c(SR, fact.level.cells)
                  ## update the mask by removing already selected cells
                  mask.env.tmp[SR] <- NA
                }
                SR <- c(SR, sampleRandom(x = mask.env.tmp,
                                         size = nb.points - length(SR),
                                         cells = T,
                                         na.rm = T)[, "cell", drop = TRUE])
                ## repeat sampling until haing the right number of points
                while(length(SR)<nb.points){
                  ## update the mask by removing already selected cells
                  mask.env.tmp[SR] <- NA
                  ## extract the missing number of points
                  SR <- c(SR, sampleRandom(x = mask.env.tmp,
                                           size = nb.points - length(SR),
                                           cells = T,
                                           na.rm = T)[, "cell", drop = TRUE])
                }
                pa.tab.tmp[,j] <- SR
              }

              # puting cells in good format
              selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
              pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              for( j in 1:ncol(pa.tab)){
                pa.tab[is.element(selected.cells, pa.tab.tmp[,j]), j] <- TRUE
              }

              # puting presences, true absences and pseudo absences together
              xy <- rbind(coordinates(sp), xyFromCell(mask.env, selected.cells))
              xy <- .add_PA_rownames(xy)
              sp <- as.numeric(unlist(c(as.vector(sp@data),
                                        rep(NA, length(selected.cells))),
                                      use.names = FALSE))
              env <- extract(env, xy)

              pa.tab <- rbind(matrix(TRUE, nrow = (nrow(xy) - length(selected.cells)),
                                     ncol = ncol(pa.tab)), pa.tab)

              return(list(xy = xy,
                          sp = sp,
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))

            }
          })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# if( !isGeneric( "random.pseudo.abs.selection" ) ) {
setGeneric( "user.defined.pseudo.abs.selection",
            def = function(sp,env, ...){
              standardGeneric( "user.defined.pseudo.abs.selection" )
            } )
# }

setMethod('user.defined.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function( sp, env, pa.table ){
            cat("\n   > User defined pseudo absences selection")

              return(list(xy = coordinates(sp),
                          sp = as.vector(sp@data),
                          env = as.data.frame(env@data),
                          pa.tab = pa.table))

          })

setMethod('user.defined.pseudo.abs.selection', signature(env="RasterStack"),
          function( sp, env, pa.table ){
#             require('raster',quietly=T)
            cat("\n   > User defined pseudo absences selection")

#             env <- as.data.frame(extract(env, coordinates(sp), method='bilinear'))
              env <- as.data.frame(extract(env, coordinates(sp)))


            return(list(xy = coordinates(sp),
                        sp = as.numeric(unlist(sp@data, use.names=FALSE)),
                        env = as.data.frame(env),
                        pa.tab = pa.table))

          })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# if( !isGeneric( "random.pseudo.abs.selection" ) ) {
  setGeneric( "sre.pseudo.abs.selection",
              def = function(sp,env, ...){
                      standardGeneric( "sre.pseudo.abs.selection" )
                      } )
# }
setMethod('sre.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function(sp, env, quant.SRE, nb.points, nb.repet){
            cat("\n   > SRE pseudo absences selection")

            # 1. calculate sre to determine availables
            mask.in <- sre(Response = sp, Explanatory = env, NewData = env@data, Quant = quant.SRE)
            ## we want to sample PA out of the SRE => have to revert the mask
            mask.in <- data.frame(mask.in = !as.logical(mask.in))
#             # mask of already sampled points (presneces/absences)
#             mask[mask == 0] <- NA
#             mask[which(as.vector(sp@data)==1),1] <- 1
#             mask[which(as.vector(sp@data)==0),1] <- 0

            # 2. Check if NA are present in sp or not to determine which dataset to use
#             if(.nb.available.pa.cells(mask) > 0 ){ # PA will be taken into response variable
              nb.cells <- .nb.available.pa.cells(mask.in$mask.in, PA.flag = TRUE)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }
              pa.tab <- matrix(FALSE, ncol=nb.repet, nrow=nrow(sp))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              # select always the presences and the true absences
              pa.tab[c(which(sp@data == 1), which(sp@data == 0)),] <- TRUE
              # and a subset of candidates cells
              cand.cells <- which(!mask.in$mask.in)
              for(j in 1:ncol(pa.tab)){
                ## force to get at least one value of each factorial variable
                fact.level.cells <- sample.factor.levels(as.data.frame(env),
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

#             }
          })



setMethod('sre.pseudo.abs.selection', signature(env="RasterStack"),
          function(sp, env, quant.SRE, nb.points, nb.repet){
            cat("\n   > SRE pseudo absences selection")

            # 1. calculate sre to determine availables
            ## mask in which we want to sample
            mask.in <- sre(Response = sp, Explanatory = env, NewData = env, Quant = quant.SRE)
            ## remove all points that are in the cpecies SRE
            mask.in[mask.in[] > 0] <- NA

            ## mask of already sampled points (presences/absences)
            mask.out <- raster::subset(env, 1)
            mask.out[] <- NA; mask.out[cellFromXY(mask.out,
                                                  coordinates(sp)[is.element(as.vector(sp@data), c(0, 1)), ])]

#             # removing cells in envelops, presences and absences
#             mask[mask == 1] <- NA
#             mask[cellFromXY(mask, coordinates(sp)[which(as.vector(sp@data) ==1 ), ])] <- NA
#             mask[cellFromXY(mask, coordinates(sp)[which(as.vector(sp@data) == 0), ])] <- NA


            # checking of nb candidates
            nb.cells <- .nb.available.pa.cells(mask.in, PA.flag = 0)
            if(nb.cells <= nb.points){
              nb.repet <- 1
              nb.points <- nb.cells
              cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
            }

            # select cells into raster
            pa.tab.tmp <- matrix(NA, ncol = nb.repet, nrow = nb.points)
            for( j in 1:ncol(pa.tab.tmp)){
              ## initialise the vector of sample cells
              SR <- NULL
              ## define a compy of the sampling mask
              mask.in.tmp <- mask.in
              ## force to get at least one value of each factorial variable
              fact.level.cells <- sample.factor.levels(env,
                                                       mask.out  = mask.out,
                                                       mask.in = mask.in)
              if(length(fact.level.cells)){
                SR <- c(SR, fact.level.cells)
                ## update the mask by removing already selected cells
                mask.in.tmp[SR] <- NA
              }
              SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                       size = nb.points - length(SR),
                                       cells = TRUE,
                                       na.rm = TRUE)[, "cell", drop = TRUE])
              ## repeat sampling until haing the right number of points
              while(length(SR) < nb.points){
                ## update the mask by removing already selected cells
                mask.env.tmp[SR] <- NA
                ## extract the missing number of points
                SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                         size = nb.points - length(SR),
                                         cells = T,
                                         na.rm = T)[, "cell", drop = TRUE])
              }
              pa.tab.tmp[, j] <- SR
            }

            # puting cells in good format
            selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
            pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
            colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
            for( j in 1:ncol(pa.tab)){
              pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
            }

            # puting presences, true absences and pseudo absences together
            xy <- rbind(coordinates(sp)[which(!is.na(as.vector(sp@data))),],
                        xyFromCell(mask.in, selected.cells))
            xy <- .add_PA_rownames(xy)
            sp <- as.numeric(unlist(c(na.omit(as.vector(sp@data)), rep(NA,length(selected.cells))), use.names=FALSE))
            env <- extract(env, xy)

            pa.tab <- rbind(matrix(TRUE,nrow=(nrow(xy)-length(selected.cells)), ncol=ncol(pa.tab)),
                           pa.tab)

            return(list(xy = xy,
                        sp = sp,
                        env = as.data.frame(env),
                        pa.tab = as.data.frame(pa.tab)))

          })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

setGeneric( "disk.pseudo.abs.selection",
              def = function(sp,env, ...){
                      standardGeneric( "disk.pseudo.abs.selection" )
                      } )

setMethod('disk.pseudo.abs.selection', signature(env="SpatialPointsDataFrame"),
          function(sp, env, distMin, distMax, nb.points, nb.repet){
            cat("\n   > Disk pseudo absences selection")

            # 1. determining selectable area
            coor <- coordinates(sp)
            pres <- which(sp@data[,1]==1)
            true.abs <- which(sp@data[,1]==0)
            tmp.abs <- which(is.na(sp@data[,1])) #(1:ncol(sp@data))[-c(pres,true.abs)]
            outside <- rep(0, length(abs))
            inside <- rep(0, length(abs))


            for(i in 1:length(pres)){
              # removing points too close from presences
              inside <- inside + ( sqrt((coor[tmp.abs,1]-coor[pres[i],1])^2 + (coor[tmp.abs,2]-coor[pres[i],2])^2) > distMin )
              # keeping points not to far from presences
              if(!is.null(distMax)){
                outside <- outside + ( sqrt((coor[tmp.abs,1]-coor[pres[i],1])^2 + (coor[tmp.abs,2]-coor[pres[i],2])^2) < distMax )
              }
            }
            if(is.null(distMax)){ # no cells are too far
              outside <- outside + 1
            }
            selected.abs <- tmp.abs[ (inside == length(pres)) & (outside > 0) ]

            # 2. adding presences and true absences and selecting randomly pseudo absences

            return(random.pseudo.abs.selection( sp[c(pres, true.abs, selected.abs),],
                                                env[c(pres, true.abs, selected.abs),],
                                                nb.points, nb.repet ))


          })


## TODO(damien): remimplement disk.pseudo.abs.selection to call random.pseudo.abs.selection
setMethod('disk.pseudo.abs.selection', signature(env="RasterStack"),
          function(sp, env, distMin, distMax, nb.points, nb.repet){
            cat("\n   > Disk pseudo absences selection")

            # 1. Check if NA are present in sp or not to determine which dataset to use
            if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
              env.tmp <- SpatialPointsDataFrame(coords = coordinates(sp),
#                                                 data = as.data.frame(extract(env,coordinates(sp),method='bilinear'))
                                                data = as.data.frame(extract(env,coordinates(sp))))

              return(disk.pseudo.abs.selection(sp, env.tmp, distMin, distMax, nb.points, nb.repet))
            } else {
              cat("\n   > Pseudo absences are selected in explanatory variables")

              # create a mask
              dist.mask <- raster::subset(env,1, drop=TRUE)
              dist.mask[] <- NA

              pres.xy <- coordinates(sp[which(sp@data[,1]==1),])
              dist.mask[cellFromXY(dist.mask,pres.xy)] <- 1

              dist.mask <- raster::distance(dist.mask)
              dist.mask <- mask(dist.mask, raster::subset(env,1, drop=TRUE))

              if(is.null(distMax)) distMax <- Inf

              mask.in <- reclassify(dist.mask, c(-Inf,distMin,NA ,distMin, distMax,-1, distMax,Inf,NA))

              # get the mask of already sampled mask
              mask.out <- mask.in; mask.out[] <- NA; mask.out[cellFromXY(mask.out, coordinates(sp))] <- 1

              # checking of nb candidates
              nb.cells <- .nb.available.pa.cells(mask.in, PA.flag = -1)
              if(nb.cells <= nb.points){
                nb.repet <- 1
                nb.points <- nb.cells
                cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
              }

              # select cells into raster
              pa.tab.tmp <- matrix(NA, ncol = nb.repet, nrow = nb.points)
              for( j in 1:ncol(pa.tab.tmp)){
                ## initialise the vector of sample cells
                SR <- NULL
                ## define a compy of the sampling mask
                mask.in.tmp <- mask.in
                ## force to get at least one value of each factorial variable
                fact.level.cells <- sample.factor.levels(env, mask.out = mask.out)
                if(length(fact.level.cells)){
                  SR <- c(SR, fact.level.cells)
                  ## update the mask by removing already selected cells
                  mask.in.tmp[SR] <- NA
                }
                SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                         size = nb.points - length(SR),
                                         cells = T,
                                         na.rm = T)[, "cell", drop = TRUE])
                ## repeat sampling until haing the right number of points
                while(length(SR)<nb.points){
                  ## update the mask by removing already selected cells
                  mask.in.tmp[SR] <- NA
                  ## extract the missing number of points
                  SR <- c(SR, sampleRandom(x = mask.in.tmp,
                                           size = nb.points - length(SR),
                                           cells = T,
                                           na.rm = T)[, "cell", drop = TRUE])
                }
                pa.tab.tmp[,j] <- SR
              }

              # puting cells in good format
              selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
              pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
              colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
              for( j in 1:ncol(pa.tab)){
                pa.tab[is.element(selected.cells, pa.tab.tmp[,j]), j] <- TRUE
              }

              # puting presences, true absences and pseudo absences together
              xy <- rbind(coordinates(sp), xyFromCell(mask.in,
                                                      selected.cells))
              xy <- .add_PA_rownames(xy)
              sp <- as.numeric(unlist(c(as.vector(sp@data), rep(NA,length(selected.cells))),
                                      use.names = FALSE))
              env <- extract(env, xy)

              pa.tab <- rbind(matrix(TRUE,
                                     nrow = (nrow(xy) - length(selected.cells)),
                                     ncol = ncol(pa.tab)), pa.tab)

              return(list(xy = xy,
                          sp = sp,
                          env = as.data.frame(env),
                          pa.tab = as.data.frame(pa.tab)))

            }
          })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# automaticaly add rownames to a data.frame
.add_PA_rownames <- function(xy){
  rn <- row.names(xy)
  missing_rn <- which(rn == "")
  if(length(missing_rn)){
    rn[missing_rn] <- paste("pa", 1:length(missing_rn), sep="")
  }
  rownames(xy) <- rn
  return(xy)
}
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# setMethod('disk.pseudo.abs.selection', signature(env="RasterStack"),
#           function(sp, env, distMin, distMax, nb.points, nb.repet){
#             cat("\n   > Disk pseudo absences selection")
#
#             # 1. Check if NA are present in sp or not to determine which dataset to use
#             if(.nb.available.pa.cells(sp) > 0 ){ # PA will be taken into response variable
#               env.tmp <- SpatialPointsDataFrame(coords = coordinates(sp),
#                                                 data = as.data.frame(extract(env,coordinates(sp),method='bilinear')))
#
#               return(disk.pseudo.abs.selection(sp, env.tmp, distMin, distMax, nb.points, nb.repet))
#             } else {
#               cat("\n   > Pseudo absences are selected in explanatory variables")
#
#               # create a mask
#               mask <- maskInside <- maskOutside <- reclassify(raster::subset(env,1), c(-Inf,Inf,0))
#               pres.xy <- coordinates(sp[which(sp@data[,1]==1),])
#
#               # to convert longitudinal degrees into metters
#               coef.conversion <- ifelse(grepl("longlat",env@crs@projargs), 111319.5, 1)
#               #               coef.conversion <- 1
#               ## progress bar
#               cat("\n")
#               pb <- txtProgressBar(min = 0, max = nrow(pres.xy), initial = 0, char = "=-",width = 20,  style = 3, file = "")
#               for(i in 1:nrow(pres.xy)){
#                 setTxtProgressBar(pb,i)
#                 if(distMin > 0){
#                   maskInside <- maskInside + (distanceFromPoints(mask, pres.xy[i,]) > (distMin * coef.conversion))
#                 }
#                 if(!is.null(distMax)){
#                   maskOutside <- maskOutside + (distanceFromPoints(mask, pres.xy[i,]) <= (distMax * coef.conversion))
#                 }
#               }
#
#               if(distMin > 0){
#                 maskInside <- maskInside == nrow(pres.xy)
#               } else { # keep all cells
#                 maskInside <- maskInside + 1
#               }
#
#               if(!is.null(distMax)){
#                 maskOutside <- maskOutside > 0
#               } else{ # keep all cells
#                 maskOutside <- maskOutside + 1
#               }
#
#
#               mask <- maskInside * maskOutside
#               mask[mask==0] <- NA
#               mask <- (-1) * mask
#
#               # remove presences and true absences from our raster
#               mask[cellFromXY(mask,coordinates(sp))] <- NA
#
#               # checking of nb candidates
#               nb.cells <- .nb.available.pa.cells(mask)
#               if(nb.cells <= nb.points){
#                 nb.repet <- 1
#                 nb.points <- nb.cells
#                 cat("\n   > All available cells have been selected (", nb.points, "pseudo absences selected )")
#               }
#
#               # select cells into raster
#               pa.tab.tmp <- matrix(NA, ncol=nb.repet, nrow=nb.points)
#               for( j in 1:ncol(pa.tab.tmp)){
#                 pa.tab.tmp[,j] <- sampleRandom(x=mask, size=nb.points, cells=T)[,"cell"]
#               }
#
#               # puting cells in good format
#               selected.cells <- sort(unique(as.vector(pa.tab.tmp)))
#               pa.tab <- matrix(FALSE, ncol = nb.repet, nrow = length(selected.cells))
#               colnames(pa.tab) <- paste("PA", 1:nb.repet, sep="")
#               for( j in 1:ncol(pa.tab)){
#                 pa.tab[selected.cells %in% pa.tab.tmp[,j], j] <- TRUE
#               }
#
#               # puting presences, true absences and pseudo absences together
#               xy <- rbind(coordinates(sp), xyFromCell(mask, selected.cells))
#               sp <- as.numeric(unlist(c(as.vector(sp@data), rep(NA,length(selected.cells))), use.names=FALSE))
#               env <- extract(env, xy)
#
#               pa.tab <- rbind(matrix(TRUE,nrow=(nrow(xy)-length(selected.cells)), ncol=ncol(pa.tab)),
#                               pa.tab)
#
#               return(list(xy = xy,
#                           sp = sp,
#                           env = as.data.frame(env),
#                           pa.tab = as.data.frame(pa.tab)))
#
#             }
#           })

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# .arranging.pa.table(pa.data, pa.tab, sp.data=NULL, xy=NULL){
#
#   # transforming sp.data into vector if it's not
#   if(!is.null(sp.data)){ # that means that PA were chosed into explanatories data
#     if(inherits(sp.data, 'SpatialPoints')){
#       xy <- coordinates(sp.data)
#       sp.data <- sp.data@data
#     }
#     if(inherits(sp.data, 'Raster')){
#       xy <- rbind(xyFromCell(sp.data, Which(sp.data >= 1), cells=TRUE), xyFromCell(sp.data, Which(sp.data == 0)))
#       sp.data.tmp <- rep(0,nrow(xy))
#       sp.data.tmp[1:length(Which(sp.data >= 1, cells=TRUE))] <- 1
#       sp.data <- sp.data.tmp
#       rm('sp.data.tmp')
#     }
#   }
#
#   # getting PA selected
#
# }

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# # additional hidden functions
# .allAvailableAbs <- function(data.biomod.species){
#   out <- data.biomod.species
#   if( sum(is.na(out)>0) )
#     out[is.na(out)] <- 0
#   return(out)
# }

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# .is.some.na.in.data <- function(sp){
#   if(is.vector(sp)){
#     if(sum(is.na(sp)) == 0){
#       cat('\nAvailable absences will be get in explanatory variables')
#       return(FALSE)
#     } else { return(TRUE) }
#   }
#
#   if(inherits(sp, 'SpatialPoints')){
#     if(sum(is.na(sp[,1])) == 0){
#       cat('\nAvailable absences will be get in explanatory variables')
#       return(FALSE)
#     } else { return(TRUE) }
#   }
#
#   if(inherits(sp, 'Raster')){
#     if(sp@data@min >= 0){
#       cat('\nAvailable absences will be get in explanatory variables')
#       return(FALSE)
#     } else { return(TRUE) }
#   }
# }
