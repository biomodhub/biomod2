setGeneric( "BIOMOD_RangeSize", 
            def = function(CurrentPred, FutureProj, ...){
              standardGeneric( "BIOMOD_RangeSize" )
            } )

# The data.frame input function......................................................................
setMethod('BIOMOD_RangeSize', signature(CurrentPred='data.frame', FutureProj='data.frame' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){
            # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
            args <- .BIOMOD_RangeSize.check.args( CurrentPred, FutureProj,  SpChange.Save )
            
            
            #Function that estimate the number of pixel gain, loss, stable (present and absent) by species
            Value <- c(-2, 0, -1, 1)
            CompteurSp <- function(Data, Value)
            {
              if(is.data.frame(Data)) {
                N <- dim(Data)[2]
                Compt <- as.data.frame(matrix(0, ncol=4, nrow=dim(Data)[2]))
                i <- 1
                while(i <= N) {
                  Compt[i, 1] <- length(Data[Data[,i] == Value[1], i])
                  Compt[i, 2] <- length(Data[Data[,i] == Value[2], i])
                  Compt[i, 3] <- length(Data[Data[,i] == Value[3], i])
                  Compt[i, 4] <- length(Data[Data[,i] == Value[4], i])
                  i <- i + 1
                }
              }
              return(Compt)
            }
            Diff.By.Pixel <- as.data.frame(FutureProj - 2 * CurrentPred)
            Compt.By.Models <- as.data.frame(CompteurSp(Diff.By.Pixel, Value))
            Compt.By.Models[, 5] <- (100 * Compt.By.Models[, 1])/(Compt.By.Models[, 1] + Compt.By.Models[,3])
            Compt.By.Models[, 6] <- (100 * Compt.By.Models[, 4])/(Compt.By.Models[, 1] + Compt.By.Models[,3])
            Compt.By.Models[, 7] <- Compt.By.Models[, 6] - Compt.By.Models[, 5]
            Compt.By.Models[, 8] <- Compt.By.Models[, 1] + Compt.By.Models[,3]
            Compt.By.Models[, 9] <- Compt.By.Models[,3]
            Compt.By.Models[, 10] <- Compt.By.Models[, 4] + Compt.By.Models[,3]
            dimnames(Compt.By.Models) <- list(colnames(CurrentPred), c("Loss","Stable0", "Stable1", "Gain", "PercLoss", "PercGain", 
                                                                       "SpeciesRangeChange", "CurrentRangeSize", 
                                                                       "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp"))
            Output <- c()
            Output <- list(Compt.By.Models=Compt.By.Models, Diff.By.Pixel=Diff.By.Pixel)
            invisible(Output)
          })

# The array input function......................................................................
setMethod('BIOMOD_RangeSize', signature(CurrentPred='array', FutureProj='array' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){
            # transform arrays into data.frame
            CurrentPred <- as.data.frame(CurrentPred)
            names(CurrentPred) <- unlist(lapply(strsplit(names(CurrentPred),".", fixed=TRUE), 
                                                function(x){
                                                  return(paste( x[3], x[2], x[1],sep="_"))
                                                }))
            
            names(CurrentPred) <- unlist(lapply(strsplit(names(CurrentPred),".", fixed=TRUE), 
                                             function(x){
                                               x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                               data.set.id <- x.rev[1]
                                               cross.valid.id <- x.rev[2]
                                               algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                               model.id <- paste(data.set.id,
                                                                 cross.valid.id,
                                                                 algo.id, sep="_")
                                               return(model.id)
                                             }))
            
            FutureProj <- as.data.frame(FutureProj)
            names(FutureProj) <- unlist(lapply(strsplit(names(FutureProj),".", fixed=TRUE), 
                                                function(x){
                                                  x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                                  data.set.id <- x.rev[1]
                                                  cross.valid.id <- x.rev[2]
                                                  algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                                  model.id <- paste(data.set.id,
                                                                    cross.valid.id,
                                                                    algo.id, sep="_")
                                                  return(model.id)
                                                }))
            return(BIOMOD_RangeSize(CurrentPred, FutureProj,  SpChange.Save))
          })

# The Raster input function......................................................................

setMethod('BIOMOD_RangeSize', signature(CurrentPred='RasterStack', FutureProj='RasterStack' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){
            # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
            args <- .BIOMOD_RangeSize.check.args( CurrentPred, FutureProj,  SpChange.Save )
            
            CBS <- matrix(ncol=10, nrow=length(CurrentPred@layers), dimnames=list(names(CurrentPred), 
                                                                                  c("Loss","Stable0", "Stable1", "Gain", "PercLoss", "PercGain", "SpeciesRangeChange", "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")))
            
            
            sp.stack <- stack() 
            for(i in 1:length(CurrentPred@layers)){
              #DiffByPixel
              Cur <- CurrentPred@layers[[i]]
              Fut <- FutureProj@layers[[i]]
              Ras <- Fut - (Cur + Cur)
              sp.stack <- addLayer(sp.stack, Ras)
              
              #ComptBySpecies
              CBS[i, 1] <- length(which(Ras[]==-2))
              CBS[i, 2] <- length(which(Ras[]== 0))
              CBS[i, 3] <- length(which(Ras[]==-1))
              CBS[i, 4] <- length(which(Ras[]== 1))
              
              CBS[i, 5] <- round(CBS[i,1] / (CBS[i,1]+CBS[i,3]) *100, digits=3)                       
              CBS[i, 6] <- round(CBS[i,4] / (CBS[i,1]+CBS[i,3]) *100, digits=3)                       
              CBS[i, 7] <- round((CBS[i,3]+CBS[i,4]) / (CBS[i,1]+CBS[i,3]) *100 -100, digits=3)       
              
              CBS[i, 8] <- CBS[i,1]+CBS[i,3]                                        
              CBS[i, 9] <- CBS[i,3]                                                 
              CBS[i, 10] <- CBS[i,3]+CBS[i,4]                                       
            }
            
            if(is.null(SpChange.Save)) SpChange.Save <- "NoName"
            #     assign(paste(SpChange.Save, "_Compt.By.Species", sep=""), CBS, pos=1)
            names(sp.stack) <- rownames(CBS)
            #     assign(paste(SpChange.Save, "_Diff.By.Pixel", sep=""), sp.stack, pos=1)
            return(list(Compt.By.Models = CBS,
                        Diff.By.Pixel = sp.stack))
          })

####
setMethod('BIOMOD_RangeSize', signature(CurrentPred='RasterLayer', FutureProj='RasterStack' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){
            
            # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
            CBS <- matrix(ncol=10, nrow=length(FutureProj@layers), dimnames=list(names(FutureProj), 
                                                                                 c("Loss","Stable0", "Stable1", "Gain", "PercLoss", "PercGain", "SpeciesRangeChange", "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")))
            sp.stack <- stack() 
            for(i in 1:length(FutureProj@layers)){
              #DiffByPixel
              Cur <- CurrentPred
              Fut <- FutureProj@layers[[i]]
              Ras <- Fut - (Cur + Cur)
              sp.stack <- addLayer(sp.stack, Ras)
              
              #ComptBySpecies
              CBS[i, 1] <- length(which(Ras[]==-2))
              CBS[i, 2] <- length(which(Ras[]== 0))
              CBS[i, 3] <- length(which(Ras[]==-1))
              CBS[i, 4] <- length(which(Ras[]== 1))
              
              CBS[i, 5] <- round(CBS[i,1] / (CBS[i,1]+CBS[i,3]) *100, digits=3)                       
              CBS[i, 6] <- round(CBS[i,4] / (CBS[i,1]+CBS[i,3]) *100, digits=3)                       
              CBS[i, 7] <- round((CBS[i,3]+CBS[i,4]) / (CBS[i,1]+CBS[i,3]) *100 -100, digits=3)       
              
              CBS[i, 8] <- CBS[i,1]+CBS[i,3]                                        
              CBS[i, 9] <- CBS[i,3]                                                 
              CBS[i, 10] <- CBS[i,3]+CBS[i,4]                                       
            }
            
            if(is.null(SpChange.Save)) SpChange.Save <- "NoName"
            #     assign(paste(SpChange.Save, "_Compt.By.Species", sep=""), CBS, pos=1)
            names(sp.stack) <- rownames(CBS)
            #     assign(paste(SpChange.Save, "_Diff.By.Pixel", sep=""), sp.stack, pos=1)
            return(list(Compt.By.Models = CBS,
                        Diff.By.Pixel = sp.stack))
          })


setMethod('BIOMOD_RangeSize', signature(CurrentPred='RasterLayer', FutureProj='RasterLayer' ),
          function(CurrentPred, FutureProj,  SpChange.Save=NULL){
            BIOMOD_RangeSize(CurrentPred = CurrentPred, FutureProj = stack(FutureProj),  SpChange.Save = SpChange.Save)
          })



####################################################################
.BIOMOD_RangeSize.check.args <- function( CurrentPred, FutureProj,  SpChange.Save ){
  # dimensions checking
  if(sum(!(dim(CurrentPred) == dim(FutureProj)) > 0)){
    stop("CurrentPred & FutureProj dimensions mismatched !")
  }
  
  # binary checking
  
  
}
