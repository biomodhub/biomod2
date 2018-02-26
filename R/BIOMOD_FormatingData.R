####################################################################################################
# BIOMOD_FormatingData
# Damien.G
# feb 2012
####################################################################################################

# AIM :
#   puting input data in the right format and doing Pseudo Absences selection if desire

# INPUT :
#   resp.var <- Response Variable (monospecific) as vector, sp.point.data.frame or rasterLayer
#               code for vector and sp.objects : 1=pres, 0=true_abs, NA=no_info
#               code for vector and sp.objects : 1=pres, 0=true_abs, -1=no_info
#   expl.var <- Explanatory Variable as matrix, data.frame, sp.point.data.frame or rasterStack
#   resp.xy <- coordiantes of reponse points (2 column matrix)
#   resp.name <- name of considered specie
#   eval.resp.var <- independent response variable for models evaluations
#   eval.expl.var <- independent explanatory variable for models evaluations
#   eval.resp.xy <- independent response variable coordinates variable for models evaluations
#   PA.nb.rep <- Nb of Pseudo Absences Run to compute
#   PA.nb.absences <- Nb of Absences selected (true absences are counted in)
#   PA.strategy <- Pseudo Absences strategy
#   PA.dist.min <- If strategy is 'disk' : Pseudo Absences minimum distance between pres and selected absences (in metters if explanatory is georeferenced or in resp.xy units in all other cases)
#   PA.dist.man <- If strategy is 'disk' : Pseudo Absences maximum distance between pres and selected absences (in metters if explanatory is georeferenced or in resp.xy units in all other cases)
# 
#   PA.sre.quant <- If strategy is 'sre' : the quantile use for sre calculation
#   PA.table <- If strategy is 'user.defined' : a boolean data.frame indiacating which points of resp.var should be sonsidered in each PA run.
#   na.rm <- if True na are automatically removed

# OUTPUT : 
#   a BIOMOD.formated.data object that will be given to BIOMOD_Modeling function

####################################################################################################

'BIOMOD_FormatingData' <- function(resp.var,
                                   expl.var,
                                   resp.xy = NULL,
                                   resp.name = NULL,
                                   eval.resp.var = NULL,
                                   eval.expl.var = NULL,
                                   eval.resp.xy = NULL,
                                   PA.nb.rep = 0,
                                   PA.nb.absences = 1000,
                                   PA.strategy = 'random',
                                   PA.dist.min = 0,
                                   PA.dist.max = NULL,
                                   PA.sre.quant = 0.025,
                                   PA.table = NULL,
                                   na.rm = TRUE){
  .bmCat(paste(resp.name, " Data Formating", sep=""))
  
  # 1 check args
  args <- .BIOMOD_FormatingData.check.args(resp.var,
                                           expl.var,
                                           resp.xy,
                                           resp.name,
                                           eval.resp.var,
                                           eval.expl.var,
                                           eval.resp.xy,
                                           PA.nb.rep,
                                           PA.nb.absences,
                                           PA.strategy,
                                           PA.dist.min,
                                           PA.dist.max,
                                           PA.sre.quant,
                                           PA.table)
  
  resp.var <- args$resp.var
  expl.var <- args$expl.var 
  resp.xy <- args$resp.xy
  resp.name <- args$resp.name
  eval.resp.var <- args$eval.resp.var
  eval.expl.var <- args$eval.expl.var 
  eval.resp.xy <- args$eval.resp.xy
  PA.nb.rep <- args$PA.nb.rep
  PA.nb.absences <- args$PA.nb.absences
  PA.strategy <- args$PA.strategy
  PA.dist.min <- args$PA.dist.min
  PA.dist.max <- args$PA.dist.max
  PA.sre.quant <- args$PA.sre.quant
  PA.table <- args$PA.table
  
  rm('args')
  gc()
  
  out <- NULL
  
  if(PA.strategy == 'none'){ # no Pseudo Absences
    out <- BIOMOD.formated.data(sp=resp.var,
                                xy=resp.xy,
                                env=expl.var,
                                sp.name=resp.name,
                                eval.sp=eval.resp.var,
                                eval.env=eval.expl.var,
                                eval.xy=eval.resp.xy,
                                na.rm=na.rm)
  } else{ # Automatic Pseudo Absences Selection
    out <- BIOMOD.formated.data.PA(sp=resp.var, xy=resp.xy, env=expl.var, sp.name=resp.name,
                                   eval.sp=eval.resp.var, eval.env=eval.expl.var, eval.xy=eval.resp.xy,
                                   PA.NbRep=PA.nb.rep, PA.strategy=PA.strategy, 
                                   PA.nb.absences = PA.nb.absences, PA.dist.min = PA.dist.min,
                                   PA.dist.max = PA.dist.max, PA.sre.quant = PA.sre.quant, PA.table=PA.table, 
                                   na.rm=na.rm)
  } 
  
  
  .bmCat("Done")
  return(out)
}

.BIOMOD_FormatingData.check.args <- function(resp.var,
                                             expl.var,
                                             resp.xy,
                                             resp.name,
                                             eval.resp.var,
                                             eval.expl.var,
                                             eval.resp.xy,
                                             PA.nb.rep,
                                             PA.nb.absences,
                                             PA.strategy,
                                             PA.dist.min,
                                             PA.dist.max,
                                             PA.sre.quant,
                                             PA.table){
  
  # 0. names checking

  
  
  ### check resp.name is available
  if(grepl('_',resp.name) | grepl(' ',resp.name)){
    resp.name <- paste(unlist(strsplit(resp.name,'_')),collapse='.')
    resp.name <- paste(unlist(strsplit(resp.name,' ')),collapse='.')
    
    cat('\n Response variable name was converted into', resp.name)
  }
  
  ### check resp.name is available
  ### Not done because no imporance
  
  # 1. Checking input params class
  available.types <- c( 'numeric', 'data.frame', 'matrix', 
                        'RasterLayer', 'RasterStack',
                        'SpatialPointsDataFrame', 'SpatialPoints')
  ###### resp.var
  if(!(class(resp.var) %in% available.types)){
    stop( paste("Response variable must be one of ", toString(available.types), sep="") )
  }
    
  ### response var raster object not supported yet
  if(inherits(resp.var, 'Raster')){
    stop("Raster response variable not supported yet ! \nPlease extract your Presences and your absences by yourself")
    #### TO DO ####
    ## extract the 0 and 1 in sp format
  }
  
  ###### expl.var
  if(!(class(expl.var) %in% available.types[-which(available.types == 'SpatialPoints')])){
    stop( paste("Explanatory variable must be one of ", toString(available.types), sep="") )
  }
    
  
  ###### resp.xy
  if(inherits(resp.var,'SpatialPoints') ){
    if(!is.null(resp.xy)){
      cat("\n      ! XY coordinates of response variable will be ignored because spatial response object is given.")
    }
    resp.xy <- data.matrix(sp::coordinates(resp.var))
    if(class(resp.var) == 'SpatialPointsDataFrame'){
      resp.var <- resp.var@data
    } else{
      cat("\n      ! Response variable is considered as a only presences one... Is it really what you want?")
      resp.var <- rep(1,nrow(resp.xy))
    }
    
  }
  
  
  ### transforming into numeric if data.frame or matrix
  if(is.matrix(resp.var) | is.data.frame(resp.var)){
    if(ncol(resp.var) > 1){
      stop("You must give a monospecific response variable (1D object)")
    } else{
      resp.var <- as.numeric(resp.var[,1])
    }
  }
  
  if(is.matrix(expl.var) | is.numeric(expl.var) ){
    expl.var <- as.data.frame(expl.var)
  }
      
  if(inherits(expl.var, 'Raster')){
    expl.var <- raster::stack(expl.var, RAT=FALSE)
  }
    
  if(inherits(expl.var, 'SpatialPoints')){
    expl.var <- as.data.frame(expl.var@data)
  }
      
  ### check of xy coordinates validity    
  if(!is.null(resp.xy)){
    if(ncol(resp.xy)!=2){
      stop("if given, resp.xy must be a 2 column matrix or data.frame")
    }
    if(nrow(resp.xy) != length(resp.var)){
      stop("Response variable and its coordinates don't match")
    }
    resp.xy <- as.data.frame(resp.xy)
  }
  
  ### convert response var into binary
  resp.var[which(resp.var>0)] <- 1
  resp.var[which(resp.var<=0)] <- 0
      
  #### At this point :
  ####  - resp.var is a numeric
  ####  - resp.xy is NULL or a data.frame
  ####  - expl.var is a data.frame or a RasterStack
  ####  - sp.name is a character
  
  ### check resp and expl var compatibility
  if(is.data.frame(expl.var)){
    if(nrow(expl.var) != length(resp.var)){
      stop("If explanatory variable is not a raster then dimentions of response variable and explanatory variable must match!")
    }
  }
  
  ### PA strategy
#   if(!is.null(PA.strategy)){ # force PA.nb.rep to be positive if PA.strategy is defined
#     PA.nb.rep = max(c(PA.nb.rep,1))
#   }
  
  if(is.null(PA.table) & PA.nb.rep < 1){
    cat("\n> No pseudo absences selection !")
    PA.strategy <- "none"
    PA.nb.rep <- 0
  }
  
  if(is.null(PA.strategy) &  PA.nb.rep > 0){
    cat("\n> Pseudo absences will be selected randomly !")
    PA.strategy <- "random"
  }
  
  
  if( !is.null(PA.table)){
    cat("\n> Pseudo absences used will be user defined ones !")
    PA.strategy <- "user.defined"
    PA.nb.rep <- 0
  }
  
  if(PA.strategy == "user.defined"){
    if(! (is.matrix(PA.table) | is.data.frame(PA.table)))
      stop("\n PA.table must be a matrix or a data.frame")
    
    if(nrow(PA.table) != length(resp.var))
      stop("\n PA.table must have as many row than the number 
           of observation of your response variable")
    
    #PA.table <- as.data.frame(sapply(PA.table,simplify=FALSE,as.logical))
    colnames(PA.table) <- paste("PA",1:ncol(PA.table),sep="")
    
  }
  
  # 2. eval.resp.var.checking
     
  if(!is.null(eval.resp.var)){
    # do the same test than previous one
    ###### eval.resp.var
    if(!(class(eval.resp.var) %in% available.types)){
      stop( paste("Response variable must be one of ", toString(available.types), sep="") )
    }
      
    ### response var raster object not supported yet
    if(inherits(eval.resp.var, 'Raster')){
      stop("Raster response variable not supported yet ! \nPlease extract your Presences and your absences by yourself")
      #### TO DO ####
      ## extract the 0 and 1 in sp format
    }
    
    ###### expl.var
    if(!is.null(eval.expl.var)){
      if(!(class(eval.expl.var) %in% available.types[-which(available.types == 'SpatialPoints')])){
        stop( paste("Explanatory variable must be one of ", toString(available.types), sep="") )
      }
    } else{
      if(!(inherits(expl.var, 'Raster'))){
        stop("If explanatory variable is not a raster and you want to consider evaluation response variable, you have to give evaluation explanatory variables")
      }
    }

    ###### resp.xy
    if(inherits(eval.resp.var,'SpatialPoints') ){
      if(!is.null(eval.resp.xy)){
        cat("\n      ! XY coordinates of response variable will be ignored because spatial response object is given.")
      }
      eval.resp.xy <- data.matrix(sp::coordinates(eval.resp.var))
      if(class(eval.resp.var) == 'SpatialPointsDataFrame'){
        eval.resp.var <- eval.resp.var@data
      } else{
        cat("\n      ! Response variable is considered as a only presences one... Is it really what you want?")
        eval.resp.var <- rep(1,nrow(eval.resp.xy))
      }
      
    }
    
    
    ### transforming into numeric if data.frame or matrix
    if(is.matrix(eval.resp.var) | is.data.frame(eval.resp.var)){
      if(ncol(eval.resp.var) > 1){
        stop("You must give a monospecific response variable (1D object)")
      } else{
        eval.resp.var <- as.numeric(eval.resp.var[,1])
      }
    }
    
    if(is.matrix(eval.expl.var) | is.numeric(eval.expl.var) ){
      eval.expl.var <- as.data.frame(eval.expl.var)
    }
        
    if(inherits(eval.expl.var, 'Raster')){
      eval.expl.var <- raster::stack(eval.expl.var)
    }
      
    if(inherits(eval.expl.var, 'SpatialPoints')){
      eval.expl.var <- as.data.frame(eval.expl.var@data)
    }
        
    ### check of xy coordinates validity    
    if(!is.null(eval.resp.xy)){
      if(ncol(eval.resp.xy)!=2){
        stop("if given, resp.xy must be a 2 column matrix or data.frame")
      }
      if(nrow(eval.resp.xy) != length(eval.resp.var)){
        stop("Response variable and its coordinates don't match")
      }
      eval.resp.xy <- as.data.frame(eval.resp.xy)
    }

    if(is.data.frame(eval.expl.var)){
      if(nrow(eval.expl.var) != length(eval.resp.var)){
        stop("If explanatory variable is not a raster then dimentions of response variable and explanatory variable must match!")
      }
    }

    ### remove NAs from evaluation data
    if( sum(is.na(eval.resp.var)) > 0 ){
      cat("\n      ! NAs have been automaticly removed from Evaluation data")
      if(!is.null(eval.resp.xy)){
        eval.resp.xy <- eval.resp.xy[-which(is.na(eval.resp.var)),]
      }
      eval.resp.var <- na.omit(eval.resp.var)
    }

    ### convert response var into binary
    eval.resp.var[which(eval.resp.var>0)] <- 1
    eval.resp.var[which(eval.resp.var<=0)] <- 0

    ### check there are both presences and absences in evaluation dataset
    if( sum(eval.resp.var == 1) < 1 | sum(eval.resp.var == 0) < 1){
      stop("Evaluation response data must have both presences and absences")
    }

  } else { 
    cat("\n      ! No data has been set aside for modeling evaluation")
    eval.expl.var <- eval.resp.xy <- NULL
  }
  
  ### PA arguments are not checked here because it will be done later... (may be will do it here later)
  
  return(list( resp.var = resp.var,
               expl.var = expl.var,
               resp.xy = resp.xy,
               resp.name = resp.name,
               eval.resp.var = eval.resp.var,
               eval.expl.var = eval.expl.var,
               eval.resp.xy = eval.resp.xy,
               PA.nb.rep = PA.nb.rep,
               PA.nb.absences = PA.nb.absences,
               PA.strategy = PA.strategy,
               PA.dist.min = PA.dist.min,
               PA.dist.max = PA.dist.max,
               PA.sre.quant = PA.sre.quant,
               PA.table = PA.table))
  
}