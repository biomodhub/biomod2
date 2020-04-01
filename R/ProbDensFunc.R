##' @name ProbDensFunc
##' @title Probability Density Function
##' @description
##' Using a variety of parameters in modelling will inevitably bring
##' variability in predictions, especially when it comes to making future
##' predictions.
##' This function enables an overall viewing of the future predictions
##' range per species and gives the likelihood of range shift estimations.
##' It will calculate the optimal way for condensing a defined proportion
##' (50, 75, 90 and 95\% per default) of the data.
##' 
##' 
##' @param initial a vector (resp. a \code{SpatialPointsDataFrame}) in a
##' binary format (ones and zeros) representing the current distribution
##' of a species which will be used as a reference for the range change
##' calculations
##' @param projections a matrix (resp; a \code{rasterStack}) grouping all
##' the predictions where each column is a single prediction. Make sure
##' you keep projections in the same order as the initial vector (
##' line1 = site1, line2 = site2, etc.).
##' @param plothist set to TRUE to plot the range change histogram
##' @param cvsn stands for "current vs new". If true, the range change
##' calculations will be of two types: the percentage of cells currently
##' occupied by the species to be lost, and the relative percentage of
##' cells currently unoccupied but projected to be, namely 'new' cells,
##' compared to current surface range.
##' @param groups an option for ungrouping the projections enabling a
##' separated visualisation of the prediction range per given group. A
##' matrix is expected where each column is a single prediction and each
##' line is giving details of one parameter (See the examples section).
##' @param resolution the step used for classes of prediction in graphics.
##' The default value is 5.
##' @param filename the name of file (with extension) where plots will be
##' stored. If not \code{NULL}, no plotting windows will be open
##' @param ... further arguments:
##' 
##' - \code{lim}: ordered numeric vector indicating the proportion of data
##' to consider for histogram representation (by default : c(0.5, 0.75,
##' 0.9, 0.95))
##' - \code{nb.points.max}: the maximum number of points to sample, 25000
##' by default (useful for huge \code{raster*} objects)
##' 
##' @details
##' The future range changes are calculated as a percentage of the
##' species' present state. For example, if a species currently occupies
##' 100 cells and is estimated by a model to cover 120 cells in the 
##' future, the range change will be + 20\%.
##' 
##' Resolution : Note that modifying the resolution will directly
##' influence the probability scale. Bigger classes will accumulate a
##' greater number of predictions and therefore represent a greater
##' fraction of the total predictions. The probability is in fact that of
##' the class and not of isolated events.
##' 
##' @return This is a plotting function, no objects are returned or
##' created.
##' 
##' @author Wilfried Thuiller, Bruno Lafourcade, Damien Georges 
##' 
##' @seealso \code{\link{BIOMOD_Projection}}, 
##' \code{\link{BIOMOD_EnsembleForecasting}}
##' 
##' @keywords optimize
##' @keywords distribution
##' 
##' @examples 
##' \dontrun{
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##'                                     package="biomod2"), row.names = 1)
##' head(DataSpecies)
##' 
##' ##' the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' ##' the presence/absences data for our species
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' ##' remove all 0 from response vector to work with
##' ##' presence only data (Pseudo Absences selections)
##' rm_id <- which(myResp==0)
##' myResp <- myResp[-rm_id]
##' 
##' 
##' ##' the XY coordinates of species data
##' myRespXY <- DataSpecies[-rm_id,c("X_WGS84","Y_WGS84")]
##' 
##' 
##' ##' Environmental variables extracted from BIOCLIM
##' myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd",
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio4.grd",
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio7.grd",
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio11.grd",
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio12.grd",
##'                              package="biomod2"))
##' 
##' ##' 1. Formatting Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName,
##'                                      PA.nb.rep=3)
##' 
##' ##' 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' ##' 3. Doing Modelisation
##' myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
##'                                      models = c('CTA','RF','GLM','GAM','ANN','MARS'),
##'                                      models.options = myBiomodOption,
##'                                      NbRunEval=5,
##'                                      DataSplit=70,
##'                                      Prevalence=0.5,
##'                                      models.eval.meth = c('TSS'),
##'                                      do.full.models = FALSE,
##'                                      rescal.all.models=T,
##'                                      modeling.id='test')
##' 
##' ##' 4. Build ensemble-models that will be taken as reference
##' myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
##'                                        chosen.models = 'all',
##'                                        em.by = 'all',
##'                                        eval.metric = c('TSS'),
##'                                        eval.metric.quality.threshold = c(0.7),
##'                                        prob.mean = TRUE,
##'                                        prob.median = TRUE)
##' 
##' ##' 5. Projection on future environmental conditions
##' 
##' ###' load future environmental conditions from biomod2 package
##' myExpl_fut <- raster::stack( system.file( "external/bioclim/future/bio3.grd",
##'                                   package="biomod2"),
##'                      system.file( "external/bioclim/future/bio4.grd",
##'                                   package="biomod2"),
##'                      system.file( "external/bioclim/future/bio7.grd",
##'                                   package="biomod2"),
##'                      system.file( "external/bioclim/future/bio11.grd",
##'                                   package="biomod2"),
##'                      system.file( "external/bioclim/future/bio12.grd",
##'                                   package="biomod2"))
##' 
##' myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
##'                                         new.env = myExpl_fut,
##'                                         proj.name = 'future',
##'                                         selected.models = 'all',
##'                                         binary.meth = 'TSS',
##'                                         compress = FALSE,
##'                                         build.clamping.mask = TRUE)
##' 
##' BIOMOD_EnsembleForecasting(projection.output=myBiomodProjection,
##'                            EM.output=myBiomodEM,
##'                            binary.meth='TSS')
##' 
##' ##' 6. load binary projections
##' consensusBin <- raster::stack('GuloGulo/proj_future/proj_future_GuloGulo_ensemble_TSSbin.grd')
##' projectionsBin <- raster::stack('GuloGulo/proj_future/proj_future_GuloGulo_TSSbin.grd')
##' 
##' ##' 7. build a ref state based on ensemble-models
##' ref <- sampleRandom(subset(consensusBin, 1, drop=T), size=5000, sp=T, na.rm=T)
##' 
##' ##' 8. autoatic creation of groups matrix
##' find_groups <- function(diff_by_pix){
##'   data.set <- sapply(names(diff_by_pix),biomod2:::.extractModelNamesInfo,info='data.set')
##'   run.eval <- sapply(names(diff_by_pix),biomod2:::.extractModelNamesInfo,info='run.eval')
##'   models <- sapply(names(diff_by_pix),biomod2:::.extractModelNamesInfo,info='models')
##'   return(rbind(data.set,run.eval,models))
##' }
##' 
##' groups <- find_groups(projectionsBin)
##' 
##' ##' 9. plot ProbDensFunct graphs
##' ProbDensFunc(initial = ref,
##'              projections = projectionsBin,
##'              plothist=TRUE,
##'              cvsn=TRUE,
##'              groups=groups,
##'              resolution=2,
##'              filename=NULL,
##'              lim=c(0.5,0.8,0.95))
##' 
##' ###' 3 plots should be produced.. Should be convenient to save it within a device
##' ###' supporting multiple plots.
##' 
##' }
##' 
ProbDensFunc <- function(
  initial,
  projections,
  groups = NULL,
  plothist = TRUE,
  cvsn = TRUE,
  resolution = 5,
  filename = NULL,
  ...
){

  args <- .ProbDensFunc.checkArgs(
    initial,
    projections,
    groups,
    plothist,
    cvsn,
    resolution,
    filename,
    ...
  )

  initial = args$initial
  projections = args$projections
  groups = args$groups
  plothist = args$plothist
  cvsn = args$cvsn
  resolution = args$resolution
  filename = args$filename
  lim = args$lim
  one_plot = args$one_plot

  rm(list='args')

  # results will be stored in out
  out <- list()

  # oppen the divice if filename
  if(length(filename)){
    switch(tools::file_ext(filename),
           pdf = pdf(filename),
           jpeg = jpeg(filename),
           tiff = tiff(filename),
           eps = postscript(filename),
           png = png(filename))
  }

  # area stores the species range change calculations
  area <- (apply(projections,2,sum, na.rm=T) / sum(initial, na.rm=T) - 1 ) * 100
  a <- round( (min(area, na.rm=TRUE)-(resolution+10))/10 ) *10
  b <- round( (max(area, na.rm=TRUE)+(resolution+10))/10 ) *10
  area_hist <- hist(area, breaks = seq(a,b,resolution), plot=FALSE)
  area_hist$density <- area_hist$counts / sum(area_hist$counts)


  # analysis of the distribution density and calculation of the probability of events
  area_sorted <- sort(area)
  nb <- round(length(area_sorted) * lim)
  lower_limit <- upper_limit <- c()

  for(i in 1:length(lim)){
     g <- rep(NA,length(area_sorted)-nb[i])
     for(j in 1:length(g)) g[j] <- diff(range(area_sorted[j:min(j+nb[i],length(area_sorted))]))
     lower_limit <- c(lower_limit,area_sorted[which.min(g)])
     upper_limit <- c(upper_limit,area_sorted[which.min(g)]+g[which.min(g)])
  }

  names(lower_limit) <- names(upper_limit) <- paste(lim*100, "%", sep="")

  out$stats <- cbind(lower_limit,upper_limit)

  if(!is.null(groups)){
    if(! (length(filename) | one_plot)) dev.new()
    par(mfrow=c(1,nrow(groups)))

    color.samp <- list()
    for(pa in 1:nrow(groups)){

      lv <- levels(as.factor(as.matrix(groups[pa,])))
      color.samp[[pa]] <- colors()[sample(c(90,417,552,616,382,11,150,468,28,31,420,476,333),length(lv))]
      g <- hist(area, breaks = seq(a,b,resolution), plot=FALSE)
      fac <- (max(g$counts) / sum(g$counts)) / max(g$density)
      g$density <- g$density * fac

      plot(g, freq=FALSE, border='grey88', main=row.names(groups)[pa], xlab="Species range change (%)", ylab="Event   occurence   probability")
      lines(density(area, width=30)$x, density(area, width=30)$y*fac)
      for(i in 1:length(lv)){
        div <- length(area) / length(area[groups[pa,]==lv[i]])
        lines(density(area[groups[pa,]==lv[i]],width=30)$x, density(area[groups[pa,]==lv[i]], width=30)$y /div*fac, col=color.samp[[pa]][i])
      }
      lv <- as.factor(as.matrix(groups[pa,]))
      leg <- list()
      for(j in 1:length(levels(lv))) leg[[j]] <- levels(lv)[j]
      legend("topright", legend=leg, bty='n',fill=color.samp[[pa]])
    }
  }

  if(cvsn){
    #calculation of the 2 axes independently (lost vs new sites)
    area2 <- (apply(projections[which(initial==1),],2,sum) / sum(initial==1) -1) * 100
    area3 <- area - area2

    if(! (length(filename) | one_plot)) dev.new()
    par(mfrow=c(1,nrow(groups)))

    for(i in 1:nrow(groups)){
      lv <- as.factor(as.matrix(groups[i,]))
      leg <- list()
      for(j in 1:length(levels(lv))) leg[[j]] <- levels(lv)[j]

      levels(lv) <- 1:length(levels(lv))
      plot(area3~area2, xlab='current', ylab='new', ylim=c(0,if(max(area3)<100){100}else{max(area3)+30}), xlim=c(-100,0), col=color.samp[[i]][lv], pch=20, main=row.names(groups)[i])
      legend("bottomleft", legend=leg, bty='n',fill=color.samp[[i]])
      abline(0,-1, col='grey80')
      abline(100,-1, col='grey80')
      text(x=-97,y=103,pos=1,label="SRC = 0", col="black", cex=0.8)
      if(max(area3)<100)text(x=-3,y=103,pos=1,label="SRC = 100", col="black", cex=0.8)
      else text(x=if(max(area3)<200){-(max(area3)+33-100)}else{-96},y=if(max(area3)<200){max(area3)+33}else{203},pos=1,label="SRC = 100", col="black", cex=0.8)
    }
  }

  #if plot of distribution plot wanted
  if(plothist){
    if(! (length(filename) | one_plot)) dev.new()
    par(mfrow=c(1,1))

    col_list <- colorRampPalette(c("dodgerblue1","steelblue1","slategray1","aliceblue"))(length(lim))

    for( l in length(lim):1){
      hist( mean(out$stat[l,], na.rm=T), breaks=out$stat[l,], col=col_list[l], xlim=c(a,b), ylim=c(0,max(area_hist$density)*1.2), xlab="",ylab="", main="", add=ifelse(l==length(lim),FALSE,TRUE))
      abline(v=out$stat[l,],col= col_list[l], lwd=1.7)
    }

    legend("topright",
          legend=rownames(out$stat), bty='n', fill=col_list, cex=0.8, title='distrib. of data' )

    par(new=TRUE)
    plot(area_hist, freq=FALSE, col="white", xlim=c(a,b), ylim=c(0,max(area_hist$density)*1.2), main="Probability density function", xlab="Species range change (%)", ylab="Event occurence probability")
  }

  if(length(filename)) dev.off()

  return(out)
}

.ProbDensFunc.checkArgs <- function(initial,
                                    projections,
                                    groups,
                                    plothist,
                                    cvsn,
                                    resolution,
                                    filename,
                                    ...){
  add.args <- list(...)
  if(is.null(add.args$nb.points.max) | !is.numeric(add.args$nb.points.max)) add.args$nb.points.max <- 25000
  if(is.null(add.args$lim) | !is.numeric(add.args$lim)) add.args$lim <- c(0.5, 0.75, 0.90, 0.95)
  if(is.null(add.args$one_plot)) add.args$one_plot <- FALSE

  # check lim arg
  if(sum(add.args$lim>1 | add.args$lim<0) > 0 ) stop("'lim' must be a numeric vector with 0 to 1 values")


  # check args types
  if(inherits(projections, 'Raster')){
    if(!inherits(initial, c('RasterLayer', 'SpatialPointsDataFrame')))
      stop("If projections is a raster object, initial should be a 'RasterLayer' or a 'SaptialPointDataFrame'")
  } else if(is.matrix(projections)){
    if(!is.numeric(initial)){
      stop("If projections is a matrix, initial should be a 'numeric'")
    }
  } else{
    stop("projections should be a 'matrix' or a 'RasterStack'")
  }

  # extract values
  if(inherits(projections, 'Raster')){
    if(inherits(initial, 'SpatialPointsDataFrame')){
      if(nrow(initial) > add.args$nb.points.max){
        initial[sort(sample(1:nrow(initial),size=add.args$nb.points.max)), drop=FALSE]
      }
    } else {
      initial <- sampleRandom(initial, size=min(add.args$nb.points.max, ncell(initial)) ,sp=TRUE, na.rm=TRUE)
    }

    projections <- extract(projections, initial, method='simple', na.rm=FALSE)
    initial <- initial@data[,1]
  } else{
    if(length(initial) != nrow(projections))
      stop("initial & projections dimensions don't match")
    if(length(initial) > add.args$nb.points.max){
      kept_rows <- sort(sample(1:length(initial), size=add.args$nb.points.max, replace=FALSE))
      initial <- initial[kept_rows]
      projections <- projections[kept_rows,,drop=FALSE]
    }
  }

  # remove NAs
  na_rows <- unique(c( which(is.na(initial)), which(is.na(projections), arr.ind=TRUE)[,1] ) )
  if(length(na_rows)){
    initial <- initial[-na_rows]
    projections <- projections[-na_rows,,drop=FALSE]
  }

  # check groups arg
  if(is.null(groups)){
    if(cvsn){
      cat("\n\t! 'cvsn' was automatically switch off because no 'groups' given")
       cvsn <- FALSE
    }
  } else {
    if(!is.matrix(groups)) stop("'groups' should be a matrix")
    if(ncol(groups)!=ncol(projections)) stop("'groups' and 'projections' do not have the same number of columns (resp. layers)")
  }

  # check saving options
  if(!is.null(filename)){
    if( ! (tools::file_ext(filename) %in% c("pdf","jpeg","tiff","eps","png"))){
      filename <- paste(tools::file_path_sans_ext(filename),".pdf",sep="")
      cat("\n\t! 'filename' extension unknown => outputs will be stored in :", filename)
    }
  }

  return(list( initial = initial,
               projections = projections,
               groups = groups,
               plothist = plothist,
               cvsn = cvsn,
               resolution = resolution,
               filename = filename,
               lim = add.args$lim,
               one_plot = add.args$one_plot))
}
