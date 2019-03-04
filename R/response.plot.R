`response.plot` <-
  function(model, Data, show.variables=seq(1:ncol(Data)), save.file="no", name="response_curve", ImageSize=480, plot=TRUE){

    cat("\n! Deprecated function, please use response.plot2 instead!")
    return(TRUE)
  }


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
`response.plot2` <- function(models,
                             Data,
                             show.variables=seq(1:ncol(Data)),
                             do.bivariate = FALSE,
                             fixed.var.metric = 'mean',
                             save.file="no",
                             name="response_curve",
                             ImageSize=480,
                             plot=TRUE,
                             ...){

  # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  add.args <- list(...)
  on_0_1000 <- add.args$on_0_1000
  col <- add.args$col
  lty <- add.args$lty
  legend <- add.args$legend
  main <- add.args$main
  display_title <- add.args$display_title
  # par <- add.args$par
  restrict_to_pres_range <- add.args$restrict_to_pres_range
  data_species <- add.args$data_species
  use.formal.names <- add.args$use.formal.names

  if(is.null(on_0_1000)) on_0_1000 <- FALSE
  if(is.null(col)) if(!do.bivariate) col <- "black" else col <- c("red","orange","green")
  if(is.null(lty)) lty <- 1
  if(is.null(legend)) legend <- FALSE
  if(is.null(display_title)) display_title <- TRUE
  if(is.null(restrict_to_pres_range)) restrict_to_pres_range <- FALSE
  if(is.null(use.formal.names)) use.formal.names <- FALSE

  formal_names <- models

  args <- .response.plot2.check.arg(models, Data, show.variables, save.file, name, ImageSize, plot, fixed.var.metric, do.bivariate, add.args)

  models <- args$models
  Data <- args$Data
  show.variables <- args$show.variables
  save.file <- args$save.file
  name <- args$name
  ImageSize <- args$ImageSize
  plot <- args$plot
  fixed.var.metric <- args$fixed.var.metric
  do.bivariate <- args$do.bivariate
  nb.pts <- args$nb.pts

  if(is.null(main)) main <- try(paste("Response curves for ", get(models[1])@resp_name, "'s ",  get(models[1])@model_class,sep=""))
  if(is.null(data_species)) data_species <- rep(1,nrow(Data)) else data_species[data_species!=1 | is.na(data_species)] <- 0


  # 2. build function outputs
  factor_id <- which(sapply(Data,is.factor))

  list.out <- list()

  # Create a ranged data table
  ref_table <- Data[1,,drop=F]
  rownames(ref_table) <- NULL

  for(i in 1:ncol(Data)){
    if(is.numeric(Data[,i])){
      ref_table[,i] <- switch(fixed.var.metric,
                              mean = mean(Data[data_species==1,i]),
                              median = median(Data[data_species==1,i]),
                              min = min(Data[data_species==1,i]),
                              max = max(Data[data_species==1,i]))
    } else{
      # return everytimes the majoritary class
      sum_level <- summary(Data[data_species==1,i], na.rm = TRUE)
      ref_table[,i] <- names(sum_level)[which.max(sum_level)]
    }
  }



  if(plot){
    # X. Open a graphic file for plotting restults
    if(save.file=="pdf") pdf(paste(name, "pdf", sep="."))
    if(save.file=="jpeg") jpeg(paste(name, "jpeg", sep="."), width=ImageSize, height=ImageSize)
    if(save.file=="tiff") tiff(paste(name, "tiff", sep="."), width=ImageSize, height=ImageSize)
    if(save.file=="postscript") postscript(paste(name, "eps", sep="."))

    # XX. parametrize our plot window

    if(!do.bivariate){
      nb.graphs <- length(show.variables)
    } else{
      nb.graphs <- length(models) *  ( (length(show.variables)-1) * length(show.variables) / 2 )
    }

    if(legend) nb.graphs <- nb.graphs + 1

    if(display_title){
      W.width <- ceiling(sqrt(nb.graphs))
      W.height <- ceiling(nb.graphs/W.width)

      mat <- matrix(c(rep(1,W.width), 1:(W.height*W.width)+1), ncol=W.width, byrow=TRUE)
      layout(mat, widths=rep(1,W.width), heights=c(0.3,rep(1,W.height)))

      par(mar = c(0.1, 0.1, 0.1, 0.1))
      plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE)
      polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
      text(x=0.5, y=0.8, pos=1, cex=1.6, labels= main ,col="#4c57eb")
      par(mar = c(2,2,3.5,1))
    } else {
      par(mfrow=.CleverCut(nb.graphs))

    }

  }


  if(!do.bivariate){
    for(vari in show.variables){
      if(plot) {
        if(on_0_1000) ylim <- c(0,1000) else ylim <- c(0,1)

        if(is.factor(Data[,vari])) xlim <- c(1,length(levels(Data[,vari]))) else xlim=c(min(Data[,vari], na.rm=T), max(Data[,vari], na.rm=T))

        plot(0,0,col="white",xlim=xlim, ylim=ylim, main=vari, ann=TRUE, bty="o",xaxs="r", xaxt="s", xaxt=ifelse(is.factor(Data[,vari]),'n','t'),
             xlab="",ylab="")
        rug(Data[ ,vari])

        if(is.factor(Data[,vari])) axis(1, at=seq(xlim[1],xlim[2],1), labels=levels(Data[,vari]))

        # define color vector
        col <- rep(col,length.out=length(models))
        lty <- rep(lty,length.out=length(models))
      }


      # creating Tmp data
      if(is.factor(Data[,vari])) pts.tmp <- as.factor(levels(Data[,vari])) else pts.tmp <- seq(min(Data[,vari]), max(Data[,vari]), length.out=nb.pts)

      Data.r.tmp <- eval(parse(text=paste("cbind(",vari,"=pts.tmp,ref_table[,-which(colnames(ref_table)==vari),drop=F])",sep="")))
      Data.r.tmp <- Data.r.tmp[,colnames(ref_table),drop=F]
      if(length(factor_id)){
        for(f in factor_id){
          Data.r.tmp[,f] <- factor(as.character(Data.r.tmp[,f]), levels=levels(Data[,f]))
        }
      }


      for(model in models){


        # 0. get model
        mod <- get(model)
        mod.name <- ifelse(use.formal.names, formal_names[which(is.element(models, model))], model)

        # cat("\n*** model = ", model, ", mod.name =  ", mod.name)


        # 2. make projections
        proj.tmp <- predict(mod, Data.r.tmp, on_0_1000=on_0_1000, do_check=FALSE)

        # 4. Ploting results
        if(plot ) {
          if(is.factor(Data[,vari])){
            points(pts.tmp[1:length(levels(Data[,vari]))], proj.tmp[1:length(levels(Data[,vari]))], col=col[which(models==model)], lty = lty[which(models==model)])
          } else{
            lines(pts.tmp[1:nb.pts], proj.tmp[1:nb.pts], col=col[which(models==model)], lty = lty[which(models==model)])
          }
        }

        # 5. Storing results
        if(length(list.out[[vari]]) == 0){ #init
          eval(parse(text=paste("list.out[['",vari,"']] <- data.frame(",vari,"=pts.tmp, ",mod.name,"=proj.tmp)",sep="")))
        } else{
          eval(parse(text=paste("list.out[['",vari,"']] <- cbind(list.out[['",vari,"']],",mod.name,"=proj.tmp)",sep="")))
        }

      }

    }
    if(legend & plot){
      plot.new()
      legend(x="center",
             legend = formal_names,
             col = col,
             lty = lty,
             bty = 'n')
    }

  } else{ ## bivariate case
    for(vari1 in show.variables[-length(show.variables)]){
      for(vari2 in show.variables[-(1:which(show.variables == vari1))]){


        # creating Tmp data
        #         if(is.factor(Data[,vari])) pts.tmp <- as.factor(levels(Data[,vari])) else pts.tmp <- seq(min(Data[,vari]), max(Data[,vari]), length.out=nb.pts)

        pts.tmp1 <- rep(seq(min(Data[,vari1]), max(Data[,vari1]), length.out=sqrt(nb.pts)),each=sqrt(nb.pts))
        pts.tmp2 <- rep(seq(min(Data[,vari2]), max(Data[,vari2]), length.out=sqrt(nb.pts)),sqrt(nb.pts))

        Data.r.tmp <- eval(parse(text=paste("cbind(",vari1,"=pts.tmp1,",vari2,"=pts.tmp2, ref_table[,-which(colnames(ref_table)%in% c(vari1,vari2)),drop=F])",sep="")))
        Data.r.tmp <- Data.r.tmp[,colnames(ref_table),drop=F]
        if(length(factor_id)){
          for(f in factor_id){
            Data.r.tmp[,f] <- factor(as.character(Data.r.tmp[,f]), levels=levels(Data[,f]))
          }
        }

        for(model in models){

          # 0. get model
          mod <- get(model)
          mod.name <- ifelse(use.formal.names, formal_names[which(is.element(models, model))], model)

          # 2. make projections
          proj.tmp <- predict(mod, Data.r.tmp, on_0_1000=on_0_1000, do_check=FALSE)

          # 4. Storing results
          vari <- paste(vari1,vari2,sep="_")
          if(length(list.out[[vari]]) == 0){ #init
            eval(parse(text=paste("list.out[['",vari,"']] <- data.frame(",vari1,"=pts.tmp1,",vari2,"=pts.tmp2, ",mod.name,"=proj.tmp)",sep="")))
          } else{
            eval(parse(text=paste("list.out[['",vari,"']] <- cbind(list.out[['",vari,"']],",mod.name,"=proj.tmp)",sep="")))
          }

          # 5. Ploting results
          if(plot) {
            # reformating results to perform a persp plot
            pts.tmp1 <- sort(unique(pts.tmp1))
            pts.tmp2 <- sort(unique(pts.tmp2))
            proj.tmp <- matrix(proj.tmp, ncol=length(pts.tmp2), byrow=FALSE)

            # build color scale
            ncz <- length(pts.tmp2)
            nrz <- length(pts.tmp1)
            # Create a function interpolating colors in the range of specified colors
            jet.colors <- colorRampPalette(col)
            # Generate the desired number of colors from this palette
            nbcol <- 50
            color <- jet.colors(nbcol)
            # Compute the z-value at the facet centres
            zfacet <- proj.tmp[-1, -1] + proj.tmp[-1, -ncz] + proj.tmp[-nrz, -1] + proj.tmp[-nrz, -ncz]
            # Recode facet z-values into color indices
            facetcol <- cut(zfacet, nbcol)

            persp(x=pts.tmp1,y=pts.tmp2,z=proj.tmp, xlab = vari1, ylab=vari2, zlab="pred", zlim=c(0,1), theta = 30, phi = 30,
                  expand = 0.5, col = color[facetcol], ltheta = 120, shade = 0.25, ticktype = "simple", main = formal_names[which(models==model)], cex.main = 0.9, cex.axis=0.7)
          }

        }
      }
    }
  }

  # XXX. Close file
  if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()

  # delete temp files if somes has been created
  if(file.exists(file.path(get(models[1])@resp_name,'RespPlotTmp'))){
    unlink(path.expand(file.path(get(models[1])@resp_name,'RespPlotTmp')), recursive=TRUE, force=TRUE)
  }

  # transform list.out into ggplot firendly shape
  if(do.bivariate){
    gg.out <- .as.ggdat.2D(list.out)
  } else {
    gg.out <- .as.ggdat.1D(list.out)
  }

  invisible(gg.out)
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.response.plot2.check.arg <- function(models, Data, show.variables, save.file, name, ImageSize, plot, fixed.var.metric, do.bivariate, add.args){

  # 1. check add.args
  #   if(sum(! (names(add.args) %in% c("nb.pts","xy"))) > 0){
  #     warning(paste(toString(names(add.args)[which(! (names(add.args) %in% c("nb.pts")))]), " are unknown arguments", sep="" ))
  #   }


  ### check of models args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(!is.character(models)){
    stop("models must be a character vector of models names")
  }

  mod_names <- NULL
  for(mod in models){
    if(!exists(mod)){
      stop("you need to load the selected models!")
    }

    if(!inherits(get(mod), 'biomod2_model')){

      # create a biomod2 modeling object
      mod_tmp <- .Construct.default.biomod2.modeling.obj(get(mod))
      assign(mod_tmp@model_name, mod_tmp, envir = parent.frame(n = 1))
      mod_names <- c(mod_names, mod_tmp@model_name)
    } else{
      mod_names <- c(mod_names, mod)
    }
  }

  models <- mod_names


  ### defining the number split in each variables range =-=-=-=-=- #
  if(!is.null(add.args$nb.pts)){
    if(do.bivariate){
      # total number of points is the square of the difined
      add.args$nb.pts <- add.args$nb.pts^2
    }
  } else{
    if(!do.bivariate){
      add.args$nb.pts <- 100
    } else{
      add.args$nb.pts <- 25^2
    }
  }

  ### check of data args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(inherits(Data,"Raster")){
    cat("\n   > Extracting raster infos..")
    DataTmp <- matrix(0,ncol=nlayers(Data), nrow=add.args$nb.pts)
    colnames(DataTmp) <- names(Data)
    maxVal <- maxValue(Data)
    minVal <- minValue(Data)
    for(i in 1:ncol(DataTmp)){
      DataTmp[,i] <- seq(minVal[i],maxVal[i],length.out=add.args$nb.pts)
    }
    Data <- DataTmp
    rm(list=c('maxVal','minVal','DataTmp'))

  }

  ### check show.variables arg -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if( ( length(show.variables) > ncol(Data) ) | (sum(!(show.variables %in% colnames(Data)))) ) stop("columns wanted in show.variables do not match the data \n")

  # remove factorial var in do.bivariate case
  if(do.bivariate){
    fact_var <- sapply(Data[,show.variables, drop=F], is.factor)
    if(sum(fact_var)>0){
      cat("\n\tFactorial variables have been automatically removed!")
      show.variables <- show.variables[!fact_var]
    }
  }

  ### check save.file arg -=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #






  # TO DO
  return(list(models = models,
              Data = Data,
              show.variables = show.variables,
              save.file = save.file,
              name = name,
              ImageSize = ImageSize,
              plot = plot,
              fixed.var.metric = fixed.var.metric,
              do.bivariate = do.bivariate,
              nb.pts = add.args$nb.pts))
}

###

.Construct.default.biomod2.modeling.obj <- function(mod){

  ## ANN ##
  if(sum(!(c("nnet") %in% class(mod))) == 0){
    return(new("ANN_biomod2_model",
               model = mod,
               model_name = paste(ifelse(is.null(mod$terms[[2]]), "species",as.character(mod$terms[[2]])),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_ANN", sep=""),
               model_class = 'ANN',
               resp_name = ifelse(is.null(mod$terms[[2]]), "species",as.character(mod$terms[[2]])),
               expl_var_names = ifelse( is.character(attr( mod$terms,"term.labels")), attr( mod$terms,"term.labels"), "") ))
  }


  ## CTA ##
  if(sum(!(c("rpart") %in% class(mod)) == 0 )){

    return(new("CTA_biomod2_model",
               model = mod,
               model_name = paste(as.character(mod$terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_CTA", sep=""),
               model_class = 'CTA',
               resp_name = as.character(mod$terms[[2]]),
               expl_var_names = attr(mod$terms,"term.labels")))
  }

  ## FDA ##
  if(sum(!(c("fda") %in% class(mod)) == 0 )){
    return(new("FDA_biomod2_model",
               model = mod,
               model_name = paste(as.character(mod$terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_FDA", sep=""),
               model_class = 'FDA',
               resp_name = as.character(mod$terms[[2]]),
               expl_var_names = attr(mod$terms,"term.labels")))
  }

  ## GAM ##
  if(sum(!(c("gam") %in% class(mod)) == 0 )){
    return(new("GAM_biomod2_model",
               model = mod,
               model_subclass = ifelse(mod$method=="glm.fit","GAM_gam","GAM_mgcv"),
               model_name = paste(as.character(mod$terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_GAM", sep=""),
               model_class = 'GAM',
               resp_name = as.character(mod$terms[[2]]),
               expl_var_names = attr(mod$terms,"term.labels")))
  }

  ## GBM ##
  if(sum(!(c("gbm") %in% class(mod)) == 0 )){
    return(new("GBM_biomod2_model",
               model = mod,
               model_name = paste(as.character(mod$Terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_GBM", sep=""),
               model_class = 'GBM',
               resp_name = as.character(mod$Terms[[2]]),
               expl_var_names = attr(mod$Terms,"term.labels")))
  }

  ## GLM ##
  if(sum(!(c("glm", "lm") %in% class(mod)) == 0 )){
    return(new("GLM_biomod2_model",
               model = mod,
               model_name = paste(as.character(mod$terms[[2]]),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_GLM", sep=""),
               model_class = 'GLM',
               resp_name = as.character(mod$terms[[2]]),
               expl_var_names = attr(mod$terms,"term.labels")))
  }

  ## MARS ##
  if(sum(!(c("mars") %in% class(mod)) == 0 )){
    return(new("MARS_biomod2_model",
               model = mod,
               model_name =paste("species_AllData_",as.character(format(Sys.time(), "%OS6")),"_MARS",sep=""),
               model_class = 'MARS',
               resp_name = "species",
               expl_var_names = as.character(colnames(mod$factor))))
  }

  ## RF ##
  if(sum(!(c("randomForest") %in% class(mod)) == 0 )){
    return(new("RF_biomod2_model",
               model = mod,
               model_name =paste(ifelse(is.null(mod$terms[[2]]), "species",as.character(mod$terms[[2]])),"_AllData_",as.character(format(Sys.time(), "%OS6")),"_RF", sep=""),
               model_class = 'RF',
               resp_name = ifelse(is.null(mod$terms[[2]]), "species",as.character(mod$terms[[2]])),
               expl_var_names = ifelse( is.character(attr( mod$terms,"term.labels")), attr( mod$terms,"term.labels"), "") ))
  }

  stop("Unknown model class")

}

.as.ggdat.1D <-
  function (rp.dat)
  {
    requireNamespace('dplyr')
    out_ <- 
      bind_rows(
        lapply(
          rp.dat, 
          function(dat_) {
            dat_$id <- rownames(dat_)
            id.col.id <- which(colnames(dat_) == "id")
            expl.dat_ <- dat_ %>% 
              dplyr::select(1, id.col.id) %>%
              tidyr::gather("expl.name", "expl.val", 1)
            pred.dat_ <- dat_ %>% 
              dplyr::select(-1, id.col.id) %>%
              tidyr::gather("pred.name", "pred.val", (1:(ncol(dat_)-2)))
            out.dat_ <- 
              dplyr::full_join(expl.dat_, pred.dat_, by = 'id') %>%
              dplyr::mutate_at(c('expl.name', 'pred.name'), as.character) %>%
              dplyr::mutate_at('expl.val', as.numeric)
            return(out.dat_)
          }
        )
      )
    
    out_ <- 
      out_ %>%
      dplyr::mutate_at('expl.name', factor)
    
    return(out_)
  }


.as.ggdat.2D <- 
  function(rp.dat){
    out_ <- 
      bind_rows(
        lapply(
          rp.dat, 
          function(dat_) {
            dat_$id <- rownames(dat_)
            #   dat_$expl.name <- as.character(dat_$expl.name)
            #   dat_$pred.name <- as.character(dat_$pred.name)
            id.col.id <- which(colnames(dat_) == "id")
            expl1.dat_ <- dat_ %>% 
              dplyr::select(1, id.col.id) %>% 
              tidyr::gather("expl1.name", "expl1.val", 1)
            expl2.dat_ <- dat_ %>% 
              dplyr::select(2, id.col.id) %>% 
              tidyr::gather("expl2.name", "expl2.val", 1)
            pred.dat_ <- dat_ %>% 
              dplyr::select(3, id.col.id) %>% 
              tidyr::gather("pred.name", "pred.val", 1)
            out.dat_  <- 
              dplyr::full_join(
                dplyr::full_join(
                  expl1.dat_, 
                  expl2.dat_,
                  by = 'id'
                ), 
                pred.dat_,
                by = 'id'
              )
            out.dat_ <- out.dat_ %>%
              dplyr::mutate_at(c('expl1.name', 'expl2.name', 'pred.name'), as.character) %>% 
              dplyr::mutate_at(c('expl1.val', 'expl2.val'), as.numeric)
            return(out.dat_)
  }))
  ## ensure that the stips are in the right order
  out_ <- 
    out_ %>%
    dplyr::mutate_at(c('expl1.name', 'expl2.name'), factor)
  return(out_)
}
