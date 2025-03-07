###################################################################################################
##' @name bm_Report
##' 
##' 
###################################################################################################


bm_Report <- function(bm.out
                      , strategy = 'report'
                      , params.ODMAP = list(O.mod.objective = NULL
                                            , O.boundary = NULL
                                            , O.obs.type = NULL
                                            , O.pred.type = NULL
                                            , D.eco.level = NULL
                                            , D.samp.design = NULL))
{
  bm.files = bm.mod = bm.ens = bm.form = NA
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_Report.check.args(bm.out, strategy, params.ODMAP)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## 1. Get modeling files ------------------------------------------------------------------------
  bm.files = list.files(path = sp.name, pattern = ".out$", recursive = TRUE, full.names = TRUE)
  if (length(bm.files) > 0) {
    bm.files = data.frame(file = bm.files
                          , path = dirname(bm.files)
                          , sp = sapply(basename(bm.files), function(x) strsplit(x, "[.]")[[1]][1])
                          , type = sapply(basename(bm.files), function(x) rev(strsplit(x, "[.]")[[1]])[2])
                          , level = sapply(basename(bm.files), function(x) strsplit(x, "[.]")[[1]][3])
                          , ID = sapply(basename(bm.files), function(x) strsplit(x, "[.]")[[1]][2])
    )
    bm.files$level = ifelse(bm.files$level == "ensemble", "ensemble", "single")
    bm.files$level = factor(bm.files$level, c("single", "ensemble"))
    bm.files = bm.files[order(bm.files$type, bm.files$level, bm.files$ID), ]
    rownames(bm.files) = NULL
    
    bm.files$refer_to = NA
    ind = which(bm.files$type == "models" & bm.files$level == "single")
    bm.files$refer_to[ind] = letters[ind]
    ind.na = which(is.na(bm.files$refer_to))
    for (i in ind.na) {
      tmp = get(load(file = bm.files$file[i]))
      # tmp.link = sub("[.]/", "", tmp@models.out@link)
      tmp.link = sub(paste0(".*", sp.name, "/"), paste0(sp.name, "/"), tmp@models.out@link)
      bm.files$refer_to[i] = bm.files$refer_to[which(bm.files$file == tmp.link)]
      suppressWarnings(rm(list = basename(bm.files$file[i])))
    }
    
    ## 
    if (is.na(name.bm.mod)) {
      ind.mod = which(bm.files$type == "models" &
                        bm.files$level == "single")
      if (length(ind.mod) > 0) {
        warning(paste0("No bm.mod selected but some available. Please select one among :"
                       , paste0(bm.files$file[ind.mod], collapse = ", ")))
      }
    } else {
      
      ## load-modeling-files
      ind.mod = which(bm.files$file == name.bm.mod)
      bm.mod = get(load(file = name.bm.mod))
      
      ## add potentially missing slots
      if (inherits(try(bm.mod@data.type), "try-error")) {
        bm.mod@data.type <- "binary"
      }
      
      ind.ens = which(bm.files$type == "models" &
                        bm.files$level == "ensemble" &
                        bm.files$refer_to == bm.files$refer_to[ind.mod])
      if (length(ind.ens) > 0) {
        name.bm.ens = bm.files$file[ind.ens]
        bm.ens = get(load(file = name.bm.ens))
        
        ## add potentially missing slots
        if (inherits(try(bm.ens@data.type), "try-error")) {
          bm.ens@data.type <- "binary"
        }
      }
      suppressWarnings(rm(list = basename(bm.files$file)))
      
      ## get-formated-data
      bm.form = get_formal_data(bm.mod)
    }
  } else {
    bm.files = NA
  }
  
  ## add potentially missing slots
  if (inherits(try(bm.form@has.filter.raster), "try-error")) {
    bm.form@has.filter.raster <- FALSE
  }
  if (inherits(try(bm.form@data.type), "try-error")) {
    bm.form@data.type <- "binary"
  }
  
  
  ## 2. Create output object ----------------------------------------------------------------------
  out <- switch(strategy,
                report = bm_Report_report(sp.name, bm.files, bm.mod, bm.ens, bm.form),
                ODMAP = bm_Report_ODMAP(sp.name, bm.files, bm.mod, bm.ens, bm.form, params.ODMAP))
  cat("\n")
  return(out)
}


###################################################################################################

.bm_Report.check.args <- function(bm.out, strategy, params.ODMAP)
{
  sp.name = name.bm.mod = bm.form = NA
  
  ## 1. Check bm.out -------------------------------------------------------
  test.format <- inherits(bm.out, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  test.models <- inherits(bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  if (test.format && !test.models) {
    sp.name <- bm.out@sp.name
    bm.form <- bm.out
  } else if (!test.format && test.models) {
    sp.name <- bm.out@sp.name
    if (inherits(bm.out, "BIOMOD.models.out")) {
      name.bm.mod <- sub(paste0(".*", sp.name, "/"), paste0(sp.name, "/"), bm.out@link)
    } else {
      name.bm.mod <- sub(paste0(".*", sp.name, "/"), paste0(sp.name, "/"), bm.out@models.out@link)
    }
  } else {
    stop("")
  }
  
  
  ## 2. Check strategy argument -----------------------------------------------
  availableReport <- c("report", "ODMAP", "code")
  if (is.null(strategy) || !(strategy %in% availableReport)) {
    strategy <- "report"
    cat("\n   ! Report was automatically selected")
  }
  
  ## 3. Check params.ODMAP argument -------------------------------------------
  if (strategy == "ODMAP") {
    .fun_testIfIn(TRUE, "names(params.ODMAP)", names(params.ODMAP)
                  , c("O.mod.objective", "O.boundary", "O.obs.type"
                      , "O.pred.type", "D.eco.level", "D.samp.design")
                  , exact = FALSE)
    if (!is.null(params.ODMAP$O.mod.objective)) {
      .fun_testIfIn(TRUE, "params.ODMAP$O.mod.objective", params.ODMAP$O.mod.objective
                    , c("Inference and explanation", "Mapping and interpolation", "Forecast and transfer"))
    } else { params.ODMAP$O.mod.objective = NA }
    if (!is.null(params.ODMAP$O.boundary)) {
      .fun_testIfIn(TRUE, "params.ODMAP$O.boundary", params.ODMAP$O.boundary
                    , c("natural", "political", "rectangle"))
    } else { params.ODMAP$O.boundary = NA }
    if (!is.null(params.ODMAP$O.obs.type)) {
      .fun_testIfIn(TRUE, "params.ODMAP$O.obs.type", params.ODMAP$O.obs.type
                    , c("citizen science", "field survey", "GPS tracking"
                        , "range map", "standardised monitoring data"))
    } else { params.ODMAP$O.obs.type = NA }
    if (!is.null(params.ODMAP$O.pred.type)) {
      .fun_testIfIn(TRUE, "params.ODMAP$O.pred.type", params.ODMAP$O.pred.type
                    , c("climatic", "edaphic", "habitat", "topographic"))
    } else { params.ODMAP$O.pred.type = NA }
    if (!is.null(params.ODMAP$D.eco.level)) {
      .fun_testIfIn(TRUE, "params.ODMAP$D.eco.level", params.ODMAP$D.eco.level
                    , c("communities", "individuals", "operational taxonomic units", "populations", "species"))
    } else { params.ODMAP$D.eco.level = NA }
    if (!is.null(params.ODMAP$D.samp.design)) {
      .fun_testIfIn(TRUE, "params.ODMAP$D.samp.design", params.ODMAP$D.samp.design
                    , c("spatial random", "spatial uniform", "spatial stratified", "temporal", "nestedness"))
    } else { params.ODMAP$D.samp.design = NA }
  }
  
  return(list(sp.name = sp.name
              , name.bm.mod = name.bm.mod
              , bm.form = bm.form
              , params.ODMAP = params.ODMAP))
}


###################################################################################################

##'
##' @rdname bm_Report
##' @export
##'

bm_Report_report <- function(sp.name, bm.files, bm.mod, bm.ens, bm.form)
{
  # rmarkdown::render(input = "biomod2_template_report.Rmd"
  rmarkdown::render(input = system.file("rmd", "biomod2_template_report.Rmd", package = "biomod2")
                    , params = list(sp.name = sp.name
                                    , bm.files = bm.files
                                    , bm.mod = bm.mod
                                    , bm.ens = bm.ens
                                    , bm.form = bm.form)
                    , output_format = "html_document"
                    , output_file = paste0("biomod2_report_", sp.name, ".html"))
}

##'
##' @rdname bm_Report
##' @export
##'

bm_Report_ODMAP <- function(sp.name, bm.files, bm.mod, bm.ens, bm.form, params.ODMAP)
{
  # rmarkdown::render(input = paste0("biomod2_template_ODMAP.Rmd")
  rmarkdown::render(input = system.file("rmd", "biomod2_template_ODMAP.Rmd", package = "biomod2")
                    , params = list(sp.name = sp.name
                                    , bm.files = bm.files
                                    , bm.mod = bm.mod
                                    , bm.ens = bm.ens
                                    , bm.form = bm.form
                                    , params.ODMAP = params.ODMAP)
                    , output_format = "html_document"
                    , output_file = paste0("biomod2_ODMAP_", sp.name, ".html"))
}