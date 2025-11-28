###################################################################################################
##' @name BIOMOD_Report
##' @author Maya Guéguen
##' 
##' @title Produce summary outputs from a simulation folder
##' 
##' @description This function allows to produce summary report or ODMAP table from a 
##' \pkg{biomod2} simulation folder.
##' 
##' 
##' @param bm.out a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function ; or a 
##' \code{\link{BIOMOD.models.out}} or \code{\link{BIOMOD.ensemble.models.out}} 
##' object that can be obtained with the \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions
##' @param strategy a \code{character} defining the type of summary file that will be produced, 
##' must be \code{report}, \code{ODMAP} or \code{code} (see Details)
##' @param params.color a \code{list} containing 3 color values to custom the reports
##' @param params.ODMAP a \code{list} containing values of some ODMAP fields to be filled in 
##' from pre-existing choices (see Details)
##' 
##' 
##' @return 
##' 
##' A standardized \code{.html} file obtained from an \code{Rmarkdown} template, and a \code{.csv} 
##' table in the case of ODMAP report.
##' 
##' 
##' @details 
##' 
##' This function gathers and formats all objects contained in one \pkg{biomod2} modeling folder 
##' to produce, based on \code{Rmarkdown} templates, standardized reports to help the user :
##' \itemize{
##'   \item summarize its modeling results
##'   \item share them through standardized informations through ODMAP protocol
##'   \item provide reproducible code \cr \cr
##' }
##' 
##' 
##' \describe{
##'   \item{Type of report}{
##'   Different data types are available, and require different values :
##'   \describe{
##'     \item{report}{\pkg{biomod2} provides functions to summarize the information, such as 
##'     \code{print}, \code{plot} or \code{summary} methods adapted to \code{BIOMOD.[...].out} 
##'     objects, as well as 
##'     \href{https://biomodhub.github.io/biomod2/reference/getters.out.html}{\code{get_[...]}} 
##'     and \code{bm_Plot[...]} functions. All these are called here and applied to objects 
##'     contained in the provided modeling folder.
##'     }
##'     \item{ODMAP}{following Zurell et al. 2020, ODMAP (Overview, Data, Model, Assessment and 
##'     Prediction) protocol aims to standardize documentation of modeling to help improve both 
##'     transparency and reproducibility of results. 
##'     \href{https://odmap.wsl.ch/}{ODMAP v1.0 website} provides an application to fill this type 
##'     of report. \pkg{biomod2} tries here to help one user to pre-fill the fields of this 
##'     protocol.
##'     }
##'     \item{code}{\code{call} slot contained within \code{\link{BIOMOD.formated.data}}, 
##'     \code{BIOMOD.models.out}, \code{BIOMOD.projection.out} and 
##'     \code{BIOMOD.ensemble.models.out} objects keep in memory the R command used to obtain 
##'     them. All these calls are gathered here in one summary file.}
##'   }
##'   }
##' }
##' 
##'   
##' @references
##' 
##' \itemize{
##'   \item Zurell D, Franklin J, König C, Bouchet PJ, Serra-Diaz JM, Dormann CF, Elith J, 
##'   Fandos Guzman G, Feng X, Guillera-Arroita G, Guisan A, Leitão PJ, Lahoz-Monfort JJ, 
##'   Park DS, Peterson AT, Rapacciuolo G, Schmatz DR, Schröder B, Thuiller W, Yates KL, 
##'   Zimmermann NE, Merow C (2020) A standard protocol for describing species distribution 
##'   models. \emph{Ecography} \bold{43}: 1261-1277. \doi{10.1111/ecog.04960}
##' }
##'
##' @keywords report ODMAP markdown html
##' 
##' 
##' @seealso \code{\link{ODMAP}}, \code{\link[base]{match.call}}
##' @family Primary functions
##' 
##' 
##' @examples
##' library(terra)
##' 
##' 
# @importFrom foreach foreach %do% %:%
##' 
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_Report <- function(bm.out
                          , strategy = 'report'
                          , params.color = list(color1 = "#eb4034"
                                                , color2 = "#e0918b"
                                                , color3 = "#658f70")
                          , params.ODMAP = list(O.mod.objective = NULL
                                                , O.boundary = NULL
                                                , O.obs.type = NULL
                                                , O.pred.type = NULL
                                                , D.eco.level = NULL
                                                , D.samp.design = NULL))
{
  .bm_cat("Do biomod2 Report")
  bm.files = bm.mod = bm.ens = bm.form = NA
  
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_Report.check.args(bm.out, strategy, params.color, params.ODMAP)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)

    
  ## 1. Get modeling files ------------------------------------------------------------------------
  cat("\n\t> Getting modeling files...")
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
        warning(paste0("No bm.mod selected but some available. Please select one among : "
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
  cat("\n\t> Building report...")
  out <- switch(strategy,
                report = BIOMOD_Report_summary(sp.name, dir.name, bm.files, bm.form, bm.mod, bm.ens, params.color),
                ODMAP = BIOMOD_Report_ODMAP(sp.name, dir.name, bm.files, bm.form, bm.mod, bm.ens, params.color, params.ODMAP),
                code = BIOMOD_Report_code(sp.name, dir.name, bm.files, bm.form, bm.mod, bm.ens, params.color))
  
  .bm_cat("Done")
}


###################################################################################################

.BIOMOD_Report.check.args <- function(bm.out, strategy, params.color, params.ODMAP)
{
  ## check namespace ----------------------------------------------------------
  if (!isNamespaceLoaded("rmarkdown")) { 
    if (!requireNamespace('rmarkdown', quietly = TRUE)) stop("Package 'rmarkdown' not found")
  }
  
  
  sp.name = dir.name = name.bm.mod = bm.form = NA
  
  ## 1. Check bm.out -------------------------------------------------------
  test.format <- inherits(bm.out, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  test.models <- inherits(bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  if (test.format && !test.models) {
    sp.name <- bm.out@sp.name
    dir.name <- bm.out@dir.name
    bm.form <- bm.out
  } else if (!test.format && test.models) {
    sp.name <- bm.out@sp.name
    dir.name <- bm.out@dir.name
    if (inherits(bm.out, "BIOMOD.models.out")) {
      name.bm.mod <- sub(paste0(".*", sp.name, "/"), paste0(sp.name, "/"), bm.out@link)
    } else {
      name.bm.mod <- sub(paste0(".*", sp.name, "/"), paste0(sp.name, "/"), bm.out@models.out@link)
    }
  } else {
    .fun_testIfInherits(TRUE, "bm.out", bm.out, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"
                                                  , "BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
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
  
  ## 4. Check params.color argument -------------------------------------------
  .fun_testIfIn(TRUE, "names(params.color)", names(params.color), c("color1", "color2", "color3"))
  
  
  return(list(sp.name = sp.name
              , dir.name = dir.name
              , name.bm.mod = name.bm.mod
              , bm.form = bm.form
              , params.color = params.color
              , params.ODMAP = params.ODMAP))
}


###################################################################################################

BIOMOD_Report_summary <- function(sp.name, dir.name, bm.files, bm.form, bm.mod, bm.ens, params.color)
{
  out <- rmarkdown::render(input = system.file("rmd", "biomod2_template_report.Rmd", package = "biomod2")
                           , params = list(sp.name = sp.name
                                           , bm.files = bm.files
                                           , bm.mod = bm.mod
                                           , bm.ens = bm.ens
                                           , bm.form = bm.form
                                           , params.color = params.color)
                           , output_format = "html_document"
                           , output_file = paste0("biomod2_report_", sp.name, ".html")
                           , output_dir = dir.name
                           , knit_root_dir = dir.name
                           , quiet = FALSE)
  cat("\n\t> Report has been created : ", out)
  cat("\n")
}

BIOMOD_Report_ODMAP <- function(sp.name, dir.name, bm.files, bm.form, bm.mod, bm.ens, params.color, params.ODMAP)
{
  out <- rmarkdown::render(input = system.file("rmd", "biomod2_template_ODMAP.Rmd", package = "biomod2")
                           , params = list(sp.name = sp.name
                                           , bm.files = bm.files
                                           , bm.mod = bm.mod
                                           , bm.ens = bm.ens
                                           , bm.form = bm.form
                                           , params.color = params.color
                                           , params.ODMAP = params.ODMAP)
                           , output_format = "html_document"
                           , output_file = paste0("biomod2_ODMAP_", sp.name, ".html")
                           , output_dir = dir.name
                           , knit_root_dir = dir.name
                           , quiet = FALSE)
  cat("\n\t> ODMAP has been created : ", out)
  cat("\n")
}

BIOMOD_Report_code <- function(sp.name, dir.name, bm.files, bm.form, bm.mod, bm.ens, params.color)
{
  out <- rmarkdown::render(input = system.file("rmd", "biomod2_template_code.Rmd", package = "biomod2")
                           , params = list(sp.name = sp.name
                                           , bm.files = bm.files
                                           , bm.mod = bm.mod
                                           , bm.ens = bm.ens
                                           , bm.form = bm.form
                                           , params.color = params.color)
                           , output_format = "html_document"
                           , output_file = paste0("biomod2_code_", sp.name, ".html")
                           , output_dir = dir.name
                           , knit_root_dir = dir.name
                           , quiet = FALSE)
  cat("\n\t> Code has been created : ", out)
  cat("\n")
}

