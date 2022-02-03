###################################################################################################
##' @name bm_PlotEvalMean
##' @author Damien Georges, Maya Gueguen
##' 
##' @title Plot mean evaluation scores
##' 
##' @description
##' 
##' This function represents mean evaluation scores of species distribution models, from 
##' \code{BIOMOD.models.out} or \code{BIOMOD.ensemble.models.out} objects that can be obtained 
##' from \code{\link{BIOMOD_Modeling}} or \code{\link{BIOMOD_EnsembleModeling}} functions. Scores 
##' are represented according to 2 different evaluation methods, and models can be grouped (see 
##' \href{bm_PlotEvalMean.html#details}{Details}).
##' 
##' @param modeling.output a \code{\link{BIOMOD.models.out}} or \code{BIOMOD.ensemble.models.out} 
##' object that can be obtained from \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions
##' @param eval.metric a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{HK}, \code{HSS}, \code{OR}, 
##' \code{ORSS}
##' @param group.by a \code{character} corresponding to the way kept models will be combined to 
##' compute mean and sd evaluation scores, must be among \code{model}, \code{algo}, \code{run}, 
##' \code{dataset}
##' @param plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the plot is to be rendered or not
##' @param \ldots some additional arguments (see \href{bm_PlotEvalMean.html#details}{Details})
##' 
##' 
##' @return  
##' 
##' A \code{ggplot} object representing mean and standard deviation of evaluation scores 
##' according to 2 different evaluation methods.
##' 
##' 
##' @details
##' 
##' \code{...} can take the following values :
##' 
##' \itemize{
##'   \item{\code{xlim}}{ : an \code{integer} corresponding to the x maximum limit to represent}
##'   \item{\code{ylim}}{ : an \code{integer} corresponding to the y maximum limit to represent}
##'   \item{\code{main}}{ : a \code{character} corresponding to the graphic title}
##'   \item{\code{col}}{ : a \code{vector} containing new color values}
##' }
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{get_evaluations}}
##' 
##' @keywords evaluation, ggplot
##' 
##' 
##' @examples
##' 
##' 
##' @importFrom ggplot2 ggplot aes_string geom_point geom_errorbarh geom_errorbar xlab ylab
##' theme element_blank element_rect coord_cartesian labs
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PlotEvalMean <- function(modeling.output, eval.metric = NULL, group.by = 'algo', plot = TRUE, ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PlotEvalMean.check.args(modeling.output, eval.metric, group.by, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  tmp = strsplit(group.by, '')[[1]]
  group.by <- paste0(toupper(tmp[1]), paste0(tmp[2:length(tmp)], collapse = ''))
  
  ## 1. Get data for graphic ----------------------------------------------------------------------
  ## Get evaluation values
  scores <- get_evaluations(modeling.output, as.data.frame = TRUE)
  
  ## Choose which dataset (calibration or validation) should be used
  eval.data <- ifelse(all(is.na(scores$Evaluating.data)), "Testing.data", "Evaluating.data")
  
  ## Compute mean and sd evaluation scores
  models_mean = tapply(X = scores[, eval.data]
                       , INDEX = list(scores$Eval.metric, scores[, group.by])
                       , FUN = mean, na.rm = TRUE)
  models_sd = tapply(X = scores[, eval.data]
                     , INDEX = list(scores$Eval.metric, scores[, group.by])
                     , FUN = sd, na.rm = TRUE)
  
  ## Prepare data table for graphic
  ggdat <- merge(data.frame(name = colnames(models_mean), t(models_mean)),
                 data.frame(name = colnames(models_sd), t(models_sd)), 
                 by = "name" )
  colnames(ggdat) <- c("name", "mean1", "mean2", "sd1", "sd2")
  
  limits1 <- aes_string(xmax = "mean1 + sd1", xmin = "mean1 - sd1", fill = NULL)
  limits2 <- aes_string(ymax = "mean2 + sd2", ymin = "mean2 - sd2", fill = NULL)
  
  ## 2. PLOT graphic ------------------------------------------------------------------------------
  gg <- ggplot(ggdat, aes_string(x = "mean1", y = "mean2", colour = "name", fill = NULL)) +
    geom_point() + ## add mean points
    geom_errorbarh(limits1, height = 0) + ## add horizontal error bars
    geom_errorbar(limits2, width = 0) + ## add vertical error bars
    xlab(eval.metric[1]) +
    ylab(eval.metric[2]) +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white"))
  
  if (length(ylim) | length(xlim)) { ## fix scale
    gg <- gg + coord_cartesian(ylim = ylim, xlim = xlim)
  }
  
  if (length(main)) { ## add title
    gg <- gg + labs(title = main)
  }
  
  if (plot){ print(gg) }
  invisible(gg)
}


###################################################################################################

.bm_PlotEvalMean.check.args <- function(modeling.output, eval.metric = NULL, group.by = 'Algo', ...)
{
  args <- list(...)
  
  ## 1. Check modeling.output argument ----------------------------------------
  .fun_testIfInherits(TRUE, "modeling.output", modeling.output, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 2. Check eval.metric argument --------------------------------------------
  scores <- get_evaluations(modeling.output, as.data.frame = TRUE)
  avail_eval.metric <- unique(scores$Eval.metric)
  if (length(avail_eval.metric) < 2) { stop("At least 2 different evaluations eval.metric needed") }
  if (is.null(eval.metric)) {
    eval.metric <- avail_eval.metric[1:2]
    warnings(toString(eval.metric), " evaluation eval.metric automatically selected")
  }
  
  ## 3. Check group.by argument -----------------------------------------------
  .fun_testIfIn(TRUE, "group.by", group.by, c('model', 'algo', 'run', 'dataset'))

  ## 4. Check extra args argument ---------------------------------------------
  .fun_testIfIn(TRUE, "names(args)", names(args), c('xlim', 'ylim', 'main'))
  
  
  return(list(eval.metric = eval.metric,
              xlim = args$xlim,
              ylim = args$ylim,
              main = args$main))
} 

