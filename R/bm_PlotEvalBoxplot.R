###################################################################################################
##' @name bm_PlotEvalBoxplot
##' @author Damien Georges, Maya Gueguen
##' 
##' @title Plot boxplot of evaluation scores
##' 
##' @description
##' 
##' This function represents boxplot of evaluation scores of species distribution models, from 
##' \code{BIOMOD.models.out} or \code{BIOMOD.ensemble.models.out} objects that can be obtained 
##' from \code{\link{BIOMOD_Modeling}} or \code{\link{BIOMOD_EnsembleModeling}} functions. Scores 
##' are represented according to 2 grouping methods (see 
##' \href{bm_PlotEvalBoxplot.html#details}{Details}).
##' 
##' @param modeling.output a \code{\link{BIOMOD.models.out}} or \code{BIOMOD.ensemble.models.out} 
##' object that can be obtained from \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions
##' @param group.by a 2-length \code{vector} containing the way kept models will be represented,
##' must be among \code{model}, \code{algo}, \code{run}, \code{dataset}
##' @param plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the plot is to be rendered or not
##' @param \ldots some additional arguments (see \href{bm_PlotEvalMean.html#details}{Details})
##' 
##' 
##' @return  
##' 
##' A \code{ggplot} object representing boxplot of evaluation scores.
##' 
##' 
##' @details
##' 
##' \code{...} can take the following values :
##' 
##' \itemize{
##'   \item{\code{main}}{ : a \code{character} corresponding to the graphic title}
##'   \item{\code{scales}}{ : a \code{character} corresponding to the \code{scales} argument of 
##'   the \code{\link[ggplot2]{facet_wrap}} function, must be either \code{fixed}, \code{free_x}, 
##'   \code{free_y} or \code{free}}
##' }
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{get_evaluations}}
##' 
##' @keywords evaluation, ggplot, boxplot
##' 
##' 
##' @examples
##' 
##' 
##' @importFrom ggplot2 ggplot aes_string geom_boxplot facet_wrap xlab ylab
##' theme element_blank element_rect labs
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PlotEvalBoxplot <- function(modeling.output, group.by = c('algo', 'run'), plot = TRUE, ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PlotEvalBoxplot.check.args(modeling.output, group.by, list(...))
  main <- args$main
  scales <- args$scales
  rm(args)
  
  for (i in 1:length(group.by)) {
    tmp = strsplit(group.by[i], '')[[1]]
    group.by[i] <- paste0(toupper(tmp[1]), paste0(tmp[2:length(tmp)], collapse = ''))
  }
  
  ## 1. Get data for graphic ----------------------------------------------------------------------
  ## Get evaluation values
  scores <- get_evaluations(modeling.output, as.data.frame = TRUE)
  
  ## Choose which dataset (calibration or validation) should be used
  eval.data <- ifelse(all(is.na(scores$Evaluating.data)), "Testing.data", "Evaluating.data")
  
  ## Prepare data table for graphic
  ggdat = scores
  
  ## 2. PLOT graphic ------------------------------------------------------------------------------
  gg <- ggplot(ggdat, aes_string(x = group.by[1], y = eval.data, fill = group.by[2])) +
    geom_boxplot() + ## add boxplot
    facet_wrap(~ Eval.metric, scales = scales) +
    xlab("") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white"))
  
  if (length(main)) { ## add title
    gg <- gg + labs(title = main)
  }
  
  if (plot){ print(gg) }
  invisible(gg)
}


###################################################################################################

.bm_PlotEvalBoxplot.check.args <- function(modeling.output, group.by = 'Algo', args)
{
  ## 1. Check modeling.output argument ----------------------------------------
  .fun_testIfInherits(TRUE, "modeling.output", modeling.output, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 3. Check group.by argument -----------------------------------------------
  if (length(group.by) != 2) { stop("2 group values needed") }
  for (i in 1:length(group.by)) {
    .fun_testIfIn(TRUE, paste0("group.by[", i, "]"), group.by[i], c('model', 'algo', 'run', 'dataset'))
  }

  ## 4. Check extra args argument ---------------------------------------------
  .fun_testIfIn(TRUE, "names(args)", names(args), c('main', 'scales'))
  if ("scales" %in% names(args)) {
    .fun_testIfIn(TRUE, "args$scales", args$scales, c('fixed', 'free_x', 'free_y', 'free'))
  } else {
    args$scales = "fixed"
  }
  
  
  return(list(main = args$main,
              scales = args$scales))
} 

