###################################################################################################
##' @name bm_PlotVarImpBoxplot
##' @author Damien Georges, Maya Gueguen
##' 
##' @title Plot boxplot of variables importance
##' 
##' @description
##' 
##' This function represents boxplot of variables importance of species distribution models, from 
##' \code{BIOMOD.models.out} or \code{BIOMOD.ensemble.models.out} objects that can be obtained 
##' from \code{\link{BIOMOD_Modeling}} or \code{\link{BIOMOD_EnsembleModeling}} functions. Scores 
##' are represented according to 3 grouping methods (see 
##' \href{bm_PlotVarImpBoxplot.html#details}{Details}).
##' 
##' @param modeling.output a \code{\link{BIOMOD.models.out}} or \code{BIOMOD.ensemble.models.out} 
##' object that can be obtained from \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions
##' @param group.by a 2-length \code{vector} containing the way kept models will be represented,
##' must be among \code{model}, \code{algo}, \code{run}, \code{dataset}, \code{expl.var}
##' @param plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the plot is to be rendered or not
##' @param \ldots some additional arguments (see \href{bm_PlotEvalMean.html#details}{Details})
##' 
##' 
##' @return  
##' 
##' A \code{ggplot} object representing boxplot of variables importance.
##' 
##' 
##' @details
##' 
##' \code{...} can take the following values :
##' 
##' \itemize{
##'   \item{\code{main}}{ : a \code{character} corresponding to the graphic title}
##' }
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{get_variables_importance}}
##' 
##' @keywords evaluation, ggplot, boxplot
##' 
##' 
##' @examples
##' 
##' 
##' @importFrom ggplot2 ggplot aes_string geom_boxplot facet_wrap xlab ylab
##' theme element_blank element_rect labs scale_y_continuous
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PlotVarImpBoxplot <- function(modeling.output, group.by = c('run', 'expl.var', 'algo'), plot = TRUE, ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PlotVarImpBoxplot.check.args(modeling.output, group.by, list(...))
  main <- args$main
  rm(args)
  
  for (i in 1:length(group.by)) {
    tmp = strsplit(group.by[i], '')[[1]]
    group.by[i] <- paste0(toupper(tmp[1]), paste0(tmp[2:length(tmp)], collapse = ''))
  }
  
  ## 1. Get data for graphic ----------------------------------------------------------------------
  ## Get variables importance values
  scores <- get_variables_importance(modeling.output, as.data.frame = TRUE)
  
  ## Prepare data table for graphic
  ggdat = scores
  
  ## 2. PLOT graphic ------------------------------------------------------------------------------
  gg <- ggplot(ggdat, aes_string(x = group.by[1], y = "Var.imp", fill = group.by[2])) +
    geom_boxplot() + ## add boxplot
    facet_wrap(group.by[3], scales = "free_x") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%")) + 
    xlab("") +
    ylab("") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white"))
  
  if (length(main)) { ## add title
    gg <- gg + labs(title = main)
  }
  
  if (plot){ print(gg) }
  invisible(gg)
}


###################################################################################################

.bm_PlotVarImpBoxplot.check.args <- function(modeling.output, group.by = 'Algo', args)
{
  ## 1. Check modeling.output argument ----------------------------------------
  .fun_testIfInherits(TRUE, "modeling.output", modeling.output, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 3. Check group.by argument -----------------------------------------------
  if (length(group.by) != 3) { stop("3 group values needed") }
  for (i in 1:length(group.by)) {
    .fun_testIfIn(TRUE, paste0("group.by[", i, "]"), group.by[i], c('model', 'algo', 'run', 'dataset', 'expl.var'))
  }

  ## 4. Check extra args argument ---------------------------------------------
  .fun_testIfIn(TRUE, "names(args)", names(args), c('main'))
  
  
  return(list(main = args$main))
} 

