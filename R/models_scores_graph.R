############################
# models_scores_graph fct
# Damien G. - 2014/10/22
############################

## Description ##
# This function is a graphic tools to represent evaluation
# scores of models produced within biomod2 according to 2
# different evalution methods. Models can be grouped in several
# ways (by algo, by cv run, ...) to highlight potential differences
# in models quality due to chosen models, cross validation sampling 
# bias,...

## Input ##
# modeling.output : an biomod2 modeling or ensemble-modeling modeling.outputect
# eval.metric : charcter vector of 2 chosen eval.metric (e.g c("ROC", "TSS"))
# by : the way models are grouped ('models', 'algos', 'cv_run' or 'data_set')
# plot : if you want to produce plot or not
# ... : several graphical options

## Ouput ##
# the ggplot2 modeling.outputect used to produce the graph is returned. That implies that 
# user should quite easily customize this plot.

##' 
##' @export
##' 

bm_PlotEvalMean <- function(modeling.output, eval.metric = NULL, group.by = 'Algo', plot = TRUE, ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PlotEvalMean.check.args(modeling.output, eval.metric, group.by, list(...))
  eval.metric <- args$eval.metric
  xlim <- args$xlim
  ylim <- args$ylim
  main <- args$main
  rm(args)
  
  ## 1. Get data for graphic ----------------------------------------------------------------------
  ## Get evaluation values
  scores <- get_evaluations(modeling.output, as.data.frame = TRUE)
  
  ## Choose which dataset (calibration or validation) should be used
  eval.data <- ifelse(all(is.na(scores$Evaluating.data)), "Testing.data", "Evaluating.data")
  
  ## Compute mean and sd evaluation scores
  models_mean = tapply(X = scores[, eval.data], INDEX = list(scores$Eval.metric, scores[, group.by]), FUN = mean, na.rm = TRUE)
  models_sd = tapply(X = scores[, eval.data], INDEX = list(scores$Eval.metric, scores[, group.by]), FUN = sd, na.rm = TRUE)
  
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
    theme(legend.title = element_blank())
  
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

.bm_PlotEvalMean.check.args <- function(modeling.output, eval.metric = NULL, group.by = 'Algo', args)
{
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
  .fun_testIfIn(TRUE, "group.by", group.by, c('Model', 'Algo', 'Run', 'Dataset'))

  ## 4. Check extra args argument ---------------------------------------------
  .fun_testIfIn(TRUE, "names(args)", names(args), c('xlim', 'ylim', 'main', 'col'))
  
  
  return(list(eval.metric = eval.metric,
              xlim = args$xlim,
              ylim = args$ylim,
              main = args$main))
} 

