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
# obj : an biomod2 modeling or ensemble-modeling object
# metrics : charcter vector of 2 chosen metrics (e.g c("ROC", "TSS"))
# by : the way models are grouped ('models', 'algos', 'cv_run' or 'data_set')
# plot : if you want to produce plot or not
# ... : several graphical options

## Ouput ##
# the ggplot2 object used to produce the graph is returned. That implies that 
# user should quite easily customize this plot.

## Main code ##
models_scores_graph <- function(obj, metrics = NULL, by = 'models', plot = TRUE, ...){
  
  ## get additional args
  args = list(...)
  
  ## check args
  ck_args <- .models_scores_graph_check_args(obj, metrics = metrics, by = by, args = args)
  
  metrics <- ck_args$metrics
  xlim <- ck_args$xlim
  ylim <- ck_args$ylim
  main <- ck_args$main
  
  ## get models scores
  scores <- get_evaluations(obj, as.data.frame = T )
  
  ## add some columns to enable different wy of grouping scores
  scores$mod = sapply(as.character(scores$Model.name), function(x) { xx <- unlist( strsplit(x,"_") ); xx[1] } )
  scores$alg = sapply(as.character(scores$Model.name), function(x) { xx <- unlist( strsplit(x,"_") ); xx[length(xx) - 2] } )
  scores$run = sapply(as.character(scores$Model.name), function(x) { xx <- unlist( strsplit(x,"_") ); xx[length(xx) - 1] } )
  scores$dat = sapply(as.character(scores$Model.name), function(x) { xx <- unlist( strsplit(x,"_") ); xx[length(xx)] } )
  
  ## extract some usefull info
  eval_data <- ifelse( all(is.na( scores$Evaluating.data )), "Testing.data", "Evaluating.data")
  
  ## calculate summaries statistics
  
  ### models mean scores calculation
  models_mean <- switch(by,
                        models = sapply( unique(scores$mod) , 
                                         function(x){ 
                                           sapply( metrics,  
                                                   function(xx, x){
                                                     mean(scores[ scores$mod == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                   },
                                                   x = x ) 
                                           }),
                        algos = sapply( unique(scores$alg) , 
                                         function(x){ 
                                           sapply( metrics,  
                                                   function(xx, x){
                                                     mean(scores[ scores$alg == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                   },
                                                   x = x ) 
                                         }),
                        cv_run = sapply( unique(scores$run) , 
                                      function(x){ 
                                        sapply( metrics,  
                                                function(xx, x){
                                                  mean(scores[ scores$run == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                },
                                                x = x ) 
                                      }),
                        data_set = sapply( unique(scores$dat) , 
                                           function(x){ 
                                             sapply( metrics,  
                                                     function(xx, x){
                                                       mean(scores[ scores$dat == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                     },
                                                     x = x ) 
                                           }) )
  
  ### sd of models scores calculation
  models_sd <- switch(by,
                        models = sapply( unique(scores$mod) , 
                                         function(x){ 
                                           sapply( metrics,  
                                                   function(xx, x){
                                                     sd(scores[ scores$mod == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                   },
                                                   x = x ) 
                                         }),
                        algos = sapply( unique(scores$alg) , 
                                       function(x){ 
                                         sapply( metrics,  
                                                 function(xx, x){
                                                   sd(scores[ scores$alg == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                 },
                                                 x = x ) 
                                       }),
                        cv_run = sapply( unique(scores$run) , 
                                      function(x){ 
                                        sapply( metrics,  
                                                function(xx, x){
                                                  sd(scores[ scores$run == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                },
                                                x = x ) 
                                      }),
                        data_set = sapply( unique(scores$dat) , 
                                           function(x){ 
                                             sapply( metrics,  
                                                     function(xx, x){
                                                       sd(scores[ scores$dat == x & scores$Eval.metric == xx,eval_data, drop=T], na.rm=T)
                                                     },
                                                     x = x ) 
                                           }) )
  
  ### merge data to fit with ggplot2 formalism
  ggdat <- merge( 
            data.frame(name = colnames(models_mean), t(models_mean) ),
            data.frame(name = colnames(models_sd), t(models_sd) ),
            by = "name" )
  
  colnames(ggdat) <- c("name", "mean1", "mean2", "sd1", "sd2")
  
  #   ggdat <- data.frame( 
  #             merge( 
  #               reshape2::melt(models_mean), 
  #               reshape2::melt(models_sd),
  #               by = c("X1","X2") ) )
  #   
  #   colnames(ggdat) <- c("metric","name","mean","sd")

  ## produce plots
  gg <- ggplot(ggdat, aes_string(x="mean1", y="mean2", colour="name", fill = NULL))
  
  ### add mean poins
  gg <- gg + geom_point()  

  
  ### add axis names and remove legend name
  gg <- gg + xlab(metrics[1]) + ylab(metrics[2]) + theme(legend.title=element_blank())
  
  ### fix scale
  if( length(ylim) | length(xlim)){
    gg <- gg + coord_cartesian(ylim=ylim, xlim=xlim)
  }
  
  ### add title
  if(length(main)){
    gg <- gg + labs(title=main)
  }
  
  ### add error bars
  limits1 <- aes_string(xmax = "mean1 + sd1", xmin= "mean1 - sd1", fill = NULL)
  limits2 <- aes_string(ymax = "mean2 + sd2", ymin= "mean2 - sd2", fill = NULL)
  gg <- gg + geom_errorbar(limits2,width=0) + geom_errorbarh(limits1, height=0)

  if(plot){
    print(gg)
  }
  
  invisible(gg)
} ## end of models_scores_graph function

.models_scores_graph_check_args <- function(obj, metrics = NULL, by = 'models', args){
  ## check obj type
  if(! ( class(obj) %in% c("BIOMOD.models.out", "BIOMOD.EnsembleModeling.out") ) ){
    stop("obj must be a 'BIOMOD.models.out' or a 'BIOMOD.EnsembleModeling.out' object")
  }
  
  ## check metrics
  scores <- get_evaluations(obj, as.data.frame = T )
  
  avail_metrics <- unique( scores$Eval.metric )
  if( length(avail_metrics)<2 ){
    stop("at least 2 different evaluations metrics needed")
  }
  if (is.null(metrics)){
    metrics <- avail_metrics[1:2]
    warnings(toString(metrics), " evaluation metrics automatically selected")
  }
  
  ## check by
  if(! (by %in% c('models', 'algos', 'cv_run', 'data_set') ) ){
    stop("by arg should be one of 'models', 'algos', 'cv_run' or 'data_set'")
  }
  
  ## check extra args
  test_args_names <- ! ( names(args) %in% c("xlim", "ylim", "main", "col"))
  if( any ( test_args_names ) ){
    stop("unknown ", toString( names(args)[ test_args_names ] )," args")
  }
  
  
  xlim <- args$xlim
  ylim <- args$ylim
  main <- args$main
  
  return(list(metrics = metrics,
              xlim = xlim,
              ylim = ylim,
              main = main))
  
} ## end of checking args function

# x11()
# models_scores_graph(myBiomodModelOut, by = 'cv_run', metrics=c("ROC","TSS"))
# models_scores_graph(myBiomodEM, by = 'algos', metrics=c("ROC","TSS"))
# 
# str(gg)
# get_evaluations(myBiomodModelOut, as.data.frame=T)
