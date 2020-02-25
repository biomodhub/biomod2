##' @name BIOMOD_ModelingOptions
##' @aliases BIOMOD_ModelingOptions
##' @title Configure the modeling options for each selected model
##'
##' @description Parametrize and/or tune biomod's single models options.
##'
##' @usage
##'   BIOMOD_ModelingOptions(GLM = NULL,
##'                          GBM = NULL,
##'                          GAM = NULL,
##'                          CTA = NULL,
##'                          ANN = NULL,
##'                          SRE = NULL,
##'                          FDA = NULL,
##'                          MARS = NULL,
##'                          RF = NULL,
##'                          MAXENT.Phillips = NULL)
##'
##' @param GLM  list, GLM options
##' @param GBM  list, GBM options
##' @param GAM  list, GAM options
##' @param CTA  list, CTA options
##' @param ANN  list, ANN options
##' @param SRE  list, SRE options
##' @param FDA  list, FDA options
##' @param MARS list, MARS options
##' @param RF list, RF options
##' @param MAXENT.Phillips  list, MAXENT.Phillips options
##'
##'
##' @details
##'   The aim of this function is to allow advanced user to change some default parameters of BIOMOD inner models.
##'   For each modeling technique, options can be set up.
##'
##'   Each argument have to be put in a list object.
##'
##'   The best way to use this function is to print defaut models options (\code{\link{Print_Default_ModelingOptions}}) or create a default 'BIOMOD.model.option object' and print it in your console. Then copy the output, change only the required parameters, and paste it as function arguments. (see example)
##'
##'   Here the detailed list of modifiable parameters. They correspond to the traditional parameters that could be setted out for each modeling technique (e.g. ?GLM)
##'
##' @section GLM (\code{\link[stats]{glm}}):
##'
##'   \itemize{
##'
##'     \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off.
##'       You can choose to either:
##'         \itemize{
##'           \item{generate automatically the GLM formula by using the type and interaction.level arguments
##'             type (default \code{'quadratic'}) : formula given to the model ('simple', 'quadratic' or 'polynomial').
##'             interaction.level (default \code{0}) : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GLM.}
##'           \item{or construct specific formula}
##'         }}
##'
##'     \item{\code{test} (default \code{'AIC'}) : Information criteria for the stepwise selection procedure: AIC for Akaike Information Criteria, and BIC for Bayesian Information Criteria ('AIC' or 'BIC'). 'none' is also a supported value which implies to concider only the full model (no stepwise selection). This can lead to convergence issu and strange results.}
##'
##'     \item{\code{family} (default \code{binomial(link = 'logit')}) : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \link{family} for details of family functions.) . BIOMOD only runs on presence-absence data so far, so binomial family by default.}
##'
##'     \item{\code{control} : a list of parameters for controlling the fitting process. For glm.fit this is passed to \code{\link{glm.control}}.}
##'
##'   }
##'
##' @section GBM (default \code{\link[gbm]{gbm}}):
##'
##'   Please refer to \code{\link[gbm]{gbm}} help file to get the meaning of this options.
##'   \itemize{
##'     \item{ \code{distribution} (default \code{'bernoulli'})}
##'     \item{ \code{n.trees} (default \code{2500})}
##'     \item{ \code{interaction.depth} (default \code{7})}
##'     \item{ \code{n.minobsinnode} (default \code{5})}
##'     \item{ \code{shrinkage} (default \code{0.001})}
##'     \item{ \code{bag.fraction} (default \code{0.5})}
##'     \item{ \code{train.fraction} (default \code{1})}
##'     \item{ \code{cv.folds} (default \code{3})}
##'     \item{ \code{keep.data} (default \code{FALSE})}
##'     \item{ \code{verbose} (default \code{FALSE})}
##'     \item{ \code{perf.method} (default \code{'cv'})}
##'     \item{ \code{n.cores} (default \code{1})}
##'   }
##'
##'
##' @section GAM (\code{\link[gam]{gam}} or \code{\link[mgcv]{gam}}):
##'   \itemize{
##'
##'     \item{algo : either "GAM_gam" (default), "GAM_mgcv" or "BAM_mgcv" defining the chosen GAM function (see \code{\link[mgcv]{gam}}, \code{\link[gam]{gam}} resp. \code{\link[mgcv]{bam}} for more details)}
##'
##'     \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off.
##'       You can choose to either:
##'         \itemize{
##'           \item{generate automatically the GAM formula by using the type and interaction.level arguments
##'             type : the smother used to generate the formula. Only "s_smoother" available at time.
##'             interaction.level : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GAM. Interaction are not considered if you choosed "GAM_gam" algo}
##'           \item{or construct specific formula}
##'         }}
##'
##'     \item{k (default \code{-1} or \code{4}): a smooth term in a formula argument to gam (see \pkg{gam} \code{\link[gam]{s}} or \pkg{mgcv} \code{\link[mgcv]{s}})}
##'     \item{family (default \code{binomial(link = 'logit')}) : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \link{family} for details of family functions.) . BIOMOD only runs on presence-absence data so far, so binomial family by default. }
##'     \item{control : see \code{\link[mgcv]{gam.control}} or \code{\link[gam]{gam.control}}}
##'
##'     \item{some extra "GAM_mgcv" specific options (ignored if algo = "GAM_gam")
##'       \itemize{
##'         \item{\code{method} (default \code{'GCV.Cp'})}
##'         \item{\code{optimizer} (default \code{c('outer','newton')})}
##'         \item{\code{select} (default \code{FALSE})}
##'         \item{\code{knots} (default \code{NULL})}
##'         \item{\code{paramPen} (default \code{NULL})}
##'       }
##'
##'     }
##'   }
##'
##'
##' @section CTA (\code{\link[rpart]{rpart}}):
##'
##'   Please refer to \code{\link[rpart]{rpart}} help file to get the meaning of the following options.
##'   \itemize{
##'     \item{\code{method} (default \code{'class'})}
##'     \item{\code{parms} (default \code{'default'}) : if \code{'default'}, default \pkg{rpart} parms value are kept}
##'     \item{\code{cost} (default \code{NULL})}
##'     \item{\code{control}: see \code{\link[rpart]{rpart.control}}}
##'   }
##'
##'   NOTE: for method and parms, you can give a 'real' value as described in the rpart help file or 'default' that implies default \code{\link[rpart]{rpart}} values.
##'
##' @section ANN (\code{\link[nnet]{nnet}}):
##'
##'   \itemize{
##'     \item{\code{NbCV} (default \code{5}) : nb of cross validation to find best size and decay parameters}
##'     \item{\code{size}} (default \code{NULL}) : number of units in the hidden layer. If \code{NULL} then size parameter will be optimised by cross validation based on model AUC (\code{NbCv} cross validation; tested size will be the following c(2,4,6, 8) ). You can also specified a vector of size you want to test. The one giving the best model AUC will be then selected.
##'     \item{\code{decay}} (default \code{NULL}) : parameter for weight decay. If \code{NULL} then decay parameter will be optimised by cross validation on model AUC (\code{NbCv} cross validation; tested decay will be the following c(0.001, 0.01, 0.05, 0.1) ). You can also specified a vector of decay you want to test. The one giving the best model AUC will be then selected.
##'     \item{\code{rang} (default \code{0.1}) : Initial random weights on [-rang, rang]}
##'     \item{\code{maxit} (default \code{200}): maximum number of iterations.}
##'   }
##'
##'
##' @section SRE (\code{\link[biomod2]{sre}}):
##'   \itemize{
##'     \item{\code{quant} (default \code{0.025}): quantile of 'extreme environmental variable' removed for selection of species envelops}
##'   }
##'
##'
##' @section FDA (\code{\link[mda]{fda}}):
##'
##'   Please refer to \code{\link[mda]{fda}} help file to get the meaning of these options.
##'   \itemize{
##'     \item{\code{method} (default \code{'mars'})}
##'     \item{\code{add_args} (default \code{NULL}) : additional arguments to \code{method} given as a list of parameters
##'       (corespond to the \ldots options of fda function) }
##'   }
##'
##' @section MARS (\code{\link[earth]{earth}}):
##'
##'   Please refer to \code{\link[earth]{earth}} help file to get the meaning of these options.
##'   \itemize{
##'     \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off.
##'       You can choose to either:
##'         \itemize{
##'           \item{generate automatically the GLM formula by using the type and interaction.level arguments
##'             type (default \code{'simple'}) : formula given to the model ('simple', 'quadratic' or 'polynomial').
##'             interaction.level (default \code{0}) : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GLM/MARS.}
##'           \item{or construct specific formula}
##'         }}
##'     %    \item{\code{degree} (default \code{2})}
##'     \item{\code{nk}} (default \code{NULL}) : an optional integer specifying the maximum number of model terms. If NULL is given then default mars function value is used ( i.e max(21, 2 * nb_expl_var + 1) )
##'     \item{\code{penalty} (default \code{2})}
##'     \item{\code{thresh} (default \code{0.001})}
##'     \item{\code{nprune} (default \code{NULL})}
##'     \item{\code{pmethod} (default \code{"backward"})}
##'   }
##'
##'
##' @section RF (\code{\link[randomForest]{randomForest}}):
##'
##'   \itemize{
##'     \item{\code{do.classif} (default \code{TRUE}) : if TRUE classification random.forest computed else regression random.forest will be done}
##'     \item{\code{ntree} (default \code{500})}
##'     \item{\code{mtry} (default \code{'default'})}
##'     \item{\code{nodesize} (default \code{5})}
##'     \item{\code{maxnodes} (default \code{NULL})}
##'   }
##'
##'   NOTE: for mtry, you can give a 'real' value as described in randomForest help file or 'default' that implies default randomForest values
##'
##' @section  [MAXENT.Phillips](https://biodiversityinformatics.amnh.org/open_source/maxent/) :
##'   \itemize{
##'     \item{\code{path_to_maxent.jar} : character, the link to \pkg{maxent.jar} file (the working directory by default) }
##'     \item{\code{memory_allocated} : integer (default \code{512}), the amount of memory (in Mo) reserved for java to run MAXENT.Phillips. should be 64, 128, 256, 512, 1024, 2048... or NULL if you want to use default java memory limitation parameter.}
##'     \item{\code{background_data_dir} : character, path to a directory where explanatory variables are stored as ASCII files (raster format).
##'       If specified MAXENT will generate it's own background data from expalantory variables rasters (as usually done in MAXENT studies). If not
##'       set, then MAXENT will use the same pseudo absences than other models (generated within biomod2 at formatting step) as background data.}
##'     \item{\code{maximumbackground} : integer, the maximum number of background data to sample. This parameter will be use only if \code{background_data_dir}
##'       option has been set to a non default value.}
##'     \item{\code{maximumiterations} : integer (default \code{200}), maximum iteration done}
##'     \item{\code{visible} : logical (default \code{FALSE}), make the Maxent user interface visible}
##'     \item{\code{linear} : logical (default \code{TRUE}), allow linear features to be used}
##'     \item{\code{quadratic} : logical (default \code{TRUE}), allow quadratic features to be used}
##'     \item{\code{product} : logical (default \code{TRUE}), allow product features to be used}
##'     \item{\code{threshold} : logical (default \code{TRUE}), allow threshold features to be used}
##'     \item{\code{hinge} : logical (default \code{TRUE}), allow hinge features to be used}
##'     \item{\code{lq2lqptthreshold} : integer (default \code{80}), number of samples at which product and threshold features start being used}
##'     \item{\code{l2lqthreshold} : integer (default \code{10}), number of samples at which quadratic features start being used}
##'     \item{\code{hingethreshold} : integer (default \code{15}), number of samples at which hinge features start being used}
##'     \item{\code{beta_threshold} : numeric (default \code{-1.0}), regularization parameter to be applied to all threshold features; negative value enables automatic setting}
##'     \item{\code{beta_categorical} : numeric (default \code{-1.0}), regularization parameter to be applied to all categorical features; negative value enables automatic setting}
##'     \item{\code{beta_lqp} : numeric (default \code{-1.0}), regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting}
##'     \item{\code{beta_hinge} : numeric (default \code{-1.0}), regularization parameter to be applied to all hinge features; negative value enables automatic setting}
##'     \item{\code{betamultiplier} : numeric (default \code{1}), multiply all automatic regularization parameters by this number. A higher number gives a more spread-out distribution.}
##'     \item{\code{defaultprevalence} : numeric (default \code{0.5}), default prevalence of the species: probability of presence at ordinary occurrence points}
##'   }
##'
##' % @section MAXENT.Tsuruoka (\code{\link[maxent]{maxent}}):
##' %
##' % \itemize{
##' %   \item{\code{l1_regularizer} (default \code{0.0}): An numeric turning on L1 regularization and setting the regularization parameter. A value of 0 will disable L1 regularization}
##' %   \item{\code{l2_regularizer} (default \code{0.0}): An numeric turning on L2 regularization and setting the regularization parameter. A value of 0 will disable L2 regularization}
##' %   \item{\code{use_sgd} (default \code{FALSE}): A logical indicating that SGD parameter estimation should be used. Defaults to FALSE}
##' %   \item{\code{set_heldout} (default \code{0}): An integer specifying the number of documents to hold out. Sets a held-out subset of your data to test against and prevent overfitting}
##' %   \item{\code{verbose} (default \code{FALSE}): A logical specifying whether to provide descriptive output about the training process}
##' % }
##'
##' % NOTE: if you use the \code{set_heldout} parameter then the data that will be held out will be taken in the
##' % calibration data pool. It can be penilizing in case of low number of occurences dataset.
##'
##'
##' @return
##'   A \code{"\link[=BIOMOD.Model.Options-class]{BIOMOD.Model.Options}"} object given to \code{\link[biomod2]{BIOMOD_Modeling}}
##'
##' @author Damien Georges, Wilfried Thuiller
##' @keywords models
##' @keywords options
##'
##' @examples
##'   ## default BIOMOD.model.option object
##'   myBiomodOptions <- BIOMOD_ModelingOptions()
##'
##'   ## print the object
##'   myBiomodOptions
##'
##'   ## you can copy a part of the print, change it and custom your options
##'   ## here we want to compute quadratic GLM and select best model with 'BIC' criterium
##'   myBiomodOptions <- BIOMOD_ModelingOptions(
##'     GLM = list( type = 'quadratic',
##'                 interaction.level = 0,
##'                 myFormula = NULL,
##'                 test = 'BIC',
##'                 family = 'binomial',
##'                 control = glm.control(epsilon = 1e-08,
##'                                       maxit = 1000,
##'                                       trace = FALSE) ))
##'
##'   ## check changes was done
##'   myBiomodOptions
##'
##'   ##' you can prefer to establish your own GLM formula
##'   myBiomodOptions <- BIOMOD_ModelingOptions(
##'     GLM = list( myFormula = formula("Sp277 ~ bio3 +
##'                     log(bio10) + poly(bio16,2) + bio19 + bio3:bio19")))
##'
##'   ## check changes was done
##'   myBiomodOptions
##'
##'   ##' you also can directly print default parameters and then follow the same processus
##'   Print_Default_ModelingOptions()
##'


####################################################################################################
'BIOMOD_ModelingOptions' <- function(
                        GLM = NULL,
                        GBM = NULL,
                        GAM = NULL,
                        CTA = NULL,
                        ANN = NULL,
                        SRE = NULL,
                        FDA = NULL,
                        MARS = NULL,
                        RF = NULL,
                        MAXENT.Phillips = NULL
                        ){
  # 1. create a defaut BIOMOD.Model.Options object
  opt <- new('BIOMOD.Model.Options')

  # 2. modify it if necessary
  if(!is.null(GLM)){
    if(!is.null(GLM$type)) { opt@GLM$type <- GLM$type }
    if(!is.null(GLM$interaction.level)) { opt@GLM$interaction.level <- GLM$interaction.level }
    if(!is.null(GLM$myFormula)) { opt@GLM$myFormula <- GLM$myFormula }
    if(!is.null(GLM$test)) { opt@GLM$test <- GLM$test }
    if(!is.null(GLM$family)) {
      fam.test <- TRUE
      if(inherits(GLM$family, 'family')){
        opt@GLM$family <- GLM$family
      } else{
        if( is.character(GLM$family)){
          if(! unlist(strsplit(GLM$family,"[/(]"))[1] %in% c('binomial', 'gaussian', 'Gamma', 'inverse.gaussian', 'poisson', 'quasi', 'quasibinomial', 'quasipoisson')){ fam.test <- FALSE}

          if(grepl(')', GLM$family)){ # check string formalisation to add () if necessary
            opt@GLM$family <- eval(parse(text=GLM$family))
          } else{
            opt@GLM$family <- eval(parse(text=paste(GLM$family,"()", sep="")))
          }
        } else{ fam.test <- FALSE }
      }
      if(!fam.test){
        cat("\n!!! invalid GLM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt@GLM$family <- binomial(link = 'logit')
      }
    }
    if(!is.null(GLM$mustart)) { opt@GLM$mustart <- GLM$mustart }
    if(!is.null(GLM$control)) { opt@GLM$control <- GLM$control }
  }

  if(!is.null(GBM)){
#     if(!is.null(GBM$type )) { opt@GBM$type <- GBM$type }
#     if(!is.null(GBM$interaction.level )) { opt@GBM$interaction.level <- GBM$interaction.level }
    if(!is.null(GBM$distribution )) { opt@GBM$distribution <- GBM$distribution }
    if(!is.null(GBM$n.trees )) { opt@GBM$n.trees <- GBM$n.trees }
    if(!is.null(GBM$interaction.depth )) { opt@GBM$interaction.depth <- GBM$interaction.depth }
    if(!is.null(GBM$n.minobsinnode )) { opt@GBM$n.minobsinnode <- GBM$n.minobsinnode }
    if(!is.null(GBM$shrinkage )) { opt@GBM$shrinkage <- GBM$shrinkage }
    if(!is.null(GBM$bag.fraction )) { opt@GBM$bag.fraction <- GBM$bag.fraction }
    if(!is.null(GBM$train.fraction )) { opt@GBM$train.fraction <- GBM$train.fraction }
    if(!is.null(GBM$cv.folds )) { opt@GBM$cv.folds <- GBM$cv.folds }
    if(!is.null(GBM$keep.data )) { opt@GBM$keep.data <- GBM$keep.data }
    if(!is.null(GBM$verbose )) { opt@GBM$verbose <- GBM$verbose }
#     if(!is.null(GBM$class.stratify.cv )) { opt@GBM$class.stratify.cv <- GBM$cv.folds }
    if(!is.null(GBM$perf.method )) { opt@GBM$perf.method <- GBM$perf.method }
    if(!is.null(GBM$n.cores)) { opt@GBM$n.cores <- GBM$n.cores } else { opt@GBM$n.cores <- NULL }
  }



  if(!is.null(GAM)){
    if(!is.null(GAM$algo )) { opt@GAM$algo <- GAM$algo }
    if(!is.null(GAM$type )) { opt@GAM$type <- GAM$type }
    if(!is.null(GAM$k )) { opt@GAM$k <- GAM$k } else{
      if(opt@GAM$algo == 'GAM_gam'){
        opt@GAM$k <- 4
      } else{
        opt@GAM$k <- -1
      }
    }
    if(!is.null(GAM$interaction.level )) { opt@GAM$interaction.level <- GAM$interaction.level }
    if(!is.null(GAM$myFormula )) { opt@GAM$myFormula <- GAM$myFormula }
    if(!is.null(GAM$family)) {
      fam.test <- TRUE
      if(inherits(GAM$family, 'family')){
        opt@GAM$family <- GAM$family
      } else{
        if( is.character(GAM$family)){
          if(! unlist(strsplit(GAM$family,"[/(]"))[1] %in% c('binomial', 'gaussian', 'Gamma', 'inverse.gaussian', 'poisson', 'quasi', 'quasibinomial', 'quasipoisson')){ fam.test <- FALSE}

          if(grepl(')', GAM$family)){ # check string formalisation to add () if necessary
            opt@GAM$family <- eval(parse(text=GAM$family))
          } else{
            opt@GAM$family <- eval(parse(text=paste(GAM$family,"()", sep="")))
          }
        } else{ fam.test <- FALSE }
      }
      if(!fam.test){
        cat("\n!!! invalid GAM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt@GAM$family <- binomial(link = 'logit')
      }
    }

    if(is.null(GAM$control )) {
      if(opt@GAM$algo == 'GAM_gam'){
        requireNamespace('gam', quietly = TRUE)
        opt@GAM$control <- gam::gam.control()
      } else{ opt@GAM$control <- mgcv::gam.control() }
    } else{
      user.control.list <- GAM$control
      if(opt@GAM$algo == 'GAM_gam'){
        default.control.list <- gam::gam.control()
      } else{
        default.control.list <- mgcv::gam.control()
      }

      control.list <- lapply(names(default.control.list), function(x){
        if(x %in% names(user.control.list)){
          return(user.control.list[[x]])
        } else {
          return(default.control.list[[x]])
        }
      })

      names(control.list) <- names(default.control.list)
      opt@GAM$control <- control.list
    }

    if(!is.null(GAM$method )) { opt@GAM$method <- GAM$method }
    if(!is.null(GAM$optimizer )) { opt@GAM$optimizer <- GAM$optimizer }
    if(!is.null(GAM$select )) { opt@GAM$select <- GAM$select }
    if(!is.null(GAM$knots )) { opt@GAM$knots <- GAM$knots }
    if(!is.null(GAM$paraPen )) { opt@GAM$paraPen <- GAM$paraPen }
  } else{
    if(opt@GAM$algo == 'GAM_gam'){
      opt@GAM$control <- gam::gam.control()
      opt@GAM$k <- 4
    } else{
      opt@GAM$control <- mgcv::gam.control()
      opt@GAM$k <- -1
    }
  }


  if(!is.null(CTA)){
#     if(!is.null(CTA$type )) { opt@CTA$type <- CTA$type }
#     if(!is.null(CTA$interaction.level )) { opt@CTA$interaction.level <- CTA$interaction.level }
    if(!is.null(CTA$method )) { opt@CTA$method <- CTA$method }
    if(!is.null(CTA$parms )) { opt@CTA$parms <- CTA$parms }
    if(!is.null(CTA$control )) { opt@CTA$control <- CTA$control }
    if(!is.null(CTA$cost )) { opt@CTA$cost <- CTA$cost }
  }

  if(!is.null(ANN)){
#     if(!is.null(ANN$type )) { opt@ANN$type <- ANN$type }
#     if(!is.null(ANN$interaction.level )) { opt@ANN$interaction.level <- ANN$interaction.level }
    if(!is.null(ANN$NbCV )) { opt@ANN$NbCV <- ANN$NbCV }
    if(!is.null(ANN$size )) { opt@ANN$size <- ANN$size }
    if(!is.null(ANN$decay )) { opt@ANN$decay <- ANN$decay }
    if(!is.null(ANN$rang )) { opt@ANN$rang <- ANN$rang }
    if(!is.null(ANN$maxit )) { opt@ANN$maxit <- ANN$maxit }
  }

  if(!is.null(SRE)){
    if(!is.null(SRE$quant )) { opt@SRE$quant <- SRE$quant }
  }

  if(!is.null(FDA)){
#     if(!is.null(FDA$type )) { opt@FDA$type <- FDA$type }
#     if(!is.null(FDA$interaction.level )) { opt@FDA$interaction.level <- FDA$interaction.level }
    if(!is.null(FDA$method )) { opt@FDA$method <- FDA$method }
    if(!is.null(FDA$add_args )) { opt@FDA$add_args <- FDA$add_args } ## additional args such as degree, nk
  }

  if(!is.null(MARS)){
    if(!is.null(MARS$type)) { opt@MARS$type <- MARS$type }
    if(!is.null(MARS$interaction.level)) { opt@MARS$interaction.level <- MARS$interaction.level }
    if(!is.null(MARS$myFormula)) { opt@MARS$myFormula <- MARS$myFormula }
#     if(!is.null(MARS$degree )) { opt@MARS$degree <- MARS$degree }
    if(!is.null(MARS$nk )) { opt@MARS$nk <- MARS$nk }
    if(!is.null(MARS$penalty )) { opt@MARS$penalty <- MARS$penalty }
    if(!is.null(MARS$thresh )) { opt@MARS$thresh <- MARS$thresh }
    if(!is.null(MARS$nprune )) { opt@MARS$nprune <- MARS$nprune }
    if(!is.null(MARS$pmethod )) { opt@MARS$pmethod <- MARS$pmethod }
  }

  if(!is.null(RF)){
    if(!is.null(RF$type )) { opt@RF$type <- RF$type }
#     if(!is.null(RF$interaction.level )) { opt@RF$interaction.level <- RF$interaction.level }
#     if(!is.null(RF$do.classif )) { opt@RF$do.classif <- RF$do.classif }
    if(!is.null(RF$ntree )) { opt@RF$ntree <- RF$ntree }
    if(!is.null(RF$mtry )) { opt@RF$mtry <- RF$mtry }
    if(!is.null(RF$nodesize )) { opt@RF$nodesize <- RF$nodesize }
    if(!is.null(RF$maxnodes )) { opt@RF$maxnodes <- RF$maxnodes }
  }

  if(!is.null(MAXENT.Phillips)){
    if(!is.null(MAXENT.Phillips$path_to_maxent.jar )) {
      opt@MAXENT.Phillips$path_to_maxent.jar <- normalizePath(sub("maxent.jar", "", MAXENT.Phillips$path_to_maxent.jar)) # ensure path format validity
      } else {opt@MAXENT.Phillips$path_to_maxent.jar <- getwd()}
    if(!is.null(MAXENT.Phillips$memory_allocated )) { opt@MAXENT.Phillips$memory_allocated <- MAXENT.Phillips$memory_allocated }
	if(!is.null(MAXENT.Phillips$background_data_dir )) { opt@MAXENT.Phillips$background_data_dir <- MAXENT.Phillips$background_data_dir }
    if(!is.null(MAXENT.Phillips$maximumbackground )) { opt@MAXENT.Phillips$maximumbackground <- MAXENT.Phillips$maximumbackground }
    if(!is.null(MAXENT.Phillips$maximumiterations )) { opt@MAXENT.Phillips$maximumiterations <- MAXENT.Phillips$maximumiterations }
    if(!is.null(MAXENT.Phillips$visible )) { opt@MAXENT.Phillips$visible <- MAXENT.Phillips$visible }
    if(!is.null(MAXENT.Phillips$linear )) { opt@MAXENT.Phillips$linear <- MAXENT.Phillips$linear }
    if(!is.null(MAXENT.Phillips$quadratic )) { opt@MAXENT.Phillips$quadratic <- MAXENT.Phillips$quadratic }
    if(!is.null(MAXENT.Phillips$product )) { opt@MAXENT.Phillips$product <- MAXENT.Phillips$product }
    if(!is.null(MAXENT.Phillips$threshold )) { opt@MAXENT.Phillips$threshold <- MAXENT.Phillips$threshold }
    if(!is.null(MAXENT.Phillips$hinge )) { opt@MAXENT.Phillips$hinge <- MAXENT.Phillips$hinge }
    if(!is.null(MAXENT.Phillips$lq2lqptthreshold )) { opt@MAXENT.Phillips$lq2lqptthreshold <- MAXENT.Phillips$lq2lqptthreshold }
    if(!is.null(MAXENT.Phillips$l2lqthreshold )) { opt@MAXENT.Phillips$l2lqthreshold <- MAXENT.Phillips$l2lqthreshold }
    if(!is.null(MAXENT.Phillips$hingethreshold )) { opt@MAXENT.Phillips$hingethreshold <- MAXENT.Phillips$hingethreshold }
    if(!is.null(MAXENT.Phillips$beta_threshold )) { opt@MAXENT.Phillips$beta_threshold <- MAXENT.Phillips$beta_threshold }
    if(!is.null(MAXENT.Phillips$beta_categorical )) { opt@MAXENT.Phillips$beta_categorical <- MAXENT.Phillips$beta_categorical }
    if(!is.null(MAXENT.Phillips$beta_lqp )) { opt@MAXENT.Phillips$beta_lqp <- MAXENT.Phillips$beta_lqp }
    if(!is.null(MAXENT.Phillips$beta_hinge )) { opt@MAXENT.Phillips$beta_hinge <- MAXENT.Phillips$beta_hinge }
	  if(!is.null(MAXENT.Phillips$betamultiplier )) { opt@MAXENT.Phillips$betamultiplier <- MAXENT.Phillips$betamultiplier }
    if(!is.null(MAXENT.Phillips$defaultprevalence )) { opt@MAXENT.Phillips$defaultprevalence <- MAXENT.Phillips$defaultprevalence }
  } else{
    opt@MAXENT.Phillips$path_to_maxent.jar <- getwd()
  }

  # if(!is.null(MAXENT.Tsuruoka)){
  #   if(!is.null(MAXENT.Tsuruoka$l1_regularizer )) { opt@MAXENT.Tsuruoka$l1_regularizer <- MAXENT.Tsuruoka$l1_regularizer }
  #   if(!is.null(MAXENT.Tsuruoka$l2_regularizer )) { opt@MAXENT.Tsuruoka$l2_regularizer <- MAXENT.Tsuruoka$l2_regularizer }
  #   if(!is.null(MAXENT.Tsuruoka$use_sgd )) { opt@MAXENT.Tsuruoka$use_sgd <- MAXENT.Tsuruoka$use_sgd }
  #   if(!is.null(MAXENT.Tsuruoka$set_heldout )) { opt@MAXENT.Tsuruoka$set_heldout <- MAXENT.Tsuruoka$set_heldout }
  #   if(!is.null(MAXENT.Tsuruoka$verbose )) { opt@MAXENT.Tsuruoka$verbose <- MAXENT.Tsuruoka$verbose }
  # }

  test <- as.logical(validObject(object = opt, test = TRUE, complete = FALSE))

  if(!test){
    cat("\n\n!!! NULL object returned because of invalid parameters given !!!")
    return(NULL)
  }

  return(opt)
}

Print_Default_ModelingOptions <- function(){
  cat('\n Defaut modeling options. copy, change what you want paste it as arg to BIOMOD_ModelingOptions\n\n')

  opt_tmp <- BIOMOD_ModelingOptions()
  print(opt_tmp)
}
