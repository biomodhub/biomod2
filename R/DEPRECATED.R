#' @title Deprecated functions in package \pkg{biomod2}.
#' @description The functions listed below are deprecated and were removed
#' in the current version. When possible, alternative functions with similar
#'   functionality are mentioned. Help pages for deprecated functions are
#'   available at \code{help("<function>-deprecated")}.
#' @name biomod2-deprecated
#' @keywords internal
#' 
#'
NULL

#' @title BIOMOD_ConvertOldRun
#' @description Deprecated function used in \pkg{biomod2} 3.5.1 to convert
#' results from older version into the current version
#' @param ... Additional arguments
#' 
#' @name BIOMOD_ConvertOldRun-deprecated
#' @usage BIOMOD_ConvertOldRun(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{BIOMOD_ConvertOldRun}:
#' \code{BIOMOD_ConvertOldRun} was removed in \pkg{biomod2} version 4.0.
#'
#' @export
BIOMOD_ConvertOldRun <- function(...) { .Deprecated() }

#' @title BIOMOD_cv
#' @description Deprecated function name for
#'   \code{\link{BIOMOD_CrossValidation}}
#' @param ... Additional arguments
#'
#' @name BIOMOD_cv-deprecated
#' @usage BIOMOD_cv(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{BIOMOD_cv}:
#' For \code{BIOMOD_cv} use \code{\link{BIOMOD_CrossValidation}}.
#'
#' @export
BIOMOD_cv <- function(...) { .Deprecated("BIOMOD_CrossValidation") }

#' @title BIOMOD_presenceonly
#' @description Deprecated function name for
#'   \code{\link{BIOMOD_PresenceOnly}}
#' @param ... Additional arguments
#'
#' @name BIOMOD_presenceonly-deprecated
#' @usage BIOMOD_presenceonly(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{BIOMOD_presenceonly}:
#' For \code{BIOMOD_presenceonly} use \code{\link{BIOMOD_PresenceOnly}}.
#'
#' @export
BIOMOD_presenceonly <- function(...) { .Deprecated("BIOMOD_PresenceOnly") }


#' @title BIOMOD_tuning
#' @description Deprecated function name for
#'   \code{\link{BIOMOD_Tuning}}
#' @param ... Additional arguments
#'
#' @name BIOMOD_tuning-deprecated
#' @usage BIOMOD_tuning(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{BIOMOD_tuning}:
#' For \code{BIOMOD_tuning} use \code{\link{BIOMOD_Tuning}}.
#'
#' @export
BIOMOD_tuning <- function(...) { .Deprecated("BIOMOD_Tuning") }

#' @title BinaryTransformation
#' @description Deprecated function name for
#'   \code{\link{bm_BinaryTransformation}}
#' @param ... Additional arguments
#'
#' @name BinaryTransformation-deprecated
#' @usage BinaryTransformation(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{BinaryTransformation}:
#' For \code{BinaryTransformation} use \code{\link{bm_BinaryTransformation}}.
#'
#' @export
BinaryTransformation <- function(...) { .Deprecated("bm_BinaryTransformation") }


#' @title calculate.stat
#' @description Deprecated function name for
#'   \code{\link{bm_CalculateStat}}
#' @param ... Additional arguments
#'
#' @name calculate.stat-deprecated
#' @usage calculate.stat(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{calculate.stat}:
#' For \code{calculate.stat} use \code{\link{bm_CalculateStat}}.
#'
#' @export
calculate.stat <- function(...) { .Deprecated("bm_CalculateStat") }


#' @title CustomIndexMaker
#' @description Deprecated function used in \pkg{biomod2} 3.5.1 to replace
#'   default html index file by a custom one if defined
#' @param ... Additional arguments
#'
#' @name CustomIndexMaker-deprecated
#' @usage CustomIndexMaker(...)
#' @seealso \code{\link{CustomIndexMaker-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{CustomIndexMaker}:
#' \code{CustomIndexMaker} was removed in \pkg{biomod2} version 4.0.
#'
#' @export
CustomIndexMaker <- function(...) { .Deprecated() }



#' @title FilteringTransformation
#' @description Deprecated function name for
#'   \code{\link{BinaryTransformation}}
#' @param ... Additional arguments
#'
#' @name FilteringTransformation-deprecated
#' @usage FilteringTransformation(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{FilteringTransformation}:
#' For \code{FilteringTransformation} use \code{\link{BinaryTransformation}}.
#'
#' @export
FilteringTransformation <- function(...) { .Deprecated("BinaryTransformation") }


#' @title Find.Optim.Stat
#' @description Deprecated function name for
#'   \code{\link{bm_FindOptimStat}}
#' @param ... Additional arguments
#'
#' @name Find.Optim.Stat-deprecated
#' @usage Find.Optim.Stat(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{Find.Optim.Stat}:
#' For \code{Find.Optim.Stat} use \code{\link{bm_FindOptimStat}}.
#'
#' @export
Find.Optim.Stat <- function(...) { .Deprecated("bm_FindOptimStat") }

#' @title getStatOptimValue
#' @description Deprecated function name for
#'   \code{\link{get_optim_value}}
#' @param ... Additional arguments
#'
#' @name getStatOptimValue-deprecated
#' @usage getStatOptimValue(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{getStatOptimValue}:
#' For \code{getStatOptimValue} use \code{\link{get_optim_value}}.
#'
#' @export
getStatOptimValue <- function(...) { .Deprecated("get_optim_value") }

#' @title makeFormula
#' @description Deprecated function name for
#'   \code{\link{bm_MakeFormula}}
#' @param ... Additional arguments
#'
#' @name makeFormula-deprecated
#' @usage makeFormula(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{makeFormula}:
#' For \code{makeFormula} use \code{\link{bm_MakeFormula}}.
#'
#' @export
makeFormula <- function(...) { .Deprecated("bm_MakeFormula") }

#' @title models_scores_graph
#' @description Deprecated function name for
#'   \code{\link{bm_PlotEvalMean}}
#' @param ... Additional arguments
#'
#' @name models_scores_graph-deprecated
#' @usage models_scores_graph(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{models_scores_graph}:
#' For \code{models_scores_graph} use \code{\link{bm_PlotEvalMean}}.
#'
#' @export
models_scores_graph <- function(...) { .Deprecated("bm_PlotEvalMean") }

#' @title Print_Default_ModelingOptions
#' @description Deprecated function name for
#'   \code{\link{bm_DefaultModelingOptions}}
#' @param ... Additional arguments
#'
#' @name Print_Default_ModelingOptions-deprecated
#' @usage Print_Default_ModelingOptions(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{Print_Default_ModelingOptions}:
#' For \code{Print_Default_ModelingOptions} use \code{\link{bm_DefaultModelingOptions}}.
#'
#' @export
Print_Default_ModelingOptions <- function(...) { .Deprecated("bm_DefaultModelingOptions") }

#' @title ProbDensFunc
#' @description Deprecated function name for
#'   \code{\link{bm_PlotRangeSize}}
#' @param ... Additional arguments
#'
#' @name ProbDensFunc-deprecated
#' @usage ProbDensFunc(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{ProbDensFunc}:
#' For \code{ProbDensFunc} use \code{\link{bm_PlotRangeSize}}.
#'
#' @export
ProbDensFunc <- function(...) { .Deprecated("bm_PlotRangeSize") }


#' @title response.plot
#' @description Deprecated function name for
#'   \code{\link{bm_PlotResponseCurves}}
#' @param ... Additional arguments
#'
#' @name response.plot-deprecated
#' @usage response.plot(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{response.plot}:
#' For \code{response.plot} use \code{\link{bm_PlotResponseCurves}}.
#'
#' @export
response.plot <- function(...) { .Deprecated("bm_PlotResponseCurves") }

#' @title response.plot2
#' @description Deprecated function name for
#'   \code{\link{bm_PlotResponseCurves}}
#' @param ... Additional arguments
#'
#' @name response.plot2-deprecated
#' @usage response.plot2(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{response.plot2}:
#' For \code{response.plot2} use \code{\link{bm_PlotResponseCurves}}.
#'
#' @export
response.plot2 <- function(...) { .Deprecated("bm_PlotResponseCurves") }

#' @title sample.factor.levels
#' @description Deprecated function name for
#'   \code{\link{bm_SampleFactorLevels}}
#' @param ... Additional arguments
#'
#' @name sample.factor.levels-deprecated
#' @usage sample.factor.levels(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{sample.factor.levels}:
#' For \code{sample.factor.levels} use \code{\link{bm_SampleFactorLevels}}.
#'
#' @export
sample.factor.levels <- function(...) { .Deprecated("bm_SampleFactorLevels") }

#' @title sre
#' @description Deprecated function name for
#'   \code{\link{bm_SRE}}
#' @param ... Additional arguments
#'
#' @name sre-deprecated
#' @usage sre(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{sre}:
#' For \code{sre} use \code{\link{bm_SRE}}.
#'
#' @export
sre <- function(...) { .Deprecated("bm_SRE") }

#' @title .transform.outputs.array
#' @description Deprecated function used in \pkg{biomod2} 3.5.1 to transform
#' outputs into array
#' @param ... Additional arguments
#'
#' @name .transform.outputs.array-deprecated
#' @usage .transform.outputs.array(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{.transform.outputs.array}:
#' \code{.transform.outputs.array} was removed in \pkg{biomod2} version 4.0.
#'
#' @export
.transform.outputs.array <- function(...) { .Deprecated() }


#' @title .transform.outputs.list
#' @description Deprecated function name for
#'   \code{.transform_outputs_list}
#' @param ... Additional arguments
#'
#' @name .transform.outputs.list-deprecated
#' @usage .transform.outputs.list(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{.transform.outputs.list}:
#' For \code{.transform.outputs.list} use \code{.transform_outputs_list}.
#'
#' @export
.transform.outputs.list <- function(...) { .Deprecated(".transform_outputs_list") }

#' @title update_objects
#' @description Deprecated function used in \pkg{biomod2} 3.5.1 to update
#'   objects construct with a old version of biomod2 to a current one
#' @param ... Additional arguments
#'
#' @name update_objects-deprecated
#' @usage update_objects(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{update_objects}:
#' \code{update_objects} was removed in \pkg{biomod2} version 4.0.
#'
#' @export
update_objects <- function(...) { .Deprecated() }


#' @title variables_importance
#' @description Deprecated function name for
#'   \code{\link{bm_VariablesImportance}}
#' @param ... Additional arguments
#'
#' @name variables_importance-deprecated
#' @usage variables_importance(...)
#' @seealso \code{\link{biomod2-deprecated}}
#' @keywords internal
NULL
#' @rdname biomod2-deprecated
#' @section \code{variables_importance}:
#' For \code{variables_importance} use \code{\link{bm_VariablesImportance}}.
#'
#' @export
variables_importance <- function(...) { .Deprecated("bm_VariablesImportance") }


# DF_to_ARRAY
# evaluate
# full_suffling
# level.plot
# multiple.plot
# randomise_data
# SampleMat2
# zzz_bm