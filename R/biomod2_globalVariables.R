
## biomod2_internal ------------
utils::globalVariables(names = c("res"))

## biomod2_classes_5 ------------
utils::globalVariables(names = c("cv"))


## get_formal_data ------------
utils::globalVariables(names = c("i"))


## BIOMOD_Tuning ------------
utils::globalVariables(names = c("i", "quant", "type", "IA"))

## BIOMOD_Projection ------------
utils::globalVariables(names = c("do.stack", "output.format", "mod.name", "on_0_1000"
                                 , "omit.na", "keep.in.memory"))

## BIOMOD_EnsembleModeling ------------
utils::globalVariables(names = c("metric.select.user"))

## .get_needed_predictions ------------
utils::globalVariables(names = c("thisPA"))

## BIOMOD_EnsembleForecasting ------------
utils::globalVariables(names = c("on_0_1000", "output.format", "keep.in.memory"))

## BIOMOD_LoadModels ------------
utils::globalVariables(names = c("full.name", "models", "run.eval", "data.set"))

## BIOMOD_PresenceOnly ------------
utils::globalVariables(names = c("ind.1"))


## bm_RunModel ------------
utils::globalVariables(names = c("expl_var_names", "resp_name", "criteria"))
## bm_RunModelsLoop ------------
utils::globalVariables(names = c("modi"))

## bm_VariablesImportance ------------
utils::globalVariables(names = c("temp_workdir", "variables", "v", "r"))


## bm_PlotEvalMean ------------
utils::globalVariables(names = c("xlim", "ylim", "main"))

## bm_PlotEvalBoxplot ------------
utils::globalVariables(names = c("scales", "main"))

## bm_PlotVarImpBoxplot ------------
utils::globalVariables(names = c("main"))

## bm_PlotResponseCurves ------------
utils::globalVariables(names = c("data_species", "models", "vari", "nb.pts", "model"
                                 , "use.formal.names", "on_0_1000", "main"))

## bm_PlotRangeSize ------------
utils::globalVariables(names = c("vali"))

