
## biomod2_internal ------------
utils::globalVariables(names = c("i.dim1",
                                 "i.dim2",
                                 "i.dim3",
                                 "sub.i"))

## biomod2_classes_0 ------------
utils::globalVariables(names = c("expected_CVnames", "OptionsBigboss"))

## biomod2_classes_1 ------------
utils::globalVariables(names = c("this_PA",
                                 "this_run",
                                 "this_dataset",
                                 "has.mask",
                                 "has.mask.eval",
                                 "y"))

## biomod2_classes_3 ------------
utils::globalVariables(names = c("i"))

## biomod2_classes_4 ------------
utils::globalVariables(names = c("thislayername"))

## biomod2_classes_5 ------------
utils::globalVariables(names = c("cv"))


## get_formal_data ------------
utils::globalVariables(names = c("i"))


## bm_Tuning ------------
utils::globalVariables(names = c("calib.i", "rep", "quant"
                                 , "i", "typ", "intlev", "fi"))

## BIOMOD_Projection ------------
utils::globalVariables(names = c("do.stack",
                                 "output.format",
                                 "mod.name",
                                 "on_0_1000",
                                 "omit.na",
                                 "keep.in.memory"))

                                 

## BIOMOD_Modeling ------------
utils::globalVariables(names = c("resp", "value", "pa"))

## BIOMOD_EnsembleModeling ------------
utils::globalVariables(names = c("eval.m",
                                 "assemb",
                                 "algo",
                                 "xx",
                                 "em.algo.long",
                                 "em.algo.class",
                                 "em.mod.assemb",
                                 "metric.select.user"))

## .get_needed_predictions ------------
utils::globalVariables(names = c("this_PA"))

## BIOMOD_EnsembleForecasting ------------
utils::globalVariables(names = c("on_0_1000",
                                 "output.format",
                                 "keep.in.memory",
                                 "em.name"))

## BIOMOD_LoadModels ------------
utils::globalVariables(names = c("full.name", "models", "run", "PA"))

## BIOMOD_RangeSize ------------
utils::globalVariables(names = c("thiscol", "pred", "proj"))


## bm_PseudoAbsences ------------
utils::globalVariables(names = c("i.abs"))

## bm_Tuning ------------
utils::globalVariables(names = c("dataset.i", "PA.i", "tuned.mod", "train.params"
                                 , "tuning.grid", "criteria.AIC"))

## bm_ModelingOptions ------------
utils::globalVariables(names = c("ModelsTable"))

## bm_CrossValidation ------------
utils::globalVariables(names = c("pa", "env", "this.colnames"))

## bm_FindOptimStat ------------
utils::globalVariables(names = c("expl.var.names"))

## bm_RunModel ------------
utils::globalVariables(names = c("expl_var_names", "resp_name", "criteria"
                                 , "weights", "data_env", "data_sp", "data_xy"))
## bm_RunModelsLoop ------------
utils::globalVariables(names = c("modi", "xx", "ii"))

## bm_VariablesImportance ------------
utils::globalVariables(names = c("temp_workdir", "variables", "v", "r"))



## bm_PlotEvalMean ------------
utils::globalVariables(names = c("xlim", "ylim", "main"))

## bm_PlotEvalBoxplot ------------
utils::globalVariables(names = c("scales", "main"))

## bm_PlotVarImpBoxplot ------------
utils::globalVariables(names = c("main"))

## bm_PlotResponseCurves ------------
utils::globalVariables(names = c("data_species",
                                 "models",
                                 "nb.pts",
                                 "use.formal.names",
                                 "on_0_1000",
                                 "main", 
                                 "vari",
                                 "vari2", 
                                 "model", 
                                 "combi", 
                                 "dat_"))


## bm_PlotRangeSize ------------
utils::globalVariables(names = c("vali"))

