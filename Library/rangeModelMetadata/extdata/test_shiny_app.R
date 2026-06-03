
# NOTE: if we're using the .Rproj file to load up, we don't need explicit paths
# # write test files to the inst dir
# if(Sys.info()['user']=='ctg') myInstDir='/Users/ctg/Dropbox/Projects/Range_Metadata/rangeModelMetadata/inst/extdata'
# if(Sys.info()['user']=='Brian') myInstDir=NULL
# if(Sys.info()['user']=='hannah') myInstDir=NULL
# if(Sys.info()['user']=='musasabi')  myInstDir=NULL

# rmm1 object 1
rmm1=rmmTemplate()
rmm1=rmmAutofillPackageCitation(rmm1,c('raster','sp'))
# rmm1AutoFillData(rmm1,species=)
raster.files=list.files(system.file("extdata/Env_Demo",package='rangeModelMetadata'),full.names = TRUE)
env=raster::stack(raster.files)
rmm1=rmmAutofillEnvironment(rmm1,env,transfer=0) # for fitting environment
rmm1=rmmAutofillEnvironment(rmm1,env,transfer=1) # for transfer environment 1 (assuming different than for fitting)
rmm1=rmmAutofillEnvironment(rmm1,env,transfer=2) # for transfer environment 2 (assuming different than for fitting)
# save as rdata (saveRDS works if you want to load the object back in after)
saveRDS(rmm1, file='inst/extdata/shiny_test_rmm1.rds')
# rmm1 <- readRDS('inst/extdata/shiny_test_rmm1.rds')
# save as csv
rmmToCSV(rmm1,file='inst/extdata/shiny_test_rmm1.csv')

# rmm1 object 2 (to test that they're the same)
rmm2=rmm1
# save as rdata
saveRDS(rmm2, file='inst/extdata/shiny_test_rmm2.rds')
# save as csv
rmmToCSV(rmm2,file='inst/extdata/shiny_test_rmm2.csv')

# alternate output type
rmms=list(rmm1,rmm2)
saveRDS(rmms, file='inst/extdata/shiny_test_rmms.rds')

# rmmCheckFinalize(rmm1)
# rmmCheckShiny()
