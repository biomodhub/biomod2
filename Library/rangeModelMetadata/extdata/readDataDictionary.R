library(googlesheets4)

# version 1
# d_title=read_sheet('https://docs.google.com/spreadsheets/d/1QNwT3O9Hp_NaO1RhS58qRQt1jAdWTdXuIS-TIFgZ3Mw/edit#gid=0')
#d=as.data.frame(gs_read(d_title),stringsAsFactors=FALSE)
d=read_sheet('https://docs.google.com/spreadsheets/d/1rEtDPigNWL8eXqIhiY8r2E2XiZWiiE3o9jW41rPPY7A/edit#gid=0',col_type='c')

#== Write out the dictionary to the package for use with building the metadata template
if(Sys.info()['user']=='ctg') userPath='/Users/ctg/Dropbox/Projects/Range_Metadata'
if(Sys.info()['user']=='Brian') userPath='C:/Users/Brian/Desktop/current_projects'
if(Sys.info()['user']=='hannah') userPath='/Users/ctg/Dropbox/Projects/Range_Metadata'
if(Sys.info()['user']=='musasabi') userPath='/Users/musasabi/Documents/github'

write.csv(d,paste0(userPath,'/rangeModelMetaData/inst/extdata/dataDictionary.csv'),row.names = F)

