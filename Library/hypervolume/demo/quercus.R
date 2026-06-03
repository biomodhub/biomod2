if (exists('doHypervolumeQuercusDemo')==TRUE)
{
  message('Demo structure and results have changed due to improvements between version 1.4.x and version 2.x of the package.\nTo replicate results seen in our 2014 GEB paper please use an older version of the package.')
  
  require(raster)
  require(maps)
  
  # load in lat/lon data
  data('quercus') 
  data_alba = subset(quercus, Species=="Quercus alba")[,c("Longitude","Latitude")]
  data_rubra = subset(quercus, Species=="Quercus rubra")[,c("Longitude","Latitude")]
  
  # get worldclim data from internet
  climatelayers <- getData('worldclim', var='bio', res=10, path=tempdir())
  
  # z-transform climate layers to make axes comparable
  climatelayers_ss = climatelayers[[c(1,4,12,15)]]
  for (i in 1:nlayers(climatelayers_ss))
  {
    climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd') 
  }
  climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60))
  
  # extract transformed climate values
  climate_alba = extract(climatelayers_ss_cropped, data_alba)
  climate_rubra = extract(climatelayers_ss_cropped, data_rubra)
  
  # compute hypervolumes with auto-bandwidth for both species
  hv_alba = hypervolume_gaussian(climate_alba,name='alba',samples.per.point=10)
  hv_rubra = hypervolume_gaussian(climate_rubra,name='rubra',samples.per.point=10)
  
  # determine intersection and unique components of the overlap
  hv_set = hypervolume_set(hv_alba, hv_rubra, check.memory=FALSE)
  
  # put all the output volumes in one convenient place
  volumes <- get_volume(hv_set)
  
  # do species distribution modeling (reduce point density by factor of 10 for demo speed)
  rubra_map = hypervolume_project(hv_rubra, climatelayers_ss_cropped,reduction.factor=0.1)
  alba_map = hypervolume_project(hv_alba, climatelayers_ss_cropped,reduction.factor=0.1)
  
  # then barplot of hypervolumes of each component
  op=par(mar=c(3,10,1,1))
  barplot(volumes,horiz=TRUE,las=2,main="Hypervolume",cex.names=0.5,col='lightblue')
  
  # then pairs plot of the set operations
  par(op)
  plot(hv_set[[c(3,5,6)]]) # only the unique components of each + intersection
  
  # plot the geographic projections of the ranges
  plot(rubra_map,col=colorRampPalette(c(rgb(1,1,1),rgb(1,0,0)))(100),legend=FALSE,main='Quercus rubra')
  map('world',add=TRUE)
  points(Latitude~Longitude,data=data_rubra,pch=3,cex=0.1)
  
  plot(alba_map,col=colorRampPalette(c(rgb(1,1,1),rgb(0,0,1)))(100),legend=FALSE,main='Quercus alba')
  map('world',add=TRUE)
  points(Latitude~Longitude,data=data_alba,pch=3,cex=0.1)
  
  
  rm(doHypervolumeQuercusDemo)
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires internet access to download 10MB of climate data and')
  message('will take approximately 3 minutes to run.')
  message('To run the demo, type')
  message('\tdoHypervolumeQuercusDemo=TRUE')
  message('\tdemo(quercus)')
  message('at the R command line prompt.')
}