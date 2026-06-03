if (exists('doHypervolumeHolesFinchDemo')==TRUE)
{
  message('Demo structure and results have changed due to improvements between version 1.4.x and version 2.x of the package.\nTo replicate results seen in our 2016 American Naturalist paper please use an older version of the package.')
  
  data(morphSnodgrassHeller)
  
  # select data for only Isabela Island
  finch_isabela <- morphSnodgrassHeller[morphSnodgrassHeller$IslandID=="Isa_Alb",]
  
  # select trait axes
  trait_axes <- c("BodyL","WingL","TailL","BeakW")
  trait_data <- finch_isabela[,trait_axes]
  # keep complete cases only
  trait_data <- na.omit(trait_data)
  # convert length units to comparable scales
  trait_data_scaled <- scale(trait_data)
  
  # the below line is commented out because of approximate two-hour runtime
  # bw_plugin <- estimate_bandwidth(trait_data_scaled, method="plug-in")
  # instead, use pre-computed output for that function
  bw_plugin <- c(0.5278828, 0.4883812, 0.5951435, 0.4480163)
  
  # get overall community hypervolume
  hv_finch <- hypervolume_gaussian(trait_data_scaled,kde.bandwidth=bw_plugin, quantile.requested=0.5, quantile.requested.type="volume")
  hv_finch@Name <- "Finches"
  
  # compute convex expectation
  # first thin the hypervolume
  hv_finch_thinned = hypervolume_thin(hv_finch, num.points=500)
  
  ec_finch <- expectation_convex(hv_finch_thinned, check.memory=FALSE, use.random=TRUE)
  ec_finch@Name <- "Convex expectation"
  # find holes
  holes_finch <- hypervolume_holes(hv_finch, ec_finch, set.check.memory=FALSE)
  holes_finch@Name <- "Holes"
  
  # extract volume statistics
  volumes <- get_volume(hypervolume_join(hv_finch, ec_finch, holes_finch))
  # plot volume fractions
  barplot(volumes)
  
  # calculate fraction of volume that is holey
  hole_volume_ratio <- volumes["Holes"] / volumes["Convex expectation"]
  print(hole_volume_ratio)
  # calculate approximate length of axis occupied
  print(length_ratio <- hole_volume_ratio ^ (1/4))
  
  # plot holes
  plot(hypervolume_join(hv_finch, holes_finch),
       col=c('purple','green'),
       names=c("Body length","Wing length","Tail length", "Beak width"),
       show.legend=FALSE,cex.names=1.5,contour.lwd=2)
  
  
  # calculate (in transformed coordinates) the centroid of the holes	
  holepos <- get_centroid(holes_finch)
  print(holepos)
  
  # calculate (in untransformed coordinates) the centroid of the holes
  cpos <- attr(trait_data_scaled, "scaled:center")
  csca <- attr(trait_data_scaled, "scaled:scale")
  hole_origcoords <- holepos * csca + cpos
  print(hole_origcoords)
  
  # determine which other species in the dataset is most similar to the hole
  # calculate species mean trait values
  speciesmeans <- as.data.frame(do.call("rbind",by(morphSnodgrassHeller[, trait_axes], morphSnodgrassHeller $TaxonOrig, colMeans,na.rm=TRUE)))
  # calculate rescaled distances
  scaled_diffs <- scale(speciesmeans - hole_origcoords)
  speciesmeans$dist <- apply(scaled_diffs, 1, function(x) { sqrt(sum(x^2))})
  # reorder species list by distance
  speciesmeans <- speciesmeans[order(speciesmeans $dist),]
  # identify 'top candidates' for filling the hole
  print(head(speciesmeans))
  
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires approximately 1 minutes to run.')  
  message('To run the demo, type')
  message('\tdoHypervolumeHolesFinchDemo=TRUE')
  message('\tdemo(holes_finch)')
  message('at the R command line prompt.')
}