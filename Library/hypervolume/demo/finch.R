if (exists('doHypervolumeFinchDemo')==TRUE)
{
  message('Demo structure and results have changed due to improvements between version 1.4.x and version 2.x of the package.\nTo replicate results seen in our 2014 GEB paper please use an older version of the package.')
  
  
  data(morphSnodgrassHeller)
  finch_isabela <- morphSnodgrassHeller[morphSnodgrassHeller$IslandID=="Isa_Alb",]

  # select trait axes
  trait_axes <- c("BodyL","WingL","TailL","BeakW")
  traitdata <- finch_isabela[,c("TaxonOrig",trait_axes)]
  # keep complete cases only
  traitdata <- na.omit(traitdata)
  
  species_list = as.character(unique(traitdata$TaxonOrig))
  num_species = length(species_list)  

  # compute hypervolumes for each species  
  hv_finches_list = new("HypervolumeList")
  hv_finches_list@HVList = vector(mode="list",length=num_species)
  for (i in 1:num_species)
  {
    # keep the trait data 
    data_this_species = traitdata[traitdata$TaxonOrig==species_list[i],trait_axes]
    # log-transform to rescale
    data_this_species_log <- log10(data_this_species)
    
    # make a hypervolume using auto-bandwidth
      hv_finches_list@HVList[[i]] <- hypervolume_gaussian(data_this_species_log,
                                          name=as.character(species_list[i]),verbose=FALSE)
  }
  
  # compute all pairwise overlaps
  overlap = matrix(NA, nrow=num_species, ncol=num_species)
  dimnames(overlap)=list(species_list, species_list)
  for (i in 1:num_species)
  {
    for (j in i:num_species)
    {
      if (i!=j)
      {
        # compute set operations on each pair
        this_set = hypervolume_set(hv_finches_list@HVList[[i]], hv_finches_list@HVList[[j]], check.memory=FALSE)
        # calculate a Sorensen overlap index (2 x shared volume / sum of |hv1| + |hv2|)
        overlap[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
      }
    }   
  }
  

  
  # show all hypervolumes
  plot(hv_finches_list)
  
  # show pairwise overlaps - note that actually very few species overlap in four dimensions
  op <- par(mar=c(10,10,1,1))
  image(x=1:nrow(overlap), y=1:nrow(overlap), z=overlap,axes=FALSE,xlab='',ylab='',col=colorRampPalette(c("lightgray","red"))(100))
  box()
  axis(side=1, at=1:(length(dimnames(overlap)[[1]])),dimnames(overlap)[[1]],las=2,cex.axis=0.75)
  axis(side=2, at=1:(length(dimnames(overlap)[[2]])),dimnames(overlap)[[2]],las=1,cex.axis=0.75)
  par(op)
  
  # reset to original state
  rm(doHypervolumeFinchDemo)
} else
{
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo requires approximately 3 minutes to run.')  
  message('To run the demo, type')
  message('\tdoHypervolumeFinchDemo=TRUE')
  message('\tdemo(finch)')
  message('at the R command line prompt.')
}