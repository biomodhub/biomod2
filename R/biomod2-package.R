### HEADER #####################################################################
##' 
##' 
## END OF HEADER ###############################################################


.onAttach <- function(libname, pkgname)
{
  if (file.exists(system.file("DESCRIPTION", package = pkgname)))
  {
    RFver <- read.dcf(file = system.file("DESCRIPTION", package = pkgname), fields = "Version")
    mess <- paste(pkgname, RFver, "loaded.\n")
    mess <- paste(mess, "/!\\ Welcome to this new version of biomod2! This is now possible to model abundance ! Take a look at the Abundance vignette on the website. \n")
    packageStartupMessage(mess)
    
    toLoad <- unique(ModelsTable$package[-which(ModelsTable$package %in% c("MAXENT", "biomod2"))])
    for (pkg in toLoad) eval(parse(text = paste0("require(", pkg, ")")))
  }
}


