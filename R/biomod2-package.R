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
    mess <- paste(mess, "/!\\ Since version 4.2 biomod2 relies on terra and may thus return SpatRaster that can easily be converted in RasterStack with `stack()`. We apologize for the trouble @('_')@")
    packageStartupMessage(mess)
  }
}


