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
    mess <- paste(mess, "/!\\ Since version 4.2 biomod2 relies on terra version >= 1.6.33. Make sure to update with `devtools::install_github('rspatial/terra')`. We apologize for the trouble >{o.o}<")
    packageStartupMessage(mess)
  }
}


