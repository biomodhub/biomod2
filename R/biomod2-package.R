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
    mess <- paste(mess, "/!\\ Update 4.2-0 now relies on terra version >= 1.6.33. Make sure to update with `devtools::install_github('rspatial/terra')`. We apologize for the trouble >{o.o}<")
    packageStartupMessage(mess)
  }
}


