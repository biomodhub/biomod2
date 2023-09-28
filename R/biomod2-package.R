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
    mess <- paste(mess, "/!\\ New set up for modeling options. We apologize for the trouble ^[*.*]^")
    packageStartupMessage(mess)
    
    toLoad <- unique(ModelsTable$package[-which(ModelsTable$package %in% c("MAXENT", "biomod2"))])
    for (pkg in toLoad) eval(parse(text = paste0("require(", pkg, ")")))
  }
}


