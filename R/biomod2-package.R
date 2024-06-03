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
    mess <- paste(mess, "/!\\ Please note that some BigBoss options have been changed between biomod2 v4.2-6 and previous versions.")
    packageStartupMessage(mess)
    
    toLoad <- unique(ModelsTable$package[-which(ModelsTable$package %in% c("MAXENT", "biomod2"))])
    for (pkg in toLoad) eval(parse(text = paste0("require(", pkg, ")")))
  }
}


