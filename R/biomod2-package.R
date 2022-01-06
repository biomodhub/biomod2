### HEADER #####################################################################
##' 
##' 
## END OF HEADER ###############################################################


.onAttach <- function(libname, pkgname)
{
  if (file.exists(system.file("DESCRIPTION", package = pkgname)))
  {
    RFver <- read.dcf(file = system.file("DESCRIPTION", package = pkgname), fields = "Version")
    
    if (CustomIndexMaker()) {
      packageStartupMessage("Custom index built!")
    }
    
    mess <- paste(pkgname, RFver, "loaded.\n\n"
                  , "Type browseVignettes(package = 'biomod2') to access directly biomod2 vignettes.")
    packageStartupMessage(mess)
  }
}


.CustomIndexMaker <- function()
{
  file.custom = paste0(system.file("", package = 'biomod2'), .Platform$file.sep, "HasBeenCustom.txt")
  if (file.exists(file.custom))
  {
    return(FALSE)
  } else {
    cat("\nbiomod2 first load...")
    
    if(file.exists(file.path(system.file("doc", package = 'biomod2'), "html", "00Index.html")))
    {
      cat("\ncustom help index setting up...")
      
      # get current and custom index files
      old.index <- system.file("html/00Index.html", package = 'biomod2')
      new.index <- system.file("doc/html/00Index.html", package = 'biomod2')
      
      # update version number
      pkg.version <- read.dcf(file = system.file("DESCRIPTION", package = 'biomod2'), fields = "Version")
      
      indexTmp <- readLines(new.index)
      indexTmp[grepl("version ", indexTmp)] <- paste0("version ", pkg.version)
      cat(indexTmp, sep = "\n", file = new.index, append = FALSE)
      
      # and replace it
      file.copy(from = old.index,
                to = paste0(tools::file_path_sans_ext(old.index), "_DEFAULT.html"))
      file.copy(from = new.index,
                to = old.index,
                overwrite = TRUE)
    }
    
    # create a file that will indicate that package customing has ever been done
    invisible(file.create(file.custom, showWarnings = FALSE))
    return(TRUE)
  }
}