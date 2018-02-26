'CustomIndexMaker' <- function(){
  # 1. test if modif has ever been done
#   return(FALSE)
  
  if(file.exists(paste(system.file("",package='biomod2'), 
                 .Platform$file.sep, "HasBeenCustom.txt", sep=""))){
    return(FALSE)
  } else{
    cat("\nbiomod2 first load...")
    
    if(file.exists(file.path(system.file("doc",package='biomod2'), "html","00Index.html"))){
      cat("\ncustom help index setting up...")
      # get currentindex files
      old.index <- system.file("html/00Index.html",package='biomod2')
      
      # get customed index file
      new.index <- system.file("doc/html/00Index.html",package='biomod2')
      
      # update version number
      pkg.version <- read.dcf(file=system.file("DESCRIPTION", package='biomod2'),
               fields="Version")
      
      indexTmp <- readLines(new.index)
      indexTmp[grepl("version ",indexTmp)] <- paste("version ", pkg.version, sep="")
      cat(indexTmp, sep="\n", file=new.index, append=FALSE)
      
      # and replace it
      file.copy(from=old.index,to=paste(tools::file_path_sans_ext(old.index),"_DEFAULT.html", sep=""))
      file.copy(from=new.index,to=old.index,overwrite=TRUE)
    }
    
#     if(!file.exists(file.path(system.file("",package='biomod2'), "Meta", "vignette.rds"))){
#       cat("\ncustom vignettes files setting up...")
#       # get vignettes files
#       vignettes_R <- list.files(path=system.file("doc",package='biomod2'),pattern=".R",ignore.case=TRUE)
#       vignettes_PDF <- list.files(path=system.file("doc",package='biomod2'),pattern=".pdf",ignore.case=TRUE)
#       
#       # create the .Rnw & vignette.rds file
#       vignettes_Rnw <- union(gsub(".R","",vignettes_R), gsub(".pdf","",vignettes_PDF))
#       vignettes_Rnw <- paste(vignettes_Rnw,".Rnw", sep="")
#       invisible(file.create(file.path(system.file("doc",package='biomod2'), vignettes_Rnw), showWarnings=F))
#       
#       vignettes_rds <- file.path(system.file("doc",package='biomod2'),"html","vignettes_rds.csv")
#       if(file.exists(vignettes_rds)){
#         vignettes_rds <- read.csv(vignettes_rds, colClasses = "character")
#         saveRDS(vignettes_rds,file=file.path(system.file("Meta",package='biomod2'),"vignette.rds"))
#       }
#     }
    # create a file that will indicate that package customing has ever been done
    invisible(file.create(file.path(system.file("",package='biomod2'), "HasBeenCustom.txt"), showWarnings=F))
    return(TRUE)
  }
}
  
  


.onAttach <- function(libname, pkgname) {
  if(file.exists(system.file("DESCRIPTION", package=pkgname))){
    RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    
    if(CustomIndexMaker()){
      packageStartupMessage("Customed index built!")
    }
    
    packageStartupMessage(paste(pkgname, RFver, "loaded.\n\nType browseVignettes(package='biomod2') to access directly biomod2 vignettes."))
  }
}