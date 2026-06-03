## ----message=FALSE------------------------------------------------------------
# Load core libraries; install these packages if you have not already
library(ridigbio)
library(tidyverse)

# Load library for making nice HTML output
library(kableExtra)

## ----echo = FALSE-------------------------------------------------------------
verify_records <- FALSE

#Test that examples will run
tryCatch({
    # Your code that might throw an error
    verify_records <- records <- idig_search_media(rq = list(genus = "acer", country = "united states"), 
  fields = c("uuid",
             "accessuri",
             "rights",
             "format",
             "records"),
           limit = 10)
}, error = function(e) {
    # Code to run if an error occurs
  cat("An error occurred during the idig_search_records call: ", e$message, "\n")
  cat("Vignettes will not be fully generated. Please try again after resolving the issue.")
    # Optionally, you can return NULL or an empty dataframe
  verify_records <- FALSE
})

## ----eval=verify_records------------------------------------------------------
# Edit the fields (e.g. `genus`) and values (e.g. "manis") in `list()` 
# to adjust your query and the fields (e.g. `uuid`) in `fields` to adjust the
# columns returned in your results; edit the number after `limit` to adjust the
# number of records you will retrieve images for
records <- idig_search_media(rq =
  list(genus = "acer",
       country = "united states"), 
       fields = c("uuid",
                  "accessuri",
                  "rights",
                  "format",
                  "records"),
                limit = 10)

records$accessuri <- if_else(grepl("^http://", records$accessuri),
  gsub("^http://", "", records$accessuri),
  records$accessuri
)
records$accessuri <- if_else(grepl("https://mam.ansp.org", records$accessuri),
  gsub("https://mam.ansp.org", "mam.ansp.org", records$accessuri),
  records$accessuri
)
records$accessuri <- if_else(grepl("https://ibss-images.calacademy.org", records$accessuri),
  gsub("https://ibss-images.calacademy.org", "ibss-images.calacademy.org", records$accessuri),
  records$accessuri
)

## ----eval=verify_records, echo = FALSE, results = 'asis'----------------------
knitr::kable(records) %>% 
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(width = "100%")

## ----eval=verify_records------------------------------------------------------
# Assemble a vector of iDigBio server download URLs from `records`
mediaurl_idigbio <- records %>% 
  mutate(mediaURL = paste("https://api.idigbio.org/v2/media/", uuid, sep = "")) %>% 
  select(mediaURL) %>% 
  pull()

# Assemble a vector of external server download URLs from `records`
mediaurl_external <- records$accessuri %>% 
  str_replace("\\?size=fullsize", "")

## ----eval=verify_records------------------------------------------------------
mediaurl_idigbio

## ----eval=verify_records------------------------------------------------------
mediaurl_external

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # Create new directories to save media files in
#  dir.create("jpgs_idigbio")
#  dir.create("jpgs_external")
#  
#  # Assemble another vector of file paths to use when saving media downloaded
#  # from iDigBio
#  mediapath_idigbio <- paste("jpgs_idigbio/", records$uuid, ".jpg", sep = "")
#  
#  # Assemble another vector of file paths to use when saving media downloaded
#  # from external servers; please note that it's probably not a great idea to
#  # assume these files are all jpgs, as we're doing here...
#  mediapath_external <- paste("jpgs_external/", records$uuid, ".jpg", sep = "")
#  
#  # Add a check to deal with URLs that are broken links
#  possibly_download.file = purrr::possibly(download.file,
#                                           otherwise = "cannot download")
#  
#  #"mode" argument (="wb") in the walk function to download.file.
#  
#  # Iterate through the action of downloading whatever file is at each
#  # iDigBio URL
#  purrr::walk2(.x = mediaurl_idigbio,
#               .y = mediapath_idigbio, possibly_download.file)
#  
#  # Iterate through the action of downloading whatever file is at each
#  # external URL
#  purrr::walk2(.x = mediaurl_external,
#               .y = mediapath_external, possibly_download.file)

## ----eval=verify_records, echo = FALSE, results = 'asis'----------------------
exampleimgpath <- paste("jpgs_idigbio/",records$uuid[1],".jpg", sep = "")

image_markdown <- "No image available."

if(file.exists(exampleimgpath))
{
  image_markdown <- paste0("![Herbarium specimen of an _Acer_ species collected in the United States](", exampleimgpath, ")")
}

image_markdown

