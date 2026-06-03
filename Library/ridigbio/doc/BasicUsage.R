## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("ridigbio")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(ridigbio)

## ----echo = FALSE-------------------------------------------------------------

verify_galax_records <- FALSE

#Test that examples will run
tryCatch({
    # Your code that might throw an error
    verify_galax_records <- idig_search_records(rq=list(scientificname="Galax urceolata"),
      limit = 10
    )
}, error = function(e) {
    # Code to run if an error occurs
    cat("An error occurred during the idig_search_records call: ", e$message, "\n")
    cat("Vignettes will not be fully generated. Please try again after resolving the issue.")
    # Optionally, you can return NULL or an empty dataframe
    verify_galax_records <- FALSE
})

## ----eval=verify_galax_records------------------------------------------------
galax_records <- idig_search_records(rq=list(scientificname="Galax urceolata"))

## ----eval=verify_galax_records------------------------------------------------
colnames(galax_records)

## ----eval=verify_galax_records------------------------------------------------
diapensiaceae_records <- idig_search_records(rq=list(family="Diapensiaceae"), limit=1000)

## ----eval=verify_galax_records------------------------------------------------
rq_input <- list("scientificname"=list("type"="exists"),
                 "family"="Diapensiaceae", 
                 geopoint=list(
                   type="geo_bounding_box",
                   top_left=list(lon = -98.16, lat = 48.92),
                   bottom_right=list(lon = -64.02, lat = 23.06)
                   )
                 )

## ----eval=verify_galax_records------------------------------------------------
diapensiaceae_records_USA <- idig_search_records(rq_input, limit=1000)

## ----eval=verify_galax_records------------------------------------------------
galax_media <- idig_search_media(rq=list(scientificname="Galax urceolata"))

## ----eval=verify_galax_records------------------------------------------------
colnames(galax_media)

## ----eval=verify_galax_records------------------------------------------------
galax_media2 <- idig_search_media(rq=list(scientificname="Galax urceolata"),
                                  mq=list("data.ac:accessURI"=list("type"="exists")))

