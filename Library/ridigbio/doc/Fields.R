## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(ridigbio)

# Load library for making nice HTML output
library(kableExtra)

## ----echo = FALSE-------------------------------------------------------------

verify_record_fields <- FALSE

#Test that examples will run
tryCatch({
    # Your code that might throw an error
    verify_record_fields <- idig_meta_fields(type = "records", subset = "raw", limit = 10)
    rfall <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(rfall) <- c("Type", "FieldNameValue")

    # Stricter schema verification
    for(i in 1:length(verify_record_fields)) {
        rf <- data.frame(
            Type = verify_record_fields[[i]]$type,
            FieldNameValue = verify_record_fields[[i]]$fieldName,
            stringsAsFactors = FALSE
        )
        rfall <- rbind(rfall, rf)
    }

    if(nrow(rfall) <= 0)
        verify_record_fields <- FALSE

}, error = function(e) {
    # Code to run if an error occurs
    cat("An error occurred during the idig_search_records call: ", e$message, "\n")
    cat("Vignettes will not be fully generated. Please try again after resolving the issue.")
    # Optionally, you can return NULL or an empty dataframe
    verify_record_fields <- FALSE
})

## ----eval=verify_record_fields------------------------------------------------
record_fields <- idig_meta_fields(type = "records", subset = "raw")
rfall <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(rfall) <- c("Type", "FieldNameValue")

for(i in 1:length(record_fields)) {
    rf <- data.frame(
        Type = record_fields[[i]]$type,
        FieldNameValue = record_fields[[i]]$fieldName,
        stringsAsFactors = FALSE
    )
    rfall <- rbind(rfall, rf)
}
colnames(rfall) <- c("type", "fieldName")
nrow(rfall)

## ----eval=verify_record_fields------------------------------------------------
record_fields_index <- idig_meta_fields(type = "records", subset = "indexed")
rfalli <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(rfalli) <- c("Type", "FieldNameValue")

for(i in 1:length(record_fields_index)) {
    rf <- data.frame(
        Type = record_fields_index[[i]]$type,
        FieldNameValue = record_fields_index[[i]]$fieldName,
        stringsAsFactors = FALSE
    )
    rfalli <- rbind(rfalli, rf)
}
colnames(rfalli) <- c("type", "fieldName")
nrow(rfalli)

## ----eval=verify_record_fields------------------------------------------------
setequal(rfall, rfalli)

## ----eval=verify_record_fields------------------------------------------------
rfall[171,]

## ----eval=verify_record_fields------------------------------------------------
rfalli[69,]

## ----eval=verify_record_fields------------------------------------------------
out <- idig_search_records(rq=list(scientificname="Galax urceolata"), 
                           fields = rfall$fieldName)

## ----eval=verify_record_fields------------------------------------------------
media_fields <- idig_meta_fields(type = "media", subset = "raw")
mfall <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(mfall) <- c("Type", "FieldNameValue")

for(i in 1:length(media_fields)) {
    mf <- data.frame(
        Type = media_fields[[i]]$type,
        FieldNameValue = media_fields[[i]]$fieldName,
        stringsAsFactors = FALSE
    )
    mfall <- rbind(mfall, mf)
}
colnames(mfall) <- c("type", "fieldName")
nrow(mfall)

## ----eval=verify_record_fields------------------------------------------------
media_fields_indexed <- idig_meta_fields(type = "media", subset = "indexed")
mfalli <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(mfalli) <- c("Type", "FieldNameValue")

for(i in 1:length(media_fields_indexed)) {
    mf <- data.frame(
        Type = media_fields_indexed[[i]]$type,
        FieldNameValue = media_fields_indexed[[i]]$fieldName,
        stringsAsFactors = FALSE
    )
    mfalli <- rbind(mfalli, mf)
}
colnames(mfalli) <- c("type", "fieldName")
nrow(mfalli)

## ----eval=verify_record_fields------------------------------------------------
mfalli$fieldName

## ----eval=verify_record_fields------------------------------------------------
out <- idig_search_media(rq=list(scientificname="Galax urceolata"), fields = mfalli$fieldName) 

out %>%
kable() %>% 
kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                font_size = 12,
                fixed_thead = T) %>% 
  scroll_box(width = "100%", height = "400px")

