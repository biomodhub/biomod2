## ----message=FALSE------------------------------------------------------------
# Load core libraries; install these packages if you have not already
library(ridigbio)
library(tidyverse)

# Load library for making nice HTML output
library(kableExtra)

## ----echo = FALSE-------------------------------------------------------------
verify_records_1A <- FALSE

#Test that examples will run
tryCatch({
    # Your code that might throw an error
  verify_records_1A <- idig_search_records(
  # `rq` is where you adjust your record query
  rq = list(genus = "shortia"),
  # `fields` is where you adjust what fields you want returned by the API
  fields = c("uuid",
             "family",
             "genus",
             "specificepithet",
             "scientificname",
             "stateprovince"),
  # `limit` is where you can set a limit on the number of records to return in
  # order to speed up your query; max is 100000
  limit = 10)
}, error = function(e) {
    # Code to run if an error occurs
    cat("An error occurred during the idig_search_records call: ", e$message, "\n")
    cat("Vignettes will not be fully generated. Please try again after resolving the issue.")
    # Optionally, you can return NULL or an empty dataframe
    verify_records_1A <- FALSE
})

## ----eval=verify_records_1A---------------------------------------------------
# Let's start with a simple search introducing the primary arguments for the
# function `idig_search_records`
records_1A <- idig_search_records(
  # `rq` is where you adjust your record query
  rq = list(genus = "shortia"),
  # `fields` is where you adjust what fields you want returned by the API
  fields = c("uuid",
             "family",
             "genus",
             "specificepithet",
             "scientificname",
             "stateprovince"),
  # `limit` is where you can set a limit on the number of records to return in
  # order to speed up your query; max is 100000
  limit = 10,
  # `sort` is where you can specify fields for sorting
  sort = c("stateprovince",
           "scientificname"))

# Display the data frame we just created above in a nice pretty table for HTML
knitr::kable(records_1A) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(height = "300px")

## ----eval=verify_records_1A---------------------------------------------------
# Now let's repeat the same search but remove all arguments other than `rq` to
# see what the defaults for the other arguments look like
records_1B <- idig_search_records(
  rq = list(genus = "shortia"))

records_1B$occurrenceid <- if_else(grepl("^http://", records_1B$occurrenceid),
  gsub("^http://", "", records_1B$occurrenceid),
  records_1B$occurrenceid
)

records_1B$occurrenceid <- if_else(grepl("data.biodiversitydata.nl/naturalis", records_1B$occurrenceid),
  gsub("data.biodiversitydata.nl/naturalis", "bioportal.naturalis.nl/nl", records_1B$occurrenceid),
  records_1B$occurrenceid
)

records_1B$occurrenceid <- if_else(grepl("https://grbio.org/cool", records_1B$occurrenceid),
  gsub("https://grbio.org/cool", "grbio.org/cool", records_1B$occurrenceid),
  records_1B$occurrenceid
)

records_1B$occurrenceid <- if_else(grepl("https://biocol.org", records_1B$occurrenceid),
  gsub("https://biocol.org", "biocol.org", records_1B$occurrenceid),
  records_1B$occurrenceid
)

# Display the data frame we just created above in a nice pretty table for HTML
knitr::kable(records_1B) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(height = "300px")

## ----eval=verify_records_1A---------------------------------------------------
# In the example above, we are only using one parameter in `rq` to define our
# query, but now let's search by multiple parameters
records_2A <- idig_search_records(
  rq = list(basisofrecord = "fossilspecimen",
            # Use `type = "exists"` to search for rows where there is a value
            # present in this field; the inverse of this is `type = "missing"`
            geopoint = list(type = "exists")),
  limit = 10)

# Display the data frame we just created above in a nice pretty table for HTML
knitr::kable(records_2A) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(height = "300px")

## ----eval=verify_records_1A---------------------------------------------------
# What if we wanted to see more fields than the default provides? Using the same
# search as above, we can retrieve all indexed fields with `fields = "all"`
records_2B <- idig_search_records(
  rq = list(basisofrecord = "fossilspecimen",
          geopoint = list(type="exists")),
  fields = "all",
  limit = 10)

records_2B$institutionid <- if_else(grepl("^http://", records_2B$institutionid),
  gsub("^http://", "https://", records_2B$institutionid),
  records_2B$institutionid
)

# Display the data frame we just created above in a nice pretty table for HTML
knitr::kable(records_2B) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(height = "300px")

## ----eval=verify_records_1A---------------------------------------------------
# But wait, there are even more fields available than just those we retrieved
# in the query above! Using the same search, we can choose exactly what fields
# to retrieve from indexed and raw data if we call the fields out by name in
# the `fields` argument; raw data fields are prefaced by "data.dwc:" and use 
# camelCase in their naming convention (vs. lowercase for iDigBio fields)
records_2C <- idig_search_records(
  rq = list(basisofrecord = "fossilspecimen",
          geopoint = list(type="exists")),
  # Here is where we are explicitly asking for specific fields
  fields = c("uuid",
             "recordset",
             "institutioncode", "data.dwc:institutionCode",
             "country", "data.dwc:country",
             "countrycode", "data.dwc:countryCode",
             "stateprovince", "data.dwc:stateProvince",
             "locality", "data.dwc:locality",
             "geopoint", "data.dwc:decimalLongitude", "data.dwc:decimalLatitude"),
  limit = 10)

# Display the data frame we just created above in a nice pretty table for HTML
knitr::kable(records_2C) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(height = "300px")

## ----eval=verify_records_1A---------------------------------------------------
# You may be curious what the difference is between indexed and raw data such as
# that we saw in the search above. Indexed data has been altered by iDigBio
# (often in an attempt to standardize and/or correct values), and raw data is
# what was provided to iDigBio by the data provider, i.e. the natural history
# collection. Here we will do a new search on a data quality flag to view
# differences between indexed and raw data
records_3A <- idig_search_records(
  # Data quality flags are a way for iDigBio to communicate how data was altered
  # during its quality control process, i.e. how the indexed and raw data differ
  rq = list(flags = "rev_geocode_lat_sign"),
  fields = c("uuid",
             "institutioncode", "data.dwc:institutionCode",
             "country", "data.dwc:country",
             "countrycode", "data.dwc:countryCode",
             "stateprovince", "data.dwc:stateProvince",
             "locality", "data.dwc:locality",
             "geopoint", "data.dwc:decimalLongitude", "data.dwc:decimalLatitude"),
  limit = 10)

# Let's format our results to be more readable by renaming and reordering columns
records_3A <- records_3A %>% 
  rename_at(vars(starts_with("data.dwc:")),
            ~str_replace(., "data.dwc:", "raw_")) %>% 
  select(uuid,
         indexed_decimalLatitude = geopoint.lat,
         raw_decimalLatitude,
         indexed_decimalLongitude = geopoint.lon,
         raw_decimalLongitude,
         everything())

# Display the data frame we just created above in a nice pretty table for HTML
knitr::kable(records_3A) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(height = "300px")

## ----eval=verify_records_1A---------------------------------------------------
# Let's test out a search using parameters we know would retrieve many records
count_1A <- idig_count_records(
  rq = list(basisofrecord = "fossilspecimen",
          geopoint = list(type="exists")))

# We can reformat our result to be more readable
count_1A <- format(count_1A, big.mark = ",")

# This number shows how many records in iDigBio have a value of "fossilspecimen"
# as well as geographic coordinate data
count_1A

## ----eval=verify_records_1A---------------------------------------------------
# Let's go back to our first simple search and see what the top values are for
# `scientificname` where the genus is "shortia"
top_1A <- idig_top_records(
  # `rq` is where you adjust your record query
  rq = list(genus = "shortia"),
  # `top_fields` is where you adjust what fields you want to see summarized
  top_fields = "scientificname",
  # `count` is where you can set a limit on the number of top values to return
  # in order to speed up your query; max is 1000
  count = 10)

# We need to convert our results from a nested list into a more readable format
top_1A <- as_tibble(top_1A$scientificname) %>% 
  pivot_longer(everything(), names_to = "scientificname", values_to = "count")

# Display the data frame we just created above in a nice pretty table for HTML
knitr::kable(top_1A) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(height = "300px")

