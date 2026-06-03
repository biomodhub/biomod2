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
    verify_records <- idig_search_records(rq = list(genus = c("manis",
                                                   "rhinolophus",
                                                   "paguma")),
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
# columns returned in your results
records <- idig_search_records(rq = list(genus = c("manis",
                                                   "rhinolophus",
                                                   "paguma")),
                       fields = c("uuid",
                                  "recordset",
                                  "institutioncode",
                                  "genus",
                                  "scientificname",
                                  "country",
                                  "data.dwc:year",
                                  "data.dwc:collectionCode",
                                  "catalognumber",
                                  "data.dwc:preparations"))

## ----eval=verify_records, echo = FALSE, results = 'asis'----------------------
knitr::kable(head(records)) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(width = "100%")

## ----eval=verify_records------------------------------------------------------
# List distinct values for the `preparation` field
prepsummary <- records %>% 
  group_by(`data.dwc:preparations`) %>% 
  tally()

# Display `prepsummary` in HTML output
knitr::kable(prepsummary) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE) %>% 
  scroll_box(height = "400px")

## ----eval=verify_records------------------------------------------------------
# Normalize values in `data.dwc:preparations` to be all lowercase; then
# filter rows that include our search terms
recordsfiltered <- records %>% 
  mutate(`data.dwc:preparations` = str_to_lower(`data.dwc:preparations`)) %>% 
  filter(grepl('freeze|froze|tissue', `data.dwc:preparations`))

## ----eval=verify_records, echo = FALSE, results = 'asis'----------------------
knitr::kable(recordsfiltered) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(height = "600px")

# If you have this code open in R, you can uncomment the line below to
# save `recordsfiltered` as a csv file to your working directory
# write_csv(recordsfiltered, "recordsfiltered.csv")

## ----eval=verify_records, echo = FALSE, results = 'asis'----------------------
# List distinct values for the `preparation` field in recordsfiltered
recordsfiltered %>% 
  group_by(`data.dwc:preparations`) %>% 
  tally() %>% 
  knitr::kable() %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE) %>% 
  scroll_box(height = "400px")

## ----eval=verify_records------------------------------------------------------
# Count how many records in the data were contributed by each recordset
recordtally <- recordsfiltered %>% 
  group_by(recordset) %>% 
  tally() %>% rename()

# Get metadata from the attributes of the `records` data frame
collections <- tibble(collection = attr(recordsfiltered, "attribution")) %>% 
  # Expand information captured in nested lists
  hoist(collection, 
        recordset_uuid = "uuid",
        recordset_name = "name",
        recordset_url= "url",
        contacts = "contacts") %>% 
  # Get rid of extraneous attribution metadata
  select(-collection) %>% 
  # Expand information captured in nested lists
  unnest_longer(contacts) %>% 
  # Expand information captured in nested lists
  unnest_wider(contacts) %>% 
  # Remove any contacts without an email address listed
  filter(!is.na(email)) %>% 
  # Get rid of duplicate contacts within the same recordset
  distinct() %>% 
  # Rename some columns
  rename(contact_role = role, contact_email = email) %>% 
  # Group first and last names together in the same column
  unite(col = "contact_name", 
        first_name, last_name, 
        sep = " ", 
        na.rm = TRUE) %>% 
  # Restructure data frame so that there is one row per recordset
  group_by(recordset_uuid) %>% 
  mutate(contact_index = row_number()) %>%
  mutate(recordset_url = if_else(grepl("^http://", recordset_url),
    gsub("^http://", "https://", recordset_url),
    recordset_url
  )) %>%
  mutate(recordset_url = if_else(grepl("www.amnh.org/our-research/vertebrate-zoology/mammalogy", recordset_url),
    gsub("www.amnh.org/our-research/vertebrate-zoology/mammalogy", "www.amnh.org/research/vertebrate-zoology/mammalogy", recordset_url),
    recordset_url
  )) %>%
  mutate(recordset_url = if_else(grepl("www.burkemuseum.org/mammalogy", recordset_url),
    gsub("www.burkemuseum.org/mammalogy", "www.burkemuseum.org/collections-and-research/biology/mammalogy", recordset_url),
    recordset_url
  )) %>%
  mutate(recordset_url = if_else(grepl("www.nhm.org", recordset_url),
    gsub("www.nhm.org", "nhm.org", recordset_url),
    recordset_url
  )) %>%
  mutate(recordset_url = if_else(
    grepl("appl003.lsu.edu/natsci/lmns.nsf/\\$Content/Mammals\\?OpenDocument", recordset_url),
    gsub("appl003.lsu.edu/natsci/lmns.nsf/\\$Content/Mammals\\?OpenDocument", "appl103.lsu.edu/natsci/Collections/natscicolsearch.nsf/OpenMainPage?OpenAgent&ID=1042", recordset_url),
    recordset_url
  )) %>%
  mutate(recordset_url = if_else(grepl("www.nsrl.ttu.edu/collections/Mammals/index.htm", recordset_url),
    gsub("www.nsrl.ttu.edu/collections/Mammals/index.htm", "www.depts.ttu.edu/nsrl/collections/mammal.php", recordset_url),
    recordset_url
  )) %>%
  mutate(recordset_url = if_else(grepl("sites01.lsu.edu/wp/mns/research-collections/genetic-resources/", recordset_url),
    gsub("sites01.lsu.edu/wp/mns/research-collections/genetic-resources/", "appl103.lsu.edu/natsci/Collections/natscicolsearch.nsf/OpenMainPage?OpenAgent&ID=1050", recordset_url),
    recordset_url
  )) %>%
  mutate(recordset_url = if_else(grepl("https://www.msb.unm.edu", recordset_url),
    gsub("www.msb.unm.edu", "www.msb.unm.edu", recordset_url),
    recordset_url
  )) %>%
  pivot_wider(names_from = contact_index,
                values_from = c(contact_name, contact_role, contact_email)) %>%
  # Include how many records in the data were contributed by each recordset
  left_join(recordtally, by = c("recordset_uuid"="recordset")) %>% 
   # Filter and remove n = 0
  filter(!is.na(n)) %>% 
  # Get rid of any rows which don't actually contribute data to `records`;
  # necessary because the attribute metadata by default includes all recordsets
  # in iDigBio that match the `idig_search_records` query, even if you filter
  # or limit those results in your own code
  filter(recordset_uuid %in% records$recordset) 
  


## ----eval=verify_records, echo = FALSE, results = 'asis'----------------------
knitr::kable(collections) %>% 
    kable_styling(bootstrap_options = 
                         c("striped", "hover", "condensed", "responsive")) %>% 
  scroll_box(height = "400px")

