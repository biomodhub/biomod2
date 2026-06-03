## ----message=FALSE------------------------------------------------------------
# Load core libraries; install these packages if you have not already
library(ridigbio)
library(tidyverse)

# Load library for making nice HTML output
library(kableExtra)

## ----echo = FALSE-------------------------------------------------------------

verify_df_names <- FALSE

#Test that examples will run
tryCatch({
    # Your code that might throw an error
    verify_df_names <- idig_search_records(rq = list(recordset = "5082e6c8-8f5b-4bf6-a930-e3e6de7bf6fb"),
                    fields = c("uuid",
                               "data.dwc:occurrenceID",
                               "data.dwc:catalogNumber",
                               "family",
                               "data.dwc:family",
                               "genus",
                               "data.dwc:genus",
                               "specificepithet",
                               "data.dwc:specificEpithet",
                               "infraspecificepithet",
                               "data.dwc:infraspecificEpithet",                             
                               "data.dwc:scientificName",
                               "flags"),
                    # Set the limit for how many records are returned by the
                    # search to a low number for the purposes of this demo
                    limit = 10)
}, error = function(e) {
    # Code to run if an error occurs
    cat("An error occurred during the idig_search_records call: ", e$message, "\n")
    cat("Vignettes will not be fully generated. Please try again after resolving the issue.")
    # Optionally, you can return NULL or an empty dataframe
    verify_df_names <- FALSE
})

## ----eval=verify_df_names-----------------------------------------------------
# Edit the value after `recordset` to search for data from a different collection
# and the fields (e.g. `uuid`) in `fields` to adjust the columns returned in
# your results
df_names <- idig_search_records(rq = list(recordset = "5082e6c8-8f5b-4bf6-a930-e3e6de7bf6fb"),
                    fields = c("uuid",
                               "data.dwc:occurrenceID",
                               "data.dwc:catalogNumber",
                               "family",
                               "data.dwc:family",
                               "genus",
                               "data.dwc:genus",
                               "specificepithet",
                               "data.dwc:specificEpithet",
                               "infraspecificepithet",
                               "data.dwc:infraspecificEpithet",                             
                               "data.dwc:scientificName",
                               "flags"),
                    # Set the limit for how many records are returned by the
                    # search to a low number for the purposes of this demo
                    limit = 1000) %>% 
  # Rename fields to more easily reflect their provenance (either from the
  # data provider directly or modified by the data aggregator)
  rename(occurrenceID = `data.dwc:occurrenceID`,
         catalogNumber = `data.dwc:catalogNumber`,
         provider_family = `data.dwc:family`,
         provider_genus = `data.dwc:genus`,
         provider_species = `data.dwc:specificEpithet`,
         provider_subspecies = `data.dwc:infraspecificEpithet`,
         provider_scientificName = `data.dwc:scientificName`,
         aggregator_family = `family`,
         aggregator_genus = `genus`,
         aggregator_species = `specificepithet`,
         aggregator_subspecies = `infraspecificepithet`) %>% 
  # Reorder columns for easier viewing
  select(uuid, occurrenceID, catalogNumber, aggregator_family, provider_family,
         aggregator_genus, aggregator_species, aggregator_subspecies, 
         provider_genus, provider_species, provider_subspecies,
         provider_scientificName, flags)

## ----eval=verify_df_names, echo = FALSE---------------------------------------
# Subset `df_names` to show example
df_names[1:50,] %>% 
  select(-flags) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                font_size = 12,
                fixed_thead = T) %>% 
  column_spec(c(4,6,7,8), color = "red") %>% 
  scroll_box(width = "100%", height = "400px")

## ----eval=verify_df_names-----------------------------------------------------
# Reformat aggregator fields to title case
df_names <- df_names %>% 
  mutate(aggregator_family = str_to_title(aggregator_family)) %>% 
  mutate(aggregator_genus = str_to_title(aggregator_genus))

# Subset `df_names` to show example
df_names[1:5,] %>% 
  select(-flags) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                font_size = 12,
                fixed_thead = T) %>% 
  column_spec(c(4,6,7,8), color = "red") %>% 
  scroll_box(width = "100%", height = "400px")

## ----eval=verify_df_names-----------------------------------------------------
# Filter for rows where genus does not match
df_names %>% 
  filter(provider_genus != aggregator_genus) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                font_size = 12,
                fixed_thead = T) %>% 
  column_spec(c(4,6,7,8), color = "red") %>% 
  scroll_box(width = "100%", height = "400px")

## ----eval=verify_df_names-----------------------------------------------------
# Summarize modifications made to genus names
df_names %>% 
  filter(provider_genus != aggregator_genus) %>% 
  # Because of the nature of scientific names, it makes sense to group data by
  # all of the primary fields that comprise a scientific name
  group_by(provider_genus, provider_species, provider_subspecies,
           aggregator_genus, aggregator_species, aggregator_subspecies,
           provider_scientificName) %>% 
  # Count how many rows are affected by this modification made to genus name
  tally() %>% 
  # Order by frequency of rows affected
  arrange(desc(n)) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                font_size = 12,
                fixed_thead = T) %>% 
  column_spec(c(4,5,6), color = "red") %>% 
  scroll_box(width = "100%", height = "400px")

## ----eval=verify_df_names-----------------------------------------------------
# Search for specimen records of an example modified genus name
df_names %>% 
  filter(provider_genus == "Glossaulax" & provider_species == "reclusiana") %>%
  select(catalogNumber)

