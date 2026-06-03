## ----message=FALSE------------------------------------------------------------
# Load core libraries; install these packages if you have not already
library(ridigbio)
library(tidyverse)

# Load library for making nice HTML output
library(kableExtra)

# Load libraries for visualizing geographic data
library(leaflet)

library(cowplot)

## ----echo = FALSE-------------------------------------------------------------

verify_df_flagCoord <- FALSE

#Test that examples will run
tryCatch({
    # Your code that might throw an error
    verify_df_flagCoord <- idig_search_records(
        rq = list(flags = "rev_geocode_corrected", institutioncode = "lacm"),
        limit = 10
    )
}, error = function(e) {
    # Code to run if an error occurs
    cat("An error occurred during the idig_search_records call: ", e$message, "\n")
    cat("Vignettes will not be fully generated. Please try again after resolving the issue.")
    # Optionally, you can return NULL or an empty dataframe
    verify_df_flagCoord <- FALSE
})

## ----eval=verify_df_flagCoord-------------------------------------------------
# Edit the fields (e.g. `flags`) and values (e.g. "rev_geocode_corrected") in
# `list()` to adjust your query and the fields (e.g. `uuid`) in `fields` to
# adjust the columns returned in your results
df_flagCoord <- idig_search_records(rq = list(flags = "rev_geocode_corrected",
                                              institutioncode = "lacm"),
                    fields = c("uuid",
                               "institutioncode",
                               "collectioncode",
                               "country",
                               "data.dwc:country",
                               "stateprovince",
                               "county",
                               "locality",
                               "geopoint",
                               "data.dwc:decimalLongitude",
                               "data.dwc:decimalLatitude",
                               "flags"),
                    limit = 100000) %>% 
  # Rename fields to more easily reflect their provenance (either from the
  # data provider directly or modified by the data aggregator)
  rename(provider_lon = `data.dwc:decimalLongitude`,
         provider_lat = `data.dwc:decimalLatitude`,
         provider_country = `data.dwc:country`,
         aggregator_lon = `geopoint.lon`,
         aggregator_lat = `geopoint.lat`,
         aggregator_country = country,
         aggregator_stateprovince = stateprovince,
         aggregator_county = county,
         aggregator_locality = locality) %>% 
  # Reorder columns for easier viewing
  select(uuid, institutioncode, collectioncode, provider_lat, aggregator_lat,
         provider_lon, aggregator_lon, provider_country, aggregator_country,
         aggregator_stateprovince, aggregator_county, aggregator_locality,
         flags)

## ----eval=verify_df_flagCoord, echo = FALSE-----------------------------------
# Subset `df_flagCoord` to show example
df_flagCoord[1:50,] %>% 
  select(-flags) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                font_size = 12,
                fixed_thead = T) %>% 
  scroll_box(width = "100%", height = "400px")

## ----eval=verify_df_flagCoord-------------------------------------------------
# Create function to allow subsetting the `df_flagCoord` dataset by other flags
# found on these same records
df_flagSubset <- function(subsetFlag) {
  df_flagCoord %>% 
  filter(grepl(subsetFlag, flags)) %>% 
  select(uuid, matches("_lat|_lon")) %>% 
  unite(provider_coords, c("provider_lat", "provider_lon"), sep = ",") %>% 
  unite(aggregator_coords, c("aggregator_lat", "aggregator_lon"), sep = ",") %>% 
  gather(key = type, value = coordinates, -uuid) %>% 
  separate(coordinates, c("lat","lon"), sep = ",") %>% 
  mutate(lat = as.numeric(lat)) %>% 
  mutate(lon = as.numeric(lon)) %>% 
  arrange(uuid, type)}

# Subset `df_flagCoord` by records flagged for having had their latitude negated
# to place point in stated country by reverse geocoding process
df_rev_geocode_lat_sign <- df_flagSubset("rev_geocode_lat_sign")

# Create map displaying a few examples of records with the
# rev_geocode_flip_lat_sign flag
pal <- leaflet::colorFactor(palette = c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba"),
                   domain = df_rev_geocode_lat_sign$uuid[1:10])

map <- df_rev_geocode_lat_sign[1:10,] %>% 
  mutate(popup = str_c(type, " = ", lat, ", ", lon, sep = "")) %>% 
  leaflet() %>%
  addTiles() %>% 
  addCircleMarkers(
    lng = ~lon,
    lat = ~lat,
    radius = 10,
    weight = 1,
    color = ~pal(uuid),
    stroke = FALSE,
    fillOpacity = 100,
    popup = ~popup) %>% 
    addLegend("bottomright", pal = pal, values = ~uuid,
    title = "Specimen Records",
    opacity = 1)

## ----eval=verify_df_flagCoord, echo = FALSE, out.width = '100%'---------------
map

## ----eval=verify_df_flagCoord-------------------------------------------------
# Summarize flagged records by collection type
spmByColl <- df_flagCoord %>% 
  group_by(collectioncode) %>% 
  tally()

# Generate graph to display counts of flagged records by collection within the
# institution
graph_spmByColl <- ggplot(spmByColl, 
                          aes(x = reorder(collectioncode, -n), 
                              y = n,
                              fill = collectioncode)) +
  geom_col() +
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(x = "collection", 
       y = "# of specimen records",
       title = "LACM records flagged with geo-coordinate data quality issues by iDigBio") +
  geom_text(aes(label = n, vjust = -0.5))

# Get count of total records published by the institution using function
# `idig_count_records`
totalInstSpm <- idig_count_records(rq = list(institutioncode = "lacm"))

# Calculate flagged records as percent of total records
percentFlagged <- sum(spmByColl$n)/totalInstSpm*100

## ----eval=verify_df_flagCoord, out.width="700px", echo = FALSE----------------
graph_spmByColl <- graph_spmByColl +
                   theme_minimal_grid() +
                   theme(
                       text = element_text(size = 22),
                       axis.text = element_text(size = 22),
                       plot.title = element_text(size = 22, face = "bold")
                   )

knitr::include_graphics(save_plot("plot.png", graph_spmByColl, base_height = 10, base_width = 24))

## ----eval=verify_df_flagCoord-------------------------------------------------
# Collate `df_flagAssoc` to describe other data quality flags that are associated
# with rev_geocode_corrected in `df_flagCoord`
df_flagAssoc <- df_flagCoord %>% 
  select(uuid, flags) %>% 
  unnest(flags) %>% 
  group_by(flags) %>% 
  tally() %>% 
  mutate("category" = case_when(str_detect(flags, "geo|country|state")
                              ~ "geography",
                      str_detect(flags, "dwc_datasetid_added|dwc_multimedia_added|datecollected_bounds")
                              ~ "other",
                      str_detect(flags, "gbif|dwc|tax")
                              ~ "taxonomy")) %>% 
  mutate("percent" = n/(nrow(df_flagCoord))*100) %>% 
  arrange(category, desc(n))

# Visualize associated data quality flags
graph_spmByColl <- ggplot(df_flagAssoc, aes(x = reorder(flags, -percent), y = percent, fill = category)) +
  geom_col() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(angle = 75, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(size = 12, face = "bold")
        ) +
  labs(x = "additional iDigBio data quality flag", 
       y = "% specimen records",
       title = "LACM records flagged for geo-coordinate issues are also flagged for...",
       fill = "flag category")

## ----eval=verify_df_flagCoord, out.width="700px", echo = FALSE----------------
graph_spmByColl <- graph_spmByColl +
                   theme_minimal_grid() +
                   theme(
                       text = element_text(size = 22),
                       axis.text = element_text(size = 22),
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                       plot.title = element_text(size = 22, face = "bold")
                   )

knitr::include_graphics(save_plot("plot2.png", graph_spmByColl, base_height = 10, base_width = 24))

