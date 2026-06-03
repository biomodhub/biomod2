## ----setup, include = FALSE---------------------------------------------------
library(BIEN)

bien_vignette_online <- !is.null(
  tryCatch(expr = BIEN:::.BIEN_sql("SELECT 1;"),
           error = function(e){NULL})
)

knitr::opts_chunk$set(echo = TRUE,
                      eval = bien_vignette_online)

if(!bien_vignette_online){
  message("BIEN database is unavailable or refusing connections; remote-query chunks were not evaluated while building this vignette.")
}

## ----load-packages, message = F,warning = FALSE, results = 'hide'-------------
library(ape) #Package for working with phylogenies in R
library(maps) #Useful for making quick maps of occurrences
library(sf) # A package for spatial data

## ----load-vignette, eval = FALSE----------------------------------------------
# vignette("BIEN")

## ----xs-occurrences-----------------------------------------------------------

Xanthium_strumarium <- BIEN_occurrence_species(species = "Xanthium strumarium")


## ----view-xs-occurrences------------------------------------------------------

str(Xanthium_strumarium)

head(Xanthium_strumarium)


## ----xs-occurrences-pt-2, eval = FALSE----------------------------------------
# 
# Xanthium_strumarium_full <- BIEN_occurrence_species(species = "Xanthium strumarium",
#                                                     cultivated = TRUE,
#                                                     all.taxonomy = TRUE,
#                                                     native.status = TRUE,
#                                                     observation.type = TRUE,
#                                                     political.boundaries = TRUE)
# 
# 

## ----map-xs, eval = FALSE-----------------------------------------------------
# # Make a quick map to plot our points on
# 
#   map('world', fill = TRUE, col= "grey", bg = "light blue")
# 
# # Plot the points from the full query in red
# 
#   points(cbind(Xanthium_strumarium_full$longitude,
#                Xanthium_strumarium_full$latitude),
#          col = "red",
#          pch = 20,
#          cex = 1)
# 
# # Plot the points from the default query in blue
# 
#   points(cbind(Xanthium_strumarium$longitude,
#                Xanthium_strumarium$latitude),
#          col = "blue",
#          pch = 20,
#          cex = 1)
# 

## ----occs-bahamas, eval = FALSE-----------------------------------------------
# 
# Bahamas <- BIEN_occurrence_country(country =  "Bahamas")
# 
# #Let's see how many species we have
# 
#   length(unique(Bahamas$scrubbed_species_binomial))
#   #About 400 species with valid occurrence records.
# 
# #Now, let's take a look at where those occurrences are:
# 
#   map(regions = "Bahamas" ,
#       fill = TRUE ,
#       col= "grey",
#       bg = "light blue")
# 
#   points(cbind(Bahamas$longitude,Bahamas$latitude),
#          col = "blue",
#          pch = 20,
#          cex = 1)
# 
# #Looks like some islands are considerably better sampled than others.
# 

## ----xs-range-----------------------------------------------------------------

Xanthium_strumarium_range <- BIEN_ranges_load_species(species = "Xanthium strumarium")


## ----xs-range-map-------------------------------------------------------------

#First, let's add a base map so that our range has some context:

  map('world', fill = TRUE ,
      col= "grey",
      bg = "light blue",
      xlim = c(-180, -20),
      ylim = c(-60, 80))

#Now, we can add the range map:
  plot(Xanthium_strumarium_range[1],
       col = "green",
       add = TRUE)




## ----xs-range-and-points------------------------------------------------------

map('world',
    fill = TRUE ,
    col = "grey",
    bg = "light blue",
    xlim = c(-180, -20),
    ylim = c(-60, 80))

plot(Xanthium_strumarium_range[1],
     col = "green",
     add = TRUE)

points(cbind(Xanthium_strumarium$longitude,Xanthium_strumarium$latitude),
       col = "blue",
       pch = 20,
       cex = 1)



## ----luq-occs, eval = FALSE---------------------------------------------------
# 
# LUQUILLO <- BIEN_plot_name(plot.name = "LUQUILLO")
# 
# head(LUQUILLO)
# 

## ----eval = FALSE-------------------------------------------------------------
# 
# LUQUILLO_full <- BIEN_plot_name(plot.name = "LUQUILLO",
#                                 cultivated = TRUE,
#                                 all.taxonomy = TRUE,
#                                 native.status = TRUE,
#                                 political.boundaries = TRUE,
#                                 all.metadata = TRUE)
# 

## ----salix-traits, eval = FALSE-----------------------------------------------
# 
# Salix_traits <- BIEN_trait_genus(genus = "Salix")
# 

## ----trait-list, eval = FALSE-------------------------------------------------
# 
# BIEN_trait_list()
# 

## ----leaf-area, eval = FALSE--------------------------------------------------
# 
# leaf_area <- BIEN_trait_trait(trait = "leaf area")
# 

## ----taxonomy-----------------------------------------------------------------

Asclepias_taxonomy <- BIEN_taxonomy_genus(genus = "Asclepias")

#We see that the genus Asclepias falls within the family Apocynaceae and the order Gentianales.

#You'll also notice that a given species may appear more than once (due to multiple circumscriptions, some of which may be illegitimate).

#If we'd just like to know all the speciess that aren't illegitimate:

Asclepias_species <- unique(Asclepias_taxonomy$scrubbed_species_binomial[Asclepias_taxonomy$scrubbed_taxonomic_status %in% c("accepted",  "no opinion")])



## ----phylogeny, eval = FALSE--------------------------------------------------
# 
# phylo <- BIEN_phylogeny_conservative()
# 
# #Let's make sure it looks alright
# 
#   plot.phylo(x = phylo, show.tip.label =  FALSE)
# 
# #If we just want to see which species are included
# 
#   phylo_species <- phylo$tip.label
# 
# 

## ----stems, eval = FALSE------------------------------------------------------
# 
# Cupressus_arizonica_stems <- BIEN_stem_species("Cupressus arizonica")
# 

## ----lists, eval = FALSE------------------------------------------------------
# 
# Bahamas_species_list <- BIEN_list_country(country = "Bahamas")
# 
# #Notice that we find many more species listed than we found occurrence records for.  What happened?  There are many records coming from the Bahamas that lack coordinates.  These records are used used in the "_list_" functions, but not the occurrence functions.
# 

## ----lists-multi-country, eval = FALSE----------------------------------------
# 
# country_vector <- c("Haiti","Dominican Republic")
# 
# Haiti_DR <- BIEN_list_country(country = country_vector)
# 

## ----lists-pol-divs-----------------------------------------------------------
#To see all of the political division names, and associated codes, we can use this function:

political_names <- BIEN_metadata_list_political_names()

#Let's take a look at what the dataframe contains:

head(political_names)

#In addition to the standardized country, state (state_province_ascii) and county (county_parish_ascii) names, we have the associated codes that can be used in BIEN functions.
#Note that 'state' refers to any primary political division (e.g. province), and 'county' refers to any secondary political division (e.g. parish).

#Looking at the political_names dataframe, we see that the Dominican Republic has country code "DO", and Haiti has country code "HT"

Haiti_DR_from_codes <- BIEN_list_country(country.code = c("HT","DO"))



## ----md-1---------------------------------------------------------------------

BIEN_metadata_database_version()


## ----selaginella-occs, eval = FALSE-------------------------------------------
# 
# Selaginella_selaginoides_occurrences <- BIEN_occurrence_species("Selaginella selaginoides", new.world = NULL)
# 

## ----selaginella-md, eval = FALSE---------------------------------------------
# 
# citation_info <- BIEN_metadata_citation(dataframe = Selaginella_selaginoides_occurrences)
# 

## ----md-2, eval = FALSE-------------------------------------------------------
# 
# temp_dir <- file.path(tempdir(), "BIEN_temp") #Set a temporary working directory
# 
# 
# citation_info <- BIEN_metadata_citation(dataframe = Selaginella_selaginoides_occurrences,
#                                      bibtex_file = file.path(temp_dir,"selaginella_selaginoides.bib"),
#                                       acknowledgement_file = file.path(temp_dir,"selaginella_selaginoides.txt"))
# 
# 

## ----md-3, eval = FALSE-------------------------------------------------------
# #First, let's get some trait data:
# selaginella_selaginoides_traits <- BIEN_trait_species(species = "Selaginella selaginoides")
# 
# #Now, we just need to modify our previous bit of code to include the trait data as well:
# 
# temp_dir <- file.path(tempdir(), "BIEN_temp")
# 
# citation_info <- BIEN_metadata_citation(dataframe = Selaginella_selaginoides_occurrences,
#                                       trait.dataframe = selaginella_selaginoides_traits,
#                                       bibtex_file = file.path(temp_dir,"selaginella_selaginoides.bib"),
#                                       acknowledgement_file = file.path(temp_dir,"selaginella_selaginoides.txt"))
# 

