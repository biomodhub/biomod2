#' Presence-Absence data to build test SDM
#'
#' A dataset covering all the continent with presence/absence data for 6 mammal
#' species
#'
#' @format A data frame with 2488 rows and 10 variables:
#' \describe{
#'   \item{X_WGS84}{Longitude}
#'   \item{Y_WGS84}{Latitude}
#'   \item{ConnochaetesGnou}{Presence (1) or Absence (0) for black wildebeest}
#'   \item{GuloGulo}{Presence (1) or Absence (0) for wolverine}
#'   \item{PantheraOnca}{Presence (1) or Absence (0) for jaguar}
#'   \item{PteropusGiganteus}{Presence (1) or Absence (0) for indian flying fox}
#'   \item{TenrecEcaudatus}{Presence (1) or Absence (0) for Tailless tenrec}
#'   \item{VulpesVulpes}{Presence (1) or Absence (0) for red fox}
#' }

"DataSpecies"

# DataSpecies  <-
#   read.csv(
#     "../biomod2_old_inst_folder/external/species/mammals_table.csv",
#     row.names = 1)
# 
# usethis::use_data(DataSpecies, overwrite = TRUE)


#' Bioclimatic variables for SDM based on current condition
#'
#' A RasterStack with 5 bioclimatic variables commonly used for SDM and
#' describing current climate. Additional information available at
#' \href{https://www.worldclim.org/data/bioclim.html}{worldclim}
#'
#' @format A RasterStack with 5 layers:
#' \describe{
#'   \item{bio3}{Isothermality}
#'   \item{bio4}{Temperature Seasonality}
#'   \item{bio7}{Temperature Annual Range}
#'   \item{bio11}{Mean Temperature of Coldest Quarter}
#'   \item{bio12}{Annual Precipitation}
#' }

"bioclim_current"

# myFiles <- paste0('../biomod2_old_inst_folder/external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
# bioclim_current <- raster::stack(myFiles)
# usethis::use_data(bioclim_current, overwrite = TRUE)

#' Bioclimatic variables for SDM based on future condition
#'
#' A RasterStack with 5 bioclimatic variables commonly used for SDM and
#' describing future climate based on  CMIP6 downscaled projections. Additional
#' information available at
#' \href{https://www.worldclim.org/data/cmip6/cmip6climate.html}{worldclim}
#'
#' @format A RasterStack with 5 layers:
#' \describe{
#'   \item{bio3}{Isothermality}
#'   \item{bio4}{Temperature Seasonality}
#'   \item{bio7}{Temperature Annual Range}
#'   \item{bio11}{Mean Temperature of Coldest Quarter}
#'   \item{bio12}{Annual Precipitation}
#' }

"bioclim_future"

# myFiles <- paste0('../biomod2_old_inst_folder/external/bioclim/future/bio', c(3, 4, 7, 11, 12), '.grd')
# bioclim_future <- raster::stack(myFiles)
# usethis::use_data(bioclim_future, overwrite = TRUE)
