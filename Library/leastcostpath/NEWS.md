# 2.0.13 (09/06/2025)
* Fixed bad parts in URL field of DESCRIPTION file for CRAN

# 2.0.12 (09/10/2023)
* Fixed the updating of extent when using crop_cs

# 2.0.11 (02/05/2023)
* modified create_distance_cs to use calculated _run_ distance when assigning distance between cells
* added check_locations argument to create_accum_cost

# 2.0.10 (17/04/2023)
* create_lcp now only creates cost column if cost_distance argument is TRUE

# 2.0.9 (09/04/2023)
* removed percentile argument in create_slope_cs and replaced with their own cost_function, e.g. "campbell 2019 50"
* Modified how extent is retrieved from a terra SpatRaster. This is now done using terra::ext()
* calculate_distance uses Pythagorean theorem when coordinate system is projected and sf::st_distance when geographic
* removed plot_cf()

# 2.0.8 (24/03/2023)
* create_accum_cost now allows for one or more supplied origins. Multiple accumulated cost surfaces will be summarised using a supplied function
* exported calculate_distance, get_coordinates, and neighbourhood functions
* Modified check_locations to stop if any locations are outside extent or not traversable. check_locations also added to functions that are suppplied locations, e.g. create_lcp
* modified create_distance_cs to return conductanceMatrix values that correspond to the supplied spatRaster resolution
* Added vignette
* leastcostpath now allows for both projected and geographic coordinate systems

# 2.0.7 (11/03/2023)
* Fixed error when using terra spatVector objects in create_FETE_lcps and create_lcps
* create_accum_cost and create_lcp_density can now take terra spatVector objects
* Removed add_local_stochasticity, calculate_rmse and calculate_slope_variance. These will be amended and re-added at a later date
* Modified create_FETE_lcps to leverage that igraph shortest.paths function is vectorised. This function is now quicker. 

# 2.0.6 (16/01/2023)
* Modified cost functions "herzog" and "llobera-sluckin" to now be ansitropic, i.e. cost uphill is different to cost downhill
* "Minetti" continues to be isotropic given that downhill slope gradient values are given negative cost values

# 2.0.5 (10/12/2022)
* create_lcp now allows for multiple destinations. If supplied least-cost paths will be calculated from a single origin to all destinations
* origin, destination, and locations arguments now accepts sf POINT and MULTIPOINT, terra spatVector, data.frame and matrix objects
* Fixed max_slope argument to 30 degrees within create_slope_cs when using 'campbell 2019' and 'campbell' cost functions
* Added vignette

# 2.0.4 (04/11/2022)
* Modified add_dem_error to now allow for different methods. See function details for more information
* Added crop_cs() function to allow conductanceMatrix to be cropped based on extent of supplied object

# 2.0.3 (22/10/2022)
* modified calculation of run in cost surface calculations to use base R rather than terra::distance. This now allows for the use of DEMs with more cells without causing memory issues

# 2.0.2 (17/10/2022)
* Ensured that Matrix::summary is made explicit rather than rely on using summary from Matrix package.

# 2.0.1 (09/10/2022)
* renamed add_stochasticity to add_global_stochasticity
* added add_local_stochasticity
* added add_global_stochasticity
* added calculate_slope_variance
* added calculate_rmse

# 2.0.0 (05/10/2022)
* From version 2.0.0 onwards the R package leastcostpath is no longer reliant on the R package gdistance. leastcostpath has been updated to work with sf and terra objects
* create_slope_cs now returns a <i>class</i> conductanceMatrix object. This object contains a record of:
- the ConductanceMatrix
- the cost function argument
- the max slope argument
- whether the slope values were exaggerated
- the critical slope argument
- the percentile argument
- the number of adjacent neighbours used in the calculation
* the cost function argument in create_slope_cs() now allows for cost functions to be stated by name (e.g. 'tobler') or as a function
* update_values() now added
* replace_values() now added
* create_cs() now added
* plot_cf() now added
* rasterise() now added