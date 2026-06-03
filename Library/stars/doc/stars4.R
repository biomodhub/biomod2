## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, dev = "png")
set.seed(13579)
knitr::opts_chunk$set(fig.height = 4.5)
knitr::opts_chunk$set(fig.width = 6)

## ----fig.width=4.5, fig.height=4----------------------------------------------
suppressPackageStartupMessages(library(stars))
m = matrix(1:20, nrow = 5, ncol = 4)
dim(m) = c(x = 5, y = 4) # named dim
(s = st_as_stars(m))

## -----------------------------------------------------------------------------
dim(s[[1]])

## ----fig.width=4.5, fig.height=4----------------------------------------------
image(s, text_values = TRUE, axes = TRUE)

## ----fig.width=4.5, fig.height=4----------------------------------------------
attr(s, "dimensions")[[2]]$delta = -1
image(s, text_values = TRUE, axes = TRUE)

## -----------------------------------------------------------------------------
tif = system.file("tif/L7_ETMs.tif", package = "stars")
st_dimensions(read_stars(tif))["y"]

## -----------------------------------------------------------------------------
str(attr(st_dimensions(s), "raster"))

## -----------------------------------------------------------------------------
attr(attr(s, "dimensions"), "raster")$affine = c(0.1, 0.1)
plot(st_as_sf(s, as_points = FALSE), axes = TRUE, nbreaks = 20)

## -----------------------------------------------------------------------------
atan2(0.1, 1) * 180 / pi

## -----------------------------------------------------------------------------
attr(attr(s, "dimensions"), "raster")$affine = c(0.1, 0.2)
plot(st_as_sf(s, as_points = FALSE), axes = TRUE, nbreaks = 20)

## -----------------------------------------------------------------------------
atan2(c(0.1, 0.2), 1) * 180 / pi

## -----------------------------------------------------------------------------
x = c(0, 0.5, 1, 2, 4, 5)  # 6 numbers: boundaries!
y = c(0.3, 0.5, 1, 2, 2.2) # 5 numbers: boundaries!
(r = st_as_stars(list(m = m), dimensions = st_dimensions(x = x, y = y)))
st_bbox(r)
image(r, axes = TRUE, col = grey((1:20)/20))

## -----------------------------------------------------------------------------
x = c(0, 0.5, 1, 2, 4)  # 5 numbers: offsets only!
y = c(0.3, 0.5, 1, 2)   # 4 numbers: offsets only!
(r = st_as_stars(list(m = m), dimensions = st_dimensions(x = x, y = y)))
st_bbox(r)

## -----------------------------------------------------------------------------
x = c(0, 1, 2, 3, 4)  # 5 numbers: offsets only!
y = c(0.5, 1, 1.5, 2)   # 4 numbers: offsets only!
(r = st_as_stars(list(m = m), dimensions = st_dimensions(x = x, y = y)))
st_bbox(r)

## -----------------------------------------------------------------------------
x = st_as_stars(matrix(1:9, 3, 3), 
                st_dimensions(x = c(1, 2, 3), y = c(2, 3, 10), cell_midpoints = TRUE))

## ----eval=FALSE---------------------------------------------------------------
# install.packages("starsdata", repos = "https://cran.uni-muenster.de/pebesma/")

## -----------------------------------------------------------------------------
(s5p = system.file("sentinel5p/S5P_NRTI_L2__NO2____20180717T120113_20180717T120613_03932_01_010002_20180717T125231.nc", package = "starsdata"))

## ----echo=FALSE---------------------------------------------------------------
EVAL = s5p != ""

## ----eval=EVAL----------------------------------------------------------------
subs = gdal_subdatasets(s5p)
subs[[6]]

## ----eval=EVAL----------------------------------------------------------------
gdal_metadata(subs[[6]], "GEOLOCATION")

## ----eval=EVAL----------------------------------------------------------------
nit.c = read_stars(subs[[6]]) 
threshold = units::set_units(9e+36, mol/m^2)
nit.c[[1]][nit.c[[1]] > threshold] = NA
nit.c

## ----eval=EVAL----------------------------------------------------------------
plot(nit.c, breaks = "equal", reset = FALSE, axes = TRUE, as_points = TRUE, 
		 pch = 16,  logz = TRUE, key.length = 1)
maps::map('world', add = TRUE, col = 'red')

## ----eval=EVAL----------------------------------------------------------------
plot(nit.c, breaks = "equal", reset = FALSE, axes = TRUE, as_points = FALSE, 
		 border = NA, logz = TRUE, key.length = 1)
maps::map('world', add = TRUE, col = 'red')

## ----eval=EVAL----------------------------------------------------------------
(nit.c_ds = stars:::st_downsample(nit.c, 8))
plot(nit.c_ds, breaks = "equal", reset = FALSE, axes = TRUE, as_points = TRUE, 
		 pch = 16, logz = TRUE, key.length = 1)
maps::map('world', add = TRUE, col = 'red')

## ----eval=EVAL----------------------------------------------------------------
plot(nit.c_ds, breaks = "equal", reset = FALSE, axes = TRUE, as_points = FALSE, 
		 border = NA, logz = TRUE, key.length = 1)
maps::map('world', add = TRUE, col = 'red')

## ----eval=EVAL----------------------------------------------------------------
w = st_warp(nit.c, crs = 4326, cellsize = 0.25)
plot(w)

