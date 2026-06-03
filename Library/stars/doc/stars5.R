## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, dev = "png")
suppressPackageStartupMessages(library(sf))
knitr::opts_chunk$set(fig.height = 4.5)
knitr::opts_chunk$set(fig.width = 6)
ev = TRUE

## -----------------------------------------------------------------------------
library(stars)
system.file("gpkg/nc.gpkg", package = "sf") |>
  read_sf() |>
  st_transform(32119) -> nc
nc$dens = nc$BIR79 / units::set_units(st_area(nc), km^2)
(nc.st = st_rasterize(nc["dens"], dx = 5000, dy = 5000))
plot(nc.st)

## -----------------------------------------------------------------------------
tif = system.file("tif/L7_ETMs.tif", package = "stars")
x = read_stars(tif)[, 1:50, 1:50, 1:2]
x[[1]] = round(x[[1]]/5)

## ----eval=ev------------------------------------------------------------------
l =  st_contour(x, contour_lines = TRUE, breaks = 11:15)
plot(l[1], key.pos = 1, pal = sf.colors, lwd = 2, key.length = 0.8)

## -----------------------------------------------------------------------------
st_as_sf(x, as_points = TRUE, merge = FALSE)

## -----------------------------------------------------------------------------
st_as_sf(x, as_points = TRUE, merge = FALSE, long = TRUE)

## -----------------------------------------------------------------------------
st_as_sf(x[1], as_points = FALSE, merge = FALSE)

## -----------------------------------------------------------------------------
p = st_as_sf(x, as_points = FALSE, merge = TRUE)

## -----------------------------------------------------------------------------
plot(p)

## -----------------------------------------------------------------------------
x.sf = st_xy2sfc(x, as_points = TRUE)
x.sf

## -----------------------------------------------------------------------------
nc.st |> st_transform("+proj=laea +lat_0=34 +lon_0=-60") -> nc.curv
nc.curv
plot(nc.curv, border = NA, graticule = TRUE)

## -----------------------------------------------------------------------------
nc |> st_transform("+proj=laea +lat_0=34 +lon_0=-60") |> st_bbox() |>
	st_as_stars() -> newgrid

## -----------------------------------------------------------------------------
nc.st |> st_warp(newgrid) -> nc.new
nc.new 
plot(nc.new)

