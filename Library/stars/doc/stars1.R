## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, dev = "png")
ev = suppressWarnings(require(starsdata, quietly = TRUE))
knitr::opts_chunk$set(fig.height = 4.5)
knitr::opts_chunk$set(fig.width = 6)

## -----------------------------------------------------------------------------
library(stars)

## -----------------------------------------------------------------------------
methods(class = "stars")

## -----------------------------------------------------------------------------
tif = system.file("tif/L7_ETMs.tif", package = "stars")
x = read_stars(tif)
plot(x, axes = TRUE)

## -----------------------------------------------------------------------------
x

## ----eval=ev------------------------------------------------------------------
library(cubelyr)
as.tbl_cube(x)

## -----------------------------------------------------------------------------
(x.spl = split(x, "band"))
merge(x.spl)

## -----------------------------------------------------------------------------
merge(x.spl) |>
  setNames(names(x)) |> 
  st_set_dimensions(3, values = paste0("band", 1:6)) |>
  st_set_dimensions(names = c("x", "y", "band"))

## -----------------------------------------------------------------------------
class(x[[1]])
dim(x[[1]])
x$two = 2 * x[[1]]
x

## -----------------------------------------------------------------------------
x["two", 1:10, , 2:4]

## -----------------------------------------------------------------------------
circle = st_sfc(st_buffer(st_point(c(293749.5, 9115745)), 400), crs = st_crs(x))
plot(x[circle][, , , 1], reset = FALSE)
plot(circle, col = NA, border = 'red', add = TRUE, lwd = 2)

## -----------------------------------------------------------------------------
x1 = read_stars(tif, options = c("OVERVIEW_LEVEL=1"))
x2 = read_stars(tif, options = c("OVERVIEW_LEVEL=2"))
x3 = read_stars(tif, options = c("OVERVIEW_LEVEL=3"))
dim(x1)
dim(x2)
dim(x3)
par(mfrow = c(1, 3), mar = rep(0.2, 4))
image(x1[,,,1])
image(x2[,,,1])
image(x3[,,,1])

## ----eval=ev------------------------------------------------------------------
system.file("nc/bcsd_obs_1999.nc", package = "stars") |>
	read_stars() -> w

## ----eval=ev------------------------------------------------------------------
w

## -----------------------------------------------------------------------------
system.file("nc/bcsd_obs_1999.nc", package = "stars") |>
    read_ncdf()

## ----eval=ev------------------------------------------------------------------
x = c(
"avhrr-only-v2.19810901.nc",
"avhrr-only-v2.19810902.nc",
"avhrr-only-v2.19810903.nc",
"avhrr-only-v2.19810904.nc",
"avhrr-only-v2.19810905.nc",
"avhrr-only-v2.19810906.nc",
"avhrr-only-v2.19810907.nc",
"avhrr-only-v2.19810908.nc",
"avhrr-only-v2.19810909.nc"
)
# see the second vignette:
# install.packages("starsdata", repos = "https://cran.uni-muenster.de/pebesma/")
file_list = system.file(paste0("netcdf/", x), package = "starsdata")
(y = read_stars(file_list, quiet = TRUE))

## ----eval=ev------------------------------------------------------------------
library(dplyr)
library(abind)
z <- y |> select(sst) |> adrop()

## ----eval=ev------------------------------------------------------------------
# convert POSIXct time to character, to please ggplot's facet_wrap()
z1 = st_set_dimensions(z, 3, values = as.character(st_get_dimension_values(z, 3)))
library(ggplot2)
library(viridis)
library(ggthemes)
ggplot() +  
  geom_stars(data = z1[1], alpha = 0.8, downsample = c(10, 10, 1)) + 
  facet_wrap("time") +
  scale_fill_viridis() +
  coord_equal() +
  theme_map() +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(2, "cm"))

## ----eval=ev------------------------------------------------------------------
write_stars(adrop(y[1]), "sst.tif")

## -----------------------------------------------------------------------------
prec_file = system.file("nc/test_stageiv_xyt.nc", package = "stars")
prec = read_ncdf(prec_file, curvilinear = c("lon", "lat"))
##plot(prec) ## gives error about unique breaks
## remove NAs, zeros, and give a large number
## of breaks (used for validating in detail)
qu_0_omit = function(x, ..., n = 22) {
  if (inherits(x, "units"))
    x = units::drop_units(na.omit(x))
  c(0, quantile(x[x > 0], seq(0, 1, length.out = n)))
}
library(dplyr) # loads slice generic
prec_slice = slice(prec, index = 17, along = "time")
plot(prec_slice, border = NA, breaks = qu_0_omit(prec_slice[[1]]), reset = FALSE)
nc = sf::read_sf(system.file("gpkg/nc.gpkg", package = "sf"), "nc.gpkg")
plot(st_geometry(nc), add = TRUE, reset = FALSE, col = NA, border = 'red')

## -----------------------------------------------------------------------------
nc = st_transform(nc, st_crs(prec_slice)) # datum transformation
plot(prec_slice[nc], border = NA, breaks = qu_0_omit(prec_slice[[1]]), reset = FALSE)
plot(st_geometry(nc), add = TRUE, reset = FALSE, col = NA, border = 'red')

## -----------------------------------------------------------------------------
nc = st_read(system.file("gpkg/nc.gpkg", package="sf")) 
to = from = st_geometry(nc) # 100 polygons: O and D regions
mode = c("car", "bike", "foot") # travel mode
day = 1:100 # arbitrary
library(units)
units(day) = as_units("days since 2015-01-01")
hour = set_units(0:23, h) # hour of day
dims = st_dimensions(origin = from, destination = to, mode = mode, day = day, hour = hour)
(n = dim(dims))
traffic = array(rpois(prod(n), 10), dim = n) # simulated traffic counts
(st = st_as_stars(list(traffic = traffic),  dimensions = dims))

## ----eval=ev------------------------------------------------------------------
st |> as.tbl_cube()

## ----eval=ev------------------------------------------------------------------
b <- st |> 
  as.tbl_cube() |> 
  filter(mode == "bike") |> 
  group_by(hour) |>
  summarise(traffic = mean(traffic)) |> 
  as.data.frame()
require(ggforce) # for plotting a units variable
ggplot() +  
  geom_line(data = b, aes(x = hour, y = traffic))

## -----------------------------------------------------------------------------
s = system.file("tif/lc.tif", package = "stars")
r = read_stars(s, proxy = FALSE) |> droplevels()
levels(r[[1]]) = abbreviate(levels(r[[1]]), 10) # shorten text labels
st_point(c(3190631, 3125)) |> st_sfc(crs = st_crs(r)) |> st_buffer(25000) -> pol1
st_point(c(3233847, 21027)) |> st_sfc(crs = st_crs(r)) |> st_buffer(10000) -> pol2
if (isTRUE(dev.capabilities()$rasterImage == "yes")) {
  plot(r, reset = FALSE, key.pos = 4)
  plot(c(pol1, pol2), col = NA, border = c('yellow', 'green'), lwd = 2, add = TRUE)
}

## -----------------------------------------------------------------------------
f = function(x) { tb = table(x); names(tb)[which.max(tb)] }

## -----------------------------------------------------------------------------
aggregate(r, c(pol1, pol2), f) |> st_as_sf()

