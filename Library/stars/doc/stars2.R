## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, dev = "png")
set.seed(13579)
knitr::opts_chunk$set(fig.height = 4.5)
knitr::opts_chunk$set(fig.width = 6)
EVAL = x = suppressWarnings(require(starsdata, quietly = TRUE))

## ----eval=FALSE---------------------------------------------------------------
# install.packages("starsdata", repos = "https://cran.uni-muenster.de/pebesma/")
# # possibly after: options(timeout = 100)
# # or from an alternative repository:
# # install.packages("starsdata", repos = "http://pebesma.staff.ifgi.de", type = "source")

## -----------------------------------------------------------------------------
library(stars)
tif = system.file("tif/L7_ETMs.tif", package = "stars")
rasterio = list(nXOff = 6, nYOff = 6, nXSize = 100, nYSize = 100, bands = c(1, 3, 4))
(x = read_stars(tif, RasterIO = rasterio))
dim(x)

## -----------------------------------------------------------------------------
st_dimensions(read_stars(tif))

## -----------------------------------------------------------------------------
rasterio = list(nXOff = 6, nYOff = 6, nXSize = 100, nYSize = 100,
                nBufXSize = 20, nBufYSize = 20, bands = c(1, 3, 4))
(x = read_stars(tif, RasterIO = rasterio))

## -----------------------------------------------------------------------------
rasterio = list(nXOff = 6, nYOff = 6, nXSize = 3, nYSize = 3,
   nBufXSize = 100, nBufYSize = 100, bands = 1)
x = read_stars(tif, RasterIO = rasterio)
dim(x)
plot(x)

## -----------------------------------------------------------------------------
rasterio = list(nXOff = 6, nYOff = 6, nXSize = 3, nYSize = 3,
   nBufXSize = 100, nBufYSize = 100, bands = 1, resample = "cubic_spline")
x = read_stars(tif, RasterIO = rasterio)
dim(x)
plot(x)

## ----eval=EVAL----------------------------------------------------------------
granule = system.file("sentinel/S2A_MSIL1C_20180220T105051_N0206_R051_T32ULE_20180221T134037.zip", package = "starsdata")
s2 = paste0("SENTINEL2_L1C:/vsizip/", granule, "/S2A_MSIL1C_20180220T105051_N0206_R051_T32ULE_20180221T134037.SAFE/MTD_MSIL1C.xml:10m:EPSG_32632")
(p = read_stars(s2, proxy = TRUE))

## ----eval=EVAL----------------------------------------------------------------
system.time(plot(p))

## ----eval = FALSE-------------------------------------------------------------
# p = read_stars(s2, proxy = FALSE)

## -----------------------------------------------------------------------------
methods(class = "stars_proxy")

## ----eval=EVAL----------------------------------------------------------------
x = c("avhrr-only-v2.19810901.nc",
"avhrr-only-v2.19810902.nc",
"avhrr-only-v2.19810903.nc",
"avhrr-only-v2.19810904.nc",
"avhrr-only-v2.19810905.nc",
"avhrr-only-v2.19810906.nc",
"avhrr-only-v2.19810907.nc",
"avhrr-only-v2.19810908.nc",
"avhrr-only-v2.19810909.nc")
file_list = system.file(paste0("netcdf/", x), package = "starsdata")
y = read_stars(file_list, quiet = TRUE, proxy = TRUE)
names(y)
y["sst"]

## ----eval=EVAL----------------------------------------------------------------
bb = st_bbox(c(xmin = 10.125, ymin = 0.125, xmax = 70.125, ymax = 70.125))
ysub = y[bb]
st_dimensions(ysub)
class(ysub) # still no data here!!
plot(ysub, reset = FALSE) # plot reads the data, at resolution that is relevant
plot(st_as_sfc(bb), add = TRUE, lwd = .5, border = 'red')

## ----eval=EVAL----------------------------------------------------------------
yy = adrop(y)
yyy = yy[,1:10,1:10,]
class(yyy) # still no data
st_dimensions(yyy) # and dimensions not adjusted
attr(yyy, "call_list") # the name of object in the call (y) is replaced with x:

## ----eval = FALSE-------------------------------------------------------------
# plot(st_apply(x, c("x", "y"), range))

## ----eval=EVAL----------------------------------------------------------------
(x = st_as_stars(yyy)) # read, adrop, subset

## ----eval=EVAL----------------------------------------------------------------
# S2 10m: band 4: near infrared, band 1: red.
#ndvi = function(x) (x[4] - x[1])/(x[4] + x[1])
ndvi = function(x1, x2, x3, x4) (x4 - x1)/(x4 + x1)
rm(x)
(s2.ndvi = st_apply(p, c("x", "y"), ndvi))
system.time(plot(s2.ndvi)) # read - compute ndvi - plot 

## -----------------------------------------------------------------------------
s1 = st_as_stars(matrix(1:16, 4))
s2 = st_as_stars(matrix(1:16, 4))
s3 = st_as_stars(matrix(1:16, 4))
attr(s1, "dimensions")$X1$offset = 0
attr(s1, "dimensions")$X2$offset = 4
attr(s2, "dimensions")$X1$offset = 0
attr(s2, "dimensions")$X2$offset = 4
attr(s3, "dimensions")$X1$offset = 0
attr(s3, "dimensions")$X2$offset = 4
attr(s1, "dimensions")$X1$delta =  1
attr(s1, "dimensions")$X2$delta = -1
attr(s2, "dimensions")$X1$delta =  2
attr(s2, "dimensions")$X2$delta = -2
attr(s3, "dimensions")$X1$delta =  3
attr(s3, "dimensions")$X2$delta = -3
plot(s1, axes = TRUE, text_values = TRUE, text_color = 'orange')
plot(s2, axes = TRUE, text_values = TRUE, text_color = 'orange')
plot(s3, axes = TRUE, text_values = TRUE, text_color = 'orange')

## ----eval=TRUE----------------------------------------------------------------
fn1 = paste0(tempdir(), .Platform$file.sep, "img1.tif")
fn2 = paste0(tempdir(), .Platform$file.sep, "img2.tif")
fn3 = paste0(tempdir(), .Platform$file.sep, "img3.tif")
write_stars(s1, fn1)
write_stars(s2, fn2)
write_stars(s3, fn3) 
(r1 = read_stars(c(fn1, fn2, fn3), proxy = TRUE))

## ----eval=TRUE----------------------------------------------------------------
st_as_stars(r1) |>
  merge() |>
  plot(breaks = "equal", text_values = TRUE, text_color = 'orange', axes = TRUE)

## ----eval=TRUE----------------------------------------------------------------
st_as_stars(r1[,2:4,2:4]) |>
  merge() |>
  plot(breaks = "equal", text_values = TRUE, text_color = 'orange', axes = TRUE)

## -----------------------------------------------------------------------------
s4 = st_as_stars(matrix(1: 16, 4))
s5 = st_as_stars(matrix(1: 64, 8))
s6 = st_as_stars(matrix(1:144,12))
attr(s4, "dimensions")$X1$offset = 0
attr(s4, "dimensions")$X2$offset = 4
attr(s5, "dimensions")$X1$offset = 0
attr(s5, "dimensions")$X2$offset = 4
attr(s6, "dimensions")$X1$offset = 0
attr(s6, "dimensions")$X2$offset = 4
attr(s4, "dimensions")$X1$delta =  1
attr(s4, "dimensions")$X2$delta = -1
attr(s5, "dimensions")$X1$delta =  1/2
attr(s5, "dimensions")$X2$delta = -1/2
attr(s6, "dimensions")$X1$delta =  1/3
attr(s6, "dimensions")$X2$delta = -1/3
plot(s4, axes = TRUE, text_values = TRUE, text_color = 'orange')
plot(s5, axes = TRUE, text_values = TRUE, text_color = 'orange')
plot(s6, axes = TRUE, text_values = TRUE, text_color = 'orange')

## ----eval=TRUE----------------------------------------------------------------
fn4 = paste0(tempdir(), .Platform$file.sep, "img4.tif")
fn5 = paste0(tempdir(), .Platform$file.sep, "img5.tif")
fn6 = paste0(tempdir(), .Platform$file.sep, "img6.tif")
write_stars(s4, fn4)
write_stars(s5, fn5)
write_stars(s6, fn6) 
(r2 = read_stars(c(fn4, fn5, fn6), proxy = TRUE))

st_as_stars(r2) |>
  merge() |>
  plot(breaks = "equal", text_values = TRUE, text_color = 'orange', axes = TRUE)
st_as_stars(r2[,2:4,2:4]) |>
  merge() |>
  plot(breaks = "equal", text_values = TRUE, text_color = 'orange', axes = TRUE)

## ----eval=TRUE----------------------------------------------------------------
(r3 = read_stars(c(fn6, fn5, fn4), proxy = TRUE))

st_as_stars(r3) |>
  merge() |>
  plot(breaks = "equal", text_values = TRUE, text_color = 'orange', axes = TRUE)
st_as_stars(r3[,2:6,3:6]) |>
  merge() |>
  plot(breaks = "equal", text_values = TRUE, text_color = 'orange', axes = TRUE)

