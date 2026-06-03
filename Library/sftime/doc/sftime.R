## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----packages-----------------------------------------------------------------
# load required packages
library(sftime)
library(sf)
library(stars)
library(spacetime)
library(ggplot2)
library(tidyr)

## ----sftime-class-1-----------------------------------------------------------
# example sfc object
x_sfc <- 
  sf::st_sfc(
    sf::st_point(1:2), 
    sf::st_point(c(1,3)), 
    sf::st_point(2:3), 
    sf::st_point(c(2,1))
  )

# create an sftime object directly from x_sfc
x_sftime1 <- sftime::st_sftime(a = 1:4, x_sfc, time = Sys.time()- 0:3 * 3600 * 24)

# first create the sf object and from this the sftime object
x_sf <- sf::st_sf(a = 1:4, x_sfc, time = x_sftime1$time)
x_sftime2 <- sftime::st_sftime(x_sf)

x_sftime3 <- sftime::st_as_sftime(x_sf) # alernative option

identical(x_sftime1, x_sftime2)
identical(x_sftime1, x_sftime3)

x_sftime1

## ----sftime-class-2-----------------------------------------------------------
methods(class = "sftime")

## ----time-column-1------------------------------------------------------------
# get the values from the time column
st_time(x_sftime1)
x_sftime1$time # alternative way

# set the values in the time column
st_time(x_sftime1) <- Sys.time()
st_time(x_sftime1)

# drop the time column to convert an sftime object to an sf object
st_drop_time(x_sftime1)
x_sftime1

# add a time column to an sf object converts it to an sftime object
st_time(x_sftime1, time_column_name = "time") <- Sys.time()
class(x_sftime1)

# These can also be used with pipes
x_sftime1 <-
  x_sftime1 |>
  st_drop_time() |>
  st_set_time(Sys.time(), time_column_name = "time")

## -----------------------------------------------------------------------------
# define the geometry column
g <- 
  st_sfc(
    st_point(c(1, 2)), 
    st_point(c(1, 3)), 
    st_point(c(2, 3)), 
    st_point(c(2, 1)), 
    st_point(c(3, 1))
  )

# crate sf object
x4_sf <- st_sf(a = 1:5, g, time = Sys.time() + 1:5)

# convert to sftime
x4_sftime <- st_as_sftime(x4_sf) 
class(x4_sftime)

## -----------------------------------------------------------------------------
# load sample data
x5_stars <- stars::read_ncdf(system.file("nc/bcsd_obs_1999.nc", package = "stars"), var = c("pr", "tas"))

# convert to sftime
x5_sftime <- st_as_sftime(x5_stars, time_column_name = "time")

## ----error = TRUE-------------------------------------------------------------
try({
# failed conversion to sftime
x5_sftime <- st_as_sftime(x5_stars, merge = TRUE, time_column_name = "time")
x5_sftime <- st_as_sftime(x5_stars, long = FALSE, time_column_name = "time")
})

## -----------------------------------------------------------------------------
# get sample data
example(STI, package = "spacetime")
class(stidf)

# conversion to sftime
x1_sftime <- st_as_sftime(stidf)

## -----------------------------------------------------------------------------
# get a sample TracksCollection
x2_TracksCollection <- trajectories::rTracksCollection(p = 2, m = 3, n = 40)

# convert to sftime
x2_TracksCollection_sftime <- st_as_sftime(x2_TracksCollection)
x2_Tracks_sftime <- st_as_sftime(x2_TracksCollection@tracksCollection[[1]])
x2_Track_sftime <- st_as_sftime(x2_TracksCollection@tracksCollection[[1]]@tracks[[1]])

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
# get a sample cubble_df object
climate_aus <- cubble::climate_aus

# convert to sftime
climate_aus_sftime <- 
  st_as_sftime(climate_aus[1:4, ])

climate_aus_sftime <- 
  st_as_sftime(cubble::face_temporal(climate_aus)[1:4, ])

## -----------------------------------------------------------------------------
st_time(x_sftime1)

## -----------------------------------------------------------------------------
# selecting rows and columns (works just as for sf objects)
x_sftime1[1, ]
x_sftime1[, 3]

# beware: the time column is not sticky. If omitted, the subset becomes an sf object
class(x_sftime1[, 1])
class(x_sftime1["a"]) # the same
x_sftime1[, 1]

# to retain the time column and an sftime object, explicitly select the time column during subsetting:
class(x_sftime1[, c(1, 3)])
class(x_sftime1[c("a", "time")]) # the same

## ----plotting-plot.sftime-1, fig.width=7--------------------------------------
coords <- matrix(runif(100), ncol = 2)
g <- sf::st_sfc(lapply(1:50, function(i) st_point(coords[i, ]) ))

x_sftime4 <- 
  st_sftime(
    a = 1:200,
    b = rnorm(200),
    id_object = as.factor(rep(1:4,each=50)),
    geometry = g, 
    time = as.POSIXct("2020-09-01 00:00:00") + 0:49 * 3600 * 6
) 

plot(x_sftime4, key.pos = 4)

## ----plotting-plot.sftime-2, fig.width=7--------------------------------------
plot(x_sftime4, number = 10, max.plot = 10, key.pos = 4)

## ----plotting-ggplot-1, fig.width=7-------------------------------------------
library(ggplot2)

ggplot() + 
  geom_sf(data = x_sftime4, aes(color = b)) + 
  facet_wrap(~ cut_number(time, n = 6)) +
  theme(
    panel.spacing.x = unit(4, "mm"), 
    panel.spacing.y = unit(4, "mm")
  )

## ----plotting-ggplot-2, fig.width=7-------------------------------------------
ggplot(x_sftime4) + 
  geom_point(aes(y = id_object, x = time, color = b))

## ----plotting-ggplot-3, fig.width=7-------------------------------------------
x_sftime4 |>
  tidyr::pivot_longer(cols = c("a", "b"), names_to = "variable", values_to = "value") |>
  ggplot() + 
  geom_path(aes(y = value, x = time, color = variable)) +
  facet_wrap(~ id_object)

## ----plotting-ggplot-4, fig.width=7-------------------------------------------
x_sftime4 |>
  tidyr::pivot_longer(cols = c("a", "b"), names_to = "variable", values_to = "value") |>
  ggplot() + 
  geom_path(aes(y = value, x = time, color = id_object)) +
  facet_wrap(~ variable, scales = "free_y")

## ----eval=TRUE----------------------------------------------------------------
(tc <- as.POSIXct("2020-09-01 08:00:00")-0:3*3600*24)

## -----------------------------------------------------------------------------
tc
order(tc)
sort(tc)

## -----------------------------------------------------------------------------
# utility functions
as.character.interval <- function(x) {
  paste0("[", x[1], ", ", x[2], "]")
}

print.interval <- function(x, ...) {
  cat("Interval:", as.character(x), "\n")
}

#'[.intervals' <- function(x, i) {
#  sx <- unclass(x)[i]
#  class(sx) <- "intervals"
#  sx
#}

## -----------------------------------------------------------------------------
# time interval definition
i1 <- c(5.3,12)
class(i1) <- "interval"
i2 <- c(3.1,6)
class(i2) <- "interval"
i3 <- c(1.4,6.9)
class(i3) <- "interval"
i4 <- c(1,21)
class(i4) <- "interval"

intrvls <- structure(list(i1, i2, i3, i4), class = "Intervals")
# provide dedicated generic to xtfrm for class intervals

## -----------------------------------------------------------------------------
xtfrm.Intervals <- function(x) sapply(x, mean)
# - sort by centre
(tc <- intrvls)
order(tc)
sort(tc)[1]

## -----------------------------------------------------------------------------
# - sort by end
xtfrm.Intervals <- function(x) sapply(x, max)
(tc <- intrvls)
order(tc)
sort(tc)[1]

## -----------------------------------------------------------------------------
# - sort by start
xtfrm.Intervals <- function(x) sapply(x, min)
tc <- intrvls
order(tc)
sort(tc)[1]

