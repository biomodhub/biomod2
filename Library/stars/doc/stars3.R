## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, dev = "png")
ev = TRUE
knitr::opts_chunk$set(fig.height = 4.5)
knitr::opts_chunk$set(fig.width = 6)

## -----------------------------------------------------------------------------
library(stars)
library(dplyr)

## -----------------------------------------------------------------------------
methods(class = "stars")

## -----------------------------------------------------------------------------
system.file("tif/L7_ETMs.tif", package = "stars") |>
	read_stars() -> x
x

## -----------------------------------------------------------------------------
x |> slice(band, 6) -> x6
x6

## -----------------------------------------------------------------------------
x |> filter(x > 289000, x < 291000, band > 3) -> x7
x7

## -----------------------------------------------------------------------------
x |> pull(1) -> x8
class(x8)
dim(x8)

## -----------------------------------------------------------------------------
x |> mutate(band2 = 2 * L7_ETMs.tif) -> x2 
x2

## -----------------------------------------------------------------------------
x2 |> select(band2) -> x9
x9

## -----------------------------------------------------------------------------
library(ggplot2)
library(viridis)
ggplot() + 
  geom_stars(data = x) +
  coord_equal() +
  facet_wrap(~band) +
  theme_void() +
  scale_fill_viridis() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

