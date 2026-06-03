## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, dev = "png")
suppressPackageStartupMessages(library(dplyr))
knitr::opts_chunk$set(fig.height = 4.5)
knitr::opts_chunk$set(fig.width = 6)
EVAL = suppressWarnings(require(starsdata, quietly = TRUE))

## -----------------------------------------------------------------------------
library(stars)
f <- system.file("nc/reduced.nc", package = "stars")
(nc <- read_ncdf(f))

## -----------------------------------------------------------------------------
old_options <- options("stars.n_proxy" = 100)
(nc <- read_ncdf(f, proxy = TRUE))
options(old_options)

## -----------------------------------------------------------------------------
(nc <- read_ncdf(f, 
				 var = "sst", 
				 ncsub = cbind(start = c(90, 45, 1 , 1), 
				 			  count = c(90, 45, 1, 1))))

plot(nc)

