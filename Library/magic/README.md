Manipulation of high-dimensional arrays in R with the magic package
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/magic.png" width = "150" align="right" />

<!-- badges: start -->

[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/magic)](https://CRAN.R-project.org/package=magic)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/magic)](https://cran.r-project.org/package=magic)
<!-- badges: end -->

# Overview

The magic package implements functionality for manipulating
high-dimensional arrays using efficient vectorised methods. The original
application was high-dimensional magic hypercubes. This README shows
some of the more useful functions in the package.

# Installation

You can install the released version of `magic` from
[CRAN](https://CRAN.R-project.org) with:

``` r
# install.packages("magic")  # uncomment this to install the package
library("magic")
```

# Package highlights

-   Function `adiag()` binds arbitrarily-dimensioned arrays
    corner-to-corner
-   Function `apad()` pads arbitrarily-dimensioned arrays
-   Function `apldrop()` is a replacement for APL’s drop
-   Function `aplus()` superimposes two arrays of different dimensions
    and returns the sum of overlapping elements
-   Function `arev()` is a multidimensional generalization of `rev()`
-   Function `arot()` is a generalization of matlab’s `rotdim`
-   Function `fnsd()` returns the first nonsingleton dimension of an
    arbitrary dimensioned array
-   Function `ashift()` shifts the origin of arbitrary dimensioned
    arrays

Much of the package functionality is vectorised in array dimension.

# Further information

For more detail, see the package vignette

`vignette("magic")`
