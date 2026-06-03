
<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![rcc](https://github.com/tidyverse/blob/workflows/rcc/badge.svg)](https://github.com/tidyverse/blob/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/blob)](https://cran.r-project.org/package=blob)
[![Coverage
Status](https://codecov.io/gh/tidyverse/blob/branch/main/graph/badge.svg)](https://app.codecov.io/gh/tidyverse/blob)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# blob

## Overview

The goal of blob is to provide a simple S3 class to represent a vector
of binary objects, aka blobs. The `blob` class is a lightweight wrapper
around a list of raw vectors, suitable for inclusion in a data frame.

In most cases you will not need to use this package explicitly: it will
be used transparently by packages that need to load BLOB columns from
databases or binary file formats.

## Installation

``` r
# The easiest way to get blob is to install the whole tidyverse:
install.packages("tidyverse")

# Alternatively, install just blob:
install.packages("blob")

# Or the development version from GitHub:
# install.packages("pak")
pak::pak("tidyverse/blob")
```

## Example

To create a blob, use `blob()`, `new_blob()` or `as_blob()`:

``` r
library(blob)

x1 <- charToRaw("Good morning")
x2 <- as.raw(c(0x48, 0x65, 0x6c, 0x6c, 0x6f))

new_blob(list(x1, x2))
#> <blob[2]>
#> [1] blob[12 B] blob[5 B]
blob(x1, x2)
#> <blob[2]>
#> [1] blob[12 B] blob[5 B]

as_blob(c("Good morning", "Good evening"))
#> <blob[2]>
#> [1] blob[12 B] blob[12 B]
```

------------------------------------------------------------------------

Please note that the ‘blob’ project is released with a [Contributor Code
of
Conduct](https://github.com/tidyverse/blob/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
