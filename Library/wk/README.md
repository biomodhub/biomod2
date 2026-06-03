
<!-- README.md is generated from README.Rmd. Please edit that file -->

# wk

<!-- badges: start -->

[![R build
status](https://github.com/paleolimbot/wk/workflows/R-CMD-check/badge.svg)](https://github.com/paleolimbot/wk/actions)
[![Codecov test
coverage](https://codecov.io/gh/paleolimbot/wk/branch/master/graph/badge.svg)](https://app.codecov.io/gh/paleolimbot/wk?branch=master)
<!-- badges: end -->

The goal of wk is to provide lightweight R, C, and C++ infrastructure
for a distributed ecosystem of packages that operate on collections of
coordinates. First, wk provides vector classes for points, circles,
rectangles, well-known text (WKT), and well-known binary (WKB). Second,
wk provides a C API and set of S3 generics for event-based iteration
over vectors of geometries.

## Installation

You can install the released version of wk from
[CRAN](https://cran.r-project.org/) with:

``` r
install.packages("wk")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("paleolimbot/wk")
```

If you can load the package, you’re good to go!

``` r
library(wk)
```

## Vector classes

Use `wkt()` to mark a character vector as containing well-known text, or
`wkb()` to mark a vector as well-known binary. Use `xy()`, `xyz()`,
`xym()`, and `xyzm()` to create vectors of points, and `rct()` to create
vectors of rectangles. These classes have full
[vctrs](https://vctrs.r-lib.org) support and `plot()`/`format()` methods
to make them as frictionless as possible working in R and RStudio.

``` r
wkt("POINT (30 10)")
#> <wk_wkt[1]>
#> [1] POINT (30 10)
as_wkb(wkt("POINT (30 10)"))
#> <wk_wkb[1]>
#> [1] <POINT (30 10)>
xy(1, 2)
#> <wk_xy[1]>
#> [1] (1 2)
rct(1, 2, 3, 4)
#> <wk_rct[1]>
#> [1] [1 2 3 4]
crc(0, 0, 1)
#> <wk_crc[1]>
#> [1] [0 0, r = 1]
```

## Generics

The wk package is made up of readers, handlers, and filters. Readers
parse the various formats supported by the wk package, handlers
calculate values based on information from the readers (e.g.,
translating a vector of geometries into another format), and filters
transform information from the readers (e.g., transforming coordinates)
on the fly. The `wk_handle()` and `wk_translate()` generics power
operations for many geometry vector formats without having to explicitly
support each one.

## C API

The distributed nature of the wk framework is powered by a [~100-line
header](https://github.com/paleolimbot/wk/blob/master/inst/include/wk-v1.h)
describing the types of information that parsers typically encounter
when reading geometries and the order in which that information is
typically organized. Detailed information is available in the [C and C++
API
article](https://paleolimbot.github.io/wk/articles/articles/programming.html).

``` r
wk_debug(
  as_wkt("LINESTRING (1 1, 2 2, 3 3)"),
  wkt_format_handler(max_coords = 2)
)
#> initialize (dirty = 0  -> 1)
#> vector_start: <Unknown type / 0>[1] <0x16d75aac0> => WK_CONTINUE
#>   feature_start (1): <0x16d75aac0>  => WK_CONTINUE
#>     geometry_start (<none>): LINESTRING[UNKNOWN] <0x16d75a950> => WK_CONTINUE
#>       coord (1): <0x16d75a950> (1.000000 1.000000)  => WK_CONTINUE
#>       coord (2): <0x16d75a950> (2.000000 2.000000)  => WK_ABORT_FEATURE
#> vector_end: <0x16d75aac0>
#> deinitialize
#> [1] "LINESTRING (1 1, 2 2..."
```

## sf support

The wk package implements a reader and writer for sfc objects so you can
use them wherever you’d use an `xy()`, `rct()`, `crc()`, `wkb()`, or
`wkt()`:

``` r
wk_debug(
  sf::st_sfc(sf::st_linestring(rbind(c(1, 1), c(2, 2), c(3, 3)))),
  wkt_format_handler(max_coords = 2)
)
#> initialize (dirty = 0  -> 1)
#> vector_start: LINESTRING B[1] <0x16d75dac8> => WK_CONTINUE
#>   feature_start (1): <0x16d75dac8>  => WK_CONTINUE
#>     geometry_start (<none>): LINESTRING[3] <0x16d75da10> => WK_CONTINUE
#>       coord (1): <0x16d75da10> (1.000000 1.000000)  => WK_CONTINUE
#>       coord (2): <0x16d75da10> (2.000000 2.000000)  => WK_ABORT_FEATURE
#> vector_end: <0x16d75dac8>
#> deinitialize
#> [1] "LINESTRING (1 1, 2 2..."
```

## Lightweight

The wk package has zero dependencies and compiles in ~10 seconds.
