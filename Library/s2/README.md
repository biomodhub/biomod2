
<!-- README.md is generated from README.Rmd. Please edit that file -->

# s2

<!-- badges: start -->

[![R-CMD-check](https://github.com/r-spatial/s2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-spatial/s2/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/r-spatial/s2/branch/main/graph/badge.svg)](https://app.codecov.io/gh/r-spatial/s2)
[![CRAN](http://www.r-pkg.org/badges/version/s2)](https://cran.r-project.org/package=s2)
[![Downloads](http://cranlogs.r-pkg.org/badges/s2?color=brightgreen)](https://www.r-pkg.org/pkg/s2)
<!-- badges: end -->

The s2 R package provides bindings to Google’s
[S2Geometry](http://s2geometry.io) library. The package exposes an API
similar to Google’s [BigQuery Geography
API](https://cloud.google.com/bigquery/docs/reference/standard-sql/geography_functions),
whose functions also operate on spherical geometries. Package
[sf](https://cran.r-project.org/package=sf) uses this package by default
for nearly all its geometrical operations on objects with ellipsoidal
(unprojected) coordinates; in cases where it doesn’t, such as
`st_relate()`, it emits a warning.

This package is a complete rewrite of an earlier CRAN package s2 with
versions up to 0.4-2, for which the sources are found
[here](https://github.com/spatstat/s2/).

## Installation

You can install the released version of s2 from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("s2")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("r-spatial/s2")
```

The S2 package requires [Abseil](https://github.com/abseil/abseil-cpp)
and OpenSSL. You can install these using a system package manager on
most platforms:

- Windows: Both OpenSSL and Abseil are available from RTools since R 4.3
- MacOS: `brew install openssl abseil`
- Debian/Ubuntu: `apt-get install libssl-dev libabsl-dev`
- Fedora: `dnf install openssl-devel abseil-cpp-devel`
- Alpine: `apk add abseil-cpp`

## Example

The s2 package provides geometry transformers and predicates similar to
those found in [GEOS](https://libgeos.org), except instead of assuming a
planar geometry, s2’s functions work in latitude and longitude and
assume a spherical geometry:

``` r
library(s2)

s2_contains(
  # polygon containing much of the northern hemisphere
  "POLYGON ((-63.5 44.6, -149.75 61.20, 116.4 40.2, 13.5 52.51, -63.5 44.6))",
  # ...should contain the north pole
  "POINT (0 90)"
)
#> [1] TRUE
```

The [sf package](https://r-spatial.github.io/sf/) uses s2 for geographic
coordinates by default (this can be confirmed by calling
`sf::sf_use_s2()`). The sf package also supports creating s2 vectors
using `as_s2_geography()`:

``` r
library(dplyr)
library(sf)

nc_s2 <- read_sf(system.file("shape/nc.shp", package = "sf")) %>%
  mutate(geometry = as_s2_geography(geometry)) %>%
  as_tibble() %>%
  select(NAME, geometry)

nc_s2
#> # A tibble: 100 × 2
#>    NAME        geometry                                                         
#>    <chr>       <s2_geography>                                                   
#>  1 Ashe        POLYGON ((-81.4528885 36.2395859, -81.4310379 36.2607193, -81.41…
#>  2 Alleghany   POLYGON ((-81.1766739 36.4154434, -81.1533661 36.4247398, -81.13…
#>  3 Surry       POLYGON ((-80.4530106 36.2570877, -80.4353104 36.5510445, -80.61…
#>  4 Currituck   MULTIPOLYGON (((-75.9419327 36.2943382, -75.9575119 36.2594528, …
#>  5 Northampton POLYGON ((-77.1419601 36.4170647, -77.1393204 36.4564781, -77.12…
#>  6 Hertford    POLYGON ((-76.7074966 36.2661324, -76.7413483 36.3151665, -76.92…
#>  7 Camden      POLYGON ((-76.0173492 36.3377304, -76.0328751 36.3359756, -76.04…
#>  8 Gates       POLYGON ((-76.46035 36.3738976, -76.5024643 36.4522858, -76.4983…
#>  9 Warren      POLYGON ((-78.1347198 36.2365837, -78.1096268 36.2135086, -78.05…
#> 10 Stokes      POLYGON ((-80.0240555 36.5450249, -80.0480957 36.5471344, -80.43…
#> # ℹ 90 more rows
```

Use accessors to extract information about geometries:

``` r
nc_s2 %>%
  mutate(
    area = s2_area(geometry),
    perimeter = s2_perimeter(geometry)
  )
#> # A tibble: 100 × 4
#>    NAME        geometry                                           area perimeter
#>    <chr>       <s2_geography>                                    <dbl>     <dbl>
#>  1 Ashe        POLYGON ((-81.4528885 36.2395859, -81.4310379 3… 1.14e9   141627.
#>  2 Alleghany   POLYGON ((-81.1766739 36.4154434, -81.1533661 3… 6.11e8   119876.
#>  3 Surry       POLYGON ((-80.4530106 36.2570877, -80.4353104 3… 1.42e9   160458.
#>  4 Currituck   MULTIPOLYGON (((-75.9419327 36.2943382, -75.957… 6.94e8   301644.
#>  5 Northampton POLYGON ((-77.1419601 36.4170647, -77.1393204 3… 1.52e9   211794.
#>  6 Hertford    POLYGON ((-76.7074966 36.2661324, -76.7413483 3… 9.68e8   160780.
#>  7 Camden      POLYGON ((-76.0173492 36.3377304, -76.0328751 3… 6.16e8   150430.
#>  8 Gates       POLYGON ((-76.46035 36.3738976, -76.5024643 36.… 9.03e8   123170.
#>  9 Warren      POLYGON ((-78.1347198 36.2365837, -78.1096268 3… 1.18e9   141073.
#> 10 Stokes      POLYGON ((-80.0240555 36.5450249, -80.0480957 3… 1.23e9   140583.
#> # ℹ 90 more rows
```

Use predicates to subset vectors:

``` r
nc_s2 %>%
  filter(s2_contains(geometry, "POINT (-80.9313 35.6196)"))
#> # A tibble: 1 × 2
#>   NAME    geometry                                                              
#>   <chr>   <s2_geography>                                                        
#> 1 Catawba POLYGON ((-80.9312744 35.6195908, -81.0035782 35.6970558, -81.0547791…
```

Use transformers to create new geometries:

``` r
nc_s2 %>%
  mutate(geometry = s2_boundary(geometry))
#> # A tibble: 100 × 2
#>    NAME        geometry                                                         
#>    <chr>       <s2_geography>                                                   
#>  1 Ashe        LINESTRING (-81.4528885 36.2395859, -81.4310379 36.2607193, -81.…
#>  2 Alleghany   LINESTRING (-81.1766739 36.4154434, -81.1533661 36.4247398, -81.…
#>  3 Surry       LINESTRING (-80.4530106 36.2570877, -80.4353104 36.5510445, -80.…
#>  4 Currituck   MULTILINESTRING ((-75.9419327 36.2943382, -75.9575119 36.2594528…
#>  5 Northampton LINESTRING (-77.1419601 36.4170647, -77.1393204 36.4564781, -77.…
#>  6 Hertford    LINESTRING (-76.7074966 36.2661324, -76.7413483 36.3151665, -76.…
#>  7 Camden      LINESTRING (-76.0173492 36.3377304, -76.0328751 36.3359756, -76.…
#>  8 Gates       LINESTRING (-76.46035 36.3738976, -76.5024643 36.4522858, -76.49…
#>  9 Warren      LINESTRING (-78.1347198 36.2365837, -78.1096268 36.2135086, -78.…
#> 10 Stokes      LINESTRING (-80.0240555 36.5450249, -80.0480957 36.5471344, -80.…
#> # ℹ 90 more rows
```

Finally, use the WKB or WKT exporters to export to sf or some other
package:

``` r
nc_s2 %>%
  mutate(geometry = st_as_sfc(s2_as_binary(geometry))) %>%
  st_as_sf()
#> Simple feature collection with 100 features and 1 field
#> Geometry type: GEOMETRY
#> Dimension:     XY
#> Bounding box:  xmin: -84.32385 ymin: 33.88199 xmax: -75.45698 ymax: 36.58965
#> CRS:           NA
#> # A tibble: 100 × 2
#>    NAME                                                                 geometry
#>    <chr>                                                              <GEOMETRY>
#>  1 Ashe        POLYGON ((-81.45289 36.23959, -81.43104 36.26072, -81.41233 36.2…
#>  2 Alleghany   POLYGON ((-81.17667 36.41544, -81.15337 36.42474, -81.1384 36.41…
#>  3 Surry       POLYGON ((-80.45301 36.25709, -80.43531 36.55104, -80.61105 36.5…
#>  4 Currituck   MULTIPOLYGON (((-75.94193 36.29434, -75.95751 36.25945, -75.9137…
#>  5 Northampton POLYGON ((-77.14196 36.41706, -77.13932 36.45648, -77.12733 36.4…
#>  6 Hertford    POLYGON ((-76.7075 36.26613, -76.74135 36.31517, -76.92408 36.39…
#>  7 Camden      POLYGON ((-76.01735 36.33773, -76.03288 36.33598, -76.04395 36.3…
#>  8 Gates       POLYGON ((-76.46035 36.3739, -76.50246 36.45229, -76.49834 36.50…
#>  9 Warren      POLYGON ((-78.13472 36.23658, -78.10963 36.21351, -78.05835 36.2…
#> 10 Stokes      POLYGON ((-80.02406 36.54502, -80.0481 36.54713, -80.43531 36.55…
#> # ℹ 90 more rows
```

## Acknowledgment

This project gratefully acknowledges financial
[support](https://r-consortium.org/) from the

<a href="https://r-consortium.org/">
<img src="man/figures/rc300.png" width="300" /> </a>
