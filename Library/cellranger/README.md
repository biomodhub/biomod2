
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![TravisCI Build Status](https://travis-ci.org/rsheets/cellranger.svg?branch=master)](https://travis-ci.org/rsheets/cellranger) <!--[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rsheets/cellranger?branch=master&svg=true)](https://ci.appveyor.com/project/rsheets/cellranger)--> [![codecov.io](https://codecov.io/github/rsheets/cellranger/coverage.svg?branch=master)](https://codecov.io/github/rsheets/cellranger?branch=master) [![DOI](https://zenodo.org/badge/16122/jennybc/cellranger.svg)](http://dx.doi.org/10.5281/zenodo.21970) [![CRAN version](http://www.r-pkg.org/badges/version/cellranger)](https://cran.r-project.org/package=cellranger) ![](http://cranlogs.r-pkg.org/badges/grand-total/cellranger)

<img src="http://i.imgur.com/RJJy15I.jpg" width="270" align="right" />

Helper package to support R scripts or packages that interact with spreadsheets.

### Installation

Option 1: Install from CRAN:

``` r
install.packages("cellranger")
```

Option 2: Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("jennybc/cellranger")
```

### What is `cellranger` for?

**Describe a rectangle of cells**. For example, what you've got is the string "D12:F15" and what you want is an R object that holds the row and column for the upper left and lower right corners of this rectangle. Read below about the `cell_limits` class. The [`googlesheets`](https://github.com/jennybc/googlesheets) and [`readODS`](https://github.com/chainsawriot/readODS) packages use `cellranger` to translate user-supplied cell range info into something more programmatically useful.

**Handle cell references found in spreadsheet formulas**. If you're parsing unevaluated spreadsheet formulas, use the `ra_ref` and `cell_addr` classes for handling absolute, relative, and mixed cell references. Classes inspired by [Spreadsheet Implementation Technology](https://mitpress.mit.edu/books/spreadsheet-implementation-technology) from Sestoft (MIT Press, 2014).

**Convert between annoying spreadsheet reference formats**. Some utility functions are exposed, such as `A1_to_R1C1()`, which converts from A1 formatted strings to R1C1, and `letter_to_num()`, which converts a Excel column ID to a number, e.g. column AQZ is more usefully known as column 1144.

### Describing rectangles via `cell_limits`

`cellranger` provides an S3 class, `cell_limits`, as the standard way to store a cell range. You can explicitly construct a `cell_limits` object by specifying the upper left and lower right cells and, optionally, the hosting worksheet:

``` r
cell_limits(ul = c(ROW_MIN, COL_MIN), lr = c(ROW_MAX, COL_MAX), sheet = "SHEET")
```

Think of it like `R3C1:R7C4` notation, but with the `R` and `C` removed.

More often you'll get a `cell_limits` object by sending diverse user input through `as.cell_limits()`. That's what's going on in calls like these from [`googlesheets`](https://github.com/jennybc/googlesheets):

``` r
library(googlesheets)
gs_read(..., range = "D12:F15")
gs_read(..., range = "raw_data!R1C12:R6C15")
gs_read(..., range = cell_limits(c(1, 1), c(6, 15)))
gs_read(..., range = cell_limits(c(2, 1), c(NA, NA)))
gs_read(..., range = cell_rows(1:100))
gs_read(..., range = cell_cols(3:8))
gs_read(..., range = cell_cols("B:MZ"))
gs_read(..., range = anchored("B4", dim = c(2, 10)))
gs_read(..., range = anchored("A1", dim = c(5, 6), col_names = TRUE))
## internal usage in functions that put data into a googlesheet
anchored(input = head(iris))
anchored(input = head(iris), col_names = FALSE)
anchored(input = head(LETTERS))
anchored(input = head(LETTERS), byrow = TRUE)
```

Read the docs for more information on some specialized helpers:

-   Row- or column-only specification: `cell_rows()`, `cell_cols()`.
-   Specification via an object you want to write and, optionally, an anchor cell: `anchored()`

``` r
library("cellranger")
(cl <- as.cell_limits("raw_data!R1C12:R6C15"))
#> <cell_limits (1, 12) x (6, 15) in 'raw_data'>
```

The `dim` method reports dimensions of the targetted cell rectangle. `as.range()` converts a `cell_limits` object back into an Excel range.

``` r
dim(cl)
#> [1] 6 4

as.range(cl)
#> [1] "raw_data!R1C12:R6C15"

as.range(cl, fo = "A1", sheet = FALSE, strict = TRUE)
#> [1] "$L$1:$O$6"
```

Use `NA` to leave a limit unspecified, i.e. describe a degenerate rectangle

``` r
cell_limits(c(3, 2), c(7, NA))
#> <cell_limits (3, 2) x (7, -)>
```

If the maximum row or column is specified but the associated minimum is not, then it is set to 1.

``` r
cell_limits(c(NA, NA), c(3, 5))
#> <cell_limits (1, 1) x (3, 5)>
```

### Utilities for spreadsheet annoyances

We've exposed utility functions which could be useful to anyone manipulating Excel-like references.

``` r
## convert character column IDs to numbers ... and vice versa
letter_to_num(c('AA', 'ZZ', 'ABD', 'ZZZ', ''))
#> [1]    27   702   732 18278    NA

num_to_letter(c(27, 702, 732, 18278, 0, -5))
#> [1] "AA"  "ZZ"  "ABD" "ZZZ" NA    NA

## convert between A1 and R1C1 cell references
A1_to_R1C1(c("$A$1", "$AZ$10"))
#> [1] "R1C1"   "R10C52"
A1_to_R1C1(c("A1", "AZ10"), strict = FALSE)
#> [1] "R1C1"   "R10C52"

R1C1_to_A1(c("R1C1", "R10C52"))
#> [1] "$A$1"   "$AZ$10"
R1C1_to_A1(c("R1C1", "R10C52"), strict = FALSE)
#> [1] "A1"   "AZ10"

## detect cell reference formats with
## is_A1() and is_R1C1()
x <- c("A1", "$A4", "$b$12", "RC1", "R[-4]C9", "R5C3")
data.frame(x, A1 = is_A1(x), R1C1 = is_R1C1(x))
#>         x    A1  R1C1
#> 1      A1  TRUE FALSE
#> 2     $A4  TRUE FALSE
#> 3   $b$12  TRUE FALSE
#> 4     RC1  TRUE  TRUE
#> 5 R[-4]C9 FALSE  TRUE
#> 6    R5C3 FALSE  TRUE

## guess format with
## guess_fo()
refs <- c("A1", "$A1", "A$1", "$A$1", "a1",
          "R1C1", "R1C[-1]", "R[-1]C1", "R[-1]C[9]")
data.frame(refs, guessed = guess_fo(refs))
#>        refs guessed
#> 1        A1      A1
#> 2       $A1      A1
#> 3       A$1      A1
#> 4      $A$1      A1
#> 5        a1      A1
#> 6      R1C1    R1C1
#> 7   R1C[-1]    R1C1
#> 8   R[-1]C1    R1C1
#> 9 R[-1]C[9]    R1C1
```
