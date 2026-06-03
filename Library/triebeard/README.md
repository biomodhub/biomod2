## triebeard

Fast key-value matching in R and Rcpp

__Author:__ Os Keyes, Drew Schmidt, Yuuki Takano<br/>
__License:__ [MIT](https://opensource.org/license/mit/)<br/>
__Status:__ Stable

[![Travis-CI Build Status](https://travis-ci.org/Ironholds/triebeard.svg?branch=master)](https://travis-ci.org/Ironholds/triebeard) ![downloads](http://cranlogs.r-pkg.org/badges/grand-total/triebeard)

### Description

Tries, or [radix trees](https://en.wikipedia.org/wiki/Radix_tree), are key-value data structures optimised for very, very fast matching of the keys against user-provided data (and then the return of the associated values!)

This is pretty useful in data cleaning and value extraction, and tries let you do it *really* efficiently. `triebeard` contains
an implementation that can be used both when writing R, and when writing Rcpp (and imported and linked against, to boot). For more information see:

1. The [vignette on Rcpp usage](https://CRAN.R-project.org/package=triebeard/vignettes/rcpp_radix.html);
2. The [vignette on R usage](https://CRAN.R-project.org/package=triebeard/vignettes/r_radix.html).

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/Ironholds/triebeard/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.

### Installation

The stable, CRAN-ready version can be retrieved with:

    install.packages("triebeard")

The latest version can be obtained via:

    devtools::install_github("ironholds/triebeard")

### Dependencies
* R.
* [Rcpp](https://cran.r-project.org/package=Rcpp)
