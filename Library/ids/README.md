# ids

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Linux Build Status](https://travis-ci.org/richfitz/ids.svg?branch=master)](https://travis-ci.org/richfitz/ids)
[![codecov.io](https://codecov.io/github/richfitz/ids/coverage.svg?branch=master)](https://codecov.io/github/richfitz/ids?branch=master)
[![](https://www.r-pkg.org/badges/version/ids)](https://cran.r-project.org/package=ids)

Generate random identifiers in a number of styles:

* random ids of any number of bytes, such as `31f6d556fe2b303c`
* UUIDs using the `uuid` package, such as `4f0efabf-0375-4a08-89ea-b8f162f07c44`
* human readable identifiers in the style `<adjective>_<animal>` (following [gfycat.com](http://gfycat.com)), such as `misanthropic_lungfish`
* human readable identifiers in the style of a sentence (following [Asana](https://blog.asana.com/2011/09/6-sad-squid-snuggle-softly), such as `33_enormous_chinchillas_tumbling_elegantly`
* [proquint](https://arxiv.org/html/0901.4016) ("PRO-nouncable QUINT-uplet") idenfiers such as `dizuz-soboz` (which can be convereted to an integer such as 40,2638,895)

These can all be tweaked with options for length, words that are used, and the case that joins words.  There is a function `ids` for creating your own human readable identifiers.

## Installation

Install this package from CRAN

```r
install.packages("ids")
```

or install the development version with

```r
remotes::install_github("richfitz/ids", upgrade = FALSE)
```
