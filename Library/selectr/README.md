# selectr

[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![GitHub Actions](https://github.com/sjp/selectr/actions/workflows/r.yml/badge.svg)](https://github.com/sjp/selectr/actions/workflows/r.yml) [![CRAN version](https://www.r-pkg.org/badges/version/selectr)](https://cran.r-project.org/package=selectr) [![codecov](https://codecov.io/gh/sjp/selectr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/sjp/selectr) ![Downloads per month](https://cranlogs.r-pkg.org/badges/last-month/selectr)

selectr is a package which makes working with HTML and XML documents easier. It does this by performing translation of CSS selectors into XPath expressions so that you can query `XML` and `xml2` documents easily.

``` r
library(selectr)
xpath <- css_to_xpath("#selectr")
xpath
#> [1] "descendant-or-self::*[@id = 'selectr']"
```

## Installation

### Install the release version from CRAN

``` r
install.packages("selectr")
```

### Install the development version from GitHub

``` r
# install.packages("devtools")
devtools::install_github("sjp/selectr")
```

## Overview

The key functions in selectr are:

* Translate a CSS selector into an XPath expression with `css_to_xpath()`.

* Query an `XML` or `xml2` document with `querySelector()` and its variants.

    * Find the first matching node with `querySelector()`.

    * Find all matching nodes with `querySelectorAll()`.

    * Find the first matching node in a namespaced document with `querySelectorNS()`.

    * Find all matching nodes in a namespaced document with `querySelectorAllNS()`.

## Examples

Here is a simple example to demonstrate how to query an `XML` or `xml2` document with `querySelector()`.

``` r
library(selectr)
xmlText <- '<foo><bar><baz id="first"/></bar><baz id="second"/></foo>'

library(XML)
doc <- xmlParse(xmlText)
querySelector(doc, "baz")
#> <baz id="first"/>
querySelectorAll(doc, "baz")
#> [[1]]
#> <baz id="first"/>
#>
#> [[2]]
#> <baz id="second"/>
#>
#> attr(,"class")
#> [1] "XMLNodeSet"

library(xml2)
doc <- read_xml(xmlText)
querySelector(doc, "baz")
#> {xml_node}
#> <baz id="first">
querySelectorAll(doc, "baz")
#> {xml_nodeset (2)}
#> [1] <baz id="first"/>
#> [2] <baz id="second"/>
```
