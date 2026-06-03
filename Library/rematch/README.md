
# rematch

> Match Regular Expressions with a Nicer ‘API’

[![R-CMD-check](https://github.com/gaborcsardi/rematch/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gaborcsardi/rematch/actions/workflows/R-CMD-check.yaml)
[![](https://www.r-pkg.org/badges/version/rematch)](https://www.r-pkg.org/pkg/rematch)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/rematch)](https://www.r-pkg.org/pkg/rematch)
[![Coverage
Status](https://img.shields.io/codecov/c/github/gaborcsardi/rematch/main.svg)](https://app.codecov.io/github/gaborcsardi/rematch?branch=main)

A small wrapper on ‘regexpr’ to extract the matches and captured groups
from the match of a regular expression to a character vector.

## Installation

``` r
source("https://install-github.me/gaborcsardi/rematch")
```

## Usage

``` r
library(rematch)
```

``` r
dates <- c("2016-04-20", "1977-08-08", "not a date", "2016",
  "76-03-02", "2012-06-30", "2015-01-21 19:58")
isodate <- "([0-9]{4})-([0-1][0-9])-([0-3][0-9])"
re_match(text = dates, pattern = isodate)
```

    #>      .match                       
    #> [1,] "2016-04-20" "2016" "04" "20"
    #> [2,] "1977-08-08" "1977" "08" "08"
    #> [3,] NA           NA     NA   NA  
    #> [4,] NA           NA     NA   NA  
    #> [5,] NA           NA     NA   NA  
    #> [6,] "2012-06-30" "2012" "06" "30"
    #> [7,] "2015-01-21" "2015" "01" "21"

``` r
isodaten <- "(?<year>[0-9]{4})-(?<month>[0-1][0-9])-(?<day>[0-3][0-9])"
re_match(text = dates, pattern = isodaten)
```

    #>      .match       year   month day 
    #> [1,] "2016-04-20" "2016" "04"  "20"
    #> [2,] "1977-08-08" "1977" "08"  "08"
    #> [3,] NA           NA     NA    NA  
    #> [4,] NA           NA     NA    NA  
    #> [5,] NA           NA     NA    NA  
    #> [6,] "2012-06-30" "2012" "06"  "30"
    #> [7,] "2015-01-21" "2015" "01"  "21"

``` r
github_repos <- c("metacran/crandb", "jeroenooms/curl@v0.9.3",
                  "jimhester/covr#47", "hadley/dplyr@*release",
                  "r-lib/remotes@550a3c7d3f9e1493a2ba",
                  "/$&@R64&3")
owner_rx <- "(?:(?<owner>[^/]+)/)?"
repo_rx <- "(?<repo>[^/@#]+)"
subdir_rx <- "(?:/(?<subdir>[^@#]*[^@#/]))?"
ref_rx <- "(?:@(?<ref>[^*].*))"
pull_rx <- "(?:#(?<pull>[0-9]+))"
release_rx <- "(?:@(?<release>[*]release))"
ref_or_pull_or_release_rx <-
  sprintf("(?:%s|%s|%s)?", ref_rx, pull_rx, release_rx)
github_rx <- sprintf("^(?:%s%s%s%s|(?<catchall>.*))$",
                     owner_rx, repo_rx, subdir_rx, ref_or_pull_or_release_rx)
out <- re_match(text = github_repos, pattern = github_rx)
out
```

    #>      .match                               owner        repo      subdir
    #> [1,] "metacran/crandb"                    "metacran"   "crandb"  ""    
    #> [2,] "jeroenooms/curl@v0.9.3"             "jeroenooms" "curl"    ""    
    #> [3,] "jimhester/covr#47"                  "jimhester"  "covr"    ""    
    #> [4,] "hadley/dplyr@*release"              "hadley"     "dplyr"   ""    
    #> [5,] "r-lib/remotes@550a3c7d3f9e1493a2ba" "r-lib"      "remotes" ""    
    #> [6,] "/$&@R64&3"                          ""           ""        ""    
    #>      ref                    pull release    catchall   
    #> [1,] ""                     ""   ""         ""         
    #> [2,] "v0.9.3"               ""   ""         ""         
    #> [3,] ""                     "47" ""         ""         
    #> [4,] ""                     ""   "*release" ""         
    #> [5,] "550a3c7d3f9e1493a2ba" ""   ""         ""         
    #> [6,] ""                     ""   ""         "/$&@R64&3"

## License

MIT © Mango Solutions; Posit Software, PBC
