


# rematch2

> Match Regular Expressions with a Nicer 'API'

[![Linux Build Status](https://travis-ci.org/r-lib/rematch2.svg?branch=master)](https://travis-ci.org/r-lib/rematch2)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/github/r-lib/rematch2?svg=true)](https://ci.appveyor.com/project/gaborcsardi/rematch2)
[![](http://www.r-pkg.org/badges/version/rematch2)](http://www.r-pkg.org/pkg/rematch2)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/rematch2)](http://www.r-pkg.org/pkg/rematch2)
[![Coverage Status](https://img.shields.io/codecov/c/github/r-lib/rematch2/master.svg)](https://codecov.io/github/r-lib/rematch2?branch=master)

A small wrapper on regular expression matching functions `regexpr`
and `gregexpr` to return the results in tidy data frames.

---

  - [Installation](#installation)
  - [Rematch vs rematch2](#rematch-vs-rematch2)
  - [Usage](#usage)
    - [First match](#first-match)
    - [All matches](#all-matches)
    - [Match positions](#match-positions)
  - [License](#license)

## Installation


```r
install.packages("rematch2")
```

## Rematch vs rematch2

Note that `rematch2` is not compatible with the original `rematch` package.
There are at least three major changes:
* The order of the arguments for the functions is different. In
  `rematch2` the `text` vector is first, and `pattern` is second.
* In the result, `.match` is the last column instead of the first.
* `rematch2` returns `tibble` data frames. See
  https://github.com/hadley/tibble.

## Usage

### First match


```r
library(rematch2)
```

With capture groups:

```r
dates <- c("2016-04-20", "1977-08-08", "not a date", "2016",
  "76-03-02", "2012-06-30", "2015-01-21 19:58")
isodate <- "([0-9]{4})-([0-1][0-9])-([0-3][0-9])"
re_match(text = dates, pattern = isodate)
```

```
#> # A tibble: 7 x 5
#>      ``    ``    ``            .text     .match
#>   <chr> <chr> <chr>            <chr>      <chr>
#> 1  2016    04    20       2016-04-20 2016-04-20
#> 2  1977    08    08       1977-08-08 1977-08-08
#> 3  <NA>  <NA>  <NA>       not a date       <NA>
#> 4  <NA>  <NA>  <NA>             2016       <NA>
#> 5  <NA>  <NA>  <NA>         76-03-02       <NA>
#> 6  2012    06    30       2012-06-30 2012-06-30
#> 7  2015    01    21 2015-01-21 19:58 2015-01-21
```

Named capture groups:

```r
isodaten <- "(?<year>[0-9]{4})-(?<month>[0-1][0-9])-(?<day>[0-3][0-9])"
re_match(text = dates, pattern = isodaten)
```

```
#> # A tibble: 7 x 5
#>    year month   day            .text     .match
#>   <chr> <chr> <chr>            <chr>      <chr>
#> 1  2016    04    20       2016-04-20 2016-04-20
#> 2  1977    08    08       1977-08-08 1977-08-08
#> 3  <NA>  <NA>  <NA>       not a date       <NA>
#> 4  <NA>  <NA>  <NA>             2016       <NA>
#> 5  <NA>  <NA>  <NA>         76-03-02       <NA>
#> 6  2012    06    30       2012-06-30 2012-06-30
#> 7  2015    01    21 2015-01-21 19:58 2015-01-21
```

A slightly more complex example:

```r
github_repos <- c(
	"metacran/crandb",
	"jeroenooms/curl@v0.9.3",
    "jimhester/covr#47",
	"hadley/dplyr@*release",
    "r-lib/remotes@550a3c7d3f9e1493a2ba",
    "/$&@R64&3"
)
owner_rx   <- "(?:(?<owner>[^/]+)/)?"
repo_rx    <- "(?<repo>[^/@#]+)"
subdir_rx  <- "(?:/(?<subdir>[^@#]*[^@#/]))?"
ref_rx     <- "(?:@(?<ref>[^*].*))"
pull_rx    <- "(?:#(?<pull>[0-9]+))"
release_rx <- "(?:@(?<release>[*]release))"

subtype_rx <- sprintf("(?:%s|%s|%s)?", ref_rx, pull_rx, release_rx)
github_rx  <- sprintf(
	"^(?:%s%s%s%s|(?<catchall>.*))$",
    owner_rx, repo_rx, subdir_rx, subtype_rx
)
re_match(text = github_repos, pattern = github_rx)
```

```
#> # A tibble: 6 x 9
#>        owner    repo subdir                  ref  pull  release  catchall
#>        <chr>   <chr>  <chr>                <chr> <chr>    <chr>     <chr>
#> 1   metacran  crandb                                                     
#> 2 jeroenooms    curl                      v0.9.3                         
#> 3  jimhester    covr                                47                   
#> 4     hadley   dplyr                                   *release          
#> 5      r-lib remotes        550a3c7d3f9e1493a2ba                         
#> 6                                                               /$&@R64&3
#> # ... with 2 more variables: .text <chr>, .match <chr>
```

### All matches

Extract all names, and also first names and last names:


```r
name_rex <- paste0(
  "(?<first>[[:upper:]][[:lower:]]+) ",
  "(?<last>[[:upper:]][[:lower:]]+)"
)
notables <- c(
  "  Ben Franklin and Jefferson Davis",
  "\tMillard Fillmore"
)
not <- re_match_all(notables, name_rex)
not
```

```
#> # A tibble: 2 x 4
#>       first      last                              .text    .match
#>      <list>    <list>                              <chr>    <list>
#> 1 <chr [2]> <chr [2]>   Ben Franklin and Jefferson Davis <chr [2]>
#> 2 <chr [1]> <chr [1]>               "\tMillard Fillmore" <chr [1]>
```


```r
not$first
```

```
#> [[1]]
#> [1] "Ben"       "Jefferson"
#> 
#> [[2]]
#> [1] "Millard"
```

```r
not$last
```

```
#> [[1]]
#> [1] "Franklin" "Davis"   
#> 
#> [[2]]
#> [1] "Fillmore"
```

```r
not$.match
```

```
#> [[1]]
#> [1] "Ben Franklin"    "Jefferson Davis"
#> 
#> [[2]]
#> [1] "Millard Fillmore"
```

### Match positions

`re_exec` and `re_exec_all` are similar to `re_match` and `re_match_all`,
but they also return match positions. These functions return match
records. A match record has three components: `match`, `start`, `end`, and
each component can be a vector. It is similar to a data frame in this
respect.


```r
pos <- re_exec(notables, name_rex)
pos
```

```
#> # A tibble: 2 x 4
#>        first       last                              .text     .match
#> *     <list>     <list>                              <chr>     <list>
#> 1 <list [3]> <list [3]>   Ben Franklin and Jefferson Davis <list [3]>
#> 2 <list [3]> <list [3]>               "\tMillard Fillmore" <list [3]>
```

Unfortunately R does not allow hierarchical data frames (i.e. a column of a
data frame cannot be another data frame), but `rematch2` defines some
special classes and an `$` operator, to make it easier to extract parts
of `re_exec` and `re_exec_all` matches. You simply query the `match`,
`start` or `end` part of a column:


```r
pos$first$match
```

```
#> [1] "Ben"     "Millard"
```

```r
pos$first$start
```

```
#> [1] 3 2
```

```r
pos$first$end
```

```
#> [1] 5 8
```

`re_exec_all` is very similar, but these queries return lists, with
arbitrary number of matches:


```r
allpos <- re_exec_all(notables, name_rex)
allpos
```

```
#> # A tibble: 2 x 4
#>        first       last                              .text     .match
#>       <list>     <list>                              <chr>     <list>
#> 1 <list [3]> <list [3]>   Ben Franklin and Jefferson Davis <list [3]>
#> 2 <list [3]> <list [3]>               "\tMillard Fillmore" <list [3]>
```


```r
allpos$first$match
```

```
#> [[1]]
#> [1] "Ben"       "Jefferson"
#> 
#> [[2]]
#> [1] "Millard"
```

```r
allpos$first$start
```

```
#> [[1]]
#> [1]  3 20
#> 
#> [[2]]
#> [1] 2
```

```r
allpos$first$end
```

```
#> [[1]]
#> [1]  5 28
#> 
#> [[2]]
#> [1] 8
```

## License

MIT © Mango Solutions, Gábor Csárdi
