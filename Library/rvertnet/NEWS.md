# rvertnet 0.8.4

* remove hardcoded figure that caused failed vignette build on 2 CRAN check runners with old pandoc versions

# rvertnet 0.8.3


### MINOR IMPROVEMENTS

* new maintainer (#71)
* fix CI
* refresh test fixtures (#71)
* remove plyr from examples
* refresh Darwin core terms and move to data object
* refresh vignette and its pre-compilation
* update Makefile
* update roxygen2 and documentation fixes
* remove outdated files

### BUG FIXES

* fix class handling in vertmap to allow tibble input
* allow vignette to build
* fix ggplot2 deprecation warning (#71)

# rvertnet 0.8.2


### MINOR IMPROVEMENTS

* vignette fix

# rvertnet 0.8.0


### NEW FEATURES

* `searchbyterm()` and `bigsearch()` reworked: both functions now have the first parameter as `...`, which accepts any valid query parameter. There were so many query parameters for these functions it was a bit overwhelming. See `?searchbyterm` docs for details  (#66)

### MINOR IMPROVEMENTS

* decode the request URL before printing to the R console so users can more easily see what request they have done (#67)
* vignette title fix (#68)

### BUG FIXES 

* `searchbyterm()` fix: booleans need to be converted to VertNet's expected `0/1` instead of `true/false` (#66)


# rvertnet 0.7.0


### BUG FIXES 

* add month and day params to `searchbyterm` (#64)


# rvertnet 0.6.2


### BUG FIXES 

* A small data source used in one function was on the web, and
was moved - that data source now within the pkg as quite small, and 
now pkg won't break when the file is moved again (#61) (#62)


# rvertnet 0.6.0


Added Code of Conduct.

### NEW FEATURES

* Now using `crul` package for HTTP requests instead of `httr` (#57)
* Note that `verbose` parameter has been replaced with `messages` throughout
the package.
* Now with function for search for trait data: `traitsearch()` (#55)

### DEFUNCT AND DEPRECATED

* All `dump` functions are now defunct. Those functions tried to help users
work with bulk Vertnet data - the setup has gotten too complex (#56)

### MINOR IMPROVEMENTS

* Improvements to documentation for `traitsearch()` function on what 
fields have given data (#58) thanks @gaurav
* `vertsearch()` and `searchbyterm()` gain new parameter `only_dwc`, which 
allows to optionally only return Darwin Core fields

### BUG FIXES 

* Small fix to `vertsummary()` (#59)


# rvertnet 0.5.0


### NEW FEATURES

* `searchbyterm()` gains new parameter `query` to allow full text search, 
much like `vertsearch()`, but with the ability to also use all the parameters
available in `searchbyterm()` (#53)

### MINOR IMPROVEMENTS

* Use `dplyr::bind_rows` instead of the deprecated `dplyr::rbind_all` (#51)
* remove personal email address from tests (#52)
* Namespace base R pkg fxn calls (`methods`/`stats`/`utils`), and removed 
some package dependencies that we didn't really need (`plyr`) (#54)

# rvertnet 0.4.4


### MINOR IMPROVEMENTS

* Updated docs to better indicate how to use the cursor feature (#49)
* Now using explicit encoding specification when using `httr::content()` (#47)

### BUG FIXES

* Fixed `externalptr` error in the internal `vert_GET()` function (#48)

# rvertnet 0.4.1


### BUG FIXES

* Fixed a bug in `bigsearch()` in which we had forgotten to do
internal conversion of logical input to 0/1 needed by the web
API (#46)

# rvertnet 0.4.0


### NEW FEATURES

* New set of functions to make working with VertNet data dumps
easier. `dump_links()` gives you links to various data dump
resources; `dump_init()` initialized a SQLite database connection;
`dump_tbl()` creates a `dplyr::tbl` object, which can then be used
in a `dplyr` query. This setup requires that the user manually
download data dumps uncompress, and load into SQLite. We hope to
make this process easier in the future. (#36)

### MINOR IMPROVEMENTS

* Fixes to `vertmap()` for new `ggplot2` version (#43)
* Added note to docs for `bigsearch()` for how to read in data
after obtaining the data (#44)

### BUG FIXES

* Fix to the `searchbyterm()` function. When the parameter `stateprovince`
was used, lead to error, as that param requires different handling than
other params. (#45)

# rvertnet 0.3.4


### NEW FEATURES

* New function `vert_id()` to get occurrence records by occurenceid,
that is, single occurrence ids. (#40)

### MINOR IMPROVEMENTS

* Explicitly import non-base R functions (#39)

### BUG FIXES

* Lowercase `occurenceID` to `occurrenceid` to simplify life (#41)

# rvertnet 0.3.0


### NEW FEATURES

* `searchbyterm()` and `bigsearch()` have some parameters that accept multiple values.
Fixed to allow this (#37)
* Internals of `searchbyterm()`, `spatialsearch()`, and `vertsearch()` reworked to
use cursor so we internally do paging for you for bigger result sets. (#25)

### MINOR IMPROVEMENTS

* Replaced `data.table` import with `dplyr`
* Using `skip_on_cran()` (#38)
* Minor vignette updates (#35)
* Metadata now returned in data requests (#33)

# rvertnet 0.2.2


Package completely reworked for the new VertNet API.

### NEW FEATURES

* The functions `vertavailablemaps()`, `vertlocations()`,
`vertoccurrence()`, `vertoccurrencecount()`, `vertproviders()`,
`verttaxa()` are now defunct. You can call these functions, but
they print an error message, saying they are defunct.
* Gained new functions `bigsearch()`, `searchbyterm()`,
`spatialsearch()`, and `vertsummary()`.
* Gained new author: Chris Ray

### MINOR IMPROVEMENTS

* `RJSONIO` replaced with `jsonlite`
* Changed from CC0 to MIT license

# rvertnet 0.0-5

### NEW FEATURES

* released to CRAN
