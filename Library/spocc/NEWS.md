spocc 1.2.4
===========

### Documentation

* Link to `ridigbio::idig_search_records()` fixed for `as.idigbio()`. 

spocc 1.2.3
===========

### Documentation

* Documented argument ‘scientificname’ in 'obis_search()' changed to 'scientificName'.
* Documented argument 'n' in 'spocc_objects()' removed

### DEFUNCT

* `datecollected` field flagged by `ridigbio` maintainers as inappropriate field--should be `eventDate`. See `scrubr::fix_names()` (#263)


spocc 1.2.2
===========

### DEFUNCT

* Package `wellknown` about to be archived. Replaced with `s2` and `wk`.

spocc 1.2.1
===========

### DEFUNCT

* `as.bison()` is now defunct. BISON API is being sunsetted,  BISON data now available
via GBIF (#252)

### Documentation

* Maintainer changed to Hannah L. Owens


spocc 1.2.0
===========

### DEFUNCT

* `as.ecoengine()` is now defunct. `ecoengine` R pkg wasn't used, we used internal code, but since `ecoengine` R is archived since Feb 2020, probably best to remove here (#239)

### MINOR IMPROVEMENTS

* using `wellknown` now instead of `wicket` for Well-known text manipulation (#235)

### BUG FIXES

* fix failing checks (#240)
* bug fixed in ecoengine bbox - but now ecoengine removed, see above (#234)


spocc 1.1.0
===========

### DEFUNCT

* `fixnames()` is now defunct. it was deprecated in a previous version. See `scrubr::fix_names()` (#231)

### NEW FEATURES

* `occ()` can now handle sf objects passed to `geometry`. `spocc` itself does not import/suggest sf, but uses some code donated by Michael Sumner to pull out well known text (WKT) needed for spatially defined queries

### MINOR IMPROVEMENTS

* refactor `occ()`: factor out functions already defined inside of `occ`, add assertions for user parameters (#228)
* package logo/sticker added (#188)

### BUG FIXES

* Fix for ALA data source in `occ()`: total records found count was always 0 because ALA changed the records found field to `totalRecords`
* Fix for Vertnet data source in `occ()`: was using an old parameter `query` passed to `rvertnet::searchbyterm()` - changed to `scientificname` instead


spocc 1.0.8
===========

### BUG FIXES

* fix tests that are failing on cran checks (#225)
* fix ecoengine data option: an `if` statement was failing because we were trying to access an element of a list that is not there sometimes, leading to `NULL` which caused the `if` statement to fail


spocc 1.0.2
===========

### BUG FIXES

* preserve exact bytes for some tests that are failing on cran checks, taxize integration tests, and identifier based search tests (#221)


spocc 1.0.0
===========

### NEW FEATURES

* `source = "inat"` can now return photos. do a query as normal for inat data, and index to the `photos` slot in the data.frame, that will give a nested list of data.frames for each record with links/metadata for the photos (#214) (#217)

### MINOR IMPROVEMENTS

* tests now using vcr (#209)
* add notes in section `iNaturalist notes` in the `?occ` manual page for help with iNaturalist pagination and rate limiting (#215)
* vignette title fix for docs pages (#200)
* `tibble::data_frame`/`tibble::as_data_frame` replaced by `tibble::tibble`/`tibble::as_tibble` throughout package

### DEPRECATED

* `fixnames()` is now deprecated; still useable here until the next version released; please move to using `scrubr::fix_names` (#196)

### BUG FIXES

* fix inat data source for `occ()` queries: change from http to https for the inat base url (#213)
* inat fixes: use rbind fill approach for combining rows of data to fill missing columns safely; work with newer version of their API; unnest lat/lon results into tidy column in resulting data.frame (#215)
* OBIS API changed; changed internals for OBIS data in line with the new API; note that pagination for OBIS has changed (see `?occ` for details); `as.obis.numeric` is gone and replaced with `as.obis.character` (#218)
* fix to `fixnames()`, coerce taxon name to character in case the name is factor class (#211)


spocc 0.9.0
===========

### NEW FEATURES

* `occ()` now attempts to collect errors from requests that fail and puts these error messages (character strings) in the `$meta$errors` spot. We can not always collect errors, and some data providers do not error well: they do not provide a meaningful error message other than that there was an error. (#189) (#207)
* `occ()` gains new parameter `throw_warnings` (logical). By default set to `TRUE` (matches previous behavior) and throws warnings about errors that occur and when no results found for a query. We now prefix each warning with the data provider so you can match up an error or warning with a data provider and (hopefully) query. If set to `FALSE`, warnings are suppressed  (#189) (#207)

### DEFUNCT

* AntWeb has been removed from `spocc`. The AntWeb API has been down for a while, and no response from maintainers (#202) (#203)

### BUG FIXES

* fixes to ebird internals - a new version of `rebird` on CRAN requires a few changes in parameters used. Importantly, ebird now wants species codes instead of full scientific names, but we internally attempt to handle this, so users still just pass scientific names (#205)

### DOCUMENTATION

* make pkgdown docs better: organize functions into meaningful sets (#193) (#197) (#199)
* in the `?spocc_duplicates` manual file for duplicate records, refer to `scrubr` and `CoordinateCleaner` packages (#198)
* in `inspect()` manual file, clarify what the function does (#194)
* now we document better when we use one or the other function for BISON data source (#204)
* `occ()` gains a `return` block with detail about what's returned from the function (#208)



spocc 0.8.0
===========

### NEW FEATURES

* `occ()` gains new parameter `date` to do date range based searches across data sources without having to know the vagaries of each data source (#181)

### BUG FIXES

* fix to idigbio geometry queries (#180)
* fix `print.occdatind` so that empty data.frame's don't throw tibble warnings (#184)
* fix to internal method for standardizing dates `stand_dates()` due to ALA giving back a timestamp now (#182) (#185)
* vertnet fixes (#179)
* fix to geometry bounding box queries (#187) thanks @timcdlucas
* fix to output of list names for gbif data source when using taxonomic IDs, was resulting in booleans, should be the taxonomic IDs  (#191)


spocc 0.7.0
===========

### NEW FEATURES

* Removed javascript and V8 package import and using
`wicket` C++ based package instead. So you no longer need V8
which should make installation easier on some platforms. (#172)

### MINOR IMPROVEMENTS

* `httr` replaced with `crul` for HTTP reqeusts (#174)
* Moved to using markdown for docs. The only thing you should 
notice that's different now is doing curl options is slightly
different - it's just `curl::curl_options()` (#176)
* All `as.*()` functions can now pass on curl options to the 
http client (#177)
* Bumped minimum versions for a number of dependencies

### BUG FIXES

* Fix to `foo_ala()` - the internal plugin for `occ()` that 
handles ALA queries: changed query from full text query using 
`q=foo bar` to `q=taxon_name="foo bar"` - in addition, improved
error handling as sometimes `occurrences` slot is returned in 
results but is empty, whereas before it seemd to always be 
absent if no results (#178)


spocc 0.6.0
===========

### NEW FEATURES

* Added a new data source: Atlas of Living Australia (ALA), under
the abbreviation `ala` (#98)
* Added a new data source: Ocean Biogeographic Information System (OBIS), 
under the abbreviation `obis` (#155)

### MINOR IMPROVEMENTS

* Added note to docs and minor tweak to internal methods to account
for max results from iDigBio of 100,000. Now when you request more than 
100K, you should get a warning saying as much (#169)

### BUG FIXES

* Made `occ2df()` more robust to varied inputs - allowing for users
that may on purpose or not have a subset of the data source slots
normally in the `occdat` class object (#171)


spocc 0.5.4
===========

### MINOR IMPROVEMENTS

* `rvertnet`, a dependency dealing with data from Vertnet, was failing 
on certain searches. `rvertnet` was fixed and a new version on CRAN now. 
No changes here other than requiring the new version of `rvertnet` (#168)
* Fix to internal INAT parsers to handle JSON data output instead of 
CSV output. And fix to internal date parsing; INAT changed field for date 
from `datetime` to `observed_on`.
* Move all `is()` to `inherits()`, and namespace all `setNames()` calls
* We are now using `rgbif::occ_data()` instead of `rgbif::occ_search()`
* We are now using `rvertnet::searchbyterm()` instead of 
`rgbif::vertsearch()`

### BUG FIXES

* Fixes to iDigBio internal plugin - we were dropping scientificname
if geometry was passed by the user. Fixed now. (#167)
* Fixed bug in GBIF internal plugin - when more than 1 result given back
(e.g., multiple searches were done, resulting in a list of objects)
we weren't parsing the output correctly. Fixed now. (#166)



spocc 0.5.0
===========

### NEW FEATURES

* `occ()` now allows queries that only pass `from` and one of the data
source opts params (e.g., `gbifopts`) - allows specifying any options
passed down to the internal functions used to do data queries without
having to use the other params in `occ` (#163)

### MINOR IMPROVEMENTS

* Now using `tibble` for representing data.frames (#164)
* Now using explicit `encoding="UTF-8"` in `httr::content()` calls 
to parse raw data from web requests (#160)
* Now using `ridigbio` as its on CRAN - was using 
internal fxns prior to this (#154)

### BUG FIXES

* There was a problem in the ebird parser where it wasn't processing 
results from ebird with no data. A problem with `has_coords` also 
fixed. (#161)


spocc 0.4.5
===========

### MINOR IMPROVEMENTS

* Using `data.table::setDF()` instead of `data.frame()` to set a `data.table` 
style table to a `data.frame`
* Added many more tests to make it less likely errors will occur
* Added `vertnet` as an option to `occ_options()` to get the options for passing
to `vertopts` in `occ()`

### BUG FIXES

* Fix to `print.occdatind()` - which in last version introduced a bug in this
print method - wasn't fatal as only applied to empty slots in the output 
of a call to `occ()`, but nonetheless, not good (#159)


spocc 0.4.4
===========

### MINOR IMPROVEMENTS

* New import `data.table` for fast list to data.frame

### BUG FIXES

* Fix to ecoengine spatial search - internally we were not making the 
bounding box correctly - fixed now (#158)


spocc 0.4.0
===========

### NEW FEATURES

* New function `as.vertnet()` to coerce various inputs (e.g., result from `occ()`, `occ2df()`, or a key itself) to occurrence data objects (#142)
* `occ()` gains two parameters `start` and `page` to facilitate paging 
through results across data sources, instead of having to page 
individually for each data source. Some sources use the `start` parameter, 
while others use the `page` parameter. See __Paging__ section in `?occ` for
details on Paging (#140)

### MINOR IMPROVEMENTS

* Added Code of Conduct

### BUG FIXES

* `wkt_vis()` now works with WKT polygons with multipe polygons, e.g., 
`spocc::wkt_vis("POLYGON((-125 38.4, -121.8 38.4, -121.8 40.9, -125 40.9, -125 38.4), (-115 22.4, -111.8 22.4, -111.8 30.9, -115 30.9, -115 22.4))")` (#147)
* Fix to `print.occdatind()` to print more helpful info when a 
geometry search is used as opposed to a taxonomy based search (#149)
* Fix to `print.occdatind()` to not fail when first element not present; 
proceeds to next slot with data (#143)
* Fixed problem where `occ()` failed when multiple `geometry` elements
passed in along with taxonomic names (#146)
* Fix to `occ2df()` for combining outputs to not fail when AntWeb 
doesn't give back dates (#144) (#145) - thanks @timcdlucas
* Fix to `occ2df()` to not fail when date field missing (#141)


spocc 0.3.2
===========

### NEW FEATURES

* Added iDigBio as a new data source in `spocc` (#136) (#124)

### MINOR IMPROVEMENTS

* Added much more detail on what parameters in child packages are being used inside of the `occ()` function. Each data source is taken care of in a separate package or set of wrapper functions, and the man file now details what API parameters are being queried (#138)

### BUG FIXES

* Fixed bug where when latitude/longitude columns missing, caused problems downstream in printing outputs, etc. Now we put in NA's when those columns missing (#139)
* Fixed bug in inat data source - `Datetime` variable changed to `datetime`
* Fixed bug in vertnet data source - `occurrenceID` variable changed to `occurrenceid`

spocc 0.3.0
===========

### NEW FEATURES

* Mapping functions all gone, and put into a new package `spoccutils` (#132)
* `occ()` gains new parameter `has_coords` - a global parameter (except for ebird and bison) to return only records with lat/long data. (#128)
* `type` (#134) and `rank` (#133) parameters dropped from `occ()` 
* When object returned by `occ()` is printed, we now include a message that total count of records found (not returned) is not completely known __if ebird is included__, because eBird does not include data on records found on their servers with requests to their API (#111)
* New functions `as.*()` (e.g., `as.gbif`) for most data sources. These functions take in occurrence keys or sets of keys, and retrieve detailed occurrence record data for each key (#112)
* New data source: VertNet (#110)
* `occ2df()` now returns more fields. This function collapses all essential fields that are easy to get in all data sources: `name`, `lat`, `long`, `prov`, `date`, `key`. The `key` field is the occurrence key for each record, which you can use to keep track of individual records, get more data on the record, etc. (#103) (#108)
* New function `inspect()` - takes output from `occ()` or individual occurrence keys and gets detailed occurrence data. 

### MINOR IMPROVEMENTS

* Now importing packages: `jsonlite`, `V8`, `utils`, and `methods`. No longer importing: `ggmap`, `maptools`, `rworldmap`, `sp`, `rgeos`, `RColorBrewer`, `rgdal`, and `leafletR`. Pkgs removed mostly due to splitting off some functionality into `spoccutils`. related issues: (#131) (#132)
* Now importing explicitly all non-base R functions that we use: now importing `methods`, `utils` (#120)
* We now attempt to standardize dates across all data sources, and return that in the output of a call to `occ2df()` (#106)
* `wkt_vis()` now only has an option to view a WKT shape in the browser.

### BUG FIXES

* Fixes to being able to pass curl options on to each data source's functions (#107)

spocc 0.2.4
===========

### MINOR IMPROVEMENTS

* Improved documentation for bounding boxes, their expected format, etc. (#96)
* Remove dependency on the following packages: `assertthat`, `plyr`, `data.table`, and `XML` (#102)
* Using package `gistr` now to post interactive geojson maps on GitHub gists (#100)
* `rgbif` now must be `v0.7.7` or greater (the latest version on CRAN).
* Removed the startup message.

### BUG FIXES

* Duplicate, but not working correctly, function `occ2sp()` removed. The function `occ_to_sp()` function is the working version. (#97)
* Fixed bug where some records returned form GBIF did not have lat/long column headers, and we internally rearranged columns, which caused complete stop when that happened. Fixed now. (#101)
* Changed all `\donttest` to `\dontrun` in examples as requested by CRAN maintainers (#99)

spocc 0.2.2
===========

### NEW FEATURES

* Added new function `occ_names()` to search only for taxonomic names. The goal here is to use ths function if there is some question about what names you want to use to search for occurrences with. (#84). Suggested by @jarioksa
* New function `occ_names_options()` to quickly get parameter options to pass to `occ_names()`.
* New `summary()` method for the `occdat` `S3` object that is output from `occ()` (#83)
* In many places in `spocc` (README, vignette, `occ()` documentation file, at package startup), we make it clear that there could be duplicate records returned in certain scenarios. And a new documentation page detailing what to watch out for: `?spocc_duplicates`. (#77)

### MINOR IMPROVEMENTS

* All latitude/longitude column headers are now changed to latitude and longitude, whereas they use to vary from `latitude`, `decimalLatitude`, `Latitude`, `lat`, and `decimal_latitude`. (#91)
* Default is 500 now for the `limit` parameter in `occ()` (#78)
* You can now pass in `limit` to each functions options parameter, and it will work. Each data source can have a different parameter internally from `limit`, but now internally within `spocc`, we allow you to use `limit` so you don't have to know what the data source specific parameter is. (#81)
* There is a now a startup message to give information on the package (#79)
* `occ_options()` gains new parameter `where` to print either in the console or to open man file in the IDE, or prints to console in command line R. 

spocc 0.2.0
===========

### NEW FEATURES

* `occ()` gains new parameter `callopts` to pass on curl debugging options to `httr::GET()` (#35)
* `wkt_vis()` now by default plots a well known text area (WKT) on an interactive mapbox map in your default browser. New parameter `which` allows you to choose the interactive map or a static ggplot2 map. (#70)
* Individual data sources `occ()` gains new class. In the previous version of this package, a `data.frame` was printed. Now the data is assigned the object `occdatind` (short for _occdat individual_).
* `occ()` now uses a print method for the `occdatind` class, adopted from `dplyr` that prints a brief `data.frame`, with columns wrapped to fit the width of your console, and additional columns not printed given at bottom with their class type. Note that the print behavior for the resulting object of an `occ()` call remains the same. (#69) (#74)

### MINOR IMPROVEMENTS

* Added `whisker` as a package import to use in the `wkt_vis()` function. (#70)
* Mapping functions now all accept the same input. Previously `mapggplot()` accepted the output of `occ()`, of class `occdat`, while the other two functions for mapping, `mapleaflet()` and `mapgist()` accepted a `data.frame`. Now all three functions accept  the output of `occ()`, an object of class `occdat`. (#75)
* The `meta` slot in each returned object (indexed by `object$meta`) contains spots for `returned` and `found`, to designate number of records returned, and number of records found. (#64)

### BUG FIXES

* Fixed bug in AntWeb output, where there was supposed to be a column titled `name`. (#71)

spocc 0.1.4
===========

### NEW FEATURES

* Can now do geometry only queries. See examples in `occ()`.
* In addition, you can pass in sp objects of SpatialPolygons or SpatialPolygonsDataFrame classes.

spocc 0.1.2
===========

### NEW FEATURES

* There were quite a few changes in one of the key packages that `spocc` depends on: `rgbif`. A number of input and output parameter names changed. A new version of `rgbif` was pushed to CRAN. (#56)
* New function `clean_spocc()` started (not finished yet) to attempt to clean data. For example, one use case is removing impossible lat/long values (i.e., longitue values greater than absolute 180). Another, not implemented yet, is to remove points that are not in the country or habitat your points are supposed to be in. (#44)
* New function `fixnames()` to trim species names with optional input parameters to make data easier to use for mapping.
* New function `wkt_vis()` to visualize a WKT (well-known text) area on a map. Uses `ggmap` to pull down a Google map so that the visualization has some geographic and natural earth context. We'll soon introduce an interactive version of this function that will bring up a small Shiny app to draw a WKT area, then return those coordinates to your R session. (#34)

### MINOR IMPROVEMENTS

* Added a CONTRIBUTING.md file to the github repo to help guide contributions (#61)
* Packages that require a certain version are forced to be X version or greater. Thes are rinat (>= 0.1.1), rbison (>= 0.3.2), rgbif (>= 0.6.2), ecoengine (>= 1.3), rebird (>= 0.1.1), AntWeb (>= 0.6.1), and leafletR (>= 0.2-0). This should help avoid problems.
* General improvement to function documentation.

spocc 0.1.0
===========

* Initial release to CRAN
