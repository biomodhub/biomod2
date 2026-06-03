oai 0.4.0
=========

### NEW FEATURES

* the requests are now made using `httr::RETRY()` rather than `httr::GET()` to
  facilitate retrying with delay control etc. Consult further the documentation
  of 'httr' package (#64)
* @mbojan takes over package maintenance from @sckott

### MINOR IMPROVEMENTS

* update packages test to catch the 'cannotDisseminateFormat' OAI-PMH-level
  error correctly
* added package hex logo

oai 0.3.2
=========

### MINOR IMPROVEMENTS

* vignette fix, use markdown in suggests
* update readme, use ropensci coc

oai 0.3.0
=========

### NEW FEATURES

* `id()` gains `as` parameter so user can ask for different outputs (parsed list/data.frame, or raw xml/text) (#54)
* most OAI functions changed their default url to `http://api.gbif.org/v1/oai-pmh/registry`, while `count_identifiers()` changed it's default url to `http://export.arxiv.org/oai2`. the previous default url for Datacite was too unreliable (was often unresponsive)

### MINOR IMPROVEMENTS

* add grant information for one author (#48) (#49)
* code of conduct urls fixed
* now using markdown supported documentation (#56)
* replace `tibble::as_data_frame` with `tibble::as_tibble` throughout package

### BUG FIXES

* fix to `update_providers()`; html page that we scrape had changed (#57)


oai 0.2.2
=========

### NEW FEATURES

* Added new parsers in `get_records()` specific to different OAI prefixes. Currently has
parsers for `oai_dc` and `oai_datacite`. For prefixes we don't have
parsers for we return raw XML so you can parse it yourself. (#45)

### MINOR IMPROVEMENTS

* Replace `xml2::xml_find_one()` with `xml2::xml_find_first()` (#39)
* Update URLs in `DESCRIPTION` file (#43)
* Using `tibble` now for compact data.frame instead of internal
methods for the same (#44)
* `as` parameter in `get_records()` now has options `parsed` or `raw`, which
replaces `df` `list`, or `raw`


oai 0.2.0
=========

### NEW FEATURES

* A set of new functions for dealing with larger data results:
`dump_raw_to_txt()`, `dump_to_rds()`, and `dump_raw_to_db()`.
They can be used with `oai` functions `list_identifiers()`, `list_sets()`,
and `list_records()` (#9) (#15) (#21) thanks @mbojan

### MINOR IMPROVEMENTS

* Sped up some tests (#19)
* Better description of OAI protocol in the `DESCRIPTION` file (#28)
* Commented about where some internal functions come from (#31)
* Import `plyr` for `rbind.fill()` (#32)
* Including now some examples using OAI-PMH with GBIF and BHL (#33)
* Better error handling! (#10) (#12) (#22) (#27) thanks @mbojan

### BUG FIXES

* Fixed bug where `list_identifiers()` threw error when no result was found (#13)
* Dealing better with bad inputs to `as` parameter - stop with informative message
now instead of failing without anything returned (#34)

oai 0.1.0
=========

* Released to CRAN.
