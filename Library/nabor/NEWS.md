# nabor 0.5.0

* add radius bounded search (#6, #7)
* remove register keyword from C++ libnabo (#9)

# nabor 0.4.7

* dev: register native routines as recommended for R 3.4.0
* dev: update travis 

# nabor 0.4.6

* make RANN test conditional
* Fix description for R 3.2.0
* roxygen2/Rcpp updates

# nabor 0.4.5

* document behaviour for invalid input points to knn
  (thanks to @cmpop1, see https://github.com/jefferis/nabor/issues/3)

# nabor 0.4.4

* fixes warning for non-portable Makefile (thanks to Brian Ripley)

# nabor 0.4.3

* fixes a UBSAN error in libnabo (thanks to Brian Ripley)

# nabor 0.4.2

* license fixes for CRAN (include GJ as copyright holder)

# nabor 0.4.1

* CRAN release
* insist on RcppEigen >= 0.3.2.2.0 (conservative but should not be restrictive)

# nabor 0.4

* corectly handle non-floating point input to knn (and note that WKNN classes 
  must be initialised with floating point input)
* standardise method names of WKNN classes so that they are identical for WKNND
  and WKNNF objects
* doc: improve main package docs
* dev: templatise WKNN C++ classes

# nabor 0.3

* preparing for CRAN
* rename package from nabo to nabor
* fix handling of non-matrix inputs by knn
* knn performs self-query when query is missing
* dev: only link to RcppEigen
* dev: remove boost/any.hpp header and insist on latest BH package (>= 1.54.0-4)
  (also bumps R requirement to 3.0.2).

# nabor 0.2

* first public version (github)
* includes knn function (with different search types, defaulting to auto)
* ... and WKNNF and WKNND classes that wrap a tree object
