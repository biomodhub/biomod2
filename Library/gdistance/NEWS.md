gdistance 1.6.5 (2025-09-24)
=========================

* Add package anchors to cross-references for CRAN

gdistance 1.6.4 (2023-06-19)
=========================

* Optimized the `transition()` function for 4 and 8 directions with user-defined functions. Run time is reduced by 10+%, and the peak memory consumption is massively reduced.
* Updated deprecated usage of the Matrix library.
* Added unit tests for the `transition()` function.
* Add workaround (fix?) for `pathInc()` code in Overview vignette.

gdistance 1.6.2 (2023-04-20)
=========================

* Maintainer change

gdistance 1.6 (2022-10-10)
=========================

gdistance 1.4 (2022-03-8)
=========================

* update the file maungawhau

gdistance 1.3-6 (2020-06-29)
=========================

* expanded Description

### BUG FIXES 
* Fix transition with Mahalanobis distance method use `x.minus.y <- raster::getValues(x)[adj[,1],] - raster::getValues(x)[adj[,2],]` instead of `x.minus.y <- xy[adj[,1],-1]-xy[adj[,2],-1]` 

gdistance 1.3-1 (2020-02-28)
=========================
### BUG FIXES 
* export method `raster,TransitionLayer-method`

gdistance 1.3-0 (2020-02-17)
=========================

### BUG FIXES

* Fix R:devel errors in running the examples

### NEW FEATURES

* gdistance builds with roxygen2
* new website and venue to report issues
* a recommended citation is added


gdistance 1.2-2 (2018-05-01)
=========================

* cBind, rBind (deprecated in Matrix package) changed to cbind, rbind

gdistance 1.2-1 (2017-02-01)
=========================

* Changes to reflect publication in JSS
* Bug fixed in AccCost (mode="out"), thanks to Henjo de Knegt

gdistance 1.1-9 (2015-07-01)
=========================

* Changes to vignette based on JSS review

gdistance 1.1-7 (2015-03-01)
=========================

* Final touches for submission to JSS

gdistance 1.1-6 (2015-01-01)
=========================

* Final touches for submission to JSS
* Correction in shorthestPath (bug reported by Sergei Petrov)
* Dependencies and imports reflect new CRAN rules

gdistance 1.1-5 (2014-02-01)
=========================

* Updated to reflect changes in igraph
* Changes in vignette
* Changes in default method pathInc (overlap)


gdistance 1.1-4 (2012-12-01)
=========================

* Class changes to BasicRaster
* Change from igraph0 to igraph 6
* Change in as(x, "transitionMatrix") to speed up transition()
* Update to the new adjacent method in raster

gdistance 1.1-3 (2012-04-01)
=========================

* Further adapting to naming changes in raster ("x" argument)
* Making the same change in gdistance ("x" is the first argument in any function)
* Bug fix in resistanceDistance: uncomment reproject
* Improved documentation of definitions of different distance functions
* Dependency on igraph0 -- to be changed to igraph 0.6 after May 31, 2012
* Removed initialize functions, using prototypes instead


gdistance 1.1-2 (2011-09-01)
=========================

* Bug fix in `transition(transitionFunction="barriers")`
* Methods added in `transition()`
* Summary functions redefined
* Adapting to changes in raster package (argument for Raster* objects is now always "x")
* adjacencyFromTransition added
