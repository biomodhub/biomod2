# nabor
[![DOI](https://zenodo.org/badge/4241/jefferis/nabor.svg)](http://dx.doi.org/10.5281/zenodo.17873) 
[![Release Version](https://img.shields.io/github/release/jefferis/nabor.svg)](https://github.com/jefferis/nabor/releases/latest) 
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nabor)](https://cran.r-project.org/package=nabor) 
[![Build Status](https://travis-ci.org/jefferis/nabor.svg)](https://travis-ci.org/jefferis/nabor) 

R package **nabor** wraps [libnabo](https://github.com/ethz-asl/libnabo), 
a fast K Nearest Neighbour library for low-dimensional spaces implemented in templated C++.
In comparison with the widely used [ANN](http://www.cs.umd.edu/~mount/ANN) library (wrapped by the
[RANN](https://cran.r-project.org/package=RANN) R package), **libnabo** is reported
to be 5% to 20% faster with more compact data structures.

## Quick start
```r
# install (see below for details)
install.packages("nabor")

# use
library(nabor)

# help
?nabor
?knn

# run examples
example(knn)
example(WKNN)

# run tests
library(testthat)
test_package("nabor")

# cite
citation("nabor")
```

## nabor vs RANN
For R users **nabor** provides a function, `knn`, that is a drop in replacement for
the `nn2` function in the [RANN](https://cran.r-project.org/package=RANN) 
R package. I have seen speedups of 2-3x fold for queries of interest (a few thousand
points in 3d, k=1) when comparing nabor::knn and RANN::nn2. See `?knn` for details.

Furthermore **nabor** provides a mechanism for reusing the k-d search tree structure for 
multiple queries. This is achieved by wrapping a libnabo k-d tree and associated points
into a C++ class. This in turn is wrapped as an R reference class (by RcppModules)
that can be used in R. See `?WKNN` for details. The `WKNNF` class has the additional
feature of using floats (4 bytes per coordinate) for the underlying storage, rather
than the doubles used by R; this may be useful for large pointsets.
## Installation
### Released version from CRAN
The current stable version of the package is available from 
[CRAN](https://cran.r-project.org/). The package requires compilation, but
installing from CRAN allows mac and windows users without the full C++ compiler toolchain to install binary packages.

```r
install.packages("nabor")
```

### Development version from github
The **nabor** package is known to compile from source with the standard C(++) 
compiler toolchains for R under MacOS X, Windows and Linux. See 
https://www.rstudio.com/products/rpackages/devtools/ for details of the
developer toolchains needed for your platform.

Once you have installed the appropriate developer toolchain mentioned above, you
can use the **devtools** package to install the development version of the package:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jefferis/nabor")
```

### Dependencies
The **nabor** package includes libnabo and all of its dependencies (boost, via 
package **BH**) and Eigen (via package **RcppEigen**) and does not
depend on any non-standard system libraries. It should therefore run out of the
box on any mac/linux/windows system.

## Acknowlegements
**libnabo** and therefore the **nabor** R package are released under the 
[BSD 3 clause license](https://www.r-project.org/Licenses/BSD_3_clause). If you
make use of **nabor** please cite the original authors:

```
> citation('nabor')

Elseberg J, Magnenat S, Siegwart R and Nüchter A (2012). “Comparison of nearest-neighbor-search
strategies and implementations for efficient shape registration.” _Journal of Software Engineering for
Robotics (JOSER)_, *3*(1), pp. 2-12. ISSN 2035-3928.

A BibTeX entry for LaTeX users is

  @Article{elsebergcomparison,
    title = {Comparison of nearest-neighbor-search strategies and implementations for efficient shape registration},
    author = {J. Elseberg and S. Magnenat and R. Siegwart and A. N{\"u}chter},
    journal = {Journal of Software Engineering for Robotics (JOSER)},
    pages = {2--12},
    volume = {3},
    number = {1},
    year = {2012},
    issn = {2035-3928},
  }

```

**nabor** also makes use of the tremendous [Rcpp](https://cran.r-project.org/package=Rcpp)
and [RcppEigen](https://cran.r-project.org/package=RcppEigen) packages –
kudos to their authors!
