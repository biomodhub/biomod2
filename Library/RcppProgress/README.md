# RcppProgress
[![Build Status](https://travis-ci.org/kforner/rcpp_progress.svg?branch=dev)](https://travis-ci.org/kforner/rcpp_progress?branch=dev)
[![codecov](https://codecov.io/github/kforner/rcpp_progress/coverage.svg?branch=dev)](https://codecov.io/github/kforner/rcpp_progress?branch=dev)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/RcppProgress)](https://cran.r-project.org/package=RcppProgress)

a R package that provides a c++ interruptible progress bar with OpenMP support for c++ code in R packages

## example
see a detailed example on Rcpp Gallery:
http://gallery.rcpp.org/articles/using-rcppprogress/

## How to build

Prerequisites:

- OpenMP support to use the multithreaded parallelized version. OpenMP is available in GCC >= 4.2

Just install it the usual way.

If you want more control, unarchive it, cd to the source directory, then type
R CMD INSTALL . in the console.


## Contribute
Send me a pull request with at least one test or example


## For developers

### tests and check

If you have all the RcppProgress dependencies (and suggests) installed:

type:
 - `make tests`: to run the tests
 - `make check`: to check the package

### docker-checker

A Dockerfile (<docker_checker/Dockerfile>) is provided to help building the
dev environment (built on rocker/r-devel) in which to develop
and test RcppProgress.

type:

 - `make docker/build`: to build the docker
 - `make docker/run`: to run a shell in the docker with the current dir mounted
 	inside
 - `make docker/check`: to check the package inside the docker
 - `make docker/tests`: to run test tests of the package inside the docker

### test on windows using rhub

```
make docker/run
make check_rhub_windows
```


