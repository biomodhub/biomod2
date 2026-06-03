---
title: "estimability"
output: html_document
date: '2022-07-03'
---

R package **estimability**: Support for determining estimability of linear functions
====

[![cran version](https://www.r-pkg.org/badges/version/estimability)](https://cran.r-project.org/package=estimability)
[![downloads](https://cranlogs.r-pkg.org/badges/estimability)](https://cranlogs.r-pkg.org/badges/estimability)
[![total downloads](https://cranlogs.r-pkg.org/badges/grand-total/estimability)](https://cranlogs.r-pkg.org/badges/grand-total/estimability)
[![Research software impact](http://depsy.org/api/package/cran/estimability/badge.svg)](http://depsy.org/package/r/estimability/)


## Features
 * A `nonest.basis()` function is provided that determines a basis for the null
   space of a matrix. This may be used in conjunction with `is.estble()` to
   determine the estimability (within a tolerance) of a given linear function of
   the regression coefficients in a linear model.
 * A set of `epredict()` methods are provided for `lm`, `glm`, and `mlm` objects.
   These work just like `predict()`, except an `NA` is returned for any cases that
   are not estimable. This is a useful alternative to the generic warning that
   "predictions from rank-deficient models are unreliable."
 * A function `estble.subspace()` that projects a set of linear functions onto an  
   estimable subspace (possibly of smaller dimension). This can be useful in
   creating a set of estimable contrasts for joint testing.
 * Package developers may wish to import this package and incorporate
   estimability checks for their `predict` methods.

## Installation
 * To install latest version from CRAN, run 
```
install.packages("estimability")
```
Release notes for the latest CRAN version are found at [https://cran.r-project.org/package=estimability/NEWS](https://cran.r-project.org/package=estimability/NEWS) -- or do `news(package = "estimability")` for notes on the version you have installed.

* To install the latest development version from Github, have the newest **devtools** package installed, then run
```
devtools::install_github("rvlenth/estimability", dependencies = TRUE)
```
For latest release notes on this development version, see the [NEWS file](https://github.com/rvlenth/estimability/blob/master/inst/NEWS)
