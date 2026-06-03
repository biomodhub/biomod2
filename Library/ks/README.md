## Introduction

Kernel smoothing for data from 1- to 6-dimensions. This package forms the basis for the practical data analysis in the book [_Multivariate Kernel Smoothing and Its Applications_](https://mvstat.net/mvksa/).

There are three main types of functions in this package:

* computing kernel estimators - these function names begin with `k`
* computing bandwidth selectors - these begin with `h` (1-d) or `H` (>1-d)
* displaying kernel estimators - these begin with `plot`.
  
The kernel used throughout is the normal (Gaussian) kernel. For 1-d data, the bandwidth h is the standard deviation of the normal kernel, whereas for multivariate data, the bandwidth matrix H is the variance matrix.
 
The main function `kde()` computes a kernel density estimate. For display, its `plot` method calls `plot.kde()`. The bandwidth choice is crucial for the performance of kernel estimators. There are several varieties of bandwidth selectors available
   
* plug-in `hpi` (1-d); `Hpi`, `Hpi.diag` (2- to 6-d) 
* least squares (or unbiased) cross validation (LSCV or UCV) `hlscv` (1-d); `Hlscv`, `Hlscv.diag` (2- to 6-d) 
* biased cross validation (BCV) `Hbcv`, `Hbcv.diag` (2- to 6-d) 
* smoothed cross validation (SCV) `hscv` (1-d); `Hscv`, `Hscv.diag` (2- to 6-d) 
* normal scale `hns` (1-d); `Hns` (2- to 6-d).

For an example with bivariate data, see `vignette("ks")`. 

The other types of kernel estimators follow a similar functionality.

## Installation

Install from CRAN:

```r
install.packages("ks") 
```

## Further reading

Chacon, J.E. & Duong, T. (2018) [_Multivariate Kernel Smoothing and Its Applications_](https://mvstat.net/mvksa/). Chapman & Hall/CRC, Boca Raton. 
  
Duong, T. (2004) [_Bandwidth Matrices for Multivariate Kernel Density Estimation_](https://mvstat.net/tduong/research/publications/duong-2005-thesis.pdf) Ph.D. Thesis, University of Western Australia. 
