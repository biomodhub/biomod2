[![CRAN Status Badge](http://www.r-pkg.org/badges/version/phylolm)](https://CRAN.R-project.org/package=phylolm)
[![Research software impact](http://depsy.org/api/package/cran/phylolm/badge.svg)](http://depsy.org/package/r/phylolm)

## phylolm: R package for Phylogenetic Linear Regression on very large trees

The package provides functions for fitting phylogenetic linear models and phylogenetic generalized linear models. 
Computations use an algorithm that is linear in the number of tips in the tree. 
The package also provides functions for simulating continuous or binary traits along the tree.
When a new version is stable, it is pushed to CRAN, so the version available here is newer than that on CRAN.

- Lam Si Tung Ho and Cécile Ané (2014). 
Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models. 
*Methods in Ecology and Evolution* 5(11):1133-1146. 
[(link)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12285)

- Lam Si Tung Ho and Cécile Ané (2014). 
A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. 
*Systematic Biology* 63(3):397-408.
[(link to pdf)](https://doi.org/10.1093/sysbio/syu005)

### Installation
Install `phylolm` from `github`:
```{r}
devtools::install_github("lamho86/phylolm")
```

### Main features

- phylogenetic signal
- phylogenetic linear, logistic and Poisson regression
- stepwise model selection (from v2.3)
- OU shift detection
- continuous and discrete trait simulators with covariates
- bootstrap-based confidence intervals for phylogenetic logistic regression (from v2.3)
  and phylogenetic regression (from v2.4.2)
- goodness-of-fit test of a population tree with the coalescent (from v.2.4)
- allowing measurement errors in phylogenetic linear regression (from v2.4.1)
- log likelihood of an one-dimensional Ornstein-Uhlenbeck model (from v2.5)
- bootstrapping can be parallelized (from v2.6)
- R-squared and Adjusted R-squared (from v2.6.1)
- phyloglmstep (from v2.6.2)
- mace (from v2.6.4)
- REML for phylolm (from v2.6.4)
- Compatible with visreg (from v2.6.5)
