
<!-- README.md is generated from README.Rmd. Please edit that file and run rmarkdown::render("README.Rmd") -->

## Robust Covariance Matrix Estimators

Model-robust standard error estimators for cross-sectional, time series,
clustered, panel, and longitudinal data. Modular object-oriented
implementation with support for many model objects, including: `lm`,
`glm`, `fixest`, `survreg`, `coxph`, `mlogit`, `polr`, `hurdle`, `zeroinfl`,
and beyond.

**Sandwich covariances for general parametric models:**

<img alt="Central limit theorem and sandwich estimator" src="man/figures/README-sandwich.svg" style="border:10px solid transparent">

**Object-oriented implementation in R:**

``` r
library("sandwich")
library("lmtest")
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)
coeftest(m, vcov = sandwich)
```

    ## t test of coefficients:
    ## 
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.0297     0.0284    1.05      0.3    
    ## x             1.0348     0.0284   36.45   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
coeftest(m, vcov = vcovCL, cluster = ~ firm)
```

    ## t test of coefficients:
    ## 
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.0297     0.0670    0.44     0.66    
    ## x             1.0348     0.0506   20.45   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
