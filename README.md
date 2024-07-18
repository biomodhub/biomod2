[![Cran Version](https://www.r-pkg.org/badges/version/biomod2?color=yellow)](https://cran.r-project.org/package=biomod2)
[![Github Version](https://img.shields.io/badge/devel%20version-4.2--6-blue.svg)](https://github.com/biomodhub/biomod2)
[![Last Commit](https://img.shields.io/github/last-commit/biomodhub/biomod2.svg)](https://github.com/biomodhub/biomod2/commits/master)
[![R-CMD-check](https://github.com/biomodhub/biomod2/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/biomodhub/biomod2/actions/workflows/R-CMD-check.yml)

<!-- [![Download](http://cranlogs.r-pkg.org/badges/grand-total/biomod2?color=yellow)](https://cran.r-project.org/package=segclust2d) -->
<!-- 
badge for github version :
badger::badge_github_version("biomodhub/biomod2", "blue") 
-->

<style>
.zoom p {
width:800px;
margin-left: auto;
margin-right: auto;
}
.zoom p:hover {
width:1500px;
position: relative;
z-index: 10;
}
</style>


<div align="center">
<p><img src="articles/pictures/LogoBiomod.png" alt="Logo biomod2" width=150px></img></p>

<b>------------------------------------------------------------<br/>
Species distribution modeling, <br/>
calibration and evaluation, <br/>
ensemble modeling <br/>
------------------------------------------------------------<br/>
</b>

https://biomodhub.github.io/biomod2/
</div>


### <i class="fas fa-tools"></i> Installation

<br/>

- **Stable version** [![v](https://www.r-pkg.org/badges/version/biomod2?color=yellow)](https://cran.r-project.org/package=biomod2) from [cran](https://CRAN.R-project.org/package=biomod2) :

```R
install.packages("biomod2", dependencies = TRUE)
```

<br/>

- **Development version** [![v](https://img.shields.io/badge/devel%20version-4.2--6-blue.svg)](https://github.com/biomodhub/biomod2) from [biomodhub](https://github.com/biomodhub/biomod2) :

```R
library(devtools)
devtools::install_github("biomodhub/biomod2", dependencies = TRUE)
```

<br/><br/>

### <i class="fas fa-envelope-open-text"></i> `biomod 4.2-6` - Improved OptionsBigBoss and new model


`/!\` Please **feel free to indicate if you notice some strange new behaviors** !

#### <i class="fas fa-exchange-alt"></i> What is changed ?

- To improve the models, we made some change in the options for [**`OptionsBigboss`**](https://biomodhub.github.io/biomod2/reference/OptionsBigboss.html). (This only affects the ANN, CTA and RF models.) You can check all your options with the `get_options()` function.
- To reduce the tuning calculation time, we update the tuning ranges for ANN, FDA and MARS models.

#### <i class="fas fa-plus-square"></i> What is new ?

- `biomod2` has a new model: **RFd**. It's a Random Forest model with a down-sampling method.
- You can now define _seed.val_ for `bm_PseudoAbsences()` and `BIOMOD_FormatingData()`.
- New _fact.aggr_ argument, for pseudo-absences selection with the random and disk methods, allows to reduce the resolution of the environment.
- Possibility to give the same options for all datasets with _"for_all_datasets"_ in `bm_ModelingOptions()`.

<br/>

<div class="zoom">
<p><img src="articles/pictures/SCHEMA_BIOMOD2_WORKFLOW_functions.png" alt="Main workflow"></img></p>
</div>

<br/><br/><br/>

<br/><br/>


### <i class="fas fa-envelope-open-text"></i> `biomod 4.2-5` - Modeling options & Tuning Update

#### <i class="fas fa-exchange-alt"></i> What is changed ?

- modeling options are now automatically retrieved from single models functions, normally allowing the use of all arguments taken into account by these functions
- tuning has been cleaned up, but keep in mind that it is still a quite long running process
- in consequence, `BIOMOD_ModelingOptions` and `BIOMOD_Tuning` functions become secundary functions (`bm_ModelingOptions` and `bm_Tuning`), and modeling options can be directly built through `BIOMOD_Modeling` function

#### <i class="fas fa-plus-square"></i> What is new ?

- `ModelsTable` and `OptionsBigboss` datasets (*note that improvement of bigboss modeling options is planned in near future*)
- 3 new vignettes have been created :
    - [data preparation](https://biomodhub.github.io/biomod2/articles/vignette_dataPreparation.html) (*questions you should ask yourself before modeling*)
    - [cross-validation](https://biomodhub.github.io/biomod2/articles/vignette_crossValidation.html) (*to prepare your own calibration / validation datasets*)
    - [modeling options](https://biomodhub.github.io/biomod2/articles/vignette_dataPreparation.html) (*to help you navigate through the new way of parameterizing single models*)

<br/>

### <i class="fas fa-envelope"></i> `biomod 4.2` - Terra Update

#### <i class="fas fa-exchange-alt"></i> What is changed ?

- `biomod2` now relies on the new [`terra`](https://github.com/rspatial/terra) package that aims at replacing `raster`and `sp`.
- `biomod2` is still compatible with old format such as `RasterStack`and `SpatialPointsDataFrame`.
- `biomod2` function will sometimes return `SpatRaster` from package `terra` that you can always convert into `RasterStack` using function `stack` in `raster`.

<br/><br/>



### <i class="fas fa-envelope"></i> `biomod 4.1` is now available

`/!\` Package fresh start... meaning some changes in function names and parameters. We apologize for the trouble `>{o.o}<` <br/>
Sorry for the inconvenience, and please **feel free to indicate if you notice some strange new behaviors** !


#### <i class="fas fa-exchange-alt"></i> What is changed ?

- all code functions have been cleaned, and old / unused functions have been removed
- function names have been standardized (`BIOMOD_` for main functions, `bm_` for secundary functions)
- parameter names have been standardized (same typo, same names for similar parameters across functions)
- all documentation and examples have been cleaned up

#### <i class="fas fa-plus-square"></i> What is new ?

- plot functions have been re-written with `ggplot2`
- [`biomod2` website](https://biomodhub.github.io/biomod2/) has been created, with proper `roxygen2` documentation and help vignettes

#### <i class="fas fa-question-circle"></i> But... why ?

- “*For every minute spent on organizing, an hour is earned.*” — Benjamin Franklin
  - better documentation, better formation, better help provided
  - new improvements to come (update of single models, implementation of abundance models, etc)

<br/><br/>




