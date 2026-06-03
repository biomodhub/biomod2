[![CRAN version](https://www.r-pkg.org/badges/version/ENMeval)](https://CRAN.R-project.org/package=ENMeval) [![downloads](https://cranlogs.r-pkg.org:443/badges/grand-total/ENMeval?color=orange)](https://cranlogs.r-pkg.org:443/badges/grand-total/ENMeval?color=orange)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/jamiemkass/ENMeval/workflows/R-CMD-check/badge.svg)](https://github.com/jamiemkass/ENMeval/actions)


# ENMeval version 2.0.5

## R package for automated tuning and evaluations of ecological niche models

[`ENMeval`](https://jamiemkass.github.io/ENMeval/index.html) is an R package that performs automated tuning and evaluations of ecological niche models / species distribution models. These models make predictions of species' niche relationships and potential geographic distributions based on presence data, environmental predictor variables, and a sample of available environmental conditions (i.e., background data). 

"Model tuning" is commonly used for machine-learning models. It means building candidate models with a range of complexity settings, evaluating the accuracy of each one (here with cross-validation), then selecting optimal settings for your data based on those of the best-performing model. This exercise is important because it is difficult to predict in advance how complex your model needs to be to make accurate and ecologically realistic predictions for your species. Too much model complexity leads to overfitting, where your model fits your data very well but it cannot predict new data accurately. Model tuning helps maximize predictive ability while avoiding model overfitting. 

The `ENMeval` package features a single function that performs model tuning based on user specifications, including methods for partitioning data for cross-validation (random, leave-one-out, spatial, custom), and evaluates models using predefined performance metrics (AUC, Continuous Boyce Index, omission rates) with the option to insert others. The package includes functionality for three models: [maxent.jar](https://doi.org/10.1016/j.ecolmodel.2005.03.026) (Java implementation of Maxent), [maxnet](https://doi.org/10.1111/ecog.03049) (R implementation of Maxent), and [BIOCLIM](https://doi.org/10.1111/ddi.12144) (climate envelope method). Users can also specify other algorithms by customizing an **ENMdetails** object (`?ENMdetails`). The package also offers comprehensive metadata output, null model evaluations, visualization tools, and an extensive  [tutorial](https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html) that walks you through a full analysis workflow. Many features in `ENMeval` were created in response to user requests -- thank you for your input! Version >=2.0.0 represents an extensive restructure and expansion of previous versions, and 2.0.5 is a big move from `raster` and `dismo` functions to those of `terra` and `predicts`. 

For a more detailed description of `ENMeval`, please reference the most recent publication:

[Kass, J. M., Muscarella, R., Galante, P. J., Bohl, C., Pinilla-Buitrago, G. E., Boria, R. A., Soley-Guardia, M., & Anderson, R. P. (2021). ENMeval 2.0: redesigned for customizable and reproducible modeling of species’ niches and distributions. Methods in Ecology and Evolution, 12: 1602-1608.](https://doi.org/10.1111/2041-210X.13628)

For the original package version, please reference this older publication:

[Muscarella, R., Galante, P. J., Soley-Guardia, M., Boria, R. A., Kass, J. M., Uriarte, M. and Anderson, R. P. (2014), ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. Methods in Ecology and Evolution, 5: 1198–1205.](https://doi.org/10.1111/2041-210X.12261)

## NOTES:

1. `ENMeval` is a work in progress, changing slowly to fix bugs when users identify them. If you find a bug, please raise an Issue in this Github repo and I will resolve it as soon as I can. The CRAN version may lag behind the Github one, so please try the development version here first if you are having any issues.
Install with: `devtools::install_github("jamiemkass/ENMeval")`

2. The vignette is not included in the CRAN version of the package due to file size constraints, but is [available](https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html) on the package's Github Pages website. 

3. Please make sure to use the most recent version of [maxent.jar](https://biodiversityinformatics.amnh.org/open_source/maxent/), as bug fixes are occasionally made.

4. Note that as of version 0.3.0, the default implementation uses the ['maxnet' R package](https://cran.r-project.org/package=maxnet). The output from this differs from that of the original Java program and so some features are not compatible (e.g., variable importance, html output).
