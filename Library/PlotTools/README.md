# PlotTools

[![Project Status: Active.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#project-statuses)
[![codecov](https://codecov.io/gh/ms609/PlotTools/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ms609/PlotTools)
[![CRAN Status Badge](https://www.r-pkg.org/badges/version/PlotTools)](https://cran.r-project.org/package=PlotTools)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/PlotTools)](https://cran.r-project.org/package=PlotTools)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7803390.svg)](https://doi.org/10.5281/zenodo.7803390)


'PlotTools' is an R package that allows the legends for continuous variables
to be added to plots using the familiar functions of the default 'graphics'
package.  It also includes utility functions to manipulate irregular polygons
and locate their centres.

Install the latest production version from CRAN using

```r
install.packages("PlotTools")
```

or install the development version from GitHub with

```r
devtools::install_github("ms609/PlotTools")
```

Please let me know of any feature requests or bugs by [opening an 
issue on GitHub](https://github.com/ms609/PlotTools/issues/).

## Usage

```r
# Select a colour palette
palette <- if (packageVersion("grDevices") > 3.6) hcl.colors else heat.colors

# Plot some example data
plot(
  cars,
  xlab = "Speed / mph",
  ylab = "Stopping distance / ft",
  col = palette(125 + 1)[cars$dist + 1], # Colour points by distance
  cex = cars$speed / 10,                 # Size points by speed
  pch = 16,                              # Use filled circle for points
  frame.plot = FALSE
)

# Display legend for colour scale
PlotTools::SpectrumLegend(
  "topleft",                             # Legend position
  palette = palette,                     # Display our chosen palette
  legend = seq(125, 0, length.out = 6),  # Annotate positions on legend
  title = "Distance",
  bty = "n"                              # Don't frame with box
)

# Display legend for plotting symbol sizes
PlotTools::SizeLegend(
  "bottomright",                         # Legend position
  palette = "darkgrey",                  # Set colour - may be continuous
  horiz = TRUE,                          # Orient horizontally
  width = c(0, 2.5), scale = "pch",      # Scale for plotting character
  legend = seq(0, 25, by = 5),           # Annotate positions on legend
  x.intersp = 0,                         # Set x spacing
  bty = "n",                             # Don't frame with box
  inset = 0.05,                          # Inset from plot edges
  title = "Speed"
)
```

![image](https://user-images.githubusercontent.com/1695515/232433317-85412267-8c01-41f9-96be-bfc71533d59b.png)

## Citation

Cite this package as:

Smith, Martin R. (2023). _PlotTools: Add continuous legends to plots._
Comprehensive R Archive Network, 
[doi:10.5281/zenodo.7803390](https://dx.doi.org/10.5281/zenodo.7803390).


## Contribute

Please note that the 'PlotTools' project is released with a
[Contributor Code of Conduct](https://ms609.github.io/packages/CODE_OF_CONDUCT.html).
By [contributing](https://ms609.github.io/packages/CONTRIBUTING.html) to this project, you agree to abide by its terms.
