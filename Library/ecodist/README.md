**ecodist R package**

Dissimilarity-based analysis functions including ordination and Mantel test functions, intended for use with spatial and community data.


** CHANGES in ecodist 2.1.2**

  - fixed a scaling error in pco that resulted in the wrong magnitude (but correct pattern) of vectors


** CHANGES in ecodist 2.1.1**

  - added relrange to standardize matrices
  - updated mgram, pmgram, xmgram to be compatible with nclass.Sturges when calculating default number of breaks
  - added equiprobable option to mgram, pmgram, xmgram to calculate equal-number bins instead of default equal-width
  - added mstdist to use a minimum spanning tree calculation to estimate distances for pairs sites with no species in common
  - fixed distance to correctly return simple difference


**CHANGES in ecodist 2.0.10**

The proxy package (as of version 0.4-27), loaded by many spatial packages including spdep, overwrites the base behavior of dim() for dist objects, which breaks ecodist functions. A fix has been added to ecodist 2.0.10 that returns dim.dist() to its base state, as long as ecodist is loaded after all other packages. The maintainer plans to remove the problematic dim.dist() eventually, since ecodist is not the only package it breaks.

Until then, load ecodist last. Since dim() is called by functions outside of ecodist, I can't simply specify ecodist::dim() in all relevant cases.

dim(dist(matrix(1:15, ncol=3))) 

should return NULL and ***not*** c(5, 5) if the correct (base) dim() is being used.


**CHANGES in ecodist 2.0**

  - fixed bug in crosstab() that affected expansion of single-row or -column tables using allrows or allcols; changed result to data frame
  - added icov argument to distance() for use with Mahalanobis distance
  - changed stress calculation in nmds() to match vegan and MASS calculations; formerly was a similar method that was monotonically related, but not identical.
  - added plot.nmds to display stress and r2 for NMDS ordinations across a range of dimensions
  - added addord method to add new data to an existing NMDS ordination. 
  - added clusterlevel to calculate Mantel tests for specified groupings
  - added logistic regression to MRM
  - added xdistance cross-distance function, and cross-dissimilarity analysis functions xmantel, xmgram
  - updated examples in help files to be more helpful
  - added a vignette listing dissimilarity-based analyses

