## Restricted permutations with R

[![CRAN version](https://www.r-pkg.org/badges/version/permute)](https://cran.r-project.org/package=permute)
[![](https://cranlogs.r-pkg.org/badges/grand-total/permute)](https://cran.r-project.org/package=permute)
[![R-CMD-check](https://github.com/gavinsimpson/permute/workflows/R-CMD-check/badge.svg)](https://github.com/gavinsimpson/cocorresp/actions)
[![codecov](https://codecov.io/gh/gavinsimpson/permute/graph/badge.svg?token=2FYEfBBSJ7)](https://app.codecov.io/gh/gavinsimpson/permute)

## What is permute?

**permute** generates permutations from a range of restricted 
permutation designs.

Permute provides an R implementation of the permutation schemes 
developed by Cajo ter Braak and made available in the Canoco software, 
version 3.1 (ter Braak, 1990). These permutation schemes draw upon 
ideas from an earlier paper by Besag & Clifford (1989).

Several types of permutation are available in **permute**:

 * Free permutation of objects
 * Time series or line transect designs, where the temporal or spatial ordering is preserved.
 * Spatial grid designs, where the spatial ordering is preserved in both coordinate directions
 * Permutation of plots or groups of samples.
 * Blocking factors which restrict permutations to within blocks. The preceding designs can be nested within blocks, allowing analysis of hierarchical designs (e.g. split plot designs)

### References

Besag, J. and Clifford, P. (1989) Generalized Monte Carlo significance tests. *Biometrika* **76**; 633&ndash;642.

ter Braak, C. J. F. (1990). *Update notes: CANOCO version 3.1*. Wageningen: Agricultural Mathematics Group. (UR).
