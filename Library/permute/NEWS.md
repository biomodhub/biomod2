# permute 0.9-10

## User visible changes

* *permute* now depends on R versions >= 3.6.0

## Bug fixes

* `numPerms()` tries harder to avoid floating point issues when computing the
  number of permutations for the current design. #41

# permute 0.9-8

* Updated reference output for examples following release of *vegan* 2.7.1.

# permute 0.9-7

### Bug fixes

* The documented behaviour of `shuffleFree()` allowed `shuffleFree(x = 1:10)`,
  with `x` passed immediately to `sample.int()`, which is incorrect and raises
  an error in the development version of R (at the time of writing).
  Reported by Brian Ripley.

# permute 0.9-6

### Bug fixes

* `allPerms()` was calling `allSeries()` with the incorrect number of
  permutations, leading to problems with `shuffleSet()` also.
  [#28 reported by Cajo ter Braak](https://github.com/gavinsimpson/permute/issues/28)

# permute 0.9-5

Changes to accommodate the new RNG in R >= 3.6.0

# permute 0.9-4

The example in `?check` was made to suppress package startup messages from
vegan.

# permute 0.9-3

This release fixed some non-canonical-form CRAN URLs.

# permute 0.9-2

This release was purely to avoid issues with CRAN as a new release of vegan had
been released and the example reference material hadn't been updated to match.

# permute 0.9-1

### General

A minor bug fix release to address a single problem.

### Bug Fixes

 * `shuffleSet()` wasn't returning a matrix if `nset = 1` *and* `allPerms` was
   invoked because of a low set of possible permutations.
   [GitHub Issue #19](https://github.com/gavinsimpson/permute/issues/19)

# permute 0.9-0

### General

This is small update to **permute**, focused mainly on ensuring the 
many combinations of restrictions on permutations allowed by the 
package work. An extensive test suite has been written which covers 
~87% of the package's codebase at the time of release.

### New features

* Permutation matrices produced by `shuffleSet()` are now printed in a 
more compact form.

* Better heuristics in `check()` allow for more reliable permutations 
(i.e. fewer duplicate permutations) when the set is small. This has 
increased the `minperms` setting. Consequently we generate all possible 
permutations up to ~500,000 more often as we now randomly sample from 
the entirely generated set rather than randomly generate permutations. 
This provides a small performance hit in some rare cases.

* `shuffleSet()` has a new argument `quietly = FALSE` which is passed 
on to `check()`.

* A number of bugs were fixed. See the Changelog and the Bug reports on 
github for details.

### Defunct

* `permControl()` and `permuplot()` are defunct and have been removed 
from the package.

# permute 0.8-0

### General

* Version 0.8-0 represents a major update of **permute**, with some 
backwards-incompatible changes to the main functions. The main addition 
is the availability of block-level restrictions on the permutations, 
which are required for whole- and split-plot designs.

### New features

* `how()`, a new function to create permutation designs. This replaces 
the deprecated function `permControl`.

* **permute** gains the addition of true blocking structures with which 
to restrict the permutations. Blocks sit as the outermost layer of the 
permutations, and can contain plots which in turn contain samples. In 
contrasts to plots, blocks are never permuted and samples are never 
shuffled between blocks. Permutation only ever happens within blocks.

    To facilitate this, plot-level strata are now specified via 
    `Plots()` instead of via the old `strata` argument of `how()`. 
    Blocks are specified via a new argument `blocks`, which takes a 
    factor variable.

* A new suite of extractor and replacement functions is provided with 
which to interact with permutation designs created by `how(). Extractor 
functions have names `getFoo()`, where `Foo()` is a component of the 
design. Replacement functions have names `setFoo`. The replacement 
function are especially for use by package authors wishing to alter 
permutation within their functions. The extractor functions are 
recommended for all users wishing to extract features of the permutation 
design. 

* As a convenience for users, the \code{update()} function will now 
work with objects of classes `"how"`, `"Plots"` or `"Within"` to allow 
quick updating of features of the permutation design. This approach is 
intended for interactive use at the top-level and not within functions, 
where the new `setFoo` replacement functions should be used.

* `shuffleSet()` is enhanced in this version. Firstly, the function now 
returns a classed object which has a `print()` method to allow for 
compact printing of the design elements used to generate the set of 
permutations. Second, `shuffleSet()` will sample `nset` permutations 
from the entire set of permutations should a small number of possible 
permutations trigger generation of the entire set. This avoids the 
generation of a set of non-unique permutations. Finally the random seed 
that generated the set is stored as an attribute.

* `allPerms()` no longer assumes that samples are in block and/or plot 
ordering.

* The package vignette is much expanded in this version with new 
sections on using **permute** within functions that will be of interest 
to package authors wishing to use **permute** in their packages.

### Deprecated

* `permControl()` is deprecated in favour of `how()`.

* `permuplot()` is broken and effectively defunct given the changes to 
the way permutation are defined and the addition of blocks. 
`permuplot()` is no longer exported from the package namespace.
