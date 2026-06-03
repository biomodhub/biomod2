# TreeTools 2.3.0 (2026-04-22) #

## Performance

- Drop `stringi` dependency.
- `as.character.Splits()` re-implemented in C++; ~3× faster on 200-tip trees.

## Usability

- `Consensus()` and `SplitFrequency()` now respond to user interrupts during
  long-running computations.
- `SL_MAX_TIPS` compile-time constant increased to 32704, allowing support
  for larger trees.


# TreeTools 2.2.0 (2026-03-18) #

## New functionality

- `EdgeRatio()` reports the ratio of external:internal edges.
- `SplitInformation()` supports `Splits` and `phylo` objects.
- `sort` and `order` support `Splits` objects.
- `SplitFrequency(reference = NULL)` returns frequency of all splits.
- `as.TreeNumber()` now supports trees with up to 51 leaves (previously 19).
  Trees with 20–51 leaves have more than 2^64 distinct topologies, so their
  `TreeNumber` is stored as a decimal character string rather than `integer64`.
  The 19-leaf limit for `integer64`-backed storage (and `as.MixedBase()`
  round-trips) is unchanged.
- `as.TreeNumber()` no longer warns for trees with 20–44 leaves.
- `inst/include/TreeTools/tree_number.h` added to support downstream packages
  (e.g. TBRDist) via `LinkingTo: TreeTools`. Provides 256-bit tree number
  encoding/decoding supporting up to 51 leaves, extended from the 44-leaf
  limit of the previous `uint64_t`-based implementation.

## Performance

- `SplitFrequency(reference = NULL)`: split normalization moved to C++;
  internal split de-duplication uses hash map instead of ordered map.
- `NodeDepth()` for unrooted trees rewritten as O(n) two-pass C++ algorithm,
  replacing iterative R while-loop.
- `duplicated.Splits()` uses hash-based O(n) de-duplication, replacing
  O(n²) pairwise comparison.
- `RenumberTips.multiPhylo()` applies tip permutation in a single C++ call,
  avoiding per-tree overhead.

## Fixes

- `PhyDatToMatrix()` no longer crashes on zero-character `phyDat` objects
  (e.g. from a star tree); returns a 0-column matrix with correct row names.
- `AddUnconstrained()` handles zero-character `phyDat` input gracefully.

## Dependencies

- `RCurl` moved from Imports to Suggests; `ReadMrBayesTrees()` uses
  `RCurl::url.exists()` when available, falling back to base R `url()` for
  URL checks.


# TreeTools 2.1.0 (2026-02-10) #

- Add method `RenumberTips.Splits()`.
- `FirstMatchingSplit(x)` provides performant alternative to `match(x)[[1]]`.
- Support logical `pole` in `PolarizeSplits()`.
- `RenumberTree()` supports numeric `tipOrder` input.
- Support larger trees in `ClusterTable` objects.
  * Retires `CT_PUSH` and `CT_POP` macros.
- Support larger trees in `Consensus()`.
  * Uses 32-bit integers, necessitating downstream changes to TreeDist.
- Remove hard limit on tree size in `SplitList`.
- `MatrixToPhyDat()` gains `tipLabels` parameter.
- Document return value for `J1Index()`.


# TreeTools 2.0.0 (2025-09-23) #

## New functionality
- `Cherries()` counts the cherries in a binary tree.
- New method `as.Splits.integer()`.
- Add methods for `NULL` objects.

## Fixes and enhancements
- Fix `RoguePlot(sort = TRUE)`
  ([Rogue#33](https://github.com/ms609/Rogue/issues/33)).
- Remove R.cache dependency: `UnrootedKeys()` now uses a native cache
  implementation.
  
## Backward incompatible changes
- Require R 3.6.
- Remove deprecated functions `AllDescendantEdges()`, `.EnforceOutgroup()`,
  `ForestSplits()`, `in.Splits()`, `PhylogeneticInfo()`, `SpectrumLegend()`,
  `SplitNumber()`, `TreeSplits()`.
- Remove deprecated C++ macro `TREETOOLS_SPLITLIST_INIT`.


# TreeTools 1.16.1 (2025-08-24) #

- Compiler-safe vector initialization, resolving M1-SAN warnings.


# TreeTools 1.16.0 (2025-08-22) #

## New functionality

- `SplitConsistent()` calculates split (dis)agreement.
- `LongBranch()` identifies long-branched taxa.
- `Treeness()` computes the treeness (=stemminess) of a tree, a proxy for its
  phylogenetic signal.
- Add `KeepTip()` methods to correspond to `DropTip()`.
- `Preorder()` gains `topologyOnly` argument.

## Enhancements

- `MakeTreeBinary()` explicitly removes edge lengths.
- Optimize `ClusterTable` class.
- Optimize `NSplits()`.
- Remove exported C++ cache objects to prefer calculation with intrinsics.
- Other performance improvements.

## Fixes

- Fix `as.ClusterTable()` when leaf order varies.
- Fix regressions in `as.ClusterTable()` that caused downstream errors.
- Fix regressions in `PhyToString()`.
- Fix handling of multiple ambiguities in `Reweight()`.


# TreeTools 1.15.0 (2025-07-16) #

- `Reweight()` sets the weight of characters in a phylogenetic dataset.
- `MatchStrings()` checks for mismatched tip labels, suggesting corrections.
- `RootTree()` gains `fallback` argument to handle outgroups that do not root a
  tree.
- Fix `MakeTreeBinary()` labelling trees as in preorder.
- Fix `as.Splits.matrix(tipLabels != NULL)`.
- Improve performance of `PhyToString()`.
- Modernize aspects of C++ code.


# TreeTools 1.14.0 (2025-05-13) #

- `AddTip(lengthBelow = NA)` adds leaf at node without adding a new edge.

- `BalancedTree()` and equivalent gain a `lengths` parameter to specify edge
  lengths.

- Fix taxa misplaced by `RoguePlot(sort = TRUE)`.

- Fix unexpected polytomies in `Consensus()`
  ([#168](https://github.com/ms609/TreeTools/issues/168)).


# TreeTools 1.13.1 (2025-04-07) #

- Support non-unique labels in `DropTip()` checks.


# TreeTools 1.13.0 (2025-01-09) #

- `%in%` and `match` methods for phylo / multiPhylo objects.

- `Decompose()` decomposes additive (ordered) phylogenetic characters by
  binary decomposition.
  
- Check for overflow in splits functions.


# TreeTools 1.12.0 (2024-07-25) #

## New methods and functions

- `TopologyOnly()` removes metadata from phylo objects.
- `J1Index()` computes the robust, universal tree balance measure of 
  Lemant _et al_. 2022 <doi:10.1093/sysbio/syac027>, incorporating [code](
  https://github.com/robjohnnoble/RUtreebalance) by Rob Noble.

## Enhancements

- Consistent sequence of list entries in phylo objects.
- `RandomTree()` returns trees for < 3 leaves.
- `root_on_node()` handles trees with < 2 leaves.
- Support larger trees in `TotalCopheneticIndex()`, fixing
  [#158](https://github.com/ms609/TreeTools/issues/158).


# TreeTools 1.11.1 (2024-06-06) #

- Set random seed and increase tolerance to avoid false negatives on tests.


# TreeTools 1.11.0 (2024-05-23) #

## New methods and functions

- `YuleTree()` generates a random tree by the Yule process.
- `DescendantTips()` complements `DescendantEdges()`, rewritten in C++,
  fixing a bug when edges were not in preorder.
- `NodeNumbers()` returns the indices of nodes within a tree.

## Enhancements

- `RandomTree(root = TRUE)` roots the tree on a random edge.
- `RoguePlot()$legendLabels` returns suggested labels for legend.
- Support node labels in `AddTip()`, `CollapseNode()`, `DropTip()`, 
  `MakeTreeBinary()`, `Renumber()`, `Reorder()`, `SortTree()`, `Subtree()`
  ([#149](https://github.com/ms609/TreeTools/issues/149)).
- `AddTip(edgeLength = NULL)` defaults to `lengthBelow`. This will become the
  default in a future release.
- An entry point to the C++ function `root_on_node()` is now exported
  (intended for expert use only).
- Fix handling of weighted trees by `root_on_node()`.
- Use `KeepTip()` internally so `SplitFrequency()` supports `Splits` objects
  as documented.

# TreeTools 1.10.0 (2023-08-18) #

## New methods and functions

- `TipTimedTree()` displays trees where leaves are associated with absolute
  ages.
- `ReadMrBayesTrees()` samples trees from posterior of MrBayes output.
- `is.TreeNumber()` method.

## Improvements

- Support zero-edge trees in `as.Splits()` and `NSplits()`.
- Support empty constraints in `AddUnconstrained()`.
- Add space between tokens in `WriteTntCharacters()` to support continuous
  characters ([#139](https://github.com/ms609/TreeTools/issues/139)).

## Deprecations and breaking changes

- Change order of parameters in `DescendantEdges()`
- Deprecate `AllDescendantEdges()`; use `DescendantEdges()` instead.
- Deprecate `EnforceOutgroup()`; use `RootTree()` instead.
- Remove `NonDuplicateRoot()` and `in.Splits()`.


# TreeTools 1.9.2 (2023-04-25) #

- Improve support for comments in `ReadNotes()`.
- Support Nexus-escaped `''`s in `ReadCharacters()`.
- Add `legend` parameter to `RoguePlot()`.
- `RoguePlot()` now returns invisibly.
- Deprecate `SpectrumLegend()` -- spun off to separate
  ["PlotTools"](https://ms609.github.io/PlotTools/) package.


# TreeTools 1.9.1 (2023-03-21) #

- `AddUnconstrained()` and `ImposeConstraint()` handle wider range of inputs.

- `PhyDatToMatrix()` can (and by default does) override levels to write
  ambiguous tokens in custom formats such as `{01}`.

- Call C functions using symbols, not strings.


# TreeTools 1.9.0 (2022-11-29) #

## New methods and functions

- `ZeroTaxonTree()` creates a `phylo` object with no leaves.

- `DropTip()` gains new methods `DropTip.list()` and `DropTip.NULL()`.

- `as.matrix.phylo()` converts a tree to a matrix representation, allowing
  a tree to be passed as a constraint to `ImposeConstraint()`.
  
- `as.matrix.Splits()` and `as.matrix.phyDat()` methods added as synonyms to
  `as.logical.Splits()` and `PhyDatToMatrix()`.

## Improvements

- Handle `TipLabels(0)` and `BalancedTree(0)`.

- Support zero-leaf trees in `as.Splits()` and `duplicated.Splits()`.

- Support non-identical tip labels in `as.Splits()`.

- Try Latin-1 encoding if `ReadCharacters()` family fail under UTF-8.


# TreeTools 1.8.0 (2022-09-15) #

## New methods and functions

- `TntOrder()` renumbers a tree's nodes to match TNT's convention.

- `head()` and `tail()` methods for Splits objects.

- Set names of splits object with `names(splits) <- ...`.

- `as.Splits()` support character vectors in the form "...***".

## Improvements

- `ReadTntTree()` reads tree tags and follows TNT node numbering conventions.

- `SpectrumLegend()` gains `title` parameter and more styling options.

- Support > 32767 trees in `Consensus()` 
  ([#127](https://github.com/ms609/TreeTools/issues/127)).

- `DropTip()` speed improved when branch lengths are present.


# TreeTools 1.7.3 (2022-07-20) #

- `ReadTntTree()` supports multi-line trees.

- `as.MixedBase()` supports larger trees (44-32767 tips).

- Add deprecation warning to `in.Splits()`.


# TreeTools 1.7.2 (2022-05-24) #

- `RenumberTips()` drops "preorder" attribute, as reordering tip labels may
  break edge ordering guarantee.

- Native implementation of `ClusterTable` class.

- Replace `throw` with `stop` in C++ scripts.


# TreeTools 1.7.1 (2022-03-25) #

- `AddTip()`: Fix bug when adding tip to root of weighted tree.


# TreeTools 1.7.0 (2022-03-21) #

## New methods and functions

- `rev.Splits()` reverses order in which splits are listed.

- `KeepTip.Splits()` is a faster alternative to `SubSplit()`.

- `%in%.Splits()` retains names when comparing small splits
  ([#40](https://github.com/ms609/TreeTools/issues/40)).

- `sort.multiPhylo()` sorts lists of trees according to their mixed base
  representation ([#84](https://github.com/ms609/TreeTools/issues/84)).
  
- Bitwise manipulation of splits with `|`, `&`, `xor`.

- `as.MixedBase()` uniquely represents binary trees as a mixed-base vector.

- `PathLengths()` describes all paths within a tree.

- `KeptVerts()` and `KeptPaths()` identify elements in reduced trees.

- `PostorderOrder()` describes a sequence of edges corresponding to a
  postorder traversal of a tree.

- `SpectrumLegend()` adds gradients to plot legends.


## Improvements

- Improve handling of zero-split trees.

- `DropTip()` no longer adds a root to unrooted trees, and retains edge lengths.

- Improve speed of `DropTip()`, by an order of magnitude in some cases.

- Support edge lengths in `Preorder()`, `RootTree()`, `UnrootTree()` and
  `Postorder()` ([#49](https://github.com/ms609/TreeTools/issues/49),
  [#89](https://github.com/ms609/TreeTools/issues/89)).

- Fix bug when tree is rooted on a discontinuous outgroup.

- `SortTree()` handles weighted and non-binary trees
  ([#25](https://github.com/ms609/TreeTools/issues/25),
  [#25](https://github.com/ms609/TreeTools/issues/49)),
  and gains option to sort by tip labels.

- `TipsInSplits(smallest = TRUE)` counts tips in smaller bipartition.

- Fix a bug with `phyDat` objects in `ArtificialExtinction()`.

- `RenumberTips()` allows `tipOrder` to contain elements not present in `tree`.

- Use lighter Rcpp headers.

- Small improvements to computational efficiency.

## Deprecations

- Remove deprecated function `PostorderEdges()`
  ([#35](https://github.com/ms609/TreeTools/issues/35)).


# TreeTools 1.6.0 (2022-01-12) #

## New functions

- `RoguePlot()` plots the positions of rogue taxa.


## Improvements

- `DropTip()` gains `check` parameter to allow slightly faster operation where
  input is guaranteed to be valid.

- `RandomTree()` gains `nodes` parameter allow the inclusion of polytomies.

- Infer `tips` parameter if missing in `StringToPhyDat()`.

- Remove dependency on "phangorn" (allowing use on R < 4.1)

- Improve parsing of information from nexus files.

- Export `DropTipPhylo()` as wrapper to `DropTip.phylo()`.


# TreeTools 1.5.1 (2021-10-06) #

- `PhyDatToMatrix()` optionally encodes ambiguous / inapplicable tokens as `NA`.

- Implement `sort.multiPhylo()`.

- Update test suite for compatibility with "testthat" > 3.0.4 (@hadley, #83).


# TreeTools 1.5.0 (2021-09-13) #

## New functions

- `ConstrainedNJ()` returns an approximation to a neighbour-joining tree
  that respects constraints.

- `PolarizeSplits()` marks a specified taxon as representing the ingroup of all
  splits.

- Add `KeepTip()` and improve performance of `DropTip()`.

- `ImposeConstraint()` makes a tree consistent with topological constraints.

- `as.phylo.Splits()` represents a `Splits` object as a tree.

- `Consensus()` is a faster C++ implementation of `ape::consensus()`.

- `ClusterTable()` C++ functionality imported from "TreeDist".


## Improved functions

- Warn when empty cells passed to `MatrixToPhyDat()`.

- Warn when `LabelSplits(labels)` lack names.

- `SplitFrequency()` drops tips from `forest` that aren't in `reference`.

- `AddTipEverywhere()` supports trees with < 3 leaves.

- Make `RootTree()` and `PhyDatToMatrix()` more robust.

- Support `encoding` option in `ReadCharacters()` function family.

- Support `CHARSTATELABELS` in `ReadCharacters()`.

- Support for more formatting quirks in `ReadNotes()`.

- Better support ambiguous tokens in `WriteTntCharacters()`.


## Optimization

- Fast matching functions from "fastmatch".

- Improve efficiency of `Preorder()` and `Postorder()`,
  and lift limit on tree size.


# TreeTools 1.4.5 (2021-06-23) #

- Correct calculation of minimum value in `TCIContext()`.
- Extract tip labels from objects in `StringToPhyDat()`.
- Support `AddTip(tree, where = "tip name")`.
- `SplitFrequency()` supports four-leaf trees.
- Add `RootTree.matrix()` method for edge matrices.
- Add `TipLabels.phyDat()` method.
- Add `NULL` methods for tree reordering functions.
- Additions and improvements to text parsing functions.


# TreeTools 1.4.4 (2021-04-23) #

- Add `NTip.phyDat()` method.
- Update `MakeTreeBinary()` docs and tests to reflect updated behaviour of 
  `ape::multi2di()` in 'ape' v5.5.


# TreeTools 1.4.3 (2021-04-12) #

 - `AddTip()` supports edge lengths.
 - `CladisticInfo()` supports `Splits` objects.
 - `as.multiPhylo()` converts trees, datasets and Splits objects into 
     `multiPhylo` objects.
 - `LabelSplits(labels = NULL)` labels each split with its associated node.
 - `PhyDatToMatrix()` supports integer-only levels.
 - `SortTree()` supports lists of trees.
 - Improvements to `ReadTntCharacters()` character block extraction
   ([#50](https://github.com/ms609/TreeTools/issues/50)).


# TreeTools 1.4.2 (2021-01-26) #

 - Support star trees in `RootTree()`.
 - Improve memory handling in `root_on_node()`.
 - Documentation linkage.


# TreeTools 1.4.1 (2020-12-09) #

 - `MSTEdges()` supports distance matrices with > 256 entries.
 - Package 'vdiffr' used conditionally.


# TreeTools 1.4.0 (2020-10-20) #

## New functions
 - `MSTLength()` reports length of minimum spanning tree.
 - `AllTipLabels()` returns all labels from all trees in a list.
 - `PairwiseDistances()` (from 'TreeDistData') computes distances between all
   pairs of trees in a list.
 - `ArtificialExtinction()` gains `replaceAll` option.
 - `WriteTntCharacters(types = ...)` writes different character types to TNT 
   file.
 - Tree characterization S3 methods: add `.default` and `.NULL`.

## Enhancements
 - `MSTEdges()` implemented in C++, improving runtime by orders of magnitude.
 - Improved parsing of TNT character files.


# TreeTools 1.3.1 #

 - Improved parsing of TNT files.
 - Fix misspecified C++ linkage.


# TreeTools 1.3.0 (2020-09-22) #

## New functions
 
 - `SisterSize()` and `RootNodeDist()` measure sister-clade size and root-node
   distance.
 - `MSTEdges()`: Edges of minimum spanning tree.
 - `SplitImbalance()`: how balanced is each split?
 - New C++ functions `root_on_node()` and `root_binary()` to root trees quickly
   and robustly.

## Enhancements

 - `TNTReadTree()` handles additional punctuation characters.
 - Import RdMacros package 'Rdpack'.
 
 - C++ implementation of `TipsInSplits()`.
 - Export C++ functions `preorder_edges_and_nodes()` and `postorder_edges()`.
 - Remove obsolete copy of C++ code from 'phangorn'.


# TreeTools 1.2.0 (2020-08-30) #

 - `ArtificialExtinction()`: Remove characters that are absent in a fossil 
   template.
 - `WriteTntCharacters()`: Write morphological dataset in TNT format.
 - Improve TNT dataset parsing.
 - Documentation improvements.


# TreeTools 1.1.0 (2020-07-07) #

## New functions

 - `RandomTree()`: Draw tree from uniform distribution, instead of via
 `ape::rtree()`.
 - `MakeTreeBinary()`: Uniform equivalent of `ape::multi2di()`.
 - `match.list()` method for lists of splits.
 - `SplitsInBinaryTree()`: How many splits occur in an _n_-leaf binary tree?
 - `vapply64()`, `sapply64()`, `replicate64()`: helper functions when a function
 returns a 64-bit integer.

## Enhancements

 - Use methods for `UnrootTree()`, `RootTree()`, `RootOnNode()` to support
 lists of trees.


# TreeTools 1.0.0 (2020-06-08) #

## New functions

- `CladisticInfo()`: Calculate the information content of a tree.
- `RootNode()`: Which node is a tree's root?
- `UnrootTree()`: Safely remove a root node.
- `NodeDepth()`: Discriminate shallow from deep nodes.
- `NodeOrder()`, `NDescendants()`: Count edges incident to each node.
- `CladeSizes()`: Count leaves / nodes descended from each node.
- `ListAncestors()`: List ancestors of a node.
- `LabelSplits()`: Label splits on plotted tree.
- `DropTip()`: Remove tip, handling weird node orders.
- `LeafLabelInterchange()`: Exchange position of _n_ tips.
- `StarTree()`: Generate unresolved tree.
- `TotalCopheneticIndex()` integrated from 'tci' package.

## Deprecations

- `PostorderEdges()`: use `Postorder()` instead.
- `NonDuplicateRoot()`: unused internal function.
- `match.Splits()`: use `match()` instead.
- `in.Splits()`: use `%in%.Splits()` instead.

## Enhancements

- Improve support for unrooted trees in `as.Splits()`.
- Use methods so `Reorder` functions can handle `multiPhylo` objects and edges.
- Handle funny node orders.
- Support continuous characters in `ReadCharacters()`.
- Improve performance of `as.logical.Splits()` and related functions.
- Fail nicely when trees are too large for memory.
- Fix memory leak in `as.Splits()`.
- Various under-the-hood improvements to functions.
- Documentation improvements.


# TreeTools 0.1.4 (2020-03-04) #

- Catch hang-inducing bugs in `RootOnNode()`.
- Update `doubleFactorials` cache to fix `as.integer()` rounding error.
- Support unrooted trees in `AddTipEverywhere()`.
- Documentation improvements.


# TreeTools 0.1.3 (2020-01-07) #

- `RootOnNode()`: Quickly root a tree on a specified node.
- Improve portability of C++ code.


# TreeTools 0.1.2 (2020-12-18) #
 
- `as.Newick`: Fast conversion to Newick format.
- `as.TreeNumber`: Tree shape enumeration.


# TreeTools 0.1.1 #
 
- Add functions to translate trees to mixed base integers.
- `RenumberTips` can extract tip order from `phylo` and `Splits` objects.
- Documentation changes to satisfy CRAN submission requirements.


# TreeTools 0.1.0 (2019-10-30) #

- Pre-release version spun out of ['TreeSearch'](https://ms609.github.io/TreeSearch/)
  package.  Some functionality is subject to change.
