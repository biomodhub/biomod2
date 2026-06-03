# pixmap 0.4-14

* In test `bugs.R` only show differences (if any) rather than making the
  test fail.

* Improve `rep()` calls in `pixmap()` to avoid replicating `NULL` vectors
  and to avoid partial argument matching.


# pixmap 0.4-13

* Achim Zeileis takes over maintenance from Friedrich Leisch.


# pixmap 0.4-12

* Fixed some `NAMESPACE` problems.


# pixmap 0.4-11

* Added dummy `NAMESPACE` file.

* Move class definitions to separate file, no `Collate` in `DESCRIPTION`.


# pixmap 0.4-10

* Explicitly mention in the help file that `read.pnm` only works for 
  files, not other conntection.

* Fixed a bug in `write.pnm` that tried to open the same connection twice.


# pixmap 0.4-9

* Fixed a minor glitch in `write.pnm`.


# pixmap 0.4-8

* Fixed a bug that prevented plotting images with only 1 row (bug
  report by Robert Esswein).


# pixmap 0.4-7

* Use `LazyLoad` instead of `SaveImage`.


# pixmap 0.4-6

* Modified one of the regression tests for changes in R 2.4: 
  `terrain.colors()` now return transparency information 
  -> do not use it in example.


# pixmap 0.4-5

* Standardized license field in `DESCRIPTION` file.


# pixmap 0.4-4

* Fixed a bug in the prototype of class `"pixmap"`.


# pixmap 0.4-3

* New example for overlaying plots in `help(pixmap)` submitted by 
  Stephan Matthiesen.


# pixmap 0.4-2

* Adjust for R 2.0.0.

* Fixed a bug in coercion from `pixmapIndexed` to `pixmapRGB`.

* There was a bug in the `methods` package of R 1.9.x which was 
  triggered by functions in `pixmap`, hence this version of the
  package depends on R >= 2.0.0. 


# pixmap 0.4-1

* Fixed a bug in `write.pnm()` that wrote grey images in PPM format.

* The channel information was not changed when converting between
  RGB and grey pixmaps.

* The `maxval` in PNM headers must be less than `65536`, not less than `256`.


# pixmap 0.4-0

* `read.pnm()`: Vectorized (and renamed)  `as.integer.bytes()` which provides
  a huge performance gain for reading "PBM" (b/w bitmaps).

* New function `addlogo()`.


# pixmap 0.3-4

* Clarified documentation of `read.pnm` (file name extensions are ignored).


# pixmap 0.3-3

* Fixed some codoc problems (missing aliases).


# pixmap 0.3-1

* `read.pnm()`: Changes made to function reading and parsing PNM file
  headers to permit comments of arbitrary length.

  
# pixmap 0.3-0

* The whole package has moved to S4 classes and methods, hence all
  classes have a new representation. This also means that the code
  and the API are not fully backwards compatible with earlier
  versions of the package.

* Added support for subsetting, see `example(pixmap)`.

* New: `addChannels()` and `getChannels()`.


# pixmap 0.2-1

* pixmap():
  - Added arguments `bbcent` and `cellres`.
  - `nrow` and `ncol` default to the respective dimensions of the data
    argument (if present). Hence, pixmap does the expected when
    given a matrix or an array.
  - `data` is rescaled to [0,1] for RGB and grey, and coerced to positive
    integers for indexed.
  - `col` can also be a function like `rainbow()`

  
# pixmap 0.1-2

* Fixed bugs in plotting, `read.pnm` and `write.pnm` which confused
  dimensions (rows versus columns), but together let plots look OK.

