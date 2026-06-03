# sftime 0.3.2

* Switch from the magrittr pipe (`%>%`) to the R-native pipe (`|>`) in examples and vignettes (#15).

# sftime 0.3.1

* Correct argument `versionCheck` the `requireNamespace` for the `cubble` package in `st_as_sftime.cubble_df()`.

# sftime 0.3.0

* Add a dedicated `tidyr::drop_na()` method for `sftime` objects. (See the same recent addition for `sf` objects [#1975](https://github.com/r-spatial/sf/pull/1975/)).

* Add a dedicated `dplyr::dplyr_reconstruct()` method for `sftime` objects. 
Relying on the method for `sf` objects caused erroneously column binding when the second object was a data frame without conflicting column names for the `sf` and time columns. In this case, a `sf` objects was returned, even though an `sftime` object should be returned. See also https://github.com/r-spatial/sf/issues/1958#issuecomment-1181982244.

* Add methods to convert `sftime` objects from:
  + Objects from the `spatstat` package classes (`ppp`, `psp`, `lpp`)
  + `sftrack` and `sftraj` objects from the `sftrack` package.
  + `cubble_df` objects from the `cubble` package.

* Bug fix in `st_time<-.sftime`:  
  + Still contained references to the old `tc`class.
  + Did not allow to give the active time column a character vector as value.

# version 0.2-0

* initial CRAN submission
