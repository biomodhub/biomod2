# gridtext 0.1.6

- Removed SystemRequirements from package DESCRIPTION to fix CRAN NOTE

# gridtext 0.1.5

- Transition to curl package and drop RCurl dependency
- Fix fontface not being processed and words spaced properly in R 4.2.0
- Maintainer changes to Brenton Wiernik
- Removed LazyData from package DESCRIPTION to fix CRAN NOTE

# gridtext 0.1.4

- Make sure tests don't fail if vdiffr is missing.

# gridtext 0.1.3

- Remove unneeded systemfonts dependency.

# gridtext 0.1.2

- Fix build for testthat 3.0.

# gridtext 0.1.1

- `richtext_grob()` and `textbox_grob()` now gracefully handle empty strings
  and NAs.

# gridtext 0.1.0

First public release. Provides the two grobs `richtext_grob()` and `textbox_grob()` for formatted text rendering without and with word wrapping, respectively.
