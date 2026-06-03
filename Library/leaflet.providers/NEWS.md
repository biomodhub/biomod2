# leaflet.providers 3.0.0

## New features

* Updated bundled leaflet-providers.js from v2.0.0 to v3.0.0 (#44).

* `get_providers()` now validates `version_num` and errors if the
  version is less than `"1.0.10"` (#44).

## Bug fixes and minor improvements

* `providers_default()` now returns an `htmltools::htmlDependency()`
  that points to the packaged file directly, rather than a temporary
  copy (#43).

# leaflet.providers 2.0.0

* Updated leaflet.providers data on 2023-10-05 from https://unpkg.com/leaflet-providers using version 2.0.0 of leaflet-providers.js

* `get_providers()` and `default_providers()` both include a stable `htmltools::htmlDependency()` in the `dep` slot of the returned object. The default dependency uses the static `leaflet-providers.js` file included in the `leaflet.providers` package, while the `get_providers()` dependency uses the `leaflet-providers.js` file from the `unpkg.com` CDN. This makes it possible for `leaflet::addProviderTiles()` to use the stable HTML dependency that can be cached by knitr (#36).

* The previous version (v1.13.0) of leaflet.providers required R >= 3.6. While versions of R older than 3.6 are not actively supported, the theoretically supported minimum R version is lower, likely R >= 2.10. We've restored the previous R version minimum requirement. (#33)

# leaflet.providers 1.13.0

* Updated leaflet.providers data on 2023-08-07 from https://unpkg.com/leaflet-providers using version 1.13.0 of leaflet-providers.js


# leaflet.providers 1.9.0
* Updated leaflet.providers data on 2019-11-06 from https://unpkg.com/leaflet-providers using version 1.9.0 of leaflet-providers.js


# leaflet.providers 1.8.0

* Added a `NEWS.md` file to track changes to the package.
* Initialized package (version number 1.8.0 matches version number of leaflet-providers.js).
