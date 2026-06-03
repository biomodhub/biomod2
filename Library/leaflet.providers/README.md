
<!-- README.md is generated from README.Rmd. Please edit that file -->

# leaflet.providers

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/leaflet.providers)](https://CRAN.R-project.org/package=leaflet.providers)
[![R-CMD-check](https://github.com/rstudio/leaflet.providers/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rstudio/leaflet.providers/actions)
[![Codecov test
coverage](https://codecov.io/gh/rstudio/leaflet.providers/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rstudio/leaflet.providers?branch=main)

<!-- badges: end -->

The goal of `leaflet.providers` is to provide regularly updated data on
the third-party tile providers supported by `leaflet`. The data is
extracted from
[leaflet-providers.js](https://github.com/leaflet-extras/leaflet-providers).

While `leaflet.providers` will be regularly released with updated
providers data, the package comes with functions `use_providers()` and
`get_providers()`, which enable users to fetch up-to-date providers
information directly from
[leaflet-providers.js](https://github.com/leaflet-extras/leaflet-providers)
between package updates and to load this providers data in `leaflet`.
Users may also fetch older versions of the providers data with the
`leaflet-providers.js` version number.

## Installation

You can install the released version of `leaflet.providers` from
[CRAN](https://CRAN.R-project.org):

``` r
# CRAN version
install.packages("leaflet.providers")
```

And the development version from [GitHub](https://github.com/) with:

``` r
pak::pak("rstudio/leaflet.providers")
```

## Example

The following are basic example of the default providers data that comes
`leaflet.providers` as well as the available methods to fetch and load
custom providers.

### Usage

In general, you do not need to load or interact with `leaflet.providers`
directly. Instead, `leaflet.providers` is used indirectly by the
[leaflet package](https://rstudio.github.io/leaflet/) in the
`leaflet::addProviderTiles()` function.

``` r
library(leaflet)

leaflet()  |>
  setView(lng = -71.0589, lat = 42.3601, zoom = 12) |>
  addTiles() |>
  addProviderTiles("OpenTopoMap")
```

### Fetch custom providers

You can choose a specific version of `leaflet-providers.js` by calling
`use_providers()` with the desired version or the providers from a
desired version:

``` r
# Shorthand method:
use_providers("1.7.0")

# Longer method:
# 1. Retrieve v 1.7.0
providers_170 <- get_providers("1.7.0")

# 2. Load custom providers data
use_providers(providers_170)

# Reset loaded providers to default providers
use_providers()
```

Now that `use_providers()` has been called with a custom
`leaflet.providers` object obtained via `get_providers()` (in this case,
a previous version of the data),
[`leaflet`](https://rstudio.github.io/leaflet/) will use the custom
providers instead of the default provider data.

> Note that the package `V8` is required for `get_providers()`.

### Default providers

``` r
library(leaflet.providers)

str(providers_default(), max.level = 2)
```

<div style="height:150px; overflow-y: scroll;">

    #> List of 5
    #>  $ version_num      : chr "3.0.0"
    #>  $ providers        :List of 224
    #>   ..$ OpenStreetMap                         : chr "OpenStreetMap"
    #>   ..$ OpenStreetMap.Mapnik                  : chr "OpenStreetMap.Mapnik"
    #>   ..$ OpenStreetMap.DE                      : chr "OpenStreetMap.DE"
    #>   ..$ OpenStreetMap.CH                      : chr "OpenStreetMap.CH"
    #>   ..$ OpenStreetMap.France                  : chr "OpenStreetMap.France"
    #>   ..$ OpenStreetMap.HOT                     : chr "OpenStreetMap.HOT"
    #>   ..$ OpenStreetMap.BZH                     : chr "OpenStreetMap.BZH"
    #>   ..$ OpenStreetMap.CAT                     : chr "OpenStreetMap.CAT"
    #>   ..$ MapTilesAPI                           : chr "MapTilesAPI"
    #>   ..$ MapTilesAPI.OSMEnglish                : chr "MapTilesAPI.OSMEnglish"
    #>   ..$ MapTilesAPI.OSMFrancais               : chr "MapTilesAPI.OSMFrancais"
    #>   ..$ MapTilesAPI.OSMEspagnol               : chr "MapTilesAPI.OSMEspagnol"
    #>   ..$ OpenSeaMap                            : chr "OpenSeaMap"
    #>   ..$ OPNVKarte                             : chr "OPNVKarte"
    #>   ..$ OpenTopoMap                           : chr "OpenTopoMap"
    #>   ..$ OpenRailwayMap                        : chr "OpenRailwayMap"
    #>   ..$ OpenFireMap                           : chr "OpenFireMap"
    #>   ..$ SafeCast                              : chr "SafeCast"
    #>   ..$ Stadia                                : chr "Stadia"
    #>   ..$ Stadia.AlidadeSmooth                  : chr "Stadia.AlidadeSmooth"
    #>   ..$ Stadia.AlidadeSmoothDark              : chr "Stadia.AlidadeSmoothDark"
    #>   ..$ Stadia.AlidadeSatellite               : chr "Stadia.AlidadeSatellite"
    #>   ..$ Stadia.OSMBright                      : chr "Stadia.OSMBright"
    #>   ..$ Stadia.Outdoors                       : chr "Stadia.Outdoors"
    #>   ..$ Stadia.StamenToner                    : chr "Stadia.StamenToner"
    #>   ..$ Stadia.StamenTonerBackground          : chr "Stadia.StamenTonerBackground"
    #>   ..$ Stadia.StamenTonerLines               : chr "Stadia.StamenTonerLines"
    #>   ..$ Stadia.StamenTonerLabels              : chr "Stadia.StamenTonerLabels"
    #>   ..$ Stadia.StamenTonerLite                : chr "Stadia.StamenTonerLite"
    #>   ..$ Stadia.StamenWatercolor               : chr "Stadia.StamenWatercolor"
    #>   ..$ Stadia.StamenTerrain                  : chr "Stadia.StamenTerrain"
    #>   ..$ Stadia.StamenTerrainBackground        : chr "Stadia.StamenTerrainBackground"
    #>   ..$ Stadia.StamenTerrainLabels            : chr "Stadia.StamenTerrainLabels"
    #>   ..$ Stadia.StamenTerrainLines             : chr "Stadia.StamenTerrainLines"
    #>   ..$ Thunderforest                         : chr "Thunderforest"
    #>   ..$ Thunderforest.OpenCycleMap            : chr "Thunderforest.OpenCycleMap"
    #>   ..$ Thunderforest.Transport               : chr "Thunderforest.Transport"
    #>   ..$ Thunderforest.TransportDark           : chr "Thunderforest.TransportDark"
    #>   ..$ Thunderforest.SpinalMap               : chr "Thunderforest.SpinalMap"
    #>   ..$ Thunderforest.Landscape               : chr "Thunderforest.Landscape"
    #>   ..$ Thunderforest.Outdoors                : chr "Thunderforest.Outdoors"
    #>   ..$ Thunderforest.Pioneer                 : chr "Thunderforest.Pioneer"
    #>   ..$ Thunderforest.MobileAtlas             : chr "Thunderforest.MobileAtlas"
    #>   ..$ Thunderforest.Neighbourhood           : chr "Thunderforest.Neighbourhood"
    #>   ..$ BaseMapDE                             : chr "BaseMapDE"
    #>   ..$ BaseMapDE.Color                       : chr "BaseMapDE.Color"
    #>   ..$ BaseMapDE.Grey                        : chr "BaseMapDE.Grey"
    #>   ..$ CyclOSM                               : chr "CyclOSM"
    #>   ..$ Jawg                                  : chr "Jawg"
    #>   ..$ Jawg.Streets                          : chr "Jawg.Streets"
    #>   ..$ Jawg.Terrain                          : chr "Jawg.Terrain"
    #>   ..$ Jawg.Lagoon                           : chr "Jawg.Lagoon"
    #>   ..$ Jawg.Sunny                            : chr "Jawg.Sunny"
    #>   ..$ Jawg.Dark                             : chr "Jawg.Dark"
    #>   ..$ Jawg.Light                            : chr "Jawg.Light"
    #>   ..$ Jawg.Matrix                           : chr "Jawg.Matrix"
    #>   ..$ MapBox                                : chr "MapBox"
    #>   ..$ MapTiler                              : chr "MapTiler"
    #>   ..$ MapTiler.Streets                      : chr "MapTiler.Streets"
    #>   ..$ MapTiler.Basic                        : chr "MapTiler.Basic"
    #>   ..$ MapTiler.Bright                       : chr "MapTiler.Bright"
    #>   ..$ MapTiler.Pastel                       : chr "MapTiler.Pastel"
    #>   ..$ MapTiler.Positron                     : chr "MapTiler.Positron"
    #>   ..$ MapTiler.Hybrid                       : chr "MapTiler.Hybrid"
    #>   ..$ MapTiler.Toner                        : chr "MapTiler.Toner"
    #>   ..$ MapTiler.Topo                         : chr "MapTiler.Topo"
    #>   ..$ MapTiler.Voyager                      : chr "MapTiler.Voyager"
    #>   ..$ MapTiler.Ocean                        : chr "MapTiler.Ocean"
    #>   ..$ MapTiler.Backdrop                     : chr "MapTiler.Backdrop"
    #>   ..$ MapTiler.Dataviz                      : chr "MapTiler.Dataviz"
    #>   ..$ MapTiler.DatavizLight                 : chr "MapTiler.DatavizLight"
    #>   ..$ MapTiler.DatavizDark                  : chr "MapTiler.DatavizDark"
    #>   ..$ MapTiler.Aquarelle                    : chr "MapTiler.Aquarelle"
    #>   ..$ MapTiler.Landscape                    : chr "MapTiler.Landscape"
    #>   ..$ MapTiler.Openstreetmap                : chr "MapTiler.Openstreetmap"
    #>   ..$ MapTiler.Outdoor                      : chr "MapTiler.Outdoor"
    #>   ..$ MapTiler.Satellite                    : chr "MapTiler.Satellite"
    #>   ..$ MapTiler.Winter                       : chr "MapTiler.Winter"
    #>   ..$ TomTom                                : chr "TomTom"
    #>   ..$ TomTom.Basic                          : chr "TomTom.Basic"
    #>   ..$ TomTom.Hybrid                         : chr "TomTom.Hybrid"
    #>   ..$ TomTom.Labels                         : chr "TomTom.Labels"
    #>   ..$ Esri                                  : chr "Esri"
    #>   ..$ Esri.WorldStreetMap                   : chr "Esri.WorldStreetMap"
    #>   ..$ Esri.WorldTopoMap                     : chr "Esri.WorldTopoMap"
    #>   ..$ Esri.WorldImagery                     : chr "Esri.WorldImagery"
    #>   ..$ Esri.WorldTerrain                     : chr "Esri.WorldTerrain"
    #>   ..$ Esri.WorldShadedRelief                : chr "Esri.WorldShadedRelief"
    #>   ..$ Esri.WorldPhysical                    : chr "Esri.WorldPhysical"
    #>   ..$ Esri.OceanBasemap                     : chr "Esri.OceanBasemap"
    #>   ..$ Esri.NatGeoWorldMap                   : chr "Esri.NatGeoWorldMap"
    #>   ..$ Esri.WorldGrayCanvas                  : chr "Esri.WorldGrayCanvas"
    #>   ..$ OpenWeatherMap                        : chr "OpenWeatherMap"
    #>   ..$ OpenWeatherMap.Clouds                 : chr "OpenWeatherMap.Clouds"
    #>   ..$ OpenWeatherMap.CloudsClassic          : chr "OpenWeatherMap.CloudsClassic"
    #>   ..$ OpenWeatherMap.Precipitation          : chr "OpenWeatherMap.Precipitation"
    #>   ..$ OpenWeatherMap.PrecipitationClassic   : chr "OpenWeatherMap.PrecipitationClassic"
    #>   ..$ OpenWeatherMap.Rain                   : chr "OpenWeatherMap.Rain"
    #>   ..$ OpenWeatherMap.RainClassic            : chr "OpenWeatherMap.RainClassic"
    #>   .. [list output truncated]
    #>  $ providers_details:List of 37
    #>   ..$ OpenStreetMap        :List of 3
    #>   ..$ MapTilesAPI          :List of 3
    #>   ..$ OpenSeaMap           :List of 2
    #>   ..$ OPNVKarte            :List of 2
    #>   ..$ OpenTopoMap          :List of 2
    #>   ..$ OpenRailwayMap       :List of 2
    #>   ..$ OpenFireMap          :List of 2
    #>   ..$ SafeCast             :List of 2
    #>   ..$ Stadia               :List of 3
    #>   ..$ Thunderforest        :List of 3
    #>   ..$ BaseMapDE            :List of 3
    #>   ..$ CyclOSM              :List of 2
    #>   ..$ Jawg                 :List of 3
    #>   ..$ MapBox               :List of 2
    #>   ..$ MapTiler             :List of 3
    #>   ..$ TomTom               :List of 3
    #>   ..$ Esri                 :List of 3
    #>   ..$ OpenWeatherMap       :List of 3
    #>   ..$ HERE                 :List of 3
    #>   ..$ FreeMapSK            :List of 2
    #>   ..$ MtbMap               :List of 2
    #>   ..$ CartoDB              :List of 3
    #>   ..$ HikeBike             :List of 3
    #>   ..$ BasemapAT            :List of 3
    #>   ..$ nlmaps               :List of 3
    #>   ..$ NASAGIBS             :List of 3
    #>   ..$ NLS                  :List of 3
    #>   ..$ JusticeMap           :List of 3
    #>   ..$ GeoportailFrance     :List of 3
    #>   ..$ OneMapSG             :List of 3
    #>   ..$ USGS                 :List of 3
    #>   ..$ WaymarkedTrails      :List of 3
    #>   ..$ OpenAIP              :List of 2
    #>   ..$ OpenSnowMap          :List of 3
    #>   ..$ AzureMaps            :List of 3
    #>   ..$ SwissFederalGeoportal:List of 3
    #>   ..$ TopPlusOpen          :List of 3
    #>  $ src              : chr "(function(root, factory) {\n\tif (typeof define === 'function' && define.amd) {\n\t\t// AMD. Register as an ano"| __truncated__
    #>  $ dep              :List of 10
    #>   ..$ name      : chr "leaflet-providers"
    #>   ..$ version   : chr "3.0.0"
    #>   ..$ src       :List of 1
    #>   ..$ meta      : NULL
    #>   ..$ script    : chr "leaflet-providers.js"
    #>   ..$ stylesheet: NULL
    #>   ..$ head      : NULL
    #>   ..$ attachment: NULL
    #>   ..$ package   : chr "leaflet.providers"
    #>   ..$ all_files : logi FALSE
    #>   ..- attr(*, "class")= chr "html_dependency"
    #>  - attr(*, "class")= chr "leaflet_providers"

</div>

### View the loaded providers data

#### The version number of the source leaflet-providers.js

``` r
providers_loaded()$version_num
#> [1] "3.0.0"
```

#### Supported tile providers

``` r
names(providers_loaded()$providers)
```

<div style="height:150px; overflow-y: scroll;">

    #>   [1] "OpenStreetMap"                         
    #>   [2] "OpenStreetMap.Mapnik"                  
    #>   [3] "OpenStreetMap.DE"                      
    #>   [4] "OpenStreetMap.CH"                      
    #>   [5] "OpenStreetMap.France"                  
    #>   [6] "OpenStreetMap.HOT"                     
    #>   [7] "OpenStreetMap.BZH"                     
    #>   [8] "OpenStreetMap.CAT"                     
    #>   [9] "MapTilesAPI"                           
    #>  [10] "MapTilesAPI.OSMEnglish"                
    #>  [11] "MapTilesAPI.OSMFrancais"               
    #>  [12] "MapTilesAPI.OSMEspagnol"               
    #>  [13] "OpenSeaMap"                            
    #>  [14] "OPNVKarte"                             
    #>  [15] "OpenTopoMap"                           
    #>  [16] "OpenRailwayMap"                        
    #>  [17] "OpenFireMap"                           
    #>  [18] "SafeCast"                              
    #>  [19] "Stadia"                                
    #>  [20] "Stadia.AlidadeSmooth"                  
    #>  [21] "Stadia.AlidadeSmoothDark"              
    #>  [22] "Stadia.AlidadeSatellite"               
    #>  [23] "Stadia.OSMBright"                      
    #>  [24] "Stadia.Outdoors"                       
    #>  [25] "Stadia.StamenToner"                    
    #>  [26] "Stadia.StamenTonerBackground"          
    #>  [27] "Stadia.StamenTonerLines"               
    #>  [28] "Stadia.StamenTonerLabels"              
    #>  [29] "Stadia.StamenTonerLite"                
    #>  [30] "Stadia.StamenWatercolor"               
    #>  [31] "Stadia.StamenTerrain"                  
    #>  [32] "Stadia.StamenTerrainBackground"        
    #>  [33] "Stadia.StamenTerrainLabels"            
    #>  [34] "Stadia.StamenTerrainLines"             
    #>  [35] "Thunderforest"                         
    #>  [36] "Thunderforest.OpenCycleMap"            
    #>  [37] "Thunderforest.Transport"               
    #>  [38] "Thunderforest.TransportDark"           
    #>  [39] "Thunderforest.SpinalMap"               
    #>  [40] "Thunderforest.Landscape"               
    #>  [41] "Thunderforest.Outdoors"                
    #>  [42] "Thunderforest.Pioneer"                 
    #>  [43] "Thunderforest.MobileAtlas"             
    #>  [44] "Thunderforest.Neighbourhood"           
    #>  [45] "BaseMapDE"                             
    #>  [46] "BaseMapDE.Color"                       
    #>  [47] "BaseMapDE.Grey"                        
    #>  [48] "CyclOSM"                               
    #>  [49] "Jawg"                                  
    #>  [50] "Jawg.Streets"                          
    #>  [51] "Jawg.Terrain"                          
    #>  [52] "Jawg.Lagoon"                           
    #>  [53] "Jawg.Sunny"                            
    #>  [54] "Jawg.Dark"                             
    #>  [55] "Jawg.Light"                            
    #>  [56] "Jawg.Matrix"                           
    #>  [57] "MapBox"                                
    #>  [58] "MapTiler"                              
    #>  [59] "MapTiler.Streets"                      
    #>  [60] "MapTiler.Basic"                        
    #>  [61] "MapTiler.Bright"                       
    #>  [62] "MapTiler.Pastel"                       
    #>  [63] "MapTiler.Positron"                     
    #>  [64] "MapTiler.Hybrid"                       
    #>  [65] "MapTiler.Toner"                        
    #>  [66] "MapTiler.Topo"                         
    #>  [67] "MapTiler.Voyager"                      
    #>  [68] "MapTiler.Ocean"                        
    #>  [69] "MapTiler.Backdrop"                     
    #>  [70] "MapTiler.Dataviz"                      
    #>  [71] "MapTiler.DatavizLight"                 
    #>  [72] "MapTiler.DatavizDark"                  
    #>  [73] "MapTiler.Aquarelle"                    
    #>  [74] "MapTiler.Landscape"                    
    #>  [75] "MapTiler.Openstreetmap"                
    #>  [76] "MapTiler.Outdoor"                      
    #>  [77] "MapTiler.Satellite"                    
    #>  [78] "MapTiler.Winter"                       
    #>  [79] "TomTom"                                
    #>  [80] "TomTom.Basic"                          
    #>  [81] "TomTom.Hybrid"                         
    #>  [82] "TomTom.Labels"                         
    #>  [83] "Esri"                                  
    #>  [84] "Esri.WorldStreetMap"                   
    #>  [85] "Esri.WorldTopoMap"                     
    #>  [86] "Esri.WorldImagery"                     
    #>  [87] "Esri.WorldTerrain"                     
    #>  [88] "Esri.WorldShadedRelief"                
    #>  [89] "Esri.WorldPhysical"                    
    #>  [90] "Esri.OceanBasemap"                     
    #>  [91] "Esri.NatGeoWorldMap"                   
    #>  [92] "Esri.WorldGrayCanvas"                  
    #>  [93] "OpenWeatherMap"                        
    #>  [94] "OpenWeatherMap.Clouds"                 
    #>  [95] "OpenWeatherMap.CloudsClassic"          
    #>  [96] "OpenWeatherMap.Precipitation"          
    #>  [97] "OpenWeatherMap.PrecipitationClassic"   
    #>  [98] "OpenWeatherMap.Rain"                   
    #>  [99] "OpenWeatherMap.RainClassic"            
    #> [100] "OpenWeatherMap.Pressure"               
    #> [101] "OpenWeatherMap.PressureContour"        
    #> [102] "OpenWeatherMap.Wind"                   
    #> [103] "OpenWeatherMap.Temperature"            
    #> [104] "OpenWeatherMap.Snow"                   
    #> [105] "HERE"                                  
    #> [106] "HERE.exploreDay"                       
    #> [107] "HERE.liteDay"                          
    #> [108] "HERE.logisticsDay"                     
    #> [109] "HERE.topoDay"                          
    #> [110] "HERE.logisticsNight"                   
    #> [111] "HERE.exploreNight"                     
    #> [112] "HERE.topoNight"                        
    #> [113] "HERE.liteNight"                        
    #> [114] "HERE.exploreSatelliteDay"              
    #> [115] "HERE.liteSatelliteDay"                 
    #> [116] "HERE.logisticsSatelliteDay"            
    #> [117] "HERE.basicMap"                         
    #> [118] "HERE.mapLabels"                        
    #> [119] "HERE.trafficFlow"                      
    #> [120] "HERE.carnavDayGrey"                    
    #> [121] "HERE.hybridDay"                        
    #> [122] "HERE.hybridDayMobile"                  
    #> [123] "HERE.hybridDayTransit"                 
    #> [124] "HERE.hybridDayGrey"                    
    #> [125] "HERE.pedestrianDay"                    
    #> [126] "HERE.pedestrianNight"                  
    #> [127] "HERE.satelliteDay"                     
    #> [128] "HERE.terrainDay"                       
    #> [129] "HERE.terrainDayMobile"                 
    #> [130] "FreeMapSK"                             
    #> [131] "MtbMap"                                
    #> [132] "CartoDB"                               
    #> [133] "CartoDB.Positron"                      
    #> [134] "CartoDB.PositronNoLabels"              
    #> [135] "CartoDB.PositronOnlyLabels"            
    #> [136] "CartoDB.DarkMatter"                    
    #> [137] "CartoDB.DarkMatterNoLabels"            
    #> [138] "CartoDB.DarkMatterOnlyLabels"          
    #> [139] "CartoDB.Voyager"                       
    #> [140] "CartoDB.VoyagerNoLabels"               
    #> [141] "CartoDB.VoyagerOnlyLabels"             
    #> [142] "CartoDB.VoyagerLabelsUnder"            
    #> [143] "HikeBike"                              
    #> [144] "HikeBike.HikeBike"                     
    #> [145] "HikeBike.HillShading"                  
    #> [146] "BasemapAT"                             
    #> [147] "BasemapAT.basemap"                     
    #> [148] "BasemapAT.grau"                        
    #> [149] "BasemapAT.overlay"                     
    #> [150] "BasemapAT.terrain"                     
    #> [151] "BasemapAT.surface"                     
    #> [152] "BasemapAT.highdpi"                     
    #> [153] "BasemapAT.orthofoto"                   
    #> [154] "nlmaps"                                
    #> [155] "nlmaps.standaard"                      
    #> [156] "nlmaps.pastel"                         
    #> [157] "nlmaps.grijs"                          
    #> [158] "nlmaps.water"                          
    #> [159] "nlmaps.luchtfoto"                      
    #> [160] "NASAGIBS"                              
    #> [161] "NASAGIBS.ModisTerraTrueColorCR"        
    #> [162] "NASAGIBS.ModisTerraBands367CR"         
    #> [163] "NASAGIBS.ViirsEarthAtNight2012"        
    #> [164] "NASAGIBS.ModisTerraLSTDay"             
    #> [165] "NASAGIBS.ModisTerraSnowCover"          
    #> [166] "NASAGIBS.ModisTerraAOD"                
    #> [167] "NASAGIBS.ModisTerraChlorophyll"        
    #> [168] "NLS"                                   
    #> [169] "NLS.osgb63k1885"                       
    #> [170] "NLS.osgb1888"                          
    #> [171] "NLS.osgb10k1888"                       
    #> [172] "NLS.osgb1919"                          
    #> [173] "NLS.osgb25k1937"                       
    #> [174] "NLS.osgb63k1955"                       
    #> [175] "NLS.oslondon1k1893"                    
    #> [176] "JusticeMap"                            
    #> [177] "JusticeMap.income"                     
    #> [178] "JusticeMap.americanIndian"             
    #> [179] "JusticeMap.asian"                      
    #> [180] "JusticeMap.black"                      
    #> [181] "JusticeMap.hispanic"                   
    #> [182] "JusticeMap.multi"                      
    #> [183] "JusticeMap.nonWhite"                   
    #> [184] "JusticeMap.white"                      
    #> [185] "JusticeMap.plurality"                  
    #> [186] "GeoportailFrance"                      
    #> [187] "GeoportailFrance.plan"                 
    #> [188] "GeoportailFrance.parcels"              
    #> [189] "GeoportailFrance.orthos"               
    #> [190] "OneMapSG"                              
    #> [191] "OneMapSG.Default"                      
    #> [192] "OneMapSG.Night"                        
    #> [193] "OneMapSG.Original"                     
    #> [194] "OneMapSG.Grey"                         
    #> [195] "OneMapSG.LandLot"                      
    #> [196] "USGS"                                  
    #> [197] "USGS.USTopo"                           
    #> [198] "USGS.USImagery"                        
    #> [199] "USGS.USImageryTopo"                    
    #> [200] "WaymarkedTrails"                       
    #> [201] "WaymarkedTrails.hiking"                
    #> [202] "WaymarkedTrails.cycling"               
    #> [203] "WaymarkedTrails.mtb"                   
    #> [204] "WaymarkedTrails.slopes"                
    #> [205] "WaymarkedTrails.riding"                
    #> [206] "WaymarkedTrails.skating"               
    #> [207] "OpenAIP"                               
    #> [208] "OpenSnowMap"                           
    #> [209] "OpenSnowMap.pistes"                    
    #> [210] "AzureMaps"                             
    #> [211] "AzureMaps.MicrosoftImagery"            
    #> [212] "AzureMaps.MicrosoftBaseDarkGrey"       
    #> [213] "AzureMaps.MicrosoftBaseRoad"           
    #> [214] "AzureMaps.MicrosoftBaseHybridRoad"     
    #> [215] "AzureMaps.MicrosoftTerraMain"          
    #> [216] "AzureMaps.MicrosoftWeatherInfraredMain"
    #> [217] "AzureMaps.MicrosoftWeatherRadarMain"   
    #> [218] "SwissFederalGeoportal"                 
    #> [219] "SwissFederalGeoportal.NationalMapColor"
    #> [220] "SwissFederalGeoportal.NationalMapGrey" 
    #> [221] "SwissFederalGeoportal.SWISSIMAGE"      
    #> [222] "TopPlusOpen"                           
    #> [223] "TopPlusOpen.Color"                     
    #> [224] "TopPlusOpen.Grey"

</div>

#### Tile providers’ details

``` r
str(providers_loaded()$providers_details)
```

<div style="height:150px; overflow-y: scroll;">

    #> List of 37
    #>  $ OpenStreetMap        :List of 3
    #>   ..$ url     : chr "https://tile.openstreetmap.org/{z}/{x}/{y}.png"
    #>   ..$ options :List of 2
    #>   .. ..$ maxZoom    : int 19
    #>   .. ..$ attribution: chr "&copy; <a href=\"https://www.openstreetmap.org/copyright\">OpenStreetMap</a> contributors"
    #>   ..$ variants:List of 7
    #>   .. ..$ Mapnik: Named list()
    #>   .. ..$ DE    :List of 2
    #>   .. .. ..$ url    : chr "https://tile.openstreetmap.de/{z}/{x}/{y}.png"
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ maxZoom: int 18
    #>   .. ..$ CH    :List of 2
    #>   .. .. ..$ url    : chr "https://tile.osm.ch/switzerland/{z}/{x}/{y}.png"
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ maxZoom: int 18
    #>   .. .. .. ..$ bounds : int [1:2, 1:2] 45 48 5 11
    #>   .. ..$ France:List of 2
    #>   .. .. ..$ url    : chr "https://{s}.tile.openstreetmap.fr/osmfr/{z}/{x}/{y}.png"
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ maxZoom    : int 20
    #>   .. .. .. ..$ attribution: chr "&copy; OpenStreetMap France | {attribution.OpenStreetMap}"
    #>   .. ..$ HOT   :List of 2
    #>   .. .. ..$ url    : chr "https://{s}.tile.openstreetmap.fr/hot/{z}/{x}/{y}.png"
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ attribution: chr "{attribution.OpenStreetMap}, Tiles style by <a href=\"https://www.hotosm.org/\" target=\"_blank\">Humanitarian "| __truncated__
    #>   .. ..$ BZH   :List of 2
    #>   .. .. ..$ url    : chr "https://tile.openstreetmap.bzh/br/{z}/{x}/{y}.png"
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ attribution: chr "{attribution.OpenStreetMap}, Tiles courtesy of <a href=\"http://www.openstreetmap.bzh/\" target=\"_blank\">Bret"| __truncated__
    #>   .. .. .. ..$ bounds     : num [1:2, 1:2] 46.2 50 -5.5 0.7
    #>   .. ..$ CAT   :List of 2
    #>   .. .. ..$ url    : chr "https://tile.openstreetmap.bzh/ca/{z}/{x}/{y}.png"
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ attribution: chr "{attribution.OpenStreetMap}, Tiles courtesy of <a href=\"https://www.openstreetmap.cat\" target=\"_blank\">Bret"| __truncated__
    #>  $ MapTilesAPI          :List of 3
    #>   ..$ url     : chr "https://maptiles.p.rapidapi.com/{variant}/{z}/{x}/{y}.png?rapidapi-key={apikey}"
    #>   ..$ options :List of 4
    #>   .. ..$ attribution: chr "&copy; <a href=\"http://www.maptilesapi.com/\">MapTiles API</a>, {attribution.OpenStreetMap}"
    #>   .. ..$ variant    : chr "en/map/v1"
    #>   .. ..$ apikey     : chr "<insert your api key here>"
    #>   .. ..$ maxZoom    : int 19
    #>   ..$ variants:List of 3
    #>   .. ..$ OSMEnglish :List of 1
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ variant: chr "en/map/v1"
    #>   .. ..$ OSMFrancais:List of 1
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ variant: chr "fr/map/v1"
    #>   .. ..$ OSMEspagnol:List of 1
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ variant: chr "es/map/v1"
    #>  $ OpenSeaMap           :List of 2
    #>   ..$ url    : chr "https://tiles.openseamap.org/seamark/{z}/{x}/{y}.png"
    #>   ..$ options:List of 1
    #>   .. ..$ attribution: chr "Map data: &copy; <a href=\"http://www.openseamap.org\">OpenSeaMap</a> contributors"
    #>  $ OPNVKarte            :List of 2
    #>   ..$ url    : chr "https://tileserver.memomaps.de/tilegen/{z}/{x}/{y}.png"
    #>   ..$ options:List of 2
    #>   .. ..$ maxZoom    : int 18
    #>   .. ..$ attribution: chr "Map <a href=\"https://memomaps.de/\">memomaps.de</a> <a href=\"http://creativecommons.org/licenses/by-sa/2.0/\""| __truncated__
    #>  $ OpenTopoMap          :List of 2
    #>   ..$ url    : chr "https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png"
    #>   ..$ options:List of 2
    #>   .. ..$ maxZoom    : int 17
    #>   .. ..$ attribution: chr "Map data: {attribution.OpenStreetMap}, <a href=\"http://viewfinderpanoramas.org\">SRTM</a> | Map style: &copy; "| __truncated__
    #>  $ OpenRailwayMap       :List of 2
    #>   ..$ url    : chr "https://{s}.tiles.openrailwaymap.org/standard/{z}/{x}/{y}.png"
    #>   ..$ options:List of 2
    #>   .. ..$ maxZoom    : int 19
    #>   .. ..$ attribution: chr "Map data: {attribution.OpenStreetMap} | Map style: &copy; <a href=\"https://www.OpenRailwayMap.org\">OpenRailwa"| __truncated__
    #>  $ OpenFireMap          :List of 2
    #>   ..$ url    : chr "http://openfiremap.org/hytiles/{z}/{x}/{y}.png"
    #>   ..$ options:List of 2
    #>   .. ..$ maxZoom    : int 19
    #>   .. ..$ attribution: chr "Map data: {attribution.OpenStreetMap} | Map style: &copy; <a href=\"http://www.openfiremap.org\">OpenFireMap</a"| __truncated__
    #>  $ SafeCast             :List of 2
    #>   ..$ url    : chr "https://s3.amazonaws.com/te512.safecast.org/{z}/{x}/{y}.png"
    #>   ..$ options:List of 2
    #>   .. ..$ maxZoom    : int 16
    #>   .. ..$ attribution: chr "Map data: {attribution.OpenStreetMap} | Map style: &copy; <a href=\"https://blog.safecast.org/about/\">SafeCast"| __truncated__
    #>  $ Stadia               :List of 3
    #>   ..$ url     : chr "https://tiles.stadiamaps.com/tiles/{variant}/{z}/{x}/{y}{r}.{ext}"
    #>   ..$ options :List of 5
    #>   .. ..$ minZoom    : int 0
    #>   .. ..$ maxZoom    : int 20
    #>   .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://openm"| __truncated__
    #>   .. ..$ variant    : chr "alidade_smooth"
    #>   .. ..$ ext        : chr "png"
    #>   ..$ variants:List of 15
    #>   .. ..$ AlidadeSmooth          : chr "alidade_smooth"
    #>   .. ..$ AlidadeSmoothDark      : chr "alidade_smooth_dark"
    #>   .. ..$ AlidadeSatellite       :List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ attribution: chr "&copy; CNES, Distribution Airbus DS, u00a9 Airbus DS, u00a9 PlanetObserver (Contains Copernicus Data) | &copy; "| __truncated__
    #>   .. .. .. ..$ variant    : chr "alidade_satellite"
    #>   .. .. .. ..$ ext        : chr "jpg"
    #>   .. ..$ OSMBright              : chr "osm_bright"
    #>   .. ..$ Outdoors               : chr "outdoors"
    #>   .. ..$ StamenToner            :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_toner"
    #>   .. ..$ StamenTonerBackground  :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_toner_background"
    #>   .. ..$ StamenTonerLines       :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_toner_lines"
    #>   .. ..$ StamenTonerLabels      :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_toner_labels"
    #>   .. ..$ StamenTonerLite        :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_toner_lite"
    #>   .. ..$ StamenWatercolor       :List of 2
    #>   .. .. ..$ url    : chr "https://tiles.stadiamaps.com/tiles/{variant}/{z}/{x}/{y}.{ext}"
    #>   .. .. ..$ options:List of 5
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_watercolor"
    #>   .. .. .. ..$ ext        : chr "jpg"
    #>   .. .. .. ..$ minZoom    : int 1
    #>   .. .. .. ..$ maxZoom    : int 16
    #>   .. ..$ StamenTerrain          :List of 1
    #>   .. .. ..$ options:List of 4
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_terrain"
    #>   .. .. .. ..$ minZoom    : int 0
    #>   .. .. .. ..$ maxZoom    : int 18
    #>   .. ..$ StamenTerrainBackground:List of 1
    #>   .. .. ..$ options:List of 4
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_terrain_background"
    #>   .. .. .. ..$ minZoom    : int 0
    #>   .. .. .. ..$ maxZoom    : int 18
    #>   .. ..$ StamenTerrainLabels    :List of 1
    #>   .. .. ..$ options:List of 4
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_terrain_labels"
    #>   .. .. .. ..$ minZoom    : int 0
    #>   .. .. .. ..$ maxZoom    : int 18
    #>   .. ..$ StamenTerrainLines     :List of 1
    #>   .. .. ..$ options:List of 4
    #>   .. .. .. ..$ attribution: chr "&copy; <a href=\"https://www.stadiamaps.com/\" target=\"_blank\">Stadia Maps</a> &copy; <a href=\"https://www.s"| __truncated__
    #>   .. .. .. ..$ variant    : chr "stamen_terrain_lines"
    #>   .. .. .. ..$ minZoom    : int 0
    #>   .. .. .. ..$ maxZoom    : int 18
    #>  $ Thunderforest        :List of 3
    #>   ..$ url     : chr "https://{s}.tile.thunderforest.com/{variant}/{z}/{x}/{y}{r}.png?apikey={apikey}"
    #>   ..$ options :List of 4
    #>   .. ..$ attribution: chr "&copy; <a href=\"http://www.thunderforest.com/\">Thunderforest</a>, {attribution.OpenStreetMap}"
    #>   .. ..$ variant    : chr "cycle"
    #>   .. ..$ apikey     : chr "<insert your api key here>"
    #>   .. ..$ maxZoom    : int 22
    #>   ..$ variants:List of 9
    #>   .. ..$ OpenCycleMap : chr "cycle"
    #>   .. ..$ Transport    :List of 1
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ variant: chr "transport"
    #>   .. ..$ TransportDark:List of 1
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ variant: chr "transport-dark"
    #>   .. ..$ SpinalMap    :List of 1
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ variant: chr "spinal-map"
    #>   .. ..$ Landscape    : chr "landscape"
    #>   .. ..$ Outdoors     : chr "outdoors"
    #>   .. ..$ Pioneer      : chr "pioneer"
    #>   .. ..$ MobileAtlas  : chr "mobile-atlas"
    #>   .. ..$ Neighbourhood: chr "neighbourhood"
    #>  $ BaseMapDE            :List of 3
    #>   ..$ url     : chr "https://sgx.geodatenzentrum.de/wmts_basemapde/tile/1.0.0/{variant}/default/GLOBAL_WEBMERCATOR/{z}/{y}/{x}.png"
    #>   ..$ options :List of 2
    #>   .. ..$ attribution: chr "Map data: &copy; <a href=\"http://www.govdata.de/dl-de/by-2-0\">dl-de/by-2-0</a>"
    #>   .. ..$ variant    : chr "de_basemapde_web_raster_farbe"
    #>   ..$ variants:List of 2
    #>   .. ..$ Color: chr "de_basemapde_web_raster_farbe"
    #>   .. ..$ Grey : chr "de_basemapde_web_raster_grau"
    #>  $ CyclOSM              :List of 2
    #>   ..$ url    : chr "https://{s}.tile-cyclosm.openstreetmap.fr/cyclosm/{z}/{x}/{y}.png"
    #>   ..$ options:List of 2
    #>   .. ..$ maxZoom    : int 20
    #>   .. ..$ attribution: chr "<a href=\"https://github.com/cyclosm/cyclosm-cartocss-style/releases\" title=\"CyclOSM - Open Bicycle render\">"| __truncated__
    #>  $ Jawg                 :List of 3
    #>   ..$ url     : chr "https://tile.jawg.io/{variant}/{z}/{x}/{y}{r}.png?access-token={accessToken}"
    #>   ..$ options :List of 5
    #>   .. ..$ attribution: chr "<a href=\"https://jawg.io\" title=\"Tiles Courtesy of Jawg Maps\" target=\"_blank\">&copy; <b>Jawg</b>Maps</a> "| __truncated__
    #>   .. ..$ minZoom    : int 0
    #>   .. ..$ maxZoom    : int 22
    #>   .. ..$ variant    : chr "jawg-streets"
    #>   .. ..$ accessToken: chr "<insert your access token here>"
    #>   ..$ variants:List of 7
    #>   .. ..$ Streets: chr "jawg-streets"
    #>   .. ..$ Terrain: chr "jawg-terrain"
    #>   .. ..$ Lagoon : chr "jawg-lagoon"
    #>   .. ..$ Sunny  : chr "jawg-sunny"
    #>   .. ..$ Dark   : chr "jawg-dark"
    #>   .. ..$ Light  : chr "jawg-light"
    #>   .. ..$ Matrix : chr "jawg-matrix"
    #>  $ MapBox               :List of 2
    #>   ..$ url    : chr "https://api.mapbox.com/styles/v1/{id}/tiles/{z}/{x}/{y}{r}?access_token={accessToken}"
    #>   ..$ options:List of 6
    #>   .. ..$ attribution: chr "&copy; <a href=\"https://www.mapbox.com/about/maps/\" target=\"_blank\">Mapbox</a> {attribution.OpenStreetMap} "| __truncated__
    #>   .. ..$ tileSize   : int 512
    #>   .. ..$ maxZoom    : int 18
    #>   .. ..$ zoomOffset : int -1
    #>   .. ..$ id         : chr "mapbox/streets-v11"
    #>   .. ..$ accessToken: chr "<insert your access token here>"
    #>  $ MapTiler             :List of 3
    #>   ..$ url     : chr "https://api.maptiler.com/maps/{variant}/{z}/{x}/{y}{r}.{ext}?key={key}"
    #>   ..$ options :List of 8
    #>   .. ..$ attribution: chr "<a href=\"https://www.maptiler.com/copyright/\" target=\"_blank\">&copy; MapTiler</a> <a href=\"https://www.ope"| __truncated__
    #>   .. ..$ variant    : chr "streets"
    #>   .. ..$ ext        : chr "png"
    #>   .. ..$ key        : chr "<insert your MapTiler Cloud API key here>"
    #>   .. ..$ tileSize   : int 512
    #>   .. ..$ zoomOffset : int -1
    #>   .. ..$ minZoom    : int 0
    #>   .. ..$ maxZoom    : int 21
    #>   ..$ variants:List of 20
    #>   .. ..$ Streets      : chr "streets-v2"
    #>   .. ..$ Basic        : chr "basic-v2"
    #>   .. ..$ Bright       : chr "bright-v2"
    #>   .. ..$ Pastel       : chr "pastel"
    #>   .. ..$ Positron     : chr "positron"
    #>   .. ..$ Hybrid       :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ variant: chr "hybrid"
    #>   .. .. .. ..$ ext    : chr "jpg"
    #>   .. ..$ Toner        : chr "toner-v2"
    #>   .. ..$ Topo         : chr "topo-v2"
    #>   .. ..$ Voyager      : chr "voyager-v2"
    #>   .. ..$ Ocean        : chr "ocean"
    #>   .. ..$ Backdrop     : chr "backdrop"
    #>   .. ..$ Dataviz      : chr "dataviz"
    #>   .. ..$ DatavizLight : chr "dataviz-light"
    #>   .. ..$ DatavizDark  : chr "dataviz-dark"
    #>   .. ..$ Aquarelle    :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ variant: chr "aquarelle"
    #>   .. .. .. ..$ ext    : chr "webp"
    #>   .. ..$ Landscape    : chr "landscape"
    #>   .. ..$ Openstreetmap:List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ variant: chr "openstreetmap"
    #>   .. .. .. ..$ ext    : chr "jpg"
    #>   .. ..$ Outdoor      : chr "outdoor-v2"
    #>   .. ..$ Satellite    :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ variant: chr "satellite"
    #>   .. .. .. ..$ ext    : chr "jpg"
    #>   .. ..$ Winter       : chr "winter-v2"
    #>  $ TomTom               :List of 3
    #>   ..$ url     : chr "https://{s}.api.tomtom.com/map/1/tile/{variant}/{style}/{z}/{x}/{y}.{ext}?key={apikey}"
    #>   ..$ options :List of 7
    #>   .. ..$ variant    : chr "basic"
    #>   .. ..$ maxZoom    : int 22
    #>   .. ..$ attribution: chr "<a href=\"https://tomtom.com\" target=\"_blank\">&copy;  1992 - 2026 TomTom.</a> "
    #>   .. ..$ subdomains : chr "abcd"
    #>   .. ..$ style      : chr "main"
    #>   .. ..$ ext        : chr "png"
    #>   .. ..$ apikey     : chr "<insert your API key here>"
    #>   ..$ variants:List of 3
    #>   .. ..$ Basic : chr "basic"
    #>   .. ..$ Hybrid: chr "hybrid"
    #>   .. ..$ Labels: chr "labels"
    #>  $ Esri                 :List of 3
    #>   ..$ url     : chr "https://server.arcgisonline.com/ArcGIS/rest/services/{variant}/MapServer/tile/{z}/{y}/{x}"
    #>   ..$ options :List of 2
    #>   .. ..$ variant    : chr "World_Street_Map"
    #>   .. ..$ attribution: chr "Tiles &copy; Esri"
    #>   ..$ variants:List of 9
    #>   .. ..$ WorldStreetMap   :List of 1
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ attribution: chr "{attribution.Esri} &mdash; Source: Esri, DeLorme, NAVTEQ, USGS, Intermap, iPC, NRCAN, Esri Japan, METI, Esri Ch"| __truncated__
    #>   .. ..$ WorldTopoMap     :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ variant    : chr "World_Topo_Map"
    #>   .. .. .. ..$ attribution: chr "{attribution.Esri} &mdash; Esri, DeLorme, NAVTEQ, TomTom, Intermap, iPC, USGS, FAO, NPS, NRCAN, GeoBase, Kadast"| __truncated__
    #>   .. ..$ WorldImagery     :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ variant    : chr "World_Imagery"
    #>   .. .. .. ..$ attribution: chr "{attribution.Esri} &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-"| __truncated__
    #>   .. ..$ WorldTerrain     :List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ variant    : chr "World_Terrain_Base"
    #>   .. .. .. ..$ maxZoom    : int 13
    #>   .. .. .. ..$ attribution: chr "{attribution.Esri} &mdash; Source: USGS, Esri, TANA, DeLorme, and NPS"
    #>   .. ..$ WorldShadedRelief:List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ variant    : chr "World_Shaded_Relief"
    #>   .. .. .. ..$ maxZoom    : int 13
    #>   .. .. .. ..$ attribution: chr "{attribution.Esri} &mdash; Source: Esri"
    #>   .. ..$ WorldPhysical    :List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ variant    : chr "World_Physical_Map"
    #>   .. .. .. ..$ maxZoom    : int 8
    #>   .. .. .. ..$ attribution: chr "{attribution.Esri} &mdash; Source: US National Park Service"
    #>   .. ..$ OceanBasemap     :List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ variant    : chr "Ocean/World_Ocean_Base"
    #>   .. .. .. ..$ maxZoom    : int 13
    #>   .. .. .. ..$ attribution: chr "{attribution.Esri} &mdash; Sources: GEBCO, NOAA, CHS, OSU, UNH, CSUMB, National Geographic, DeLorme, NAVTEQ, and Esri"
    #>   .. ..$ NatGeoWorldMap   :List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ variant    : chr "NatGeo_World_Map"
    #>   .. .. .. ..$ maxZoom    : int 16
    #>   .. .. .. ..$ attribution: chr "{attribution.Esri} &mdash; National Geographic, Esri, DeLorme, NAVTEQ, UNEP-WCMC, USGS, NASA, ESA, METI, NRCAN,"| __truncated__
    #>   .. ..$ WorldGrayCanvas  :List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ variant    : chr "Canvas/World_Light_Gray_Base"
    #>   .. .. .. ..$ maxZoom    : int 16
    #>   .. .. .. ..$ attribution: chr "{attribution.Esri} &mdash; Esri, DeLorme, NAVTEQ"
    #>  $ OpenWeatherMap       :List of 3
    #>   ..$ url     : chr "http://{s}.tile.openweathermap.org/map/{variant}/{z}/{x}/{y}.png?appid={apiKey}"
    #>   ..$ options :List of 4
    #>   .. ..$ maxZoom    : int 19
    #>   .. ..$ attribution: chr "Map data &copy; <a href=\"http://openweathermap.org\">OpenWeatherMap</a>"
    #>   .. ..$ apiKey     : chr "<insert your api key here>"
    #>   .. ..$ opacity    : num 0.5
    #>   ..$ variants:List of 11
    #>   .. ..$ Clouds              : chr "clouds"
    #>   .. ..$ CloudsClassic       : chr "clouds_cls"
    #>   .. ..$ Precipitation       : chr "precipitation"
    #>   .. ..$ PrecipitationClassic: chr "precipitation_cls"
    #>   .. ..$ Rain                : chr "rain"
    #>   .. ..$ RainClassic         : chr "rain_cls"
    #>   .. ..$ Pressure            : chr "pressure"
    #>   .. ..$ PressureContour     : chr "pressure_cntr"
    #>   .. ..$ Wind                : chr "wind"
    #>   .. ..$ Temperature         : chr "temp"
    #>   .. ..$ Snow                : chr "snow"
    #>  $ HERE                 :List of 3
    #>   ..$ url     : chr "https://maps.hereapi.com/v3/base/mc/{z}/{x}/{y}/{format}?style={variant}&size={size}&apiKey={app_id}&lg={language}"
    #>   ..$ options :List of 11
    #>   .. ..$ attribution: chr "Map &copy; 1987-2026 <a href=\"http://platform.here.com\">HERE</a>"
    #>   .. ..$ subdomains : chr "1234"
    #>   .. ..$ mapID      : chr "newest"
    #>   .. ..$ apiKey     : chr "<insert your apiKey here>"
    #>   .. ..$ base       : chr "base"
    #>   .. ..$ variant    : chr "explore.day"
    #>   .. ..$ maxZoom    : int 20
    #>   .. ..$ type       : chr "maptile"
    #>   .. ..$ language   : chr "eng"
    #>   .. ..$ format     : chr "png8"
    #>   .. ..$ size       : chr "256"
    #>   ..$ variants:List of 24
    #>   .. ..$ exploreDay           : chr "explore.day"
    #>   .. ..$ liteDay              : chr "lite.day"
    #>   .. ..$ logisticsDay         : chr "logistics.day"
    #>   .. ..$ topoDay              : chr "topo.day"
    #>   .. ..$ logisticsNight       : chr "logistics.night"
    #>   .. ..$ exploreNight         : chr "explore.night"
    #>   .. ..$ topoNight            : chr "topo.night"
    #>   .. ..$ liteNight            : chr "lite.night"
    #>   .. ..$ exploreSatelliteDay  : chr "explore.satellite.day"
    #>   .. ..$ liteSatelliteDay     : chr "lite.satellite.day"
    #>   .. ..$ logisticsSatelliteDay: chr "logistics.satellite.day"
    #>   .. ..$ basicMap             :List of 1
    #>   .. .. ..$ options:List of 1
    #>   .. .. .. ..$ type: chr "basetile"
    #>   .. ..$ mapLabels            :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ type  : chr "labeltile"
    #>   .. .. .. ..$ format: chr "png"
    #>   .. ..$ trafficFlow          :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ base: chr "traffic"
    #>   .. .. .. ..$ type: chr "flowtile"
    #>   .. ..$ carnavDayGrey        : chr "carnav.day.grey"
    #>   .. ..$ hybridDay            :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ base   : chr "aerial"
    #>   .. .. .. ..$ variant: chr "hybrid.day"
    #>   .. ..$ hybridDayMobile      :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ base   : chr "aerial"
    #>   .. .. .. ..$ variant: chr "hybrid.day.mobile"
    #>   .. ..$ hybridDayTransit     :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ base   : chr "aerial"
    #>   .. .. .. ..$ variant: chr "hybrid.day.transit"
    #>   .. ..$ hybridDayGrey        :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ base   : chr "aerial"
    #>   .. .. .. ..$ variant: chr "hybrid.grey.day"
    #>   .. ..$ pedestrianDay        : chr "pedestrian.day"
    #>   .. ..$ pedestrianNight      : chr "pedestrian.night"
    #>   .. ..$ satelliteDay         :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ base   : chr "aerial"
    #>   .. .. .. ..$ variant: chr "satellite.day"
    #>   .. ..$ terrainDay           :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ base   : chr "aerial"
    #>   .. .. .. ..$ variant: chr "terrain.day"
    #>   .. ..$ terrainDayMobile     :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ base   : chr "aerial"
    #>   .. .. .. ..$ variant: chr "terrain.day.mobile"
    #>  $ FreeMapSK            :List of 2
    #>   ..$ url    : chr "https://{s}.freemap.sk/T/{z}/{x}/{y}.jpeg"
    #>   ..$ options:List of 5
    #>   .. ..$ minZoom    : int 8
    #>   .. ..$ maxZoom    : int 16
    #>   .. ..$ subdomains : chr "abcd"
    #>   .. ..$ bounds     : num [1:2, 1:2] 47.2 49.8 16 22.6
    #>   .. ..$ attribution: chr "{attribution.OpenStreetMap}, visualization CC-By-SA 2.0 <a href=\"http://freemap.sk\">Freemap.sk</a>"
    #>  $ MtbMap               :List of 2
    #>   ..$ url    : chr "http://tile.mtbmap.cz/mtbmap_tiles/{z}/{x}/{y}.png"
    #>   ..$ options:List of 1
    #>   .. ..$ attribution: chr "{attribution.OpenStreetMap} &amp; USGS"
    #>  $ CartoDB              :List of 3
    #>   ..$ url     : chr "https://{s}.basemaps.cartocdn.com/{variant}/{z}/{x}/{y}{r}.png"
    #>   ..$ options :List of 4
    #>   .. ..$ attribution: chr "{attribution.OpenStreetMap} &copy; <a href=\"https://carto.com/attributions\">CARTO</a>"
    #>   .. ..$ subdomains : chr "abcd"
    #>   .. ..$ maxZoom    : int 20
    #>   .. ..$ variant    : chr "light_all"
    #>   ..$ variants:List of 10
    #>   .. ..$ Positron            : chr "light_all"
    #>   .. ..$ PositronNoLabels    : chr "light_nolabels"
    #>   .. ..$ PositronOnlyLabels  : chr "light_only_labels"
    #>   .. ..$ DarkMatter          : chr "dark_all"
    #>   .. ..$ DarkMatterNoLabels  : chr "dark_nolabels"
    #>   .. ..$ DarkMatterOnlyLabels: chr "dark_only_labels"
    #>   .. ..$ Voyager             : chr "rastertiles/voyager"
    #>   .. ..$ VoyagerNoLabels     : chr "rastertiles/voyager_nolabels"
    #>   .. ..$ VoyagerOnlyLabels   : chr "rastertiles/voyager_only_labels"
    #>   .. ..$ VoyagerLabelsUnder  : chr "rastertiles/voyager_labels_under"
    #>  $ HikeBike             :List of 3
    #>   ..$ url     : chr "https://tiles.wmflabs.org/{variant}/{z}/{x}/{y}.png"
    #>   ..$ options :List of 3
    #>   .. ..$ maxZoom    : int 19
    #>   .. ..$ attribution: chr "{attribution.OpenStreetMap}"
    #>   .. ..$ variant    : chr "hikebike"
    #>   ..$ variants:List of 2
    #>   .. ..$ HikeBike   : Named list()
    #>   .. ..$ HillShading:List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ maxZoom: int 15
    #>   .. .. .. ..$ variant: chr "hillshading"
    #>  $ BasemapAT            :List of 3
    #>   ..$ url     : chr "https://mapsneu.wien.gv.at/basemap/{variant}/{type}/google3857/{z}/{y}/{x}.{format}"
    #>   ..$ options :List of 6
    #>   .. ..$ maxZoom    : int 19
    #>   .. ..$ attribution: chr "Datenquelle: <a href=\"https://www.basemap.at\">basemap.at</a>"
    #>   .. ..$ type       : chr "normal"
    #>   .. ..$ format     : chr "png"
    #>   .. ..$ bounds     : num [1:2, 1:2] 46.36 49.04 8.78 17.19
    #>   .. ..$ variant    : chr "geolandbasemap"
    #>   ..$ variants:List of 7
    #>   .. ..$ basemap  :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ maxZoom: int 20
    #>   .. .. .. ..$ variant: chr "geolandbasemap"
    #>   .. ..$ grau     : chr "bmapgrau"
    #>   .. ..$ overlay  : chr "bmapoverlay"
    #>   .. ..$ terrain  :List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ variant: chr "bmapgelaende"
    #>   .. .. .. ..$ type   : chr "grau"
    #>   .. .. .. ..$ format : chr "jpeg"
    #>   .. ..$ surface  :List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ variant: chr "bmapoberflaeche"
    #>   .. .. .. ..$ type   : chr "grau"
    #>   .. .. .. ..$ format : chr "jpeg"
    #>   .. ..$ highdpi  :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ variant: chr "bmaphidpi"
    #>   .. .. .. ..$ format : chr "jpeg"
    #>   .. ..$ orthofoto:List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ maxZoom: int 20
    #>   .. .. .. ..$ variant: chr "bmaporthofoto30cm"
    #>   .. .. .. ..$ format : chr "jpeg"
    #>  $ nlmaps               :List of 3
    #>   ..$ url     : chr "https://service.pdok.nl/brt/achtergrondkaart/wmts/v2_0/{variant}/EPSG:3857/{z}/{x}/{y}.png"
    #>   ..$ options :List of 4
    #>   .. ..$ minZoom    : int 6
    #>   .. ..$ maxZoom    : int 19
    #>   .. ..$ bounds     : num [1:2, 1:2] 50.5 54 3.25 7.6
    #>   .. ..$ attribution: chr "Kaartgegevens &copy; <a href=\"https://www.kadaster.nl\">Kadaster</a>"
    #>   ..$ variants:List of 5
    #>   .. ..$ standaard: chr "standaard"
    #>   .. ..$ pastel   : chr "pastel"
    #>   .. ..$ grijs    : chr "grijs"
    #>   .. ..$ water    : chr "water"
    #>   .. ..$ luchtfoto:List of 1
    #>   .. .. ..$ url: chr "https://service.pdok.nl/hwh/luchtfotorgb/wmts/v1_0/Actueel_ortho25/EPSG:3857/{z}/{x}/{y}.jpeg"
    #>  $ NASAGIBS             :List of 3
    #>   ..$ url     : chr "https://map1.vis.earthdata.nasa.gov/wmts-webmerc/{variant}/default/{time}/{tilematrixset}{maxZoom}/{z}/{y}/{x}.{format}"
    #>   ..$ options :List of 7
    #>   .. ..$ attribution  : chr "Imagery provided by services from the Global Imagery Browse Services (GIBS), operated by the NASA/GSFC/Earth Sc"| __truncated__
    #>   .. ..$ bounds       : num [1:2, 1:2] -85.1 85.1 -180 180
    #>   .. ..$ minZoom      : int 1
    #>   .. ..$ maxZoom      : int 9
    #>   .. ..$ format       : chr "jpg"
    #>   .. ..$ time         : chr ""
    #>   .. ..$ tilematrixset: chr "GoogleMapsCompatible_Level"
    #>   ..$ variants:List of 7
    #>   .. ..$ ModisTerraTrueColorCR: chr "MODIS_Terra_CorrectedReflectance_TrueColor"
    #>   .. ..$ ModisTerraBands367CR : chr "MODIS_Terra_CorrectedReflectance_Bands367"
    #>   .. ..$ ViirsEarthAtNight2012:List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ variant: chr "VIIRS_CityLights_2012"
    #>   .. .. .. ..$ maxZoom: int 8
    #>   .. ..$ ModisTerraLSTDay     :List of 1
    #>   .. .. ..$ options:List of 4
    #>   .. .. .. ..$ variant: chr "MODIS_Terra_Land_Surface_Temp_Day"
    #>   .. .. .. ..$ format : chr "png"
    #>   .. .. .. ..$ maxZoom: int 7
    #>   .. .. .. ..$ opacity: num 0.75
    #>   .. ..$ ModisTerraSnowCover  :List of 1
    #>   .. .. ..$ options:List of 4
    #>   .. .. .. ..$ variant: chr "MODIS_Terra_NDSI_Snow_Cover"
    #>   .. .. .. ..$ format : chr "png"
    #>   .. .. .. ..$ maxZoom: int 8
    #>   .. .. .. ..$ opacity: num 0.75
    #>   .. ..$ ModisTerraAOD        :List of 1
    #>   .. .. ..$ options:List of 4
    #>   .. .. .. ..$ variant: chr "MODIS_Terra_Aerosol"
    #>   .. .. .. ..$ format : chr "png"
    #>   .. .. .. ..$ maxZoom: int 6
    #>   .. .. .. ..$ opacity: num 0.75
    #>   .. ..$ ModisTerraChlorophyll:List of 1
    #>   .. .. ..$ options:List of 4
    #>   .. .. .. ..$ variant: chr "MODIS_Terra_L2_Chlorophyll_A"
    #>   .. .. .. ..$ format : chr "png"
    #>   .. .. .. ..$ maxZoom: int 7
    #>   .. .. .. ..$ opacity: num 0.75
    #>  $ NLS                  :List of 3
    #>   ..$ url     : chr "https://api.maptiler.com/tiles/{variant}/{z}/{x}/{y}.jpg?key={apikey}"
    #>   ..$ options :List of 5
    #>   .. ..$ attribution: chr "<a href=\"http://maps.nls.uk/projects/subscription-api\">National Library of Scotland Historic Maps</a>"
    #>   .. ..$ bounds     : num [1:2, 1:2] 49.6 61.7 -12 3
    #>   .. ..$ minZoom    : int 1
    #>   .. ..$ maxZoom    : int 18
    #>   .. ..$ apikey     : chr "<insert your api key here>"
    #>   ..$ variants:List of 7
    #>   .. ..$ osgb63k1885   : chr "uk-osgb63k1885"
    #>   .. ..$ osgb1888      : chr "uk-osgb1888"
    #>   .. ..$ osgb10k1888   : chr "uk-osgb10k1888"
    #>   .. ..$ osgb1919      : chr "uk-osgb1919"
    #>   .. ..$ osgb25k1937   : chr "uk-osgb25k1937"
    #>   .. ..$ osgb63k1955   : chr "uk-osgb63k1955"
    #>   .. ..$ oslondon1k1893: chr "uk-oslondon1k1893"
    #>  $ JusticeMap           :List of 3
    #>   ..$ url     : chr "https://www.justicemap.org/tile/{size}/{variant}/{z}/{x}/{y}.png"
    #>   ..$ options :List of 3
    #>   .. ..$ attribution: chr "<a href=\"http://www.justicemap.org/terms.php\">Justice Map</a>"
    #>   .. ..$ size       : chr "county"
    #>   .. ..$ bounds     : int [1:2, 1:2] 14 72 -180 -56
    #>   ..$ variants:List of 9
    #>   .. ..$ income        : chr "income"
    #>   .. ..$ americanIndian: chr "indian"
    #>   .. ..$ asian         : chr "asian"
    #>   .. ..$ black         : chr "black"
    #>   .. ..$ hispanic      : chr "hispanic"
    #>   .. ..$ multi         : chr "multi"
    #>   .. ..$ nonWhite      : chr "nonwhite"
    #>   .. ..$ white         : chr "white"
    #>   .. ..$ plurality     : chr "plural"
    #>  $ GeoportailFrance     :List of 3
    #>   ..$ url     : chr "https://data.geopf.fr/wmts?REQUEST=GetTile&SERVICE=WMTS&VERSION=1.0.0&STYLE={style}&TILEMATRIXSET=PM&FORMAT={fo"| __truncated__
    #>   ..$ options :List of 7
    #>   .. ..$ attribution: chr "<a target=\"_blank\" href=\"https://www.geoportail.gouv.fr/\">Geoportail France</a>"
    #>   .. ..$ bounds     : int [1:2, 1:2] -75 81 -180 180
    #>   .. ..$ minZoom    : int 2
    #>   .. ..$ maxZoom    : int 18
    #>   .. ..$ format     : chr "image/png"
    #>   .. ..$ style      : chr "normal"
    #>   .. ..$ variant    : chr "GEOGRAPHICALGRIDSYSTEMS.PLANIGNV2"
    #>   ..$ variants:List of 3
    #>   .. ..$ plan   : chr "GEOGRAPHICALGRIDSYSTEMS.PLANIGNV2"
    #>   .. ..$ parcels:List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ variant: chr "CADASTRALPARCELS.PARCELLAIRE_EXPRESS"
    #>   .. .. .. ..$ style  : chr "PCI vecteur"
    #>   .. .. .. ..$ maxZoom: int 20
    #>   .. ..$ orthos :List of 1
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ maxZoom: int 19
    #>   .. .. .. ..$ format : chr "image/jpeg"
    #>   .. .. .. ..$ variant: chr "ORTHOIMAGERY.ORTHOPHOTOS"
    #>  $ OneMapSG             :List of 3
    #>   ..$ url     : chr "https://maps-{s}.onemap.sg/v3/{variant}/{z}/{x}/{y}.png"
    #>   ..$ options :List of 5
    #>   .. ..$ variant    : chr "Default"
    #>   .. ..$ minZoom    : int 11
    #>   .. ..$ maxZoom    : int 18
    #>   .. ..$ bounds     : num [1:2, 1:2] 1.56 1.16 104.11 103.5
    #>   .. ..$ attribution: chr "<img src=\"https://docs.onemap.sg/maps/images/oneMap64-01.png\" style=\"height:20px;width:20px;\"/> New OneMap "| __truncated__
    #>   ..$ variants:List of 5
    #>   .. ..$ Default : chr "Default"
    #>   .. ..$ Night   : chr "Night"
    #>   .. ..$ Original: chr "Original"
    #>   .. ..$ Grey    : chr "Grey"
    #>   .. ..$ LandLot : chr "LandLot"
    #>  $ USGS                 :List of 3
    #>   ..$ url     : chr "https://basemap.nationalmap.gov/arcgis/rest/services/USGSTopo/MapServer/tile/{z}/{y}/{x}"
    #>   ..$ options :List of 2
    #>   .. ..$ maxZoom    : int 20
    #>   .. ..$ attribution: chr "Tiles courtesy of the <a href=\"https://usgs.gov/\">U.S. Geological Survey</a>"
    #>   ..$ variants:List of 3
    #>   .. ..$ USTopo       : Named list()
    #>   .. ..$ USImagery    :List of 1
    #>   .. .. ..$ url: chr "https://basemap.nationalmap.gov/arcgis/rest/services/USGSImageryOnly/MapServer/tile/{z}/{y}/{x}"
    #>   .. ..$ USImageryTopo:List of 1
    #>   .. .. ..$ url: chr "https://basemap.nationalmap.gov/arcgis/rest/services/USGSImageryTopo/MapServer/tile/{z}/{y}/{x}"
    #>  $ WaymarkedTrails      :List of 3
    #>   ..$ url     : chr "https://tile.waymarkedtrails.org/{variant}/{z}/{x}/{y}.png"
    #>   ..$ options :List of 2
    #>   .. ..$ maxZoom    : int 18
    #>   .. ..$ attribution: chr "Map data: {attribution.OpenStreetMap} | Map style: &copy; <a href=\"https://waymarkedtrails.org\">waymarkedtrai"| __truncated__
    #>   ..$ variants:List of 6
    #>   .. ..$ hiking : chr "hiking"
    #>   .. ..$ cycling: chr "cycling"
    #>   .. ..$ mtb    : chr "mtb"
    #>   .. ..$ slopes : chr "slopes"
    #>   .. ..$ riding : chr "riding"
    #>   .. ..$ skating: chr "skating"
    #>  $ OpenAIP              :List of 2
    #>   ..$ url    : chr "https://{s}.tile.maps.openaip.net/geowebcache/service/tms/1.0.0/openaip_basemap@EPSG%3A900913@png/{z}/{x}/{y}.{ext}"
    #>   ..$ options:List of 7
    #>   .. ..$ attribution : chr "<a href=\"https://www.openaip.net/\">openAIP Data</a> (<a href=\"https://creativecommons.org/licenses/by-sa/3.0"| __truncated__
    #>   .. ..$ ext         : chr "png"
    #>   .. ..$ minZoom     : int 4
    #>   .. ..$ maxZoom     : int 14
    #>   .. ..$ tms         : logi TRUE
    #>   .. ..$ detectRetina: logi TRUE
    #>   .. ..$ subdomains  : chr "12"
    #>  $ OpenSnowMap          :List of 3
    #>   ..$ url     : chr "https://tiles.opensnowmap.org/{variant}/{z}/{x}/{y}.png"
    #>   ..$ options :List of 3
    #>   .. ..$ minZoom    : int 9
    #>   .. ..$ maxZoom    : int 18
    #>   .. ..$ attribution: chr "Map data: {attribution.OpenStreetMap} & ODbL, &copy; <a href=\"https://www.opensnowmap.org/iframes/data.html\">"| __truncated__
    #>   ..$ variants:List of 1
    #>   .. ..$ pistes: chr "pistes"
    #>  $ AzureMaps            :List of 3
    #>   ..$ url     : chr "https://atlas.microsoft.com/map/tile?api-version={apiVersion}&tilesetId={variant}&x={x}&y={y}&zoom={z}&language"| __truncated__
    #>   ..$ options :List of 5
    #>   .. ..$ attribution    : chr "See https://docs.microsoft.com/en-us/rest/api/maps/render-v2/get-map-tile for details."
    #>   .. ..$ apiVersion     : chr "2.0"
    #>   .. ..$ variant        : chr "microsoft.imagery"
    #>   .. ..$ subscriptionKey: chr "<insert your subscription key here>"
    #>   .. ..$ language       : chr "en-US"
    #>   ..$ variants:List of 7
    #>   .. ..$ MicrosoftImagery            : chr "microsoft.imagery"
    #>   .. ..$ MicrosoftBaseDarkGrey       : chr "microsoft.base.darkgrey"
    #>   .. ..$ MicrosoftBaseRoad           : chr "microsoft.base.road"
    #>   .. ..$ MicrosoftBaseHybridRoad     : chr "microsoft.base.hybrid.road"
    #>   .. ..$ MicrosoftTerraMain          : chr "microsoft.terra.main"
    #>   .. ..$ MicrosoftWeatherInfraredMain:List of 2
    #>   .. .. ..$ url    : chr "https://atlas.microsoft.com/map/tile?api-version={apiVersion}&tilesetId={variant}&x={x}&y={y}&zoom={z}&timeStam"| __truncated__
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ timeStamp  : chr "2021-05-08T09:03:00Z"
    #>   .. .. .. ..$ attribution: chr "See https://docs.microsoft.com/en-us/rest/api/maps/render-v2/get-map-tile#uri-parameters for details."
    #>   .. .. .. ..$ variant    : chr "microsoft.weather.infrared.main"
    #>   .. ..$ MicrosoftWeatherRadarMain   :List of 2
    #>   .. .. ..$ url    : chr "https://atlas.microsoft.com/map/tile?api-version={apiVersion}&tilesetId={variant}&x={x}&y={y}&zoom={z}&timeStam"| __truncated__
    #>   .. .. ..$ options:List of 3
    #>   .. .. .. ..$ timeStamp  : chr "2021-05-08T09:03:00Z"
    #>   .. .. .. ..$ attribution: chr "See https://docs.microsoft.com/en-us/rest/api/maps/render-v2/get-map-tile#uri-parameters for details."
    #>   .. .. .. ..$ variant    : chr "microsoft.weather.radar.main"
    #>  $ SwissFederalGeoportal:List of 3
    #>   ..$ url     : chr "https://wmts.geo.admin.ch/1.0.0/{variant}/default/current/3857/{z}/{x}/{y}.jpeg"
    #>   ..$ options :List of 4
    #>   .. ..$ attribution: chr "&copy; <a href=\"https://www.swisstopo.admin.ch/\">swisstopo</a>"
    #>   .. ..$ minZoom    : int 2
    #>   .. ..$ maxZoom    : int 18
    #>   .. ..$ bounds     : num [1:2, 1:2] 45.4 48.23 5.14 11.48
    #>   ..$ variants:List of 3
    #>   .. ..$ NationalMapColor: chr "ch.swisstopo.pixelkarte-farbe"
    #>   .. ..$ NationalMapGrey : chr "ch.swisstopo.pixelkarte-grau"
    #>   .. ..$ SWISSIMAGE      :List of 1
    #>   .. .. ..$ options:List of 2
    #>   .. .. .. ..$ variant: chr "ch.swisstopo.swissimage"
    #>   .. .. .. ..$ maxZoom: int 19
    #>  $ TopPlusOpen          :List of 3
    #>   ..$ url     : chr "http://sgx.geodatenzentrum.de/wmts_topplus_open/tile/1.0.0/{variant}/default/WEBMERCATOR/{z}/{y}/{x}.png"
    #>   ..$ options :List of 3
    #>   .. ..$ maxZoom    : int 18
    #>   .. ..$ attribution: chr "Map data: &copy; <a href=\"http://www.govdata.de/dl-de/by-2-0\">dl-de/by-2-0</a>"
    #>   .. ..$ variant    : chr "web"
    #>   ..$ variants:List of 2
    #>   .. ..$ Color: chr "web"
    #>   .. ..$ Grey : chr "web_grau"

</div>

## Developer

### Updating to a new leaflet-providers version

To update this package to a new upstream release of
[leaflet-providers.js](https://github.com/leaflet-extras/leaflet-providers),
run the following [Claude
Code](https://platform.claude.com/docs/en/docs/claude-code) skill:

``` bash
/update-leaflet-providers
```

This walks through the full update workflow: fetching the new JS,
regenerating package data, bumping the version, running checks, and
creating a PR.
