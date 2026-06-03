
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rebird: wrapper to the eBird API

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![Build
Status](https://api.travis-ci.org/ropensci/rebird.png)](https://travis-ci.org/ropensci/rebird)
[![Build
status](https://ci.appveyor.com/api/projects/status/s3dobn991c20t2kg?svg=true)](https://ci.appveyor.com/project/sckott/rebird)
[![cran
checks](https://cranchecks.info/badges/worst/rebird)](https://cranchecks.info/pkgs/rebird)
[![Coverage
Status](https://coveralls.io/repos/ropensci/rebird/badge.svg)](https://coveralls.io/github/ropensci/rebird)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/rebird)](https://github.com/r-hub/cranlogs.app)
[![cran
version](https://www.r-pkg.org/badges/version/rebird)](https://cran.r-project.org/package=rebird/)

`rebird` is a package to interface with the eBird webservices.

eBird is a real-time, online bird checklist program. For more
information, visit their website: <https://ebird.org/home>

The API for the eBird webservices can be accessed here:
<https://documenter.getpostman.com/view/664302/S1ENwy59?version=latest>

## Install

You can install the stable version from CRAN

``` r
install.packages("rebird")
```

Or the development version from Github

``` r
install.packages("devtools")
devtools::install_github("ropensci/rebird")
```

# Direct use of `rebird`

Load the package:

``` r
library("rebird")
```

The [eBird API
server](https://documenter.getpostman.com/view/664302/S1ENwy59?version=latest)
has been updated and thus there are a couple major changes in the way
`rebird` works. API requests to eBird now require users to provide an
API key, which is linked to your eBird user account. You can pass it to
the ‘key’ argument in `rebird` functions, but we highly recommend
storing it as an environment variable called EBIRD\_KEY in your
.Renviron file. If you don’t have a key, you can obtain one from
<https://ebird.org/api/keygen>.

You can keep your .Renviron file in your global R home directory
(`R.home()`), your user’s home directory (`Sys.getenv("HOME")`), or your
current working directory (`getwd()`). Remember that .Renviron is loaded
once when you start R, so if you add your API key to the file you will
have to restart your R session. See `?Startup` for more information on
R’s startup files.

Furthermore, functions now use species codes, rather than scientific
names, for species-specific requests. We’ve made the switch easy by
providing the `species_code` function, which converts a scientific name
to its species code:

``` r
species_code('sula variegata')
#> Peruvian Booby (Sula variegata): perboo1
#> [1] "perboo1"
```

The `species_code` function can be called within other `rebird`
functions, or the species code can be specified directly.

## eBird Taxonomy

The eBird taxonomy is internally stored in `rebird` and can be called
using

``` r
rebird:::tax
#> # A tibble: 16,513 x 14
#>    sciName comName speciesCode category taxonOrder bandingCodes comNameCodes
#>    <chr>   <chr>   <chr>       <chr>         <dbl> <chr>        <chr>       
#>  1 Struth… Common… ostric2     species           1 <NA>         COOS        
#>  2 Struth… Somali… ostric3     species           6 <NA>         SOOS        
#>  3 Struth… Common… y00934      slash             7 <NA>         SOOS,COOS   
#>  4 Rhea a… Greate… grerhe1     species           8 <NA>         GRRH        
#>  5 Rhea p… Lesser… lesrhe2     species          14 <NA>         LERH        
#>  6 Rhea p… Lesser… lesrhe4     issf             15 <NA>         LERH        
#>  7 Rhea p… Lesser… lesrhe3     issf             18 <NA>         LERH        
#>  8 Nothoc… Tawny-… tabtin1     species          19 <NA>         TBTI        
#>  9 Nothoc… Highla… higtin1     species          20 HITI         <NA>        
#> 10 Nothoc… Highla… higtin2     issf             21 <NA>         HITI        
#> # … with 16,503 more rows, and 7 more variables: sciNameCodes <chr>,
#> #   order <chr>, familyComName <chr>, familySciName <chr>, reportAs <chr>,
#> #   extinct <lgl>, extinctYear <int>
```

While the internal taxonomy is kept up to date with each package
release, it could be outdated if a new taxonomy is made available before
the package is updated. You can obtain the latest eBird taxonomy by

``` r
new_tax <- ebirdtaxonomy()
```

## Sightings at location determined by latitude/longitude

Search for bird occurrences by latitude and longitude point

``` r
ebirdgeo(species = species_code('spinus tristis'), lat = 42, lng = -76)
#> American Goldfinch (Spinus tristis): amegfi
#> # A tibble: 17 x 13
#>    speciesCode comName sciName locId locName obsDt howMany   lat   lng obsValid
#>    <chr>       <chr>   <chr>   <chr> <chr>   <chr>   <int> <dbl> <dbl> <lgl>   
#>  1 amegfi      Americ… Spinus… L133… "545 R… 2021…       8  42.0 -76.1 TRUE    
#>  2 amegfi      Americ… Spinus… L197… "esthe… 2021…      12  42.1 -75.9 TRUE    
#>  3 amegfi      Americ… Spinus… L275… "Home " 2021…       8  42.1 -76.0 TRUE    
#>  4 amegfi      Americ… Spinus… L109… "Hillc… 2021…       1  42.2 -75.9 TRUE    
#>  5 amegfi      Americ… Spinus… L186… "Otsin… 2021…       1  42.1 -75.9 TRUE    
#>  6 amegfi      Americ… Spinus… L895… "Nowla… 2021…       2  42.1 -75.9 TRUE    
#>  7 amegfi      Americ… Spinus… L207… "Workw… 2021…       4  42.1 -75.9 TRUE    
#>  8 amegfi      Americ… Spinus… L133… "4457 … 2021…       2  42.0 -75.9 TRUE    
#>  9 amegfi      Americ… Spinus… L870… "325 D… 2021…       1  42.2 -76.0 TRUE    
#> 10 amegfi      Americ… Spinus… L121… "1312 … 2021…       3  42.1 -76.0 TRUE    
#> 11 amegfi      Americ… Spinus… L133… "216 W… 2021…       1  42.1 -76.0 TRUE    
#> 12 amegfi      Americ… Spinus… L524… "Victo… 2021…       4  42.1 -76.0 TRUE    
#> 13 amegfi      Americ… Spinus… L505… "Bolan… 2021…       1  42.2 -75.9 TRUE    
#> 14 amegfi      Americ… Spinus… L850… "Sandy… 2021…      20  42.1 -75.9 TRUE    
#> 15 amegfi      Americ… Spinus… L351… "Anson… 2021…      16  42.1 -76.1 TRUE    
#> 16 amegfi      Americ… Spinus… L270… "Gripp… 2021…       1  42.1 -76.1 TRUE    
#> 17 amegfi      Americ… Spinus… L564… "Kinne… 2021…       1  42.1 -76.2 TRUE    
#> # … with 3 more variables: obsReviewed <lgl>, locationPrivate <lgl>,
#> #   subId <chr>
```

## Recent observations at a region

Search for bird occurrences by region and species name

``` r
ebirdregion(loc = 'US', species = 'btbwar')
#> # A tibble: 81 x 13
#>    speciesCode comName sciName locId locName obsDt howMany   lat   lng obsValid
#>    <chr>       <chr>   <chr>   <chr> <chr>   <chr>   <int> <dbl> <dbl> <lgl>   
#>  1 btbwar      Black-… Setoph… L577… Merrit… 2021…       1  28.6 -80.7 TRUE    
#>  2 btbwar      Black-… Setoph… L863… 104 7t… 2021…       1  32.0 -80.8 TRUE    
#>  3 btbwar      Black-… Setoph… L407… Thomps… 2021…       1  38.9 -77.1 TRUE    
#>  4 btbwar      Black-… Setoph… L193… Rye     2021…       1  43.0 -70.8 TRUE    
#>  5 btbwar      Black-… Setoph… L195… 1 My H… 2021…       1  27.0 -80.1 TRUE    
#>  6 btbwar      Black-… Setoph… L104… Feathe… 2021…       1  25.6 -80.3 TRUE    
#>  7 btbwar      Black-… Setoph… L324… Wither… 2021…       1  31.0 -82.9 TRUE    
#>  8 btbwar      Black-… Setoph… L128… Zoo Mi… 2021…       1  25.6 -80.4 TRUE    
#>  9 btbwar      Black-… Setoph… L992… Kendal… 2021…       1  25.7 -80.4 TRUE    
#> 10 btbwar      Black-… Setoph… L133… 603 S … 2021…       1  26.2 -98.2 TRUE    
#> # … with 71 more rows, and 3 more variables: obsReviewed <lgl>,
#> #   locationPrivate <lgl>, subId <chr>
```

## Recent observations at hotspots

Search for bird occurrences by a given hotspot

``` r
ebirdregion(loc = 'L99381')
#> # A tibble: 38 x 13
#>    speciesCode comName sciName locId locName obsDt howMany   lat   lng obsValid
#>    <chr>       <chr>   <chr>   <chr> <chr>   <chr>   <int> <dbl> <dbl> <lgl>   
#>  1 cangoo      Canada… Branta… L993… Stewar… 2021…     300  42.5 -76.5 TRUE    
#>  2 mallar3     Mallard Anas p… L993… Stewar… 2021…      20  42.5 -76.5 TRUE    
#>  3 commer      Common… Mergus… L993… Stewar… 2021…       2  42.5 -76.5 TRUE    
#>  4 ribgul      Ring-b… Larus … L993… Stewar… 2021…      NA  42.5 -76.5 TRUE    
#>  5 hergul      Herrin… Larus … L993… Stewar… 2021…      NA  42.5 -76.5 TRUE    
#>  6 gbbgul      Great … Larus … L993… Stewar… 2021…       3  42.5 -76.5 TRUE    
#>  7 baleag      Bald E… Haliae… L993… Stewar… 2021…       1  42.5 -76.5 TRUE    
#>  8 eursta      Europe… Sturnu… L993… Stewar… 2021…      40  42.5 -76.5 TRUE    
#>  9 doccor      Double… Phalac… L993… Stewar… 2021…       5  42.5 -76.5 TRUE    
#> 10 ambduc      Americ… Anas r… L993… Stewar… 2021…       1  42.5 -76.5 TRUE    
#> # … with 28 more rows, and 3 more variables: obsReviewed <lgl>,
#> #   locationPrivate <lgl>, subId <chr>
```

## Nearest observations of a species

Search for a species’ occurrences near a given latitude and longitude

``` r
nearestobs(species_code('branta canadensis'), 42, -76)
#> Canada Goose (Branta canadensis): cangoo
#> # A tibble: 25 x 13
#>    speciesCode comName sciName locId locName obsDt howMany   lat   lng obsValid
#>    <chr>       <chr>   <chr>   <chr> <chr>   <chr>   <int> <dbl> <dbl> <lgl>   
#>  1 cangoo      Canada… Branta… L109… Hillcr… 2021…      34  42.2 -75.9 TRUE    
#>  2 cangoo      Canada… Branta… L186… Otsini… 2021…     117  42.1 -75.9 TRUE    
#>  3 cangoo      Canada… Branta… L186… Cheri … 2021…       4  42.1 -75.9 TRUE    
#>  4 cangoo      Canada… Branta… L809… Port D… 2021…      74  42.1 -75.9 TRUE    
#>  5 cangoo      Canada… Branta… L527… R Tee … 2021…     100  42.2 -75.9 TRUE    
#>  6 cangoo      Canada… Branta… L133… I-81 N… 2021…      45  42.1 -75.9 TRUE    
#>  7 cangoo      Canada… Branta… L245… Water … 2021…       3  42.1 -75.9 TRUE    
#>  8 cangoo      Canada… Branta… L116… Homest… 2021…     230  42.1 -76.0 TRUE    
#>  9 cangoo      Canada… Branta… L106… IBM CC… 2021…       1  42.1 -76.0 TRUE    
#> 10 cangoo      Canada… Branta… L273… Schnur… 2021…       2  42.1 -75.8 TRUE    
#> # … with 15 more rows, and 3 more variables: obsReviewed <lgl>,
#> #   locationPrivate <lgl>, subId <chr>
```

<!--
## Frequency of observations at hotspots or regions

Obtain historical frequencies of bird occurrences by hotspot or region


```r
# ebirdfreq(loctype = 'hotspots', loc = 'L196159')
```
-->

## Recent notable sightings

Search for notable sightings at a given latitude and longitude

``` r
ebirdnotable(lat = 42, lng = -70)
#> # A tibble: 3,578 x 13
#>    speciesCode comName sciName locId locName obsDt howMany   lat   lng obsValid
#>    <chr>       <chr>   <chr>   <chr> <chr>   <chr>   <int> <dbl> <dbl> <lgl>   
#>  1 foxsp1      Fox Sp… Passer… L382… Yard    2021…       1  42.2 -71.3 FALSE   
#>  2 foxsp1      Fox Sp… Passer… L276… Standi… 2021…       1  42.3 -71.3 FALSE   
#>  3 reshaw      Red-sh… Buteo … L575… The 20… 2021…       1  44.2 -69.4 FALSE   
#>  4 gadwal      Gadwall Mareca… L143… Holyok… 2021…       2  42.2 -72.6 FALSE   
#>  5 lbbgul      Lesser… Larus … L106… Goulds… 2021…       1  42.2 -71.4 FALSE   
#>  6 pinwar      Pine W… Setoph… L919… Northb… 2021…       1  42.1 -71.7 FALSE   
#>  7 bnhcow      Brown-… Moloth… L825… Westwo… 2021…       1  42.2 -71.2 FALSE   
#>  8 comred      Common… Acanth… L358… Fort H… 2021…       2  41.8 -70.0 FALSE   
#>  9 redcro10    Red Cr… Loxia … L480… West B… 2021…       2  41.7 -70.4 FALSE   
#> 10 foxsp1      Fox Sp… Passer… L276… Standi… 2021…       1  42.3 -71.3 FALSE   
#> # … with 3,568 more rows, and 3 more variables: obsReviewed <lgl>,
#> #   locationPrivate <lgl>, subId <chr>
```

or a region

``` r
ebirdnotable(locID = 'US-NY-109')
#> # A tibble: 81 x 13
#>    speciesCode comName sciName locId locName obsDt howMany   lat   lng obsValid
#>    <chr>       <chr>   <chr>   <chr> <chr>   <chr>   <int> <dbl> <dbl> <lgl>   
#>  1 blkvul      Black … Coragy… L212… Steven… 2021…       1  42.4 -76.4 FALSE   
#>  2 redcro      Red Cr… Loxia … L550… Cornel… 2021…       1  42.5 -76.5 FALSE   
#>  3 yerwar      Yellow… Setoph… L351… Fuerte… 2021…       1  42.5 -76.5 FALSE   
#>  4 redcro      Red Cr… Loxia … L123… Boyer … 2021…       5  42.3 -76.3 FALSE   
#>  5 whwcro      White-… Loxia … L550… Cornel… 2021…       2  42.5 -76.5 FALSE   
#>  6 x00684      Canvas… Aythya… L140… East S… 2021…       1  42.5 -76.5 FALSE   
#>  7 x00684      Canvas… Aythya… L140… East S… 2021…       1  42.5 -76.5 FALSE   
#>  8 blksco2     Black … Melani… L353… Salt P… 2021…       1  42.5 -76.5 FALSE   
#>  9 evegro      Evenin… Coccot… L133… 571 So… 2021…      17  42.3 -76.4 FALSE   
#> 10 hoared2     Hoary … Acanth… L686… George… 2021…       1  42.5 -76.3 FALSE   
#> # … with 71 more rows, and 3 more variables: obsReviewed <lgl>,
#> #   locationPrivate <lgl>, subId <chr>
```

## Historic Observations

Obtain a list of species reported on a specific date in a given region

``` r
ebirdhistorical(loc = 'US-VA-003', date = '2019-02-14',max = 10)
#> # A tibble: 10 x 13
#>    speciesCode comName sciName locId locName obsDt howMany   lat   lng obsValid
#>    <chr>       <chr>   <chr>   <chr> <chr>   <chr>   <int> <dbl> <dbl> <lgl>   
#>  1 cangoo      Canada… Branta… L139… Lickin… 2019…      30  38.1 -78.7 TRUE    
#>  2 mallar3     Mallard Anas p… L139… Lickin… 2019…       5  38.1 -78.7 TRUE    
#>  3 gnwtea      Green-… Anas c… L139… Lickin… 2019…       8  38.1 -78.7 TRUE    
#>  4 killde      Killde… Charad… L139… Lickin… 2019…       1  38.1 -78.7 TRUE    
#>  5 baleag      Bald E… Haliae… L139… Lickin… 2019…       1  38.1 -78.7 TRUE    
#>  6 belkin1     Belted… Megace… L139… Lickin… 2019…       1  38.1 -78.7 TRUE    
#>  7 carwre      Caroli… Thryot… L139… Lickin… 2019…       1  38.1 -78.7 TRUE    
#>  8 whtspa      White-… Zonotr… L139… Lickin… 2019…       2  38.1 -78.7 TRUE    
#>  9 norcar      Northe… Cardin… L139… Lickin… 2019…       1  38.1 -78.7 TRUE    
#> 10 canvas      Canvas… Aythya… L331… Montic… 2019…      19  38.0 -78.5 TRUE    
#> # … with 3 more variables: obsReviewed <lgl>, locationPrivate <lgl>,
#> #   subId <chr>
```

or a hotspot

``` r
ebirdhistorical(loc = 'L196159', date = '2019-02-14', fieldSet = 'full')
#> # A tibble: 14 x 27
#>    speciesCode comName sciName locId locName obsDt howMany   lat   lng obsValid
#>    <chr>       <chr>   <chr>   <chr> <chr>   <chr>   <int> <dbl> <dbl> <lgl>   
#>  1 annhum      Anna's… Calypt… L196… Vancou… 2019…       4  49.3 -123. TRUE    
#>  2 ribgul      Ring-b… Larus … L196… Vancou… 2019…       4  49.3 -123. TRUE    
#>  3 glwgul      Glauco… Larus … L196… Vancou… 2019…      29  49.3 -123. TRUE    
#>  4 norcro      Northw… Corvus… L196… Vancou… 2019…     100  49.3 -123. TRUE    
#>  5 bkcchi      Black-… Poecil… L196… Vancou… 2019…      16  49.3 -123. TRUE    
#>  6 bushti      Bushtit Psaltr… L196… Vancou… 2019…      20  49.3 -123. TRUE    
#>  7 pacwre1     Pacifi… Troglo… L196… Vancou… 2019…       1  49.3 -123. TRUE    
#>  8 houfin      House … Haemor… L196… Vancou… 2019…       2  49.3 -123. TRUE    
#>  9 purfin      Purple… Haemor… L196… Vancou… 2019…       3  49.3 -123. TRUE    
#> 10 amegfi      Americ… Spinus… L196… Vancou… 2019…      15  49.3 -123. TRUE    
#> 11 daejun      Dark-e… Junco … L196… Vancou… 2019…      37  49.3 -123. TRUE    
#> 12 sonspa      Song S… Melosp… L196… Vancou… 2019…      12  49.3 -123. TRUE    
#> 13 spotow      Spotte… Pipilo… L196… Vancou… 2019…       1  49.3 -123. TRUE    
#> 14 rewbla      Red-wi… Agelai… L196… Vancou… 2019…       6  49.3 -123. TRUE    
#> # … with 17 more variables: obsReviewed <lgl>, locationPrivate <lgl>,
#> #   subId <chr>, subnational2Code <chr>, subnational2Name <chr>,
#> #   subnational1Code <chr>, subnational1Name <chr>, countryCode <chr>,
#> #   countryName <chr>, userDisplayName <chr>, obsId <chr>, checklistId <chr>,
#> #   presenceNoted <lgl>, hasComments <lgl>, firstName <chr>, lastName <chr>,
#> #   hasRichMedia <lgl>
```

## Information on a given region or hotspot

Obtain detailed information on any valid eBird region

``` r
ebirdregioninfo("CA-BC-GV")
#> # A tibble: 1 x 5
#>   region                                     minX  maxX  minY  maxY
#>   <chr>                                     <dbl> <dbl> <dbl> <dbl>
#> 1 Metro Vancouver, British Columbia, Canada -123. -122.  49.0  49.6
```

or hotspot

``` r
ebirdregioninfo("L196159")
#> # A tibble: 1 x 16
#>   locId name  latitude longitude countryCode countryName subnational1Name
#>   <chr> <chr>    <dbl>     <dbl> <chr>       <chr>       <chr>           
#> 1 L196… Vanc…     49.3     -123. CA          Canada      British Columbia
#> # … with 9 more variables: subnational1Code <chr>, subnational2Code <chr>,
#> #   subnational2Name <chr>, isHotspot <lgl>, locID <chr>, locName <chr>,
#> #   lat <dbl>, lng <dbl>, hierarchicalName <chr>
```

Obtain a list of eBird species codes for all species recorded in a
region

``` r
ebirdregionspecies("GB-ENG-LND")
#> # A tibble: 304 x 1
#>    speciesCode
#>    <chr>      
#>  1 bahgoo     
#>  2 snogoo     
#>  3 gragoo     
#>  4 gwfgoo     
#>  5 tunbeg1    
#>  6 pifgoo     
#>  7 brant      
#>  8 bargoo     
#>  9 cangoo     
#> 10 rebgoo1    
#> # … with 294 more rows
```

or a hotspot

``` r
ebirdregionspecies("L5803024")
#> # A tibble: 156 x 1
#>    speciesCode
#>    <chr>      
#>  1 gragoo     
#>  2 gwfgoo     
#>  3 bargoo     
#>  4 cangoo     
#>  5 mutswa     
#>  6 egygoo     
#>  7 comshe     
#>  8 manduc     
#>  9 gargan     
#> 10 norsho     
#> # … with 146 more rows
```

Obtain a list of all subregions within an eBird region

``` r
ebirdsubregionlist("subnational1","US")
#> # A tibble: 51 x 2
#>    code  name                
#>    <chr> <chr>               
#>  1 US-AL Alabama             
#>  2 US-AK Alaska              
#>  3 US-AZ Arizona             
#>  4 US-AR Arkansas            
#>  5 US-CA California          
#>  6 US-CO Colorado            
#>  7 US-CT Connecticut         
#>  8 US-DE Delaware            
#>  9 US-DC District of Columbia
#> 10 US-FL Florida             
#> # … with 41 more rows
```

## Checklist Feed

Obtain a list of checklists submitted on a given date at a region or
hotspot

``` r
ebirdchecklistfeed(loc = "L207391", date = "2020-03-24", max = 5)
#> # A tibble: 5 x 8
#>   locId  subId  userDisplayName  numSpecies obsDt obsTime subID loc             
#>   <chr>  <chr>  <chr>                 <int> <chr> <chr>   <chr> <chr>           
#> 1 L2073… S6617… David Wood               10 24 M… 14:47   S661… L207391,Mt. Aub…
#> 2 L2073… S6617… Sofia Prado-Irw…         15 24 M… 14:31   S661… L207391,Mt. Aub…
#> 3 L2073… S6619… Jeffrey Gantz            19 24 M… 13:30   S661… L207391,Mt. Aub…
#> 4 L2073… S6617… Ann Gurka                21 24 M… 13:00   S661… L207391,Mt. Aub…
#> 5 L2073… S7098… Barbara Olson            20 24 M… 10:30   S709… L207391,Mt. Aub…
```

## Hotspots in a region or nearby coordinates

Obtain a list of hotspots within a region

``` r
ebirdhotspotlist("CA-NS-HL")
#> # A tibble: 220 x 9
#>    locId locName countryCode subnational1Code subnational2Code   lat   lng
#>    <chr> <chr>   <chr>       <chr>            <chr>            <dbl> <dbl>
#>  1 L233… Abraha… CA          CA-NS            CA-NS-HL          45.2 -62.6
#>  2 L700… Admira… CA          CA-NS            CA-NS-HL          44.7 -63.7
#>  3 L176… Admira… CA          CA-NS            CA-NS-HL          44.8 -63.1
#>  4 L584… Albro … CA          CA-NS            CA-NS-HL          44.7 -63.6
#>  5 L437… Aldern… CA          CA-NS            CA-NS-HL          44.7 -63.6
#>  6 L122… Armdal… CA          CA-NS            CA-NS-HL          44.6 -63.6
#>  7 L624… Atlant… CA          CA-NS            CA-NS-HL          44.7 -63.3
#>  8 L239… Bald R… CA          CA-NS            CA-NS-HL          44.5 -63.6
#>  9 L759… Bayers… CA          CA-NS            CA-NS-HL          44.6 -63.7
#> 10 L642… Beaufo… CA          CA-NS            CA-NS-HL          44.7 -63.5
#> # … with 210 more rows, and 2 more variables: latestObsDt <chr>,
#> #   numSpeciesAllTime <int>
```

or within a radius of up to 50 kilometers, from a given set of
coordinates.

``` r
ebirdhotspotlist(lat = 30, lng = -90, dist = 10)
#> No region code provided, locating hotspots using lat/lng
#> # A tibble: 52 x 9
#>    locId locName countryCode subnational1Code subnational2Code   lat   lng
#>    <chr> <chr>   <chr>       <chr>            <chr>            <dbl> <dbl>
#>  1 L602… Algier… US          US-LA            US-LA-071         30.0 -90.1
#>  2 L388… Armstr… US          US-LA            US-LA-071         30.0 -90.1
#>  3 L727… Audubo… US          US-LA            US-LA-071         30.0 -90.0
#>  4 L666… BAEA N… US          US-LA            US-LA-087         30.0 -90.0
#>  5 L666… BAEA N… US          US-LA            US-LA-071         29.9 -90.0
#>  6 L242… Bayou … US          US-LA            US-LA-071         30.0 -90.0
#>  7 L725… Bayou … US          US-LA            US-LA-071         30.1 -89.9
#>  8 L727… Chalme… US          US-LA            US-LA-087         29.9 -90.0
#>  9 L453… City P… US          US-LA            US-LA-071         30.0 -90.1
#> 10 L522… City P… US          US-LA            US-LA-071         30.0 -90.1
#> # … with 42 more rows, and 2 more variables: latestObsDt <chr>,
#> #   numSpeciesAllTime <int>
```

## `rebird` and other packages

### How to use `rebird`

This package is part of a richer suite called [spocc - Species
Occurrence Data](https://github.com/ropensci/spocc), along with several
other packages, that provide access to occurrence records from multiple
databases. We recommend using `spocc` as the primary R interface to
`rebird` unless your needs are limited to this single source.

### `auk` vs. `rebird`

Those interested in eBird data may also want to consider
[`auk`](https://github.com/CornellLabofOrnithology/auk), an R package
that helps extracting and processing the whole eBird dataset. The
functions in `rebird` are faster but mostly limited to accessing recent
(i.e. within the last 30 days) observations, although `ebirdfreq()` does
provide historical frequency of observation data. In contrast, `auk`
gives access to the full set of \~ 500 million eBird observations. For
most ecological applications, users will require `auk`; however, for
some use cases, e.g. building tools for birders, `rebird` provides a
quicker and easier way to access data. `rebird` and `auk` are both part
of the rOpenSci project.

## API requests covered by `rebird`

The 2.0 APIs have considerably been expanded from the previous version,
and `rebird` only covers some of them. The webservices covered are
listed below; if you’d like to contribute wrappers to APIs not yet
covered by this package, feel free to submit a pull request\!

### data/obs

  - [x] Recent observations in a region: `ebirdregion()`
  - [x] Recent notable observations in a region: `ebirdnotable()`
  - [x] Recent observations of a species in a region: `ebirdregion()`
  - [x] Recent nearby observations: `ebirdgeo()`
  - [x] Recent nearby observations of a species: `ebirdgeo()`
  - [x] Nearest observations of a species: `nearestobs()`
  - [x] Recent nearby notable observations: `ebirdnotable()`
  - [ ] Recent checklists feed
  - [x] Historic observations on a date: `ebirdhistorical()`

### product

  - [ ] Top 100
  - [x] Checklist feed on a date: `ebirdchecklistfeed()`
  - [ ] Regional statistics on a date
  - [x] Species list for a region: `ebirdregionspecies()`
  - [ ] View Checklist BETA

### ref/geo

  - [ ] Adjacent Regions

### ref/hotspot

  - [x] Hotspots in a region: `ebirdhotspotlist()`
  - [x] Nearby hotspots: `ebirdhotspotlist()`
  - [x] Hotspot Info: `ebirdregioninfo()`

### ref/taxonomy

  - [x] eBird Taxonomy: `ebirdtaxonomy()`
  - [ ] Taxonomic Forms
  - [ ] Taxonomy Versions
  - [ ] Taxonomic Groups

### ref/region

  - [x] Region Info: `ebirdregioninfo()`
  - [x] Sub Region List `ebirdsubregionlist()`

## Meta

  - Please [report any issues or
    bugs](https://github.com/ropensci/rebird/issues).
  - License: MIT
  - Get citation information for `rebird` in R doing `citation(package =
    'rebird')`
  - Please note that the ‘rebird’ project is released with a
    [Contributor Code of
    Conduct](https://github.com/ropensci/rebird/blob/master/CODE_OF_CONDUCT.md).
    By contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
