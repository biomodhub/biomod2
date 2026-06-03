
# `oai`: General Purpose ‘Oai-PMH’ Services Client <img src="man/figures/logo.png" align="right" width="20%" />

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/oai/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ropensci/oai/actions/workflows/R-CMD-check.yaml)
[![cran
checks](https://cranchecks.info/badges/worst/oai)](https://cranchecks.info/pkgs/oai)
[![codecov.io](https://codecov.io/github/ropensci/oai/coverage.svg?branch=master)](https://codecov.io/github/ropensci/oai?branch=master)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/oai?color=2ED968)](https://github.com/r-hub/cranlogs.app)
[![cran
version](https://www.r-pkg.org/badges/version/oai)](https://cran.r-project.org/package=oai)
[![](https://badges.ropensci.org/19_status.svg)](https://github.com/ropensci/software-review/issues/19)

`oai` is an R client to work with OAI-PMH (Open Archives Initiative
Protocol for Metadata Harvesting) services, a protocol developed by the
Open Archives Initiative
(<https://en.wikipedia.org/wiki/Open_Archives_Initiative>). OAI-PMH uses
XML data format transported over HTTP.

OAI-PMH Info:

-   Wikipedia
    (<https://en.wikipedia.org/wiki/Open_Archives_Initiative_Protocol_for_Metadata_Harvesting>)
-   OAI V2 specification
    (<http://www.openarchives.org/OAI/openarchivesprotocol.html>)

`oai` is built on `xml2` and `httr`. In addition, we give back
data.frame’s whenever possible to make data comprehension, manipulation,
and visualization easier. We also have functions to fetch a large
directory of OAI-PMH services - it isn’t exhaustive, but does contain a
lot.

OAI-PMH instead of paging with e.g., `page` and `per_page` parameters,
uses (optionally) `resumptionTokens`, optionally with an expiration
date. These tokens can be used to continue on to the next chunk of data,
if the first request did not get to the end. Often, OAI-PMH services
limit each request to 50 records, but this may vary by provider, I don’t
know for sure. The API of this package is such that we `while` loop for
you internally until we get all records. We may in the future expose
e.g., a `limit` parameter so you can say how many records you want, but
we haven’t done this yet.

## Install

Install from CRAN

``` r
install.packages("oai")
```

Development version

``` r
devtools::install_github("ropensci/oai")
```

``` r
library("oai")
```

## Identify

``` r
id("http://oai.datacite.org/oai")
#>   repositoryName                      baseURL protocolVersion
#> 1       DataCite https://oai.datacite.org/oai             2.0
#>             adminEmail    earliestDatestamp deletedRecord          granularity
#> 1 support@datacite.org 2011-01-01T00:00:00Z    persistent YYYY-MM-DDThh:mm:ssZ
#>   compression compression.1                                    description
#> 1        gzip       deflate oaioai.datacite.org:oai:oai.datacite.org:12425
```

## ListIdentifiers

``` r
list_identifiers(from = '2018-05-01T', until = '2018-06-01T')
#> # A tibble: 75 × 5
#>    identifier                           datestamp        setSpec setSp…¹ setSp…²
#>    <chr>                                <chr>            <chr>   <chr>   <chr>  
#>  1 4b64d1f2-31c2-40c9-80aa-bb7ddb424684 2018-05-30T13:5… instal… datase… countr…
#>  2 884378d6-d591-4760-bb70-7b4851784d96 2018-05-29T19:1… instal… datase… countr…
#>  3 18799ce9-1a66-40fc-ad18-5ac54cd3417b 2018-05-14T12:1… instal… datase… countr…
#>  4 7e91aacb-c994-41ee-a7b7-bd23c02cd5bf 2018-05-21T10:5… instal… datase… countr…
#>  5 f83746ee-4cf2-4e60-a720-dd508b559794 2018-05-08T09:4… instal… datase… countr…
#>  6 a3533a61-6f88-443e-89ae-37611ea88267 2018-05-08T13:5… instal… datase… countr…
#>  7 ba9b66a3-2d11-4193-922e-ace4d5909239 2018-05-05T23:5… instal… datase… countr…
#>  8 78b696d9-8f0d-41ab-9c23-1c3547da411d 2018-05-05T23:0… instal… datase… countr…
#>  9 c791b255-a184-4600-b828-ef9d4092a212 2018-05-05T14:2… instal… datase… countr…
#> 10 b929ccda-03b1-4166-9e5b-34588339d61d 2018-05-09T02:5… instal… datase… countr…
#> # … with 65 more rows, and abbreviated variable names ¹​setSpec.1, ²​setSpec.2
```

## Count Identifiers

``` r
count_identifiers()
#>                            url   count
#> 1 http://export.arxiv.org/oai2 2158148
```

## ListRecords

``` r
list_records(from = '2018-05-01T', until = '2018-05-15T')
#> # A tibble: 41 × 26
#>    identi…¹ dates…² setSpec setSp…³ setSp…⁴ title publi…⁵ ident…⁶ subject source
#>    <chr>    <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>   <chr> 
#>  1 18799ce… 2018-0… instal… datase… countr… Bird… Sokoin… https:… Occurr… ""    
#>  2 f83746e… 2018-0… instal… datase… countr… NDFF… Dutch … https:… Metada… "http…
#>  3 a3533a6… 2018-0… instal… datase… countr… EDP … EDP - … https:… Occurr… ""    
#>  4 ba9b66a… 2018-0… instal… datase… countr… Ende… Sokoin… https:… Occurr… ""    
#>  5 78b696d… 2018-0… instal… datase… countr… Ende… Sokoin… https:… Occurr… ""    
#>  6 c791b25… 2018-0… instal… datase… countr… Ende… Sokoin… https:… Occurr… ""    
#>  7 b929ccd… 2018-0… instal… datase… countr… List… Sokoin… https:… Occurr… ""    
#>  8 da285c2… 2018-0… instal… datase… countr… Moni… Corpor… https:… seguim… ""    
#>  9 8737287… 2018-0… instal… datase… countr… Moni… Corpor… https:… seguim… ""    
#> 10 ed7d4c2… 2018-0… instal… datase… countr… Samo… Minist… https:… Occurr… ""    
#> # … with 31 more rows, 16 more variables: description <chr>,
#> #   description.1 <chr>, type <chr>, creator <chr>, date <chr>, language <chr>,
#> #   coverage <chr>, coverage.1 <chr>, format <chr>, source.1 <chr>,
#> #   subject.1 <chr>, creator.1 <chr>, coverage.2 <chr>, description.2 <chr>,
#> #   creator.2 <chr>, subject.2 <chr>, and abbreviated variable names
#> #   ¹​identifier, ²​datestamp, ³​setSpec.1, ⁴​setSpec.2, ⁵​publisher, ⁶​identifier.1
```

## GetRecords

``` r
ids <- c("87832186-00ea-44dd-a6bf-c2896c4d09b4", "d981c07d-bc43-40a2-be1f-e786e25106ac")
get_records(ids)
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$header
#> # A tibble: 1 × 3
#>   identifier                           datestamp            setSpec             
#>   <chr>                                <chr>                <chr>               
#> 1 87832186-00ea-44dd-a6bf-c2896c4d09b4 2018-06-29T12:08:17Z installation:729a73…
#> 
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$metadata
#> # A tibble: 0 × 0
#> 
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$header
#> # A tibble: 1 × 3
#>   identifier                           datestamp            setSpec             
#>   <chr>                                <chr>                <chr>               
#> 1 d981c07d-bc43-40a2-be1f-e786e25106ac 2021-09-28T13:58:57Z installation:804b8d…
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$metadata
#> # A tibble: 1 × 12
#>   title       publi…¹ ident…² subject source descr…³ type  creator date  langu…⁴
#>   <chr>       <chr>   <chr>   <chr>   <chr>  <chr>   <chr> <chr>   <chr> <chr>  
#> 1 Peces de l… Instit… https:… Occurr… http:… Caract… Data… Fernan… 2021… es     
#> # … with 2 more variables: coverage <chr>, format <chr>, and abbreviated
#> #   variable names ¹​publisher, ²​identifier, ³​description, ⁴​language
```

## List MetadataFormats

``` r
list_metadataformats(id = "87832186-00ea-44dd-a6bf-c2896c4d09b4")
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`
#>   metadataPrefix                                                   schema
#> 1         oai_dc           http://www.openarchives.org/OAI/2.0/oai_dc.xsd
#> 2            eml http://rs.gbif.org/schema/eml-gbif-profile/1.0.2/eml.xsd
#>                             metadataNamespace
#> 1 http://www.openarchives.org/OAI/2.0/oai_dc/
#> 2          eml://ecoinformatics.org/eml-2.1.1
```

## List Sets

``` r
list_sets("http://api.gbif.org/v1/oai-pmh/registry")
#> # A tibble: 621 × 2
#>    setSpec                     setName         
#>    <chr>                       <chr>           
#>  1 dataset_type                per dataset type
#>  2 dataset_type:OCCURRENCE     occurrence      
#>  3 dataset_type:CHECKLIST      checklist       
#>  4 dataset_type:METADATA       metadata        
#>  5 dataset_type:SAMPLING_EVENT sampling_event  
#>  6 country                     per country     
#>  7 country:AD                  Andorra         
#>  8 country:AM                  Armenia         
#>  9 country:AO                  Angola          
#> 10 country:AQ                  Antarctica      
#> # … with 611 more rows
```

## Examples of other OAI providers

### Biodiversity Heritage Library

Identify

``` r
id("http://www.biodiversitylibrary.org/oai")
#>                                 repositoryName
#> 1 Biodiversity Heritage Library OAI Repository
#>                                   baseURL protocolVersion
#> 1 https://www.biodiversitylibrary.org/oai             2.0
#>                    adminEmail earliestDatestamp deletedRecord granularity
#> 1 oai@biodiversitylibrary.org        2006-01-01            no  YYYY-MM-DD
#>                                                        description
#> 1 oaibiodiversitylibrary.org:oai:biodiversitylibrary.org:item/1000
```

Get records

``` r
get_records(c("oai:biodiversitylibrary.org:item/7", "oai:biodiversitylibrary.org:item/9"),
            url = "http://www.biodiversitylibrary.org/oai")
#> $`oai:biodiversitylibrary.org:item/7`
#> $`oai:biodiversitylibrary.org:item/7`$header
#> # A tibble: 1 × 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/7 2016-01-26T06:05:19Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/7`$metadata
#> # A tibble: 1 × 11
#>   title    creator subject descr…¹ publi…² contr…³ type  ident…⁴ langu…⁵ relat…⁶
#>   <chr>    <chr>   <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>  
#> 1 Die Mus… Fleisc… Bogor;… pt.5:v… Leiden… Missou… text… https:… Dutch   https:…
#> # … with 1 more variable: rights <chr>, and abbreviated variable names
#> #   ¹​description, ²​publisher, ³​contributor, ⁴​identifier, ⁵​language, ⁶​relation
#> 
#> 
#> $`oai:biodiversitylibrary.org:item/9`
#> $`oai:biodiversitylibrary.org:item/9`$header
#> # A tibble: 1 × 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/9 2016-01-26T06:05:19Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/9`$metadata
#> # A tibble: 1 × 11
#>   title    creator subject descr…¹ publi…² contr…³ type  ident…⁴ langu…⁵ relat…⁶
#>   <chr>    <chr>   <chr>   <chr>   <chr>   <chr>   <chr> <chr>   <chr>   <chr>  
#> 1 Die Mus… Fleisc… Bogor;… pt.5:v… Leiden… Missou… text… https:… Dutch   https:…
#> # … with 1 more variable: rights <chr>, and abbreviated variable names
#> #   ¹​description, ²​publisher, ³​contributor, ⁴​identifier, ⁵​language, ⁶​relation
```

## Acknowledgements

Michał Bojanowski thanks National Science Centre for support through
grant 2012/07/D/HS6/01971.

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/oai/issues).
-   License: MIT
-   Get citation information for `oai` in R doing
    `citation(package = 'oai')`
-   Please note that this project is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By participating
    in this project you agree to abide by its terms.
