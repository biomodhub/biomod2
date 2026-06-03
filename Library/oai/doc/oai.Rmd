<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{oai introduction}
%\VignetteEncoding{UTF-8}
-->



oai introduction
================

A general purpose client to work with any 'OAI-PMH' service. 

The 'OAI-PMH' protocol is described at <http://www.openarchives.org/OAI/openarchivesprotocol.html>.

The main functions follow the OAI-PMH verbs:

* `GetRecord`
* `Identify`
* `ListIdentifiers`
* `ListMetadataFormats`
* `ListRecords`
* `ListSets`

## Get oai

Install from CRAN


```r
install.packages("oai")
```

Or install the development version from GitHub


```r
devtools::install_github("ropensci/oai")
```

Load `oai`


```r
library("oai")
```

## Identify


```r
id("http://oai.datacite.org/oai")
#>   repositoryName                      baseURL protocolVersion
#> 1   DataCite MDS https://oai.datacite.org/oai             2.0
#>             adminEmail    earliestDatestamp deletedRecord
#> 1 support@datacite.org 2011-01-01T00:00:00Z    persistent
#>            granularity compression compression.1
#> 1 YYYY-MM-DDThh:mm:ssZ        gzip       deflate
#>                                      description
#> 1 oaioai.datacite.org:oai:oai.datacite.org:12425
```

## ListIdentifiers


```r
list_identifiers(from = '2018-05-01T', until = '2018-09-01T')
#> # A tibble: 2,904 x 5
#>    identifier        datestamp    setSpec           setSpec.1     setSpec.2
#>    <chr>             <chr>        <chr>             <chr>         <chr>    
#>  1 cbe4cdda-2045-41… 2018-08-30T… installation:842… dataset_type… country:…
#>  2 3fb7ddd8-07c0-49… 2018-08-28T… installation:842… dataset_type… country:…
#>  3 34585f24-1ffe-47… 2018-08-28T… installation:842… dataset_type… country:…
#>  4 8ac24ec8-1e7b-40… 2018-08-28T… installation:842… dataset_type… country:…
#>  5 e381b970-9b62-46… 2018-08-28T… installation:842… dataset_type… country:…
#>  6 4496b1b1-1ea9-4e… 2018-08-28T… installation:842… dataset_type… country:…
#>  7 fc5d68a6-b511-4c… 2018-08-27T… installation:73e… dataset_type… country:…
#>  8 5e0073bc-de2a-4b… 2018-08-27T… installation:688… dataset_type… country:…
#>  9 ab2684cf-62e5-4f… 2018-08-27T… installation:804… dataset_type… country:…
#> 10 34fbcf59-d9bb-47… 2018-08-27T… installation:7c5… dataset_type… country:…
#> # … with 2,894 more rows
```

## Count Identifiers


```r
count_identifiers()
```

## ListRecords


```r
list_records(from = '2018-05-01T', until = '2018-05-15T')
#> # A tibble: 44 x 26
#>    identifier datestamp setSpec setSpec.1 setSpec.2 title publisher
#>    <chr>      <chr>     <chr>   <chr>     <chr>     <chr> <chr>    
#>  1 18799ce9-… 2018-05-… instal… dataset_… country:… Bird… Sokoine …
#>  2 79f51633-… 2018-05-… instal… dataset_… country:… Impl… Aïgos SAS
#>  3 f83746ee-… 2018-05-… instal… dataset_… country:… NDFF… Dutch Na…
#>  4 a3533a61-… 2018-05-… instal… dataset_… country:… EDP … EDP - En…
#>  5 ba9b66a3-… 2018-05-… instal… dataset_… country:… Ende… Sokoine …
#>  6 78b696d9-… 2018-05-… instal… dataset_… country:… Ende… Sokoine …
#>  7 c791b255-… 2018-05-… instal… dataset_… country:… Ende… Sokoine …
#>  8 b929ccda-… 2018-05-… instal… dataset_… country:… List… Sokoine …
#>  9 da285c2a-… 2018-05-… instal… dataset_… country:… Moni… Corporac…
#> 10 87372877-… 2018-05-… instal… dataset_… country:… Moni… Corporac…
#> # … with 34 more rows, and 19 more variables: identifier.1 <chr>,
#> #   subject <chr>, source <chr>, description <chr>, description.1 <chr>,
#> #   type <chr>, creator <chr>, date <chr>, language <chr>, coverage <chr>,
#> #   coverage.1 <chr>, format <chr>, source.1 <chr>, subject.1 <chr>,
#> #   coverage.2 <chr>, creator.1 <chr>, description.2 <chr>,
#> #   creator.2 <chr>, subject.2 <chr>
```

## GetRecords


```r
ids <- c("87832186-00ea-44dd-a6bf-c2896c4d09b4", "d981c07d-bc43-40a2-be1f-e786e25106ac")
get_records(ids)
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$header
#> # A tibble: 1 x 3
#>   identifier             datestamp      setSpec                            
#>   <chr>                  <chr>          <chr>                              
#> 1 87832186-00ea-44dd-a6… 2018-06-29T12… installation:729a7375-b120-4e4f-bb…
#> 
#> $`87832186-00ea-44dd-a6bf-c2896c4d09b4`$metadata
#> # A tibble: 0 x 0
#> 
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$header
#> # A tibble: 1 x 3
#>   identifier             datestamp      setSpec                            
#>   <chr>                  <chr>          <chr>                              
#> 1 d981c07d-bc43-40a2-be… 2018-01-21T21… installation:804b8dd0-07ac-4a30-bf…
#> 
#> $`d981c07d-bc43-40a2-be1f-e786e25106ac`$metadata
#> # A tibble: 1 x 12
#>   title publisher identifier subject source description type  creator date 
#>   <chr> <chr>     <chr>      <chr>   <chr>  <chr>       <chr> <chr>   <chr>
#> 1 Pece… Institut… https://w… peces,… http:… Caracteriz… Data… Fernan… 2018…
#> # … with 3 more variables: language <chr>, coverage <chr>, format <chr>
```

## List MetadataFormats


```r
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


```r
list_sets("http://api.gbif.org/v1/oai-pmh/registry")
#> # A tibble: 572 x 2
#>    setSpec                     setName         
#>    <chr>                       <chr>           
#>  1 dataset_type                per dataset type
#>  2 dataset_type:OCCURRENCE     occurrence      
#>  3 dataset_type:CHECKLIST      checklist       
#>  4 dataset_type:METADATA       metadata        
#>  5 dataset_type:SAMPLING_EVENT sampling_event  
#>  6 country                     per country     
#>  7 country:AD                  Andorra         
#>  8 country:AO                  Angola          
#>  9 country:AR                  Argentina       
#> 10 country:AT                  Austria         
#> # … with 562 more rows
```

### Biodiversity Heritage Library

Identify


```r
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


```r
get_records(c("oai:biodiversitylibrary.org:item/7", "oai:biodiversitylibrary.org:item/9"),
            url = "http://www.biodiversitylibrary.org/oai")
#> $`oai:biodiversitylibrary.org:item/7`
#> $`oai:biodiversitylibrary.org:item/7`$header
#> # A tibble: 1 x 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/7 2016-07-13T09:13:41Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/7`$metadata
#> # A tibble: 1 x 11
#>   title creator subject description publisher contributor date  type 
#>   <chr> <chr>   <chr>   <chr>       <chr>     <chr>       <chr> <chr>
#> 1 Die … Fleisc… Bogor;… pt.5:v.1 (… Leiden :… Missouri B… 1900… text…
#> # … with 3 more variables: identifier <chr>, language <chr>, rights <chr>
#> 
#> 
#> $`oai:biodiversitylibrary.org:item/9`
#> $`oai:biodiversitylibrary.org:item/9`$header
#> # A tibble: 1 x 3
#>   identifier                         datestamp            setSpec
#>   <chr>                              <chr>                <chr>  
#> 1 oai:biodiversitylibrary.org:item/9 2016-07-13T09:13:41Z item   
#> 
#> $`oai:biodiversitylibrary.org:item/9`$metadata
#> # A tibble: 1 x 11
#>   title creator subject description publisher contributor date  type 
#>   <chr> <chr>   <chr>   <chr>       <chr>     <chr>       <chr> <chr>
#> 1 Die … Fleisc… Bogor;… pt.5:v.3 (… Leiden :… Missouri B… 1906… text…
#> # … with 3 more variables: identifier <chr>, language <chr>, rights <chr>
```
