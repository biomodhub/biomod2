# ridigbio <img src="man/figures/ridigbioLogo4.png" align="right" alt="" width="200">

[![Build Status](https://api.travis-ci.com/iDigBio/ridigbio.svg?branch=master)](https://app.travis-ci.com/github/iDigBio/ridigbio)


**ridigbio** is an R package to access [iDigBio](https://www.idigbio.org/) (Integrated Digitized Biocollections). 

## Installation

ridigbio is avaliable via [CRAN](https://cran.r-project.org/). 

```r
install.packages("ridigbio")
```

### Error Messages 	

If R says the package is unavailable, you may not have set a CRAN mirror. You can do so with:

```r
chooseCRANmirror()
```

If R says that a binary package is not available, your version of R may be too old. Please 
review the versions of R that CRAN has built packages for on the [CRAN ridigbio package page.]( https://cran.r-project.org/package=ridigbio)

You can download the source package and install manually if there is no package built for 
your version of R. You may also need to install any dependencies.

```r
install.packages("ridigbio", type="source")
```

On Linux, you may encounter an error during the installation process if you do not have `libcurl` installed. The method for installing libcurl will vary between distributions, but on Ubuntu you can install the latest version via:

```
sudo apt install libcurl4
```
    
## Getting Started
There are several articles that can help get you started:

* [Introduction to ridigbio](https://idigbio.github.io/ridigbio/articles/BasicUsage.html)
* [Record API Demo](https://idigbio.github.io/ridigbio/articles/RecordAPIDemo.html)
* [Media API Demo](https://idigbio.github.io/ridigbio/articles/MediaAPIDemo.html)
* [Fields in ridigibio](https://idigbio.github.io/ridigbio/articles/Fields.html)
* [Tissue Samples Locator Demo](https://idigbio.github.io/ridigbio/articles/FindTissue.html)
* [Identification of Modified Data](https://idigbio.github.io/ridigbio/articles/ModifiedDataID.html)
* [Identification of Suspicious Coordinates](https://idigbio.github.io/ridigbio/articles/BadCoordinateID.html)
* [Identification of Data Flags](https://idigbio.github.io/ridigbio/articles/IDDataFlags.html)


Most iDigBio users are interested in downloading occurrence records:

```r
library("ridigbio")
idig_search_records(rq=list(genus="galax"))
idig_search_records(rq=list(family="diapensiaceae"), limit=1000)
```
# Meta

* [Please report any bugs or issues.](https://github.com/iDigBio/ridigbio/issues)
* License: MIT
