# Proto 

[![Travis-CI Build Status](https://travis-ci.org/hadley/proto.svg?branch=master)](https://travis-ci.org/hadley/proto)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/proto)](https://CRAN.R-project.org/package=proto)
[![Coverage Status](https://img.shields.io/codecov/c/github/hadley/proto/master.svg)](https://codecov.io/github/hadley/proto?branch=master)

Proto is an R package that facilitates prototype
programming, a type of object-oriented programming that
does not use classes as an atomic concept (but is powerful
enough to encompass them).

The package is lightweight providing a thin layer on top of
R environments.  Unlike other packages which grow over time
proto has become smaller over time as it was successively
polished to reduce it to its essentials.  Despite its small
size prototype-based systems can be more powerful than more 
complex class-based systems.  

## Ease of use

The proto package is easy to use because:

1. few names. There is only one new function name to learn
among the core functions.  The 'proto' function constructs
new proto objects.  Remaining core functions include various
as.proto methods, assignment and retrieval via $ and
is.proto.

2. documentation. A 15 page report, a reference card, a demo
and help files are provided.

3.  consistency.  Proto objects form an subclass of the
environment class.  They are and work like environments.  
One can leverage everything one knows about environments 
to use proto.  The package is highly consistent with R and
works the way R works.

4. concise implementation.  The source code, excluding
dot.proto, is about one page of code making it possible to
rapidly understand it in its entirety not only from the
documentation but also by reading its source. (This should
not be necessary but its there for those who wish.)

5. tested.  The package has been independently tested by
multiple people and has no known bugs.  (If you find any
please let the developers know!)

## Example

The proto package is used like this:

```R
library(proto)

# new object with variable a and method addtwice
oo <- proto(a = 1, addtwice = function(., x) .$a <- .$a + 2*x)

oo$addtwice(3)  # add twice 3 to 1
oo$ls()         # "a" "addtwice"
oo$a            # 7

# create child object overriding a
ooc <- oo$proto(a = 10)

ooc$addtwice(1) # inherit addtwice from parent oo
ooc$ls()        # "a"
ooc$a           # 12 - addtwice from oo used "a" in ooc!
```

### Documentation

More information is available via:

```R
demo(proto)           # a self running demo

vignette("proto")     # 15 page report
vignette("protoref")  # reference card

?proto                # constructs proto objects
?dot.proto            # visualize a proto project
```
