### R code from vignette source 'tri.Rnw'

###################################################
### code chunk number 1: init
###################################################
set.seed(42)
options(width=80)
options(continue=" ")
options(SweaveHooks=list(fig=function()
    par(mar=c(5.1, 4.1, 1.1, 2.1))))
library(interp)


###################################################
### code chunk number 2: tri.mesh
###################################################
data(tritest)
tr <- tri.mesh(tritest)
tr


###################################################
### code chunk number 3: triangles
###################################################
triangles(tr)


###################################################
### code chunk number 4: plottri
###################################################
MASS::eqscplot(tritest)
plot(tr, do.circumcircles=TRUE, add=TRUE)


###################################################
### code chunk number 5: tri.Rnw:302-303
###################################################
getOption("SweaveHooks")[["fig"]]()
MASS::eqscplot(tritest)
plot(tr, do.circumcircles=TRUE, add=TRUE)


###################################################
### code chunk number 6: vm
###################################################
vm <- voronoi.mosaic(tr)
vm


###################################################
### code chunk number 7: plotvm
###################################################
MASS::eqscplot(tritest)
plot(vm, add=TRUE)
plot(tr, add=TRUE)


###################################################
### code chunk number 8: tri.Rnw:357-358
###################################################
getOption("SweaveHooks")[["fig"]]()
MASS::eqscplot(tritest)
plot(vm, add=TRUE)
plot(tr, add=TRUE)


