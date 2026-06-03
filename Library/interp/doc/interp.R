### R code from vignette source 'interp.Rnw'

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
### code chunk number 2: akima
###################################################
data(akima)
library(scatterplot3d)
scatterplot3d(akima, type="h", angle=60, asp=0.2, lab=c(4,4,0))


###################################################
### code chunk number 3: interp.Rnw:278-279
###################################################
getOption("SweaveHooks")[["fig"]]()
data(akima)
library(scatterplot3d)
scatterplot3d(akima, type="h", angle=60, asp=0.2, lab=c(4,4,0))


###################################################
### code chunk number 4: lininterp
###################################################
li <- interp(akima$x, akima$y, akima$z, nx=150, ny=150)
MASS::eqscplot(akima$x, akima$y)
contour(li, nlevels=30, add=TRUE)
plot(tri.mesh(akima), add=TRUE)


###################################################
### code chunk number 5: interp.Rnw:297-298
###################################################
getOption("SweaveHooks")[["fig"]]()
li <- interp(akima$x, akima$y, akima$z, nx=150, ny=150)
MASS::eqscplot(akima$x, akima$y)
contour(li, nlevels=30, add=TRUE)
plot(tri.mesh(akima), add=TRUE)


###################################################
### code chunk number 6: splinterp
###################################################
si <- interp(akima$x, akima$y, akima$z, method="akima", nx=150, ny=150)
MASS::eqscplot(akima$x, akima$y)
contour(si, nlevels=30, add=TRUE)
plot(tri.mesh(akima), add=TRUE)


###################################################
### code chunk number 7: interp.Rnw:410-411
###################################################
getOption("SweaveHooks")[["fig"]]()
si <- interp(akima$x, akima$y, akima$z, method="akima", nx=150, ny=150)
MASS::eqscplot(akima$x, akima$y)
contour(si, nlevels=30, add=TRUE)
plot(tri.mesh(akima), add=TRUE)


###################################################
### code chunk number 8: splinterpnobw
###################################################
si.nobw <- interp(akima$x, akima$y, akima$z, method="akima", nx=150, ny=150,
                  baryweight=FALSE)
MASS::eqscplot(akima$x, akima$y)
contour(si.nobw, nlevels=30, add=TRUE)
plot(tri.mesh(akima), add=TRUE)


###################################################
### code chunk number 9: interp.Rnw:552-553
###################################################
getOption("SweaveHooks")[["fig"]]()
si.nobw <- interp(akima$x, akima$y, akima$z, method="akima", nx=150, ny=150,
                  baryweight=FALSE)
MASS::eqscplot(akima$x, akima$y)
contour(si.nobw, nlevels=30, add=TRUE)
plot(tri.mesh(akima), add=TRUE)


###################################################
### code chunk number 10: bilinear
###################################################
nx <- 8; ny <- 8
xg<-seq(0,1,length=nx)
yg<-seq(0,1,length=ny)
xyg<-expand.grid(xg,yg)
fg <- outer(xg,yg,function(x,y)franke.fn(x,y,1))
# not yet implemented this way:
# bil <- interp(xg,yg,fg,input="grid",output="grid",method="bilinear")
bil <- bilinear.grid(xg, yg, fg, dx=0.01, dy=0.01)
MASS::eqscplot(xyg[,1], xyg[,2])
contour(bil, add=TRUE)


###################################################
### code chunk number 11: interp.Rnw:617-618
###################################################
getOption("SweaveHooks")[["fig"]]()
nx <- 8; ny <- 8
xg<-seq(0,1,length=nx)
yg<-seq(0,1,length=ny)
xyg<-expand.grid(xg,yg)
fg <- outer(xg,yg,function(x,y)franke.fn(x,y,1))
# not yet implemented this way:
# bil <- interp(xg,yg,fg,input="grid",output="grid",method="bilinear")
bil <- bilinear.grid(xg, yg, fg, dx=0.01, dy=0.01)
MASS::eqscplot(xyg[,1], xyg[,2])
contour(bil, add=TRUE)


###################################################
### code chunk number 12: aspline
###################################################
x <- c(-3, -2, -1, 0,  1,  2, 2.5, 3)
y <- c( 0,  0,  0, 0, -1, -1, 0,   2)
MASS::eqscplot(x, y, ylim=c(-2, 3))
lines(aspline(x, y, n=200, method="original"), col="red")
lines(aspline(x, y, n=200, method="improved"), col="black", lty="dotted")
lines(aspline(x, y, n=200, method="improved", degree=10), col="green", lty="dashed")


###################################################
### code chunk number 13: interp.Rnw:654-655
###################################################
getOption("SweaveHooks")[["fig"]]()
x <- c(-3, -2, -1, 0,  1,  2, 2.5, 3)
y <- c( 0,  0,  0, 0, -1, -1, 0,   2)
MASS::eqscplot(x, y, ylim=c(-2, 3))
lines(aspline(x, y, n=200, method="original"), col="red")
lines(aspline(x, y, n=200, method="improved"), col="black", lty="dotted")
lines(aspline(x, y, n=200, method="improved", degree=10), col="green", lty="dashed")


