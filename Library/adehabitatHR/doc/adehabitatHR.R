### R code from vignette source 'adehabitatHR.Rnw'

###################################################
### code chunk number 1: adehabitatHR.Rnw:28-35
###################################################
owidth <- getOption("width")
options("width"=80)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0
wi <- 600
pt <- 25


###################################################
### code chunk number 2: afig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
##           .PngNo, ".png", sep="")
## png(file=file, width = wi, height = wi, pointsize = pt)


###################################################
### code chunk number 3: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 4: zfigg (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 5: adehabitatHR.Rnw:103-104
###################################################
library(adehabitatHR)


###################################################
### code chunk number 6: adehabitatHR.Rnw:107-109
###################################################
suppressWarnings(RNGversion("3.5.0"))
set.seed(13431)


###################################################
### code chunk number 7: adehabitatHR.Rnw:172-174
###################################################
xy <- matrix(runif(60), ncol=2)
head(xy)


###################################################
### code chunk number 8: adehabitatHR.Rnw:180-181
###################################################
xysp <- SpatialPoints(xy)


###################################################
### code chunk number 9: adehabitatHR.Rnw:189-191
###################################################
clu <- clusthr(xysp)
class(clu)


###################################################
### code chunk number 10: fig1 (eval = FALSE)
###################################################
## plot(clu)


###################################################
### code chunk number 11: adehabitatHR.Rnw:200-201 (eval = FALSE)
###################################################
## plot(clu)


###################################################
### code chunk number 12: adehabitatHR.Rnw:205-208
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(clu)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 13: adehabitatHR.Rnw:215-216
###################################################
xy2 <- matrix(runif(60), ncol=2)


###################################################
### code chunk number 14: adehabitatHR.Rnw:221-222
###################################################
xyt <- rbind(xy, xy2)


###################################################
### code chunk number 15: adehabitatHR.Rnw:228-229
###################################################
id <- gl(2,30)


###################################################
### code chunk number 16: adehabitatHR.Rnw:235-238
###################################################
idsp <- data.frame(id)
coordinates(idsp) <- xyt
class(idsp)


###################################################
### code chunk number 17: adehabitatHR.Rnw:243-246
###################################################
clu2 <- clusthr(idsp)
class(clu2)
clu2


###################################################
### code chunk number 18: adehabitatHR.Rnw:253-256
###################################################
length(clu2)
class(clu2[[1]])
class(clu2[[2]])


###################################################
### code chunk number 19: fig2 (eval = FALSE)
###################################################
## plot(clu2)


###################################################
### code chunk number 20: adehabitatHR.Rnw:268-269 (eval = FALSE)
###################################################
## plot(clu2)


###################################################
### code chunk number 21: adehabitatHR.Rnw:273-276
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(clu2)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 22: adehabitatHR.Rnw:295-297
###################################################
data(puechabonsp)
names(puechabonsp)


###################################################
### code chunk number 23: adehabitatHR.Rnw:305-306
###################################################
head(as.data.frame(puechabonsp$relocs))


###################################################
### code chunk number 24: fig3 (eval = FALSE)
###################################################
## ## Map of the elevation
## image(puechabonsp$map, col=grey(c(1:10)/10))
## ## map of the relocations
## plot(puechabonsp$relocs, add=TRUE,
##      col=as.data.frame(puechabonsp$relocs)[,1])


###################################################
### code chunk number 25: adehabitatHR.Rnw:325-326 (eval = FALSE)
###################################################
## ## Map of the elevation
## image(puechabonsp$map, col=grey(c(1:10)/10))
## ## map of the relocations
## plot(puechabonsp$relocs, add=TRUE,
##      col=as.data.frame(puechabonsp$relocs)[,1])


###################################################
### code chunk number 26: adehabitatHR.Rnw:330-333
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
## Map of the elevation
image(puechabonsp$map, col=grey(c(1:10)/10))
## map of the relocations
plot(puechabonsp$relocs, add=TRUE,
     col=as.data.frame(puechabonsp$relocs)[,1])
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 27: adehabitatHR.Rnw:367-369
###################################################
data(puechabonsp)
cp <- mcp(puechabonsp$relocs[,1], percent=95)


###################################################
### code chunk number 28: adehabitatHR.Rnw:375-376
###################################################
class(cp)


###################################################
### code chunk number 29: fig4 (eval = FALSE)
###################################################
## plot(cp)
## plot(puechabonsp$relocs, add=TRUE)


###################################################
### code chunk number 30: adehabitatHR.Rnw:389-390 (eval = FALSE)
###################################################
## plot(cp)
## plot(puechabonsp$relocs, add=TRUE)


###################################################
### code chunk number 31: adehabitatHR.Rnw:394-397
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(cp)
plot(puechabonsp$relocs, add=TRUE)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 32: adehabitatHR.Rnw:406-407
###################################################
library(sf)


###################################################
### code chunk number 33: adehabitatHR.Rnw:410-411 (eval = FALSE)
###################################################
## st_write(st_as_sf(cp), "homerange.shp")


###################################################
### code chunk number 34: adehabitatHR.Rnw:425-426
###################################################
head(as.data.frame(puechabonsp$map))


###################################################
### code chunk number 35: adehabitatHR.Rnw:435-438
###################################################
enc <- over(puechabonsp$map, as(cp[2,],"SpatialPolygons"), fn=function(x) return(1))
enc <- as.data.frame(enc)
head(enc)


###################################################
### code chunk number 36: adehabitatHR.Rnw:447-449
###################################################
coordinates(enc) <- coordinates(puechabonsp$map)
gridded(enc) <- TRUE


###################################################
### code chunk number 37: fig5 (eval = FALSE)
###################################################
## image(puechabonsp$map)
## image(enc, add=TRUE, col="black", useRasterImage=FALSE)


###################################################
### code chunk number 38: adehabitatHR.Rnw:459-460 (eval = FALSE)
###################################################
## image(puechabonsp$map)
## image(enc, add=TRUE, col="black", useRasterImage=FALSE)


###################################################
### code chunk number 39: adehabitatHR.Rnw:464-467
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(puechabonsp$map)
image(enc, add=TRUE, col="black", useRasterImage=FALSE)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 40: adehabitatHR.Rnw:475-476
###################################################
cprast <- hr.rast(cp, puechabonsp$map)


###################################################
### code chunk number 41: sldkslsd (eval = FALSE)
###################################################
## par(mar=c(0,0,0,0))
## par(mfrow=c(2,2))
## for (i in 1:4) {
##     image(cprast[,i], useRasterImage=FALSE)
##     box()
## }


###################################################
### code chunk number 42: adehabitatHR.Rnw:491-492 (eval = FALSE)
###################################################
## par(mar=c(0,0,0,0))
## par(mfrow=c(2,2))
## for (i in 1:4) {
##     image(cprast[,i], useRasterImage=FALSE)
##     box()
## }


###################################################
### code chunk number 43: adehabitatHR.Rnw:496-499
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
par(mar=c(0,0,0,0))
par(mfrow=c(2,2))
for (i in 1:4) {
    image(cprast[,i], useRasterImage=FALSE)
    box()
}
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 44: adehabitatHR.Rnw:510-511
###################################################
as.data.frame(cp)


###################################################
### code chunk number 45: fig6 (eval = FALSE)
###################################################
## hrs <- mcp.area(puechabonsp$relocs[,1], percent=seq(50, 100, by = 5))


###################################################
### code chunk number 46: adehabitatHR.Rnw:532-533 (eval = FALSE)
###################################################
## hrs <- mcp.area(puechabonsp$relocs[,1], percent=seq(50, 100, by = 5))


###################################################
### code chunk number 47: adehabitatHR.Rnw:537-540
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
hrs <- mcp.area(puechabonsp$relocs[,1], percent=seq(50, 100, by = 5))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 48: adehabitatHR.Rnw:557-558
###################################################
hrs


###################################################
### code chunk number 49: fig7 (eval = FALSE)
###################################################
## xy1 <- matrix(rnorm(200), ncol=2)
## xy2 <- matrix(rnorm(200, mean=3), ncol=2)
## xy <- rbind(xy1,xy2)
## xy <- SpatialPoints(xy)
## kud <- kernelUD(xy)
## par(mfrow=c(1,2))
## plot(xy)
## title("Input: Relocations")
## par(mar=c(0,0,2,0))
## xxyz <- as.image.SpatialGridDataFrame(kud)
## persp(xxyz, theta=30, phi=45, box=FALSE, shade=0.5, border=NA, ltheta=300)
## title("Output: the UD")


###################################################
### code chunk number 50: adehabitatHR.Rnw:608-614
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = 2*wi, height = wi, pointsize = pt)
xy1 <- matrix(rnorm(200), ncol=2)
xy2 <- matrix(rnorm(200, mean=3), ncol=2)
xy <- rbind(xy1,xy2)
xy <- SpatialPoints(xy)
kud <- kernelUD(xy)
par(mfrow=c(1,2))
plot(xy)
title("Input: Relocations")
par(mar=c(0,0,2,0))
xxyz <- as.image.SpatialGridDataFrame(kud)
persp(xxyz, theta=30, phi=45, box=FALSE, shade=0.5, border=NA, ltheta=300)
title("Output: the UD")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 51: kdksks (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## par(mar=c(0,0,2,0))
## image(kernelUD(xy, h=0.2))
## title("h = 0.2")
## par(mar=c(0,0,2,0))
## image(kernelUD(xy, h=0.5))
## title("h = 0.5")
## par(mar=c(0,0,2,0))
## image(kernelUD(xy, h=1))
## title("h = 1")
## par(mar=c(0,0,2,0))
## image(kernelUD(xy, h=2))
## title("h = 2")


###################################################
### code chunk number 52: adehabitatHR.Rnw:699-702
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
par(mfrow=c(2,2))
par(mar=c(0,0,2,0))
image(kernelUD(xy, h=0.2))
title("h = 0.2")
par(mar=c(0,0,2,0))
image(kernelUD(xy, h=0.5))
title("h = 0.5")
par(mar=c(0,0,2,0))
image(kernelUD(xy, h=1))
title("h = 1")
par(mar=c(0,0,2,0))
image(kernelUD(xy, h=2))
title("h = 2")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 53: adehabitatHR.Rnw:769-772
###################################################
data(puechabonsp)
kud <- kernelUD(puechabonsp$relocs[,1], h="href")
kud


###################################################
### code chunk number 54: fig8 (eval = FALSE)
###################################################
## image(kud)


###################################################
### code chunk number 55: adehabitatHR.Rnw:784-785 (eval = FALSE)
###################################################
## image(kud)


###################################################
### code chunk number 56: adehabitatHR.Rnw:790-793
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(kud)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 57: adehabitatHR.Rnw:801-802
###################################################
kud[[1]]@h


###################################################
### code chunk number 58: fig9 (eval = FALSE)
###################################################
## kudl <- kernelUD(puechabonsp$relocs[,1], h="LSCV")
## image(kudl)


###################################################
### code chunk number 59: adehabitatHR.Rnw:822-823 (eval = FALSE)
###################################################
## kudl <- kernelUD(puechabonsp$relocs[,1], h="LSCV")
## image(kudl)


###################################################
### code chunk number 60: adehabitatHR.Rnw:828-831
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
kudl <- kernelUD(puechabonsp$relocs[,1], h="LSCV")
image(kudl)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 61: fig10 (eval = FALSE)
###################################################
## plotLSCV(kudl)


###################################################
### code chunk number 62: adehabitatHR.Rnw:848-849 (eval = FALSE)
###################################################
## plotLSCV(kudl)


###################################################
### code chunk number 63: adehabitatHR.Rnw:854-857
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotLSCV(kudl)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 64: kskdkkds (eval = FALSE)
###################################################
## ## The relocations of "Brock"
## locs <- puechabonsp$relocs
## firs <- locs[as.data.frame(locs)[,1]=="Brock",]
## 
## ## Graphical parameters
## par(mar=c(0,0,2,0))
## par(mfrow=c(2,2))
## 
## ## Estimation of the UD with grid=20 and extent=0.2
## image(kernelUD(firs, grid=20, extent=0.2))
## title(main="grid=20, extent=0.2")
## 
## ## Estimation of the UD with grid=20 and extent=0.2
## image(kernelUD(firs, grid=80, extent=0.2))
## title(main="grid=80, extent=0.2")
## 
## ## Estimation of the UD with grid=20 and extent=0.2
## image(kernelUD(firs, grid=20, extent=3))
## title(main="grid=20, extent=3")
## 
## ## Estimation of the UD with grid=20 and extent=0.2
## image(kernelUD(firs, grid=80, extent=3))
## title(main="grid=80, extent=3")


###################################################
### code chunk number 65: adehabitatHR.Rnw:947-948 (eval = FALSE)
###################################################
## ## The relocations of "Brock"
## locs <- puechabonsp$relocs
## firs <- locs[as.data.frame(locs)[,1]=="Brock",]
## 
## ## Graphical parameters
## par(mar=c(0,0,2,0))
## par(mfrow=c(2,2))
## 
## ## Estimation of the UD with grid=20 and extent=0.2
## image(kernelUD(firs, grid=20, extent=0.2))
## title(main="grid=20, extent=0.2")
## 
## ## Estimation of the UD with grid=20 and extent=0.2
## image(kernelUD(firs, grid=80, extent=0.2))
## title(main="grid=80, extent=0.2")
## 
## ## Estimation of the UD with grid=20 and extent=0.2
## image(kernelUD(firs, grid=20, extent=3))
## title(main="grid=20, extent=3")
## 
## ## Estimation of the UD with grid=20 and extent=0.2
## image(kernelUD(firs, grid=80, extent=3))
## title(main="grid=80, extent=3")


###################################################
### code chunk number 66: adehabitatHR.Rnw:952-955
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
## The relocations of "Brock"
locs <- puechabonsp$relocs
firs <- locs[as.data.frame(locs)[,1]=="Brock",]

## Graphical parameters
par(mar=c(0,0,2,0))
par(mfrow=c(2,2))

## Estimation of the UD with grid=20 and extent=0.2
image(kernelUD(firs, grid=20, extent=0.2))
title(main="grid=20, extent=0.2")

## Estimation of the UD with grid=20 and extent=0.2
image(kernelUD(firs, grid=80, extent=0.2))
title(main="grid=80, extent=0.2")

## Estimation of the UD with grid=20 and extent=0.2
image(kernelUD(firs, grid=20, extent=3))
title(main="grid=20, extent=3")

## Estimation of the UD with grid=20 and extent=0.2
image(kernelUD(firs, grid=80, extent=3))
title(main="grid=80, extent=3")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 67: kskkssqq (eval = FALSE)
###################################################
## kus <- kernelUD(puechabonsp$relocs[,1], same4all=TRUE)
## image(kus)


###################################################
### code chunk number 68: adehabitatHR.Rnw:969-970 (eval = FALSE)
###################################################
## kus <- kernelUD(puechabonsp$relocs[,1], same4all=TRUE)
## image(kus)


###################################################
### code chunk number 69: adehabitatHR.Rnw:974-977
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
kus <- kernelUD(puechabonsp$relocs[,1], same4all=TRUE)
image(kus)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 70: adehabitatHR.Rnw:984-986
###################################################
ii <- estUDm2spixdf(kus)
class(ii)


###################################################
### code chunk number 71: adehabitatHR.Rnw:1006-1007
###################################################
kudm <- kernelUD(puechabonsp$relocs[,1], grid=puechabonsp$map)


###################################################
### code chunk number 72: adehabitatHR.Rnw:1034-1036
###################################################
homerange <- getverticeshr(kudl)
class(homerange)


###################################################
### code chunk number 73: fig11 (eval = FALSE)
###################################################
## plot(homerange, col=1:4)


###################################################
### code chunk number 74: adehabitatHR.Rnw:1051-1052 (eval = FALSE)
###################################################
## plot(homerange, col=1:4)


###################################################
### code chunk number 75: adehabitatHR.Rnw:1056-1059
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(homerange, col=1:4)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 76: adehabitatHR.Rnw:1076-1078
###################################################
vud <- getvolumeUD(kudl)
vud


###################################################
### code chunk number 77: fig12 (eval = FALSE)
###################################################
## 
## ## Set up graphical parameters
## par(mfrow=c(2,1))
## par(mar=c(0,0,2,0))
## 
## ## The output of kernelUD for the first animal
## image(kudl[[1]])
## title("Output of kernelUD")
## 
## ## Convert into a suitable data structure for
## ## the use of contour
## xyz <- as.image.SpatialGridDataFrame(kudl[[1]])
## contour(xyz, add=TRUE)
## 
## 
## ## and similarly for the output of getvolumeUD
## par(mar=c(0,0,2,0))
## image(vud[[1]])
## title("Output of getvolumeUD")
## xyzv <- as.image.SpatialGridDataFrame(vud[[1]])
## contour(xyzv, add=TRUE)


###################################################
### code chunk number 78: adehabitatHR.Rnw:1109-1110 (eval = FALSE)
###################################################
## 
## ## Set up graphical parameters
## par(mfrow=c(2,1))
## par(mar=c(0,0,2,0))
## 
## ## The output of kernelUD for the first animal
## image(kudl[[1]])
## title("Output of kernelUD")
## 
## ## Convert into a suitable data structure for
## ## the use of contour
## xyz <- as.image.SpatialGridDataFrame(kudl[[1]])
## contour(xyz, add=TRUE)
## 
## 
## ## and similarly for the output of getvolumeUD
## par(mar=c(0,0,2,0))
## image(vud[[1]])
## title("Output of getvolumeUD")
## xyzv <- as.image.SpatialGridDataFrame(vud[[1]])
## contour(xyzv, add=TRUE)


###################################################
### code chunk number 79: adehabitatHR.Rnw:1114-1117
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)

## Set up graphical parameters
par(mfrow=c(2,1))
par(mar=c(0,0,2,0))

## The output of kernelUD for the first animal
image(kudl[[1]])
title("Output of kernelUD")

## Convert into a suitable data structure for
## the use of contour
xyz <- as.image.SpatialGridDataFrame(kudl[[1]])
contour(xyz, add=TRUE)


## and similarly for the output of getvolumeUD
par(mar=c(0,0,2,0))
image(vud[[1]])
title("Output of getvolumeUD")
xyzv <- as.image.SpatialGridDataFrame(vud[[1]])
contour(xyzv, add=TRUE)
dev.null <- dev.off()
cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 80: kdkdkdk (eval = FALSE)
###################################################
## ## store the volume under the UD (as computed by getvolumeUD)
## ## of the first animal in fud
## fud <- vud[[1]]
## 
## ## store the value of the volume under UD in a vector hr95
## hr95 <- as.data.frame(fud)[,1]
## 
## ## if hr95 is <= 95 then the pixel belongs to the home range
## ## (takes the value 1, 0 otherwise)
## hr95 <- as.numeric(hr95 <= 95)
## 
## ## Converts into a data frame
## hr95 <- data.frame(hr95)
## 
## ## Converts to a SpatialPixelsDataFrame
## coordinates(hr95) <- coordinates(vud[[1]])
## gridded(hr95) <- TRUE
## 
## ## display the results
## image(hr95)


###################################################
### code chunk number 81: adehabitatHR.Rnw:1151-1152 (eval = FALSE)
###################################################
## ## store the volume under the UD (as computed by getvolumeUD)
## ## of the first animal in fud
## fud <- vud[[1]]
## 
## ## store the value of the volume under UD in a vector hr95
## hr95 <- as.data.frame(fud)[,1]
## 
## ## if hr95 is <= 95 then the pixel belongs to the home range
## ## (takes the value 1, 0 otherwise)
## hr95 <- as.numeric(hr95 <= 95)
## 
## ## Converts into a data frame
## hr95 <- data.frame(hr95)
## 
## ## Converts to a SpatialPixelsDataFrame
## coordinates(hr95) <- coordinates(vud[[1]])
## gridded(hr95) <- TRUE
## 
## ## display the results
## image(hr95)


###################################################
### code chunk number 82: adehabitatHR.Rnw:1156-1159
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
## store the volume under the UD (as computed by getvolumeUD)
## of the first animal in fud
fud <- vud[[1]]

## store the value of the volume under UD in a vector hr95
hr95 <- as.data.frame(fud)[,1]

## if hr95 is <= 95 then the pixel belongs to the home range
## (takes the value 1, 0 otherwise)
hr95 <- as.numeric(hr95 <= 95)

## Converts into a data frame
hr95 <- data.frame(hr95)

## Converts to a SpatialPixelsDataFrame
coordinates(hr95) <- coordinates(vud[[1]])
gridded(hr95) <- TRUE

## display the results
image(hr95)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 83: adehabitatHR.Rnw:1172-1173
###################################################
as.data.frame(homerange)


###################################################
### code chunk number 84: adehabitatHR.Rnw:1185-1187
###################################################
ii <- kernel.area(kudl, percent=seq(50, 95, by=5))
ii


###################################################
### code chunk number 85: figggf (eval = FALSE)
###################################################
## par(mar=c(0,0,0,0))
## plot(c(0,1),asp=1, ty="n")
## bound <- structure(list(x = c(1.11651014886660, 1.31884344750132, 1.49071796999748,
## 1.56686491034388, 1.65171435815844, 1.627782462621, 1.59079680588132,
## 1.63648497008916), y = c(0.819277593991856, 0.797467606139553,
## 0.740761637723568, 0.609901710609753, 0.417973817509493, 0.263122903758146,
## 0.106090991221569, -0.00732094561040263)), .Names = c("x", "y"
## ))
## pts <- structure(list(x = c(1.14479329813812, 1.20135959668116, 1.18395458174484,
## 1.30143843256500, 1.29056029822980, 1.28620904449572, 1.31449219376724,
## 1.46896170132708, 1.39063913411364, 1.50377173119972, 1.627782462621,
## 1.57121616407796, 1.37758537291140, 1.21006210414932, 1.32537032810244,
## 1.519001119269, 1.46678607446004, 1.32754595496948, 1.17960332801076,
## 1.12303702946772, 1.19700834294708, 1.57121616407796, 1.57774304467908,
## 1.30578968629908, 1.60167494021652, 1.62995808948804, 1.5625136566098,
## 1.46896170132708, 1.43632729832148, 1.31449219376724), y = c(0.784381613428172,
## 0.732037642582647, 0.605539713039293, 0.559738738549458, 0.67315067538143,
## 0.749485632864488, 0.76257162557587, 0.714589652300805, 0.633892697247286,
## 0.330733866100284, 0.417973817509493, 0.527023756771005, 0.548833744623307,
## 0.378715839375349, 0.156253963281865, 0.287113890395679, 0.199873938986470,
## 0.354724852737816, 0.39180183208673, 0.311104877033211, 0.20423593655693,
## 0.106090991221569, 0.212959931697851, 0.780019615857712, 0.367810845449197,
## 0.424516813865184, 0.518299761630084, 0.601177715468833, 0.712408653515574,
## 0.703684658374653)), .Names = c("x", "y"))
## lines(bound)
## bound <- do.call("cbind",bound)
## Slo1 <- Line(bound)
## Sli1 <- Lines(list(Slo1), ID="frontier1")
## barrier <- SpatialLines(list(Sli1))
## 
## points(pts, pch=16)
## pts <- as.data.frame(pts)
## coordinates(pts) <- c("x","y")
## ii <- adehabitatHR:::.boundaryk(pts, barrier, 0.04)
## points(ii, pch=3)


###################################################
### code chunk number 86: adehabitatHR.Rnw:1250-1253
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
par(mar=c(0,0,0,0))
plot(c(0,1),asp=1, ty="n")
bound <- structure(list(x = c(1.11651014886660, 1.31884344750132, 1.49071796999748,
1.56686491034388, 1.65171435815844, 1.627782462621, 1.59079680588132,
1.63648497008916), y = c(0.819277593991856, 0.797467606139553,
0.740761637723568, 0.609901710609753, 0.417973817509493, 0.263122903758146,
0.106090991221569, -0.00732094561040263)), .Names = c("x", "y"
))
pts <- structure(list(x = c(1.14479329813812, 1.20135959668116, 1.18395458174484,
1.30143843256500, 1.29056029822980, 1.28620904449572, 1.31449219376724,
1.46896170132708, 1.39063913411364, 1.50377173119972, 1.627782462621,
1.57121616407796, 1.37758537291140, 1.21006210414932, 1.32537032810244,
1.519001119269, 1.46678607446004, 1.32754595496948, 1.17960332801076,
1.12303702946772, 1.19700834294708, 1.57121616407796, 1.57774304467908,
1.30578968629908, 1.60167494021652, 1.62995808948804, 1.5625136566098,
1.46896170132708, 1.43632729832148, 1.31449219376724), y = c(0.784381613428172,
0.732037642582647, 0.605539713039293, 0.559738738549458, 0.67315067538143,
0.749485632864488, 0.76257162557587, 0.714589652300805, 0.633892697247286,
0.330733866100284, 0.417973817509493, 0.527023756771005, 0.548833744623307,
0.378715839375349, 0.156253963281865, 0.287113890395679, 0.199873938986470,
0.354724852737816, 0.39180183208673, 0.311104877033211, 0.20423593655693,
0.106090991221569, 0.212959931697851, 0.780019615857712, 0.367810845449197,
0.424516813865184, 0.518299761630084, 0.601177715468833, 0.712408653515574,
0.703684658374653)), .Names = c("x", "y"))
lines(bound)
bound <- do.call("cbind",bound)
Slo1 <- Line(bound)
Sli1 <- Lines(list(Slo1), ID="frontier1")
barrier <- SpatialLines(list(Sli1))

points(pts, pch=16)
pts <- as.data.frame(pts)
coordinates(pts) <- c("x","y")
ii <- adehabitatHR:::.boundaryk(pts, barrier, 0.04)
points(ii, pch=3)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 87: mapelev (eval = FALSE)
###################################################
## image(puechabonsp$map)


###################################################
### code chunk number 88: adehabitatHR.Rnw:1296-1297 (eval = FALSE)
###################################################
## image(puechabonsp$map)


###################################################
### code chunk number 89: adehabitatHR.Rnw:1301-1304
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(puechabonsp$map)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 90: figboundary (eval = FALSE)
###################################################
## bound <- structure(list(x = c(701751.385381925, 701019.24105475,
##                         700739.303517889,
##                         700071.760160759, 699522.651915378,
##                         698887.40904327, 698510.570051342,
##                         698262.932999504, 697843.026694212,
##                         698058.363261028),
##                         y = c(3161824.03387414,
##                         3161824.03387414, 3161446.96718494,
##                         3161770.16720425, 3161479.28718687,
##                         3161231.50050539, 3161037.5804938,
##                         3160294.22044937, 3159389.26039528,
##                         3157482.3802813)), .Names = c("x", "y"))
## image(puechabonsp$map)
## lines(bound, lwd=3)


###################################################
### code chunk number 91: adehabitatHR.Rnw:1333-1334 (eval = FALSE)
###################################################
## bound <- structure(list(x = c(701751.385381925, 701019.24105475,
##                         700739.303517889,
##                         700071.760160759, 699522.651915378,
##                         698887.40904327, 698510.570051342,
##                         698262.932999504, 697843.026694212,
##                         698058.363261028),
##                         y = c(3161824.03387414,
##                         3161824.03387414, 3161446.96718494,
##                         3161770.16720425, 3161479.28718687,
##                         3161231.50050539, 3161037.5804938,
##                         3160294.22044937, 3159389.26039528,
##                         3157482.3802813)), .Names = c("x", "y"))
## image(puechabonsp$map)
## lines(bound, lwd=3)


###################################################
### code chunk number 92: adehabitatHR.Rnw:1338-1341
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
bound <- structure(list(x = c(701751.385381925, 701019.24105475,
                        700739.303517889,
                        700071.760160759, 699522.651915378,
                        698887.40904327, 698510.570051342,
                        698262.932999504, 697843.026694212,
                        698058.363261028),
                        y = c(3161824.03387414,
                        3161824.03387414, 3161446.96718494,
                        3161770.16720425, 3161479.28718687,
                        3161231.50050539, 3161037.5804938,
                        3160294.22044937, 3159389.26039528,
                        3157482.3802813)), .Names = c("x", "y"))
image(puechabonsp$map)
lines(bound, lwd=3)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 93: adehabitatHR.Rnw:1348-1353
###################################################
## We convert bound to SpatialLines:
bound <- do.call("cbind",bound)
Slo1 <- Line(bound)
Sli1 <- Lines(list(Slo1), ID="frontier1")
barrier <- SpatialLines(list(Sli1))


###################################################
### code chunk number 94: imagekudbound (eval = FALSE)
###################################################
## kud <- kernelUD(puechabonsp$relocs[,1], h=100, grid=100, boundary=barrier)
## image(kud)


###################################################
### code chunk number 95: adehabitatHR.Rnw:1364-1365 (eval = FALSE)
###################################################
## kud <- kernelUD(puechabonsp$relocs[,1], h=100, grid=100, boundary=barrier)
## image(kud)


###################################################
### code chunk number 96: adehabitatHR.Rnw:1369-1372
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
kud <- kernelUD(puechabonsp$relocs[,1], h=100, grid=100, boundary=barrier)
image(kud)
dev.null <- dev.off()
cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 97: brownbrid (eval = FALSE)
###################################################
## par(mar=c(0,0,0,0))
## xx <- c(0,1)
## yy <- c(0,1)
## date <- c(0,1)
## class(date) <- c("POSIXt", "POSIXct")
## tr <- as.ltraj(data.frame(x = xx,y = yy), date, id="a")
## kba <- kernelbb(tr, 0.6, 0.2, extent=0.8, grid=50)
## kba <- as(kba, "SpatialPixelsDataFrame")
## fullgrid(kba) <- TRUE
## uu <- slot(kba, "data")[,1]
## uu <- matrix(uu, nrow=50)
## persp(uu, theta = 135, phi = 30, scale = FALSE,
##       ltheta = -120, shade = 0.75, border = NA, box = FALSE)


###################################################
### code chunk number 98: adehabitatHR.Rnw:1413-1416
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
par(mar=c(0,0,0,0))
xx <- c(0,1)
yy <- c(0,1)
date <- c(0,1)
class(date) <- c("POSIXt", "POSIXct")
tr <- as.ltraj(data.frame(x = xx,y = yy), date, id="a")
kba <- kernelbb(tr, 0.6, 0.2, extent=0.8, grid=50)
kba <- as(kba, "SpatialPixelsDataFrame")
fullgrid(kba) <- TRUE
uu <- slot(kba, "data")[,1]
uu <- matrix(uu, nrow=50)
persp(uu, theta = 135, phi = 30, scale = FALSE,
      ltheta = -120, shade = 0.75, border = NA, box = FALSE)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 99: adehabitatHR.Rnw:1426-1427
###################################################
pt <- 16


###################################################
### code chunk number 100: tutorkbb (eval = FALSE)
###################################################
## suppressWarnings(RNGversion("3.5.0"))
## set.seed(2098)
## pts1 <- data.frame(x = rnorm(25, mean = 4.5, sd = 0.05),
##                    y = rnorm(25, mean = 4.5, sd = 0.05))
## pts1b <- data.frame(x = rnorm(25, mean = 4.5, sd = 0.05),
##                     y = rnorm(25, mean = 4.5, sd = 0.05))
## pts2 <- data.frame(x = rnorm(25, mean = 4, sd = 0.05),
##                    y = rnorm(25, mean = 4, sd = 0.05))
## pts3 <- data.frame(x = rnorm(25, mean = 5, sd = 0.05),
##                    y = rnorm(25, mean = 4, sd = 0.05))
## pts3b <- data.frame(x = rnorm(25, mean = 5, sd = 0.05),
##                     y = rnorm(25, mean = 4, sd = 0.05))
## pts2b <- data.frame(x = rnorm(25, mean = 4, sd = 0.05),
##                     y = rnorm(25, mean = 4, sd = 0.05))
## pts <- do.call("rbind", lapply(1:25, function(i) {
##     rbind(pts1[i,], pts1b[i,], pts2[i,], pts3[i,],
##           pts3b[i,], pts2b[i,])
## }))
## dat <- 1:150
## class(dat) <- c("POSIXct","POSIXt")
## x <- as.ltraj(pts, date=dat, id = rep("A", 150))
## ## Now, we suppose that there is a precision of 0.05
## ## on the relocations
## sig2 <- 0.05
## ## and that sig1=0.1
## sig1 <- 0.1
## ## Now fits the brownian bridge home range
## kbb <- kernelbb(x, sig1 = sig1,
##                  sig2 = sig2, grid=60)
## ## Now fits the classical kernel home range
## coordinates(pts) <- c("x","y")
## kud <- kernelUD(pts)
## ###### The results
## opar <- par(mfrow=c(2,2), mar=c(0.1,0.1,2,0.1))
## plot(pts, pch=16)
## title(main="The relocation pattern")
## box()
## plot(x, axes=FALSE, main="The trajectory")
## box()
## image(kud)
## title(main="Classical kernel home range")
## plot(getverticeshr(kud, 95), add=TRUE)
## box()
## image(kbb)
## title(main="Brownian bridge kernel home range")
## plot(getverticeshr(kbb, 95), add=TRUE)
## box()


###################################################
### code chunk number 101: adehabitatHR.Rnw:1482-1485
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
suppressWarnings(RNGversion("3.5.0"))
set.seed(2098)
pts1 <- data.frame(x = rnorm(25, mean = 4.5, sd = 0.05),
                   y = rnorm(25, mean = 4.5, sd = 0.05))
pts1b <- data.frame(x = rnorm(25, mean = 4.5, sd = 0.05),
                    y = rnorm(25, mean = 4.5, sd = 0.05))
pts2 <- data.frame(x = rnorm(25, mean = 4, sd = 0.05),
                   y = rnorm(25, mean = 4, sd = 0.05))
pts3 <- data.frame(x = rnorm(25, mean = 5, sd = 0.05),
                   y = rnorm(25, mean = 4, sd = 0.05))
pts3b <- data.frame(x = rnorm(25, mean = 5, sd = 0.05),
                    y = rnorm(25, mean = 4, sd = 0.05))
pts2b <- data.frame(x = rnorm(25, mean = 4, sd = 0.05),
                    y = rnorm(25, mean = 4, sd = 0.05))
pts <- do.call("rbind", lapply(1:25, function(i) {
    rbind(pts1[i,], pts1b[i,], pts2[i,], pts3[i,],
          pts3b[i,], pts2b[i,])
}))
dat <- 1:150
class(dat) <- c("POSIXct","POSIXt")
x <- as.ltraj(pts, date=dat, id = rep("A", 150))
## Now, we suppose that there is a precision of 0.05
## on the relocations
sig2 <- 0.05
## and that sig1=0.1
sig1 <- 0.1
## Now fits the brownian bridge home range
kbb <- kernelbb(x, sig1 = sig1,
                 sig2 = sig2, grid=60)
## Now fits the classical kernel home range
coordinates(pts) <- c("x","y")
kud <- kernelUD(pts)
###### The results
opar <- par(mfrow=c(2,2), mar=c(0.1,0.1,2,0.1))
plot(pts, pch=16)
title(main="The relocation pattern")
box()
plot(x, axes=FALSE, main="The trajectory")
box()
image(kud)
title(main="Classical kernel home range")
plot(getverticeshr(kud, 95), add=TRUE)
box()
image(kbb)
title(main="Brownian bridge kernel home range")
plot(getverticeshr(kbb, 95), add=TRUE)
box()
dev.null <- dev.off()
cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 102: adehabitatHR.Rnw:1489-1490
###################################################
pt <- 25


###################################################
### code chunk number 103: adehabitatHR.Rnw:1532-1534
###################################################
wi <- 800
pt <- 18


###################################################
### code chunk number 104: ploth (eval = FALSE)
###################################################
## xx <- c(0,1)
## yy <- c(0,1)
## date <- c(0,1)
## class(date) <- c("POSIXt", "POSIXct")
## tr <- as.ltraj(data.frame(x = xx,y = yy), date, id="a")
## 
## ## Use of different smoothing parameters
## sig1 <- c(0.05, 0.1, 0.2, 0.4, 0.6)
## sig2 <- c(0.05, 0.1, 0.2, 0.5, 0.7)
## 
## y <- list()
## for (i in 1:5) {
##     for (j in 1:5) {
##         k <- paste("s1=", sig1[i], ", s2=", sig2[j], sep = "")
##         y[[k]]<-kernelbb(tr, sig1[i], sig2[j])
##     }
## }
## 
## ## Displays the results
## opar <- par(mar = c(0,0,2,0), mfrow = c(5,5))
## foo <- function(x)
## {
##     image(y[[x]])
##     title(main = names(y)[x])
##     points(tr[[1]][,c("x","y")], pch = 16)
## }
## tmp <- lapply(1:length(y), foo)
## 


###################################################
### code chunk number 105: adehabitatHR.Rnw:1570-1573
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
xx <- c(0,1)
yy <- c(0,1)
date <- c(0,1)
class(date) <- c("POSIXt", "POSIXct")
tr <- as.ltraj(data.frame(x = xx,y = yy), date, id="a")

## Use of different smoothing parameters
sig1 <- c(0.05, 0.1, 0.2, 0.4, 0.6)
sig2 <- c(0.05, 0.1, 0.2, 0.5, 0.7)

y <- list()
for (i in 1:5) {
    for (j in 1:5) {
        k <- paste("s1=", sig1[i], ", s2=", sig2[j], sep = "")
        y[[k]]<-kernelbb(tr, sig1[i], sig2[j])
    }
}

## Displays the results
opar <- par(mar = c(0,0,2,0), mfrow = c(5,5))
foo <- function(x)
{
    image(y[[x]])
    title(main = names(y)[x])
    points(tr[[1]][,c("x","y")], pch = 16)
}
tmp <- lapply(1:length(y), foo)

dev.null <- dev.off()
cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 106: adehabitatHR.Rnw:1578-1580
###################################################
wi <- 600
pt <- 25


###################################################
### code chunk number 107: adehabitatHR.Rnw:1597-1600
###################################################
data(puechcirc)
x <- puechcirc[1]
x


###################################################
### code chunk number 108: plotpuechcirc (eval = FALSE)
###################################################
## plot(x)


###################################################
### code chunk number 109: adehabitatHR.Rnw:1609-1610 (eval = FALSE)
###################################################
## plot(x)


###################################################
### code chunk number 110: adehabitatHR.Rnw:1614-1617
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(x)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 111: liker1 (eval = FALSE)
###################################################
## lik <- liker(x, sig2 = 58, rangesig1 = c(10, 100))


###################################################
### code chunk number 112: adehabitatHR.Rnw:1634-1635 (eval = FALSE)
###################################################
## lik <- liker(x, sig2 = 58, rangesig1 = c(10, 100))


###################################################
### code chunk number 113: adehabitatHR.Rnw:1639-1642
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
lik <- liker(x, sig2 = 58, rangesig1 = c(10, 100))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 114: liker2 (eval = FALSE)
###################################################
## lik2 <- liker(x, sig2 = 58, rangesig1 = c(1, 10))


###################################################
### code chunk number 115: adehabitatHR.Rnw:1653-1654 (eval = FALSE)
###################################################
## lik2 <- liker(x, sig2 = 58, rangesig1 = c(1, 10))


###################################################
### code chunk number 116: adehabitatHR.Rnw:1658-1661
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
lik2 <- liker(x, sig2 = 58, rangesig1 = c(1, 10))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 117: adehabitatHR.Rnw:1667-1668
###################################################
lik2


###################################################
### code chunk number 118: adehabitatHR.Rnw:1674-1676
###################################################
tata <- kernelbb(x, sig1 = 6.23, sig2 = 58, grid = 50)
tata


###################################################
### code chunk number 119: UDkbb (eval = FALSE)
###################################################
## image(tata)
## plot(getverticeshr(tata, 95), add=TRUE, lwd=2)


###################################################
### code chunk number 120: adehabitatHR.Rnw:1688-1689 (eval = FALSE)
###################################################
## image(tata)
## plot(getverticeshr(tata, 95), add=TRUE, lwd=2)


###################################################
### code chunk number 121: adehabitatHR.Rnw:1693-1696
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(tata)
plot(getverticeshr(tata, 95), add=TRUE, lwd=2)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 122: figexe1 (eval = FALSE)
###################################################
## par(mar=c(0,0,0,0))
## plot(c(0,1,2,3),c(0,1,0,1), pch=16, cex=4, xlim=c(-0.5, 3.5),
## ylim=c(-0.5, 1.5))
## segments(0,0,1,1)
## segments(2,0,3,1)
## points(seq(2,3,length=8), seq(0,1,length=8), pch=16, cex=2)
## arrows(1,0.5, 2, 0.5, lwd=5)


###################################################
### code chunk number 123: adehabitatHR.Rnw:1723-1726
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
par(mar=c(0,0,0,0))
plot(c(0,1,2,3),c(0,1,0,1), pch=16, cex=4, xlim=c(-0.5, 3.5),
ylim=c(-0.5, 1.5))
segments(0,0,1,1)
segments(2,0,3,1)
points(seq(2,3,length=8), seq(0,1,length=8), pch=16, cex=2)
arrows(1,0.5, 2, 0.5, lwd=5)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 124: figillustrbb (eval = FALSE)
###################################################
## x1 <- c(0,0)
## x2 <- c(1,1)
## D <- 0.01
## tt <- 0.1
## foo <- function(x)
## {
##     m <- c(tt,tt)
##     xc <- x-m
##     (1/(4*pi*D*(1-tt)))*exp(-(sum(xc*xc/D))/(4*pi*tt*(1-tt)))
## }
## xy <- as.matrix(expand.grid(seq(-0.1, 1.1, length=100), seq(-0.1, 1.1, length=100)))
## den <- sapply(1:nrow(xy), function(i) foo(xy[i,]))
## df <- data.frame(x=xy[,1], y=xy[,2], den=den)
## coordinates(df) <- c("x","y")
## gridded(df) <- TRUE
## 
## par(mfrow=c(1,2), mar=c(0.1,0.1,0.1,0.1))
## image(df)
## points(c(0,1), c(0,1), pch=16, cex=4)
## lines(c(0,1), c(0,1))
## points(tt,tt, pch=16)
## box()
## tt <- 0.5
## den <- sapply(1:nrow(xy), function(i) foo(xy[i,]))
## df <- data.frame(x=xy[,1], y=xy[,2], den=den)
## coordinates(df) <- c("x","y")
## gridded(df) <- TRUE
## image(df)
## points(c(0,1), c(0,1), pch=16, cex=4)
## lines(c(0,1), c(0,1))
## points(tt,tt, pch=16)
## box()


###################################################
### code chunk number 125: adehabitatHR.Rnw:1821-1824
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
x1 <- c(0,0)
x2 <- c(1,1)
D <- 0.01
tt <- 0.1
foo <- function(x)
{
    m <- c(tt,tt)
    xc <- x-m
    (1/(4*pi*D*(1-tt)))*exp(-(sum(xc*xc/D))/(4*pi*tt*(1-tt)))
}
xy <- as.matrix(expand.grid(seq(-0.1, 1.1, length=100), seq(-0.1, 1.1, length=100)))
den <- sapply(1:nrow(xy), function(i) foo(xy[i,]))
df <- data.frame(x=xy[,1], y=xy[,2], den=den)
coordinates(df) <- c("x","y")
gridded(df) <- TRUE

par(mfrow=c(1,2), mar=c(0.1,0.1,0.1,0.1))
image(df)
points(c(0,1), c(0,1), pch=16, cex=4)
lines(c(0,1), c(0,1))
points(tt,tt, pch=16)
box()
tt <- 0.5
den <- sapply(1:nrow(xy), function(i) foo(xy[i,]))
df <- data.frame(x=xy[,1], y=xy[,2], den=den)
coordinates(df) <- c("x","y")
gridded(df) <- TRUE
image(df)
points(c(0,1), c(0,1), pch=16, cex=4)
lines(c(0,1), c(0,1))
points(tt,tt, pch=16)
box()
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 126: adehabitatHR.Rnw:1890-1893
###################################################
data(buffalo)
tr <- buffalo$traj
tr


###################################################
### code chunk number 127: trajbuff (eval = FALSE)
###################################################
## plot(tr, spixdf=buffalo$habitat)


###################################################
### code chunk number 128: adehabitatHR.Rnw:1904-1905 (eval = FALSE)
###################################################
## plot(tr, spixdf=buffalo$habitat)


###################################################
### code chunk number 129: adehabitatHR.Rnw:1909-1912
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(tr, spixdf=buffalo$habitat)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 130: adehabitatHR.Rnw:1926-1928
###################################################
vv <- BRB.D(tr, Tmax=180*60, Lmin=50, habitat=buffalo$habitat, activity="act")
vv


###################################################
### code chunk number 131: adehabitatHR.Rnw:1943-1947
###################################################
ud <- BRB(tr, D = vv, Tmax = 180*60, tau = 300, Lmin = 50, hmin=100,
          habitat = buffalo$habitat, activity = "act", grid = 100, b=0,
          same4all=FALSE, extent=0.01)
ud


###################################################
### code chunk number 132: imagebrb (eval = FALSE)
###################################################
## image(ud)


###################################################
### code chunk number 133: adehabitatHR.Rnw:1961-1962 (eval = FALSE)
###################################################
## image(ud)


###################################################
### code chunk number 134: adehabitatHR.Rnw:1966-1969
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(ud)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 135: figidrd (eval = FALSE)
###################################################
## id <- BRB(buffalo$traj, D = 7.36, Tmax = 3*3600, Lmin = 50, type = "ID",
##           hmin=100, radius = 300, maxt = 2*3600, tau = 300,
##           grid = 100, extent=0.1)
## rd <- BRB(buffalo$traj, D = 7.36, Tmax = 3*3600, Lmin = 50, type = "RD",
##           hmin=100, radius = 300, maxt = 2*3600, tau = 300,
##           grid = 100, extent=0.1)
## ud <- BRB(buffalo$traj, D = 7.36, Tmax = 3*3600, Lmin = 50, type = "UD",
##           hmin=100, radius = 300, maxt = 2*3600, tau = 300,
##           grid = 100, extent=0.1)
## 
## 
## par(mfrow = c(2,2), mar=c(0,0,2,0))
## vid <- getvolumeUD(id)
## image(vid)
## contour(vid, add=TRUE, nlevels=1, levels=30,
##         lwd=3, drawlabels=FALSE)
## title("ID")
## 
## vrd <- getvolumeUD(rd)
## image(vrd)
## contour(vrd, add=TRUE, nlevels=1, levels=30,
##         lwd=3, drawlabels=FALSE)
## title("RD")
## 
## vud <- getvolumeUD(ud)
## image(vud)
## contour(vud, add=TRUE, nlevels=1, levels=95,
##         lwd=3, drawlabels=FALSE)
## title("UD")


###################################################
### code chunk number 136: adehabitatHR.Rnw:2069-2070 (eval = FALSE)
###################################################
## id <- BRB(buffalo$traj, D = 7.36, Tmax = 3*3600, Lmin = 50, type = "ID",
##           hmin=100, radius = 300, maxt = 2*3600, tau = 300,
##           grid = 100, extent=0.1)
## rd <- BRB(buffalo$traj, D = 7.36, Tmax = 3*3600, Lmin = 50, type = "RD",
##           hmin=100, radius = 300, maxt = 2*3600, tau = 300,
##           grid = 100, extent=0.1)
## ud <- BRB(buffalo$traj, D = 7.36, Tmax = 3*3600, Lmin = 50, type = "UD",
##           hmin=100, radius = 300, maxt = 2*3600, tau = 300,
##           grid = 100, extent=0.1)
## 
## 
## par(mfrow = c(2,2), mar=c(0,0,2,0))
## vid <- getvolumeUD(id)
## image(vid)
## contour(vid, add=TRUE, nlevels=1, levels=30,
##         lwd=3, drawlabels=FALSE)
## title("ID")
## 
## vrd <- getvolumeUD(rd)
## image(vrd)
## contour(vrd, add=TRUE, nlevels=1, levels=30,
##         lwd=3, drawlabels=FALSE)
## title("RD")
## 
## vud <- getvolumeUD(ud)
## image(vud)
## contour(vud, add=TRUE, nlevels=1, levels=95,
##         lwd=3, drawlabels=FALSE)
## title("UD")


###################################################
### code chunk number 137: adehabitatHR.Rnw:2074-2077
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
id <- BRB(buffalo$traj, D = 7.36, Tmax = 3*3600, Lmin = 50, type = "ID",
          hmin=100, radius = 300, maxt = 2*3600, tau = 300,
          grid = 100, extent=0.1)
rd <- BRB(buffalo$traj, D = 7.36, Tmax = 3*3600, Lmin = 50, type = "RD",
          hmin=100, radius = 300, maxt = 2*3600, tau = 300,
          grid = 100, extent=0.1)
ud <- BRB(buffalo$traj, D = 7.36, Tmax = 3*3600, Lmin = 50, type = "UD",
          hmin=100, radius = 300, maxt = 2*3600, tau = 300,
          grid = 100, extent=0.1)


par(mfrow = c(2,2), mar=c(0,0,2,0))
vid <- getvolumeUD(id)
image(vid)
contour(vid, add=TRUE, nlevels=1, levels=30,
        lwd=3, drawlabels=FALSE)
title("ID")

vrd <- getvolumeUD(rd)
image(vrd)
contour(vrd, add=TRUE, nlevels=1, levels=30,
        lwd=3, drawlabels=FALSE)
title("RD")

vud <- getvolumeUD(ud)
image(vud)
contour(vud, add=TRUE, nlevels=1, levels=95,
        lwd=3, drawlabels=FALSE)
title("UD")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 138: adehabitatHR.Rnw:2155-2157
###################################################
data(bear)
bear


###################################################
### code chunk number 139: adehabitatHR.Rnw:2166-2169
###################################################
uu <- kernelkc(bear, h = c(1000,1000,72*3600),
               tcalc= as.POSIXct("2004-05-01"), grid=50)
uu


###################################################
### code chunk number 140: plotuu (eval = FALSE)
###################################################
## plot(bear, spixdf=uu)


###################################################
### code chunk number 141: adehabitatHR.Rnw:2180-2181 (eval = FALSE)
###################################################
## plot(bear, spixdf=uu)


###################################################
### code chunk number 142: adehabitatHR.Rnw:2185-2188
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(bear, spixdf=uu)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 143: adehabitatHR.Rnw:2196-2197
###################################################
sum(slot(uu,"data")[,1])*(gridparameters(uu)[1,2]^2)


###################################################
### code chunk number 144: plotver (eval = FALSE)
###################################################
## ver <- getverticeshr(uu, percent=95, standardize=TRUE)
## image(uu)
## plot(ver, add=TRUE)


###################################################
### code chunk number 145: adehabitatHR.Rnw:2221-2222 (eval = FALSE)
###################################################
## ver <- getverticeshr(uu, percent=95, standardize=TRUE)
## image(uu)
## plot(ver, add=TRUE)


###################################################
### code chunk number 146: adehabitatHR.Rnw:2226-2229
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
ver <- getverticeshr(uu, percent=95, standardize=TRUE)
image(uu)
plot(ver, add=TRUE)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 147: adehabitatHR.Rnw:2241-2273 (eval = FALSE)
###################################################
## ## compute the sequence of dates at which the UD is to be
## ## estimated
## vv <- seq(min(bear[[1]]$date), max(bear[[1]]$date), length=50)
## head(vv)
## 
## ## estimates the UD at each time point
## re <- lapply(1:length(vv), function(i) {
## 
##     ## estimate the UD. We choose a smoothing parameter of
##     ## 1000 meters for X and Y coordinates, and of 72 hours
##     ## for the time (after a visual exploration)
##     uu <- kernelkc(bear, h = c(1000,1000,72*3600),
##                    tcalc= vv[i], grid=50)
## 
##     ## To save the result in a file, type
##     jpeg(paste("UD", i, ".jpg", sep=""))
## 
## 
##     ## draw the image
##     image(uu, col=grey(seq(1,0,length=10)))
##     title(main=vv[i])
## 
##     ## highlight the 95 percent home range
##     ## we set standardize = TRUE because we want to estimate
##     ## the home range in space from a UD estimated in space and
##     ## time
##     plot(getverticeshr(uu, 95, standardize=TRUE), lwd=2,
##          border="red", add=TRUE)
## 
##     ## close the device
##     dev.off()
## })


###################################################
### code chunk number 148: adehabitatHR.Rnw:2297-2298
###################################################
t0 <- as.POSIXct("2012-12-25 00:00")


###################################################
### code chunk number 149: exwc (eval = FALSE)
###################################################
## exwc(0.2)


###################################################
### code chunk number 150: adehabitatHR.Rnw:2309-2310 (eval = FALSE)
###################################################
## exwc(0.2)


###################################################
### code chunk number 151: adehabitatHR.Rnw:2314-2317
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
exwc(0.2)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 152: imagekernelkc1 (eval = FALSE)
###################################################
## uu <- kernelkc(bear, h=c(1000,1000,0.2), cycle=24*3600,
##                tcalc=3*3600, t0=t0, circular=TRUE)
## image(uu)


###################################################
### code chunk number 153: adehabitatHR.Rnw:2333-2334 (eval = FALSE)
###################################################
## uu <- kernelkc(bear, h=c(1000,1000,0.2), cycle=24*3600,
##                tcalc=3*3600, t0=t0, circular=TRUE)
## image(uu)


###################################################
### code chunk number 154: adehabitatHR.Rnw:2338-2341
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
uu <- kernelkc(bear, h=c(1000,1000,0.2), cycle=24*3600,
               tcalc=3*3600, t0=t0, circular=TRUE)
image(uu)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 155: adehabitatHR.Rnw:2403-2405
###################################################
uu <- clusthr(puechabonsp$relocs[,1])
class(uu)


###################################################
### code chunk number 156: kkkssssm (eval = FALSE)
###################################################
## plot(uu)


###################################################
### code chunk number 157: adehabitatHR.Rnw:2418-2419 (eval = FALSE)
###################################################
## plot(uu)


###################################################
### code chunk number 158: adehabitatHR.Rnw:2423-2426
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(uu)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 159: adehabitatHR.Rnw:2433-2434
###################################################
as.data.frame(uu[[3]])


###################################################
### code chunk number 160: ksqkskq (eval = FALSE)
###################################################
## plot(uu[[3]][26,])


###################################################
### code chunk number 161: adehabitatHR.Rnw:2446-2447 (eval = FALSE)
###################################################
## plot(uu[[3]][26,])


###################################################
### code chunk number 162: adehabitatHR.Rnw:2451-2454
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(uu[[3]][26,])
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 163: sldlsl (eval = FALSE)
###################################################
## uur <- MCHu.rast(uu, puechabonsp$map, percent=90)
## image(uur, 3)


###################################################
### code chunk number 164: adehabitatHR.Rnw:2471-2472 (eval = FALSE)
###################################################
## uur <- MCHu.rast(uu, puechabonsp$map, percent=90)
## image(uur, 3)


###################################################
### code chunk number 165: adehabitatHR.Rnw:2476-2479
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
uur <- MCHu.rast(uu, puechabonsp$map, percent=90)
image(uur, 3)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 166: adehabitatHR.Rnw:2490-2491
###################################################
head(as.data.frame(uu[[1]]))


###################################################
### code chunk number 167: slqllqq (eval = FALSE)
###################################################
## ii <- MCHu2hrsize(uu, percent=seq(50, 100, by=5))


###################################################
### code chunk number 168: adehabitatHR.Rnw:2505-2506 (eval = FALSE)
###################################################
## ii <- MCHu2hrsize(uu, percent=seq(50, 100, by=5))


###################################################
### code chunk number 169: adehabitatHR.Rnw:2510-2513
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
ii <- MCHu2hrsize(uu, percent=seq(50, 100, by=5))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 170: adehabitatHR.Rnw:2568-2569 (eval = FALSE)
###################################################
## res <- LoCoH.k(puechabonsp$relocs[,1], 10)


###################################################
### code chunk number 171: adehabitatHR.Rnw:2579-2580 (eval = FALSE)
###################################################
## plot(res)


###################################################
### code chunk number 172: adehabitatHR.Rnw:2587-2589 (eval = FALSE)
###################################################
## LoCoH.k.area(puechabonsp$relocs[,1], krange=seq(4, 15, length=5),
##              percent=100)


###################################################
### code chunk number 173: adehabitatHR.Rnw:2611-2613
###################################################
res <- CharHull(puechabonsp$relocs[,1])
class(res)


###################################################
### code chunk number 174: plotcharhul (eval = FALSE)
###################################################
## plot(res)


###################################################
### code chunk number 175: adehabitatHR.Rnw:2629-2630 (eval = FALSE)
###################################################
## plot(res)


###################################################
### code chunk number 176: adehabitatHR.Rnw:2634-2637
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(res)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


