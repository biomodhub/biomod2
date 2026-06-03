### R code from vignette source 'adehabitatLT.Rnw'

###################################################
### code chunk number 1: adehabitatLT.Rnw:29-36
###################################################
owidth <- getOption("width")
options("width"=80)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0
wi <- 600
pt <- 15


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
## cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 4: zfigg (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 5: adehabitatLT.Rnw:105-106
###################################################
library(adehabitatLT)


###################################################
### code chunk number 6: adehabitatLT.Rnw:109-111
###################################################
suppressWarnings(RNGversion("3.5.0"))
set.seed(13431)


###################################################
### code chunk number 7: ploltraj (eval = FALSE)
###################################################
## par(mar=c(0.1,0.1,0.1,0.1))
## plot(c(0,1), c(0,1), ty="n", axes=FALSE)
## x <- c(0.2, 0.3, 0.5, 0.7, 0.4)
## y <- c(0.5, 0.3, 0.35, 0.6, 0.9)
## points(x,y, pch=16, cex=2)
## points(x[1],y[1], pch=16, col="green")
## lines(x,y, lty=1)
## arrows(0.4, 0.9, 0.1, 0.8, lty=2, length=0.1)
## lines(c(0.5, 0.7),c(0.35, 0.6), lwd=6)
## lines(c(0.5, 0.7),c(0.35, 0.6), lwd=2, col="red")
## 
## ## dx
## arrows(0.5, 0.32, 0.7, 0.32, code=3, length=0.15)
## text(0.6, 0.3, "dx")
## 
## ## dy
## arrows(0.75, 0.35, 0.75, 0.6, code=3, length=0.15)
## text(0.77, 0.5, "dy")
## 
## ## abs.angle
## ang <- atan2(0.25, 0.2)
## ang <- seq(0,ang, length=10)
## xa <- 0.05*cos(ang) + 0.5
## ya <- 0.05*sin(ang) + 0.35
## lines(c(0.5, 0.7),c(0.35, 0.35), lty=3, col="red")
## lines(xa, ya, col="red", lwd=2)
## text(0.56, 0.38, expression(alpha), col="red")
## 
## ## rel.angle
## lines(c(0.3, 0.5),c(0.3, 0.35), col="blue", lwd=2)
## lines(c(0.5, 0.7),c(0.35, 0.4), lty=3, col="blue")
## 
## ang1 <- atan2(0.25, 0.2)
## ang0 <- atan2(0.05, 0.2)
## ang <- ang0+seq(0, ang1-ang0, length=10)
## xa <- 0.1*cos(ang) + 0.5
## ya <- 0.1*sin(ang) + 0.35
## lines(xa, ya, col="blue", lwd=2)
## xa <- 0.107*cos(ang) + 0.5
## ya <- 0.107*sin(ang) + 0.35
## lines(xa, ya, col="blue", lwd=2)
## text(0.61, 0.43,expression(beta), col="blue")
## 
## ## dist
## arrows(0.48, 0.37, 0.68, 0.62, code=3, length=0.1)
## text(0.53, 0.5, "dist, dt")
## 
## ## R2n
## arrows(0.21, 0.5, 0.49, 0.35, col="darkgreen", code=3, length=0.1)
## text(0.35, 0.49, expression(R[n]^2))
## 
## ## The legend:
## text(0.2, 0.2, expression(paste(alpha, ": abs.angle")), col="red")
## text(0.2, 0.17, expression(paste(beta, ": rel.angle")), col="blue")
## 
## ## x0, y0, t0
## text(0.14, 0.5, expression(paste(x[0], ", ", y[0], ", ", t[0])),
##      col="darkgreen")
## 
## box()


###################################################
### code chunk number 8: adehabitatLT.Rnw:239-242
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
par(mar=c(0.1,0.1,0.1,0.1))
plot(c(0,1), c(0,1), ty="n", axes=FALSE)
x <- c(0.2, 0.3, 0.5, 0.7, 0.4)
y <- c(0.5, 0.3, 0.35, 0.6, 0.9)
points(x,y, pch=16, cex=2)
points(x[1],y[1], pch=16, col="green")
lines(x,y, lty=1)
arrows(0.4, 0.9, 0.1, 0.8, lty=2, length=0.1)
lines(c(0.5, 0.7),c(0.35, 0.6), lwd=6)
lines(c(0.5, 0.7),c(0.35, 0.6), lwd=2, col="red")

## dx
arrows(0.5, 0.32, 0.7, 0.32, code=3, length=0.15)
text(0.6, 0.3, "dx")

## dy
arrows(0.75, 0.35, 0.75, 0.6, code=3, length=0.15)
text(0.77, 0.5, "dy")

## abs.angle
ang <- atan2(0.25, 0.2)
ang <- seq(0,ang, length=10)
xa <- 0.05*cos(ang) + 0.5
ya <- 0.05*sin(ang) + 0.35
lines(c(0.5, 0.7),c(0.35, 0.35), lty=3, col="red")
lines(xa, ya, col="red", lwd=2)
text(0.56, 0.38, expression(alpha), col="red")

## rel.angle
lines(c(0.3, 0.5),c(0.3, 0.35), col="blue", lwd=2)
lines(c(0.5, 0.7),c(0.35, 0.4), lty=3, col="blue")

ang1 <- atan2(0.25, 0.2)
ang0 <- atan2(0.05, 0.2)
ang <- ang0+seq(0, ang1-ang0, length=10)
xa <- 0.1*cos(ang) + 0.5
ya <- 0.1*sin(ang) + 0.35
lines(xa, ya, col="blue", lwd=2)
xa <- 0.107*cos(ang) + 0.5
ya <- 0.107*sin(ang) + 0.35
lines(xa, ya, col="blue", lwd=2)
text(0.61, 0.43,expression(beta), col="blue")

## dist
arrows(0.48, 0.37, 0.68, 0.62, code=3, length=0.1)
text(0.53, 0.5, "dist, dt")

## R2n
arrows(0.21, 0.5, 0.49, 0.35, col="darkgreen", code=3, length=0.1)
text(0.35, 0.49, expression(R[n]^2))

## The legend:
text(0.2, 0.2, expression(paste(alpha, ": abs.angle")), col="red")
text(0.2, 0.17, expression(paste(beta, ": rel.angle")), col="blue")

## x0, y0, t0
text(0.14, 0.5, expression(paste(x[0], ", ", y[0], ", ", t[0])),
     col="darkgreen")

box()
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 9: adehabitatLT.Rnw:291-295
###################################################
data(puechabonsp)
locs <- puechabonsp$relocs
locs <- as.data.frame(locs)
head(locs)


###################################################
### code chunk number 10: adehabitatLT.Rnw:324-327
###################################################
da <- as.character(locs$Date)
head(da)
da <- as.POSIXct(strptime(as.character(locs$Date),"%y%m%d", tz="Europe/Paris"))


###################################################
### code chunk number 11: adehabitatLT.Rnw:333-335
###################################################
puech <- as.ltraj(xy = locs[,c("X","Y")], date = da, id = locs$Name)
puech


###################################################
### code chunk number 12: adehabitatLT.Rnw:346-347
###################################################
head(puech[[1]])


###################################################
### code chunk number 13: plottrajjj (eval = FALSE)
###################################################
## plot(puech)


###################################################
### code chunk number 14: adehabitatLT.Rnw:362-363 (eval = FALSE)
###################################################
## plot(puech)


###################################################
### code chunk number 15: adehabitatLT.Rnw:367-370
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(puech)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 16: adehabitatLT.Rnw:403-404
###################################################
head(puech[[1]])


###################################################
### code chunk number 17: adehabitatLT.Rnw:410-413
###################################################
puech2 <- puech
puech2[[1]][2,1] <- 700146
head(puech2[[1]])


###################################################
### code chunk number 18: adehabitatLT.Rnw:420-421
###################################################
head(rec(puech2)[[1]])


###################################################
### code chunk number 19: adehabitatLT.Rnw:439-441
###################################################
puech2 <- ld(puech)
head(puech2)


###################################################
### code chunk number 20: adehabitatLT.Rnw:449-450
###################################################
dl(puech2)


###################################################
### code chunk number 21: adehabitatLT.Rnw:466-467
###################################################
is.regular(puech)


###################################################
### code chunk number 22: plotdtltr (eval = FALSE)
###################################################
## plotltr(puech, "dt/3600/24")


###################################################
### code chunk number 23: adehabitatLT.Rnw:485-486 (eval = FALSE)
###################################################
## plotltr(puech, "dt/3600/24")


###################################################
### code chunk number 24: adehabitatLT.Rnw:490-493
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotltr(puech, "dt/3600/24")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 25: adehabitatLT.Rnw:504-507
###################################################
foo <- function(dt) {
    return(dt> (100*3600*24))
}


###################################################
### code chunk number 26: adehabitatLT.Rnw:515-517
###################################################
puech2 <- cutltraj(puech, "foo(dt)", nextr = TRUE)
puech2


###################################################
### code chunk number 27: adehabitatLT.Rnw:525-527
###################################################
burst(puech2)[3:4] <- c("Chou.1992", "Chou.1993")
puech2


###################################################
### code chunk number 28: adehabitatLT.Rnw:539-540
###################################################
puech2


###################################################
### code chunk number 29: adehabitatLT.Rnw:547-549
###################################################
puech2b <- puech2[c(1,2,5)]
puech2b


###################################################
### code chunk number 30: adehabitatLT.Rnw:555-557
###################################################
puech2c <- c(puech2b, puech2[4])
puech2c


###################################################
### code chunk number 31: adehabitatLT.Rnw:569-571
###################################################
bu <- which.ltraj(puech2, "dist>2000")
bu


###################################################
### code chunk number 32: adehabitatLT.Rnw:578-579
###################################################
puech2[burst(puech2)%in%bu$burst]


###################################################
### code chunk number 33: plotltr2bu (eval = FALSE)
###################################################
## plotltr(puech2, "dt/3600/24")


###################################################
### code chunk number 34: adehabitatLT.Rnw:592-593 (eval = FALSE)
###################################################
## plotltr(puech2, "dt/3600/24")


###################################################
### code chunk number 35: adehabitatLT.Rnw:597-600
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotltr(puech2, "dt/3600/24")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 36: adehabitatLT.Rnw:610-612
###################################################
refda <- strptime("00:00", "%H:%M", tz="Europe/Paris")
refda


###################################################
### code chunk number 37: adehabitatLT.Rnw:621-623
###################################################
puech3 <- setNA(puech2, refda, 1, units = "day")
puech3


###################################################
### code chunk number 38: adehabitatLT.Rnw:640-642
###################################################
data(ibexraw)
ibexraw


###################################################
### code chunk number 39: plotibex (eval = FALSE)
###################################################
## plotltr(ibexraw, "dt/3600")


###################################################
### code chunk number 40: adehabitatLT.Rnw:652-653 (eval = FALSE)
###################################################
## plotltr(ibexraw, "dt/3600")


###################################################
### code chunk number 41: adehabitatLT.Rnw:657-660
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotltr(ibexraw, "dt/3600")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 42: adehabitatLT.Rnw:669-672
###################################################
refda <- strptime("2003-06-01 00:00", "%Y-%m-%d %H:%M", tz="Europe/Paris")
ib2 <- setNA(ibexraw, refda, 4, units = "hour")
ib2


###################################################
### code chunk number 43: plotltrib2 (eval = FALSE)
###################################################
## plotltr(ib2, "dt/3600")


###################################################
### code chunk number 44: adehabitatLT.Rnw:683-684 (eval = FALSE)
###################################################
## plotltr(ib2, "dt/3600")


###################################################
### code chunk number 45: adehabitatLT.Rnw:688-691
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotltr(ib2, "dt/3600")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 46: adehabitatLT.Rnw:699-701
###################################################
ib3 <- sett0(ib2, refda, 4, units = "hour")
ib3


###################################################
### code chunk number 47: adehabitatLT.Rnw:724-725
###################################################
is.sd(ib3)


###################################################
### code chunk number 48: adehabitatLT.Rnw:733-734
###################################################
ib3


###################################################
### code chunk number 49: adehabitatLT.Rnw:744-748
###################################################
ib4 <- set.limits(ib3, begin = "2003-06-01 00:00",
                  dur = 14, units = "day", pattern = "%Y-%m-%d %H:%M",
                  tz="Europe/Paris")
ib4


###################################################
### code chunk number 50: adehabitatLT.Rnw:753-754
###################################################
is.sd(ib4)


###################################################
### code chunk number 51: adehabitatLT.Rnw:777-779
###################################################
di <- sd2df(ib4, "dist")
head(di)


###################################################
### code chunk number 52: adehabitatLT.Rnw:802-804
###################################################
data(capreochiz)
head(capreochiz)


###################################################
### code chunk number 53: adehabitatLT.Rnw:816-820
###################################################
capreo <- as.ltraj(xy = capreochiz[,c("x","y")], date = capreochiz$date,
                   id = "Roe.Deer",
                   infolocs = capreochiz[,4:8])
capreo


###################################################
### code chunk number 54: adehabitatLT.Rnw:827-829
###################################################
inf <- infolocs(capreo)
head(inf[[1]])


###################################################
### code chunk number 55: plotdop (eval = FALSE)
###################################################
## plotltr(capreo, "log(Dop)")


###################################################
### code chunk number 56: adehabitatLT.Rnw:857-858 (eval = FALSE)
###################################################
## plotltr(capreo, "log(Dop)")


###################################################
### code chunk number 57: adehabitatLT.Rnw:862-865
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotltr(capreo, "log(Dop)")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 58: NAtest (eval = FALSE)
###################################################
## runsNAltraj(ib4)


###################################################
### code chunk number 59: adehabitatLT.Rnw:898-899 (eval = FALSE)
###################################################
## runsNAltraj(ib4)


###################################################
### code chunk number 60: adehabitatLT.Rnw:904-907
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
runsNAltraj(ib4)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 61: adehabitatLT.Rnw:933-935
###################################################
data(bear)
bear


###################################################
### code chunk number 62: runsNAbear (eval = FALSE)
###################################################
## runsNAltraj(bear)


###################################################
### code chunk number 63: adehabitatLT.Rnw:947-948 (eval = FALSE)
###################################################
## runsNAltraj(bear)


###################################################
### code chunk number 64: adehabitatLT.Rnw:952-955
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
runsNAltraj(bear)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 65: plotNAbear (eval = FALSE)
###################################################
## plotNAltraj(bear)


###################################################
### code chunk number 66: adehabitatLT.Rnw:966-967 (eval = FALSE)
###################################################
## plotNAltraj(bear)


###################################################
### code chunk number 67: adehabitatLT.Rnw:971-974
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotNAltraj(bear)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 68: adehabitatLT.Rnw:983-985
###################################################
missval <- as.numeric(is.na(bear[[1]]$x))
head(missval)


###################################################
### code chunk number 69: adehabitatLT.Rnw:1011-1013
###################################################
bearI <- typeII2typeI(bear)
bearI


###################################################
### code chunk number 70: trajetours (eval = FALSE)
###################################################
## plot(bearI)


###################################################
### code chunk number 71: adehabitatLT.Rnw:1035-1036 (eval = FALSE)
###################################################
## plot(bearI)


###################################################
### code chunk number 72: adehabitatLT.Rnw:1040-1043
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(bearI)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 73: adehabitatLT.Rnw:1050-1052
###################################################
bearIr <- redisltraj(bearI, 500)
bearIr


###################################################
### code chunk number 74: plotltrajredisbear (eval = FALSE)
###################################################
## plot(bearIr)


###################################################
### code chunk number 75: adehabitatLT.Rnw:1062-1063 (eval = FALSE)
###################################################
## plot(bearIr)


###################################################
### code chunk number 76: adehabitatLT.Rnw:1067-1070
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(bearIr)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 77: sliwinltr (eval = FALSE)
###################################################
## sliwinltr(bearIr, function(x) mean(cos(x$rel.angle)), type="locs", step=30)


###################################################
### code chunk number 78: adehabitatLT.Rnw:1083-1084 (eval = FALSE)
###################################################
## sliwinltr(bearIr, function(x) mean(cos(x$rel.angle)), type="locs", step=30)


###################################################
### code chunk number 79: adehabitatLT.Rnw:1088-1091
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
sliwinltr(bearIr, function(x) mean(cos(x$rel.angle)), type="locs", step=30)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 80: adehabitatLT.Rnw:1101-1103
###################################################
cosrelangle <- cos(bearIr[[1]]$rel.angle)
head(cosrelangle)


###################################################
### code chunk number 81: adehabitatLT.Rnw:1133-1135
###################################################
data(porpoise)
porpoise


###################################################
### code chunk number 82: adehabitatLT.Rnw:1145-1146
###################################################
(porpoise2 <- redisltraj(na.omit(porpoise[1:3]), 86400, type="time"))


###################################################
### code chunk number 83: adehabitatLT.Rnw:1163-1164 (eval = FALSE)
###################################################
## trajdyn(ib4)


###################################################
### code chunk number 84: adehabitatLT.Rnw:1202-1203
###################################################
wawotest(bear)


###################################################
### code chunk number 85: plotltrbeardist (eval = FALSE)
###################################################
## plotltr(bear, "dist")


###################################################
### code chunk number 86: adehabitatLT.Rnw:1218-1219 (eval = FALSE)
###################################################
## plotltr(bear, "dist")


###################################################
### code chunk number 87: adehabitatLT.Rnw:1223-1226
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotltr(bear, "dist")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 88: adehabitatLT.Rnw:1233-1235 (eval = FALSE)
###################################################
## sliwinltr(bear, function(x) mean(na.omit(x$dist)),
##           5*48, type="locs")


###################################################
### code chunk number 89: acfdistltra (eval = FALSE)
###################################################
## acfdist.ltraj(bear, lag=5, which="dist")


###################################################
### code chunk number 90: adehabitatLT.Rnw:1263-1264 (eval = FALSE)
###################################################
## acfdist.ltraj(bear, lag=5, which="dist")


###################################################
### code chunk number 91: adehabitatLT.Rnw:1268-1271
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
acfdist.ltraj(bear, lag=5, which="dist")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 92: critereangle (eval = FALSE)
###################################################
## opar <- par(mar = c(0,0,4,0))
## plot(0,0, asp=1, xlim=c(-1, 1), ylim=c(-1, 1), ty="n", axes=FALSE,
##      main="Criteria f for the measure of independence between successive
##      angles at time i-1 and i")
## box()
## symbols(0,0,circle=1, inches=FALSE, lwd=2, add=TRUE)
## abline(h=0, v=0)
## x <- c( cos(pi/3), cos(pi/2 + pi/4))
## y <- c( sin(pi/3), sin(pi/2 + pi/4))
## arrows(c(0,0), c(0,0), x, y)
## lines(x,y, lwd=2, col="red")
## text(0, 0.9, expression(f^2 == 2*sum((1 - cos(alpha[i]-alpha[i-1])),
##     i==1, n-1)), col="red")
## foo <- function(t, alpha)
## {
##     xa <- sapply(seq(0, alpha, length=20), function(x) t*cos(x))
##     ya <- sapply(seq(0, alpha, length=20), function(x) t*sin(x))
##     lines(xa, ya)
## }
## foo(0.3, pi/3)
## foo(0.1, pi/2 + pi/4)
## foo(0.11, pi/2 + pi/4)
## text(0.34,0.18,expression(alpha[i]), cex=1.5)
## text(0.15,0.11,expression(alpha[i-1]), cex=1.5)
## par(opar)


###################################################
### code chunk number 93: adehabitatLT.Rnw:1317-1320
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
opar <- par(mar = c(0,0,4,0))
plot(0,0, asp=1, xlim=c(-1, 1), ylim=c(-1, 1), ty="n", axes=FALSE,
     main="Criteria f for the measure of independence between successive
     angles at time i-1 and i")
box()
symbols(0,0,circle=1, inches=FALSE, lwd=2, add=TRUE)
abline(h=0, v=0)
x <- c( cos(pi/3), cos(pi/2 + pi/4))
y <- c( sin(pi/3), sin(pi/2 + pi/4))
arrows(c(0,0), c(0,0), x, y)
lines(x,y, lwd=2, col="red")
text(0, 0.9, expression(f^2 == 2*sum((1 - cos(alpha[i]-alpha[i-1])),
    i==1, n-1)), col="red")
foo <- function(t, alpha)
{
    xa <- sapply(seq(0, alpha, length=20), function(x) t*cos(x))
    ya <- sapply(seq(0, alpha, length=20), function(x) t*sin(x))
    lines(xa, ya)
}
foo(0.3, pi/3)
foo(0.1, pi/2 + pi/4)
foo(0.11, pi/2 + pi/4)
text(0.34,0.18,expression(alpha[i]), cex=1.5)
text(0.15,0.11,expression(alpha[i-1]), cex=1.5)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 94: adehabitatLT.Rnw:1329-1330
###################################################
testang.ltraj(bear, "relative")


###################################################
### code chunk number 95: acfangbear (eval = FALSE)
###################################################
## acfang.ltraj(bear, lag=5)


###################################################
### code chunk number 96: adehabitatLT.Rnw:1344-1345 (eval = FALSE)
###################################################
## acfang.ltraj(bear, lag=5)


###################################################
### code chunk number 97: adehabitatLT.Rnw:1349-1352
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
acfang.ltraj(bear, lag=5)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 98: adehabitatLT.Rnw:1379-1382
###################################################
data(porpoise)
gus <- porpoise[1]
gus


###################################################
### code chunk number 99: pltgus (eval = FALSE)
###################################################
## plot(gus)


###################################################
### code chunk number 100: adehabitatLT.Rnw:1392-1393 (eval = FALSE)
###################################################
## plot(gus)


###################################################
### code chunk number 101: adehabitatLT.Rnw:1397-1400
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(gus)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 102: acfgus (eval = FALSE)
###################################################
## acfdist.ltraj(gus, "dist", lag=20)


###################################################
### code chunk number 103: adehabitatLT.Rnw:1417-1418 (eval = FALSE)
###################################################
## acfdist.ltraj(gus, "dist", lag=20)


###################################################
### code chunk number 104: adehabitatLT.Rnw:1422-1425
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
acfdist.ltraj(gus, "dist", lag=20)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 105: plodisgus (eval = FALSE)
###################################################
## plotltr(gus, "dist")


###################################################
### code chunk number 106: adehabitatLT.Rnw:1437-1438 (eval = FALSE)
###################################################
## plotltr(gus, "dist")


###################################################
### code chunk number 107: adehabitatLT.Rnw:1442-1445
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotltr(gus, "dist")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 108: adehabitatLT.Rnw:1455-1456
###################################################
(tested.means <- round(seq(0, 130000, length = 10), 0))


###################################################
### code chunk number 109: adehabitatLT.Rnw:1464-1467
###################################################
(limod <- as.list(paste("dnorm(dist, mean =",
                        tested.means,
                        ", sd = 5000)")))


###################################################
### code chunk number 110: adehabitatLT.Rnw:1483-1485
###################################################
mod <- modpartltraj(gus, limod)
mod


###################################################
### code chunk number 111: bpmgus (eval = FALSE)
###################################################
## bestpartmod(mod)


###################################################
### code chunk number 112: adehabitatLT.Rnw:1497-1498 (eval = FALSE)
###################################################
## bestpartmod(mod)


###################################################
### code chunk number 113: adehabitatLT.Rnw:1502-1505
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
bestpartmod(mod)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 114: adehabitatLT.Rnw:1543-1544
###################################################
(pm <- partmod.ltraj(gus, 4, mod))


###################################################
### code chunk number 115: plotpartitiongus (eval = FALSE)
###################################################
## plot(pm)


###################################################
### code chunk number 116: adehabitatLT.Rnw:1554-1555 (eval = FALSE)
###################################################
## plot(pm)


###################################################
### code chunk number 117: adehabitatLT.Rnw:1559-1562
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(pm)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 118: plotpartltr (eval = FALSE)
###################################################
## ## Shows the segmentation on the distances:
## plotltr(gus, "dist")
## tmp <- lapply(1:length(pm$ltraj), function(i) {
##     coul <- c("red","green","blue")[as.numeric(factor(pm$stats$mod))[i]]
##     lines(pm$ltraj[[i]]$date, rep(tested.means[pm$stats$mod[i]],
##                                   nrow(pm$ltraj[[i]])),
##           col=coul, lwd=2)
## })


###################################################
### code chunk number 119: adehabitatLT.Rnw:1587-1588 (eval = FALSE)
###################################################
## ## Shows the segmentation on the distances:
## plotltr(gus, "dist")
## tmp <- lapply(1:length(pm$ltraj), function(i) {
##     coul <- c("red","green","blue")[as.numeric(factor(pm$stats$mod))[i]]
##     lines(pm$ltraj[[i]]$date, rep(tested.means[pm$stats$mod[i]],
##                                   nrow(pm$ltraj[[i]])),
##           col=coul, lwd=2)
## })


###################################################
### code chunk number 120: adehabitatLT.Rnw:1592-1595
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
## Shows the segmentation on the distances:
plotltr(gus, "dist")
tmp <- lapply(1:length(pm$ltraj), function(i) {
    coul <- c("red","green","blue")[as.numeric(factor(pm$stats$mod))[i]]
    lines(pm$ltraj[[i]]$date, rep(tested.means[pm$stats$mod[i]],
                                  nrow(pm$ltraj[[i]])),
          col=coul, lwd=2)
})
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 121: plotrespart (eval = FALSE)
###################################################
## ## Computes the residuals
## res <- unlist(lapply(1:length(pm$ltraj), function(i) {
##     pm$ltraj[[i]]$dist - rep(tested.means[pm$stats$mod[i]],
##                              nrow(pm$ltraj[[i]]))
## }))
## plot(res, ty = "l")


###################################################
### code chunk number 122: adehabitatLT.Rnw:1613-1614 (eval = FALSE)
###################################################
## ## Computes the residuals
## res <- unlist(lapply(1:length(pm$ltraj), function(i) {
##     pm$ltraj[[i]]$dist - rep(tested.means[pm$stats$mod[i]],
##                              nrow(pm$ltraj[[i]]))
## }))
## plot(res, ty = "l")


###################################################
### code chunk number 123: adehabitatLT.Rnw:1618-1621
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
## Computes the residuals
res <- unlist(lapply(1:length(pm$ltraj), function(i) {
    pm$ltraj[[i]]$dist - rep(tested.means[pm$stats$mod[i]],
                             nrow(pm$ltraj[[i]]))
}))
plot(res, ty = "l")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 124: adehabitatLT.Rnw:1628-1629
###################################################
wawotest(res)


###################################################
### code chunk number 125: adehabitatLT.Rnw:1734-1757
###################################################
par(mar=c(0,0,0,0))
plot(c(0,1), xlim=c(-0.3,1.3), ylim=c(-0.3,1.3), axes=FALSE)
ii <- seq(0.05, 0.95, by=0.1)
tmp <- lapply(ii, function(i) {
    lapply(ii, function(j) {
        if (j-(1-i)> -0.01) {
            polygon(c(i-0.05, i+0.05, i+0.05, i-0.05),
                    c(j-0.05, j-0.05, j+0.05, j+0.05), col="grey")
        } else {
            polygon(c(i-0.05, i+0.05, i+0.05, i-0.05),
                    c(j-0.05, j-0.05, j+0.05, j+0.05))
        }
    })})
polygon(c(0,1,1,0), c(0,0,1,1), lwd=2)
text(rep(-0.05,10), ii, 10:1)
text(ii, rep(1.05,10), 1:10)
i <- 0.65
j <- 0.85
polygon(c(i-0.05, i+0.05, i+0.05, i-0.05),
        c(j-0.05, j-0.05, j+0.05, j+0.05), col="red")
segments(0.4, 1.18, i, j)
points(i,j, pch=16)
text(0.4, 1.25, "Contrast function measured on the segment\n beginning at observation 2 and ending at observation 7")


###################################################
### code chunk number 126: adehabitatLT.Rnw:1765-1785
###################################################
par(mar=c(0,0,0,0))
plot(c(0,1), xlim=c(-0.3,1.3), ylim=c(-0.3,1.3), axes=FALSE)
ii <- seq(0.05, 0.95, by=0.1)
tmp <- lapply(ii, function(i) {
    lapply(ii, function(j) {
        if (j-(1-i)> -0.01) {
            polygon(c(i-0.05, i+0.05, i+0.05, i-0.05),
                    c(j-0.05, j-0.05, j+0.05, j+0.05), col="grey")
        } else {
            polygon(c(i-0.05, i+0.05, i+0.05, i-0.05),
                    c(j-0.05, j-0.05, j+0.05, j+0.05))
        }
    })})
polygon(c(0,1,1,0), c(0,0,1,1), lwd=2)
text(rep(-0.05,10), ii, 10:1)
text(ii, rep(1.05,10), 1:10)
lines(c(0.00, 0.45, 0.50, 0.75, 0.80, 0.95),
      c(0.95, 0.95, 0.45, 0.45, 0.15, 0.15), col="red", lwd=3)
points(c(0.45, 0.75, 0.95),
       c(0.95, 0.45, 0.15), pch=16, cex=2)


###################################################
### code chunk number 127: figseri1 (eval = FALSE)
###################################################
## suppressWarnings(RNGversion("3.5.0"))
## set.seed(129)
## seri <- c(rnorm(100), rnorm(100, mean=2),
##           rnorm(100), rnorm(100, mean=-3),
##           rnorm(100), rnorm(100, mean=2))
## plot(seri, ty="l", xlab="time", ylab="Series")


###################################################
### code chunk number 128: adehabitatLT.Rnw:1861-1862 (eval = FALSE)
###################################################
## suppressWarnings(RNGversion("3.5.0"))
## set.seed(129)
## seri <- c(rnorm(100), rnorm(100, mean=2),
##           rnorm(100), rnorm(100, mean=-3),
##           rnorm(100), rnorm(100, mean=2))
## plot(seri, ty="l", xlab="time", ylab="Series")


###################################################
### code chunk number 129: adehabitatLT.Rnw:1866-1869
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
suppressWarnings(RNGversion("3.5.0"))
set.seed(129)
seri <- c(rnorm(100), rnorm(100, mean=2),
          rnorm(100), rnorm(100, mean=-3),
          rnorm(100), rnorm(100, mean=2))
plot(seri, ty="l", xlab="time", ylab="Series")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 130: adehabitatLT.Rnw:1878-1879
###################################################
(l <- lavielle(seri, Lmin=10, Kmax=20, type="mean"))


###################################################
### code chunk number 131: adehabitatLT.Rnw:1890-1891
###################################################
chooseseg(l)


###################################################
### code chunk number 132: figchoose1 (eval = FALSE)
###################################################
## fp <- findpath(l, 6)


###################################################
### code chunk number 133: adehabitatLT.Rnw:1910-1911 (eval = FALSE)
###################################################
## fp <- findpath(l, 6)


###################################################
### code chunk number 134: adehabitatLT.Rnw:1915-1918
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
fp <- findpath(l, 6)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 135: adehabitatLT.Rnw:1926-1927
###################################################
fp


###################################################
### code chunk number 136: adehabitatLT.Rnw:1939-1940
###################################################
plotltr(gus, "dist")


###################################################
### code chunk number 137: adehabitatLT.Rnw:1950-1952
###################################################
lav <- lavielle(gus, Lmin=2, Kmax=20)
chooseseg(lav)


###################################################
### code chunk number 138: adehabitatLT.Rnw:1959-1961
###################################################
kk <- findpath(lav, 4)
kk


###################################################
### code chunk number 139: adehabitatLT.Rnw:1968-1969
###################################################
plot(kk)


###################################################
### code chunk number 140: imagerast1 (eval = FALSE)
###################################################
## data(puechcirc)
## data(puechabonsp)
## mimage(puechabonsp$map)


###################################################
### code chunk number 141: adehabitatLT.Rnw:2020-2021 (eval = FALSE)
###################################################
## data(puechcirc)
## data(puechabonsp)
## mimage(puechabonsp$map)


###################################################
### code chunk number 142: adehabitatLT.Rnw:2025-2028
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
data(puechcirc)
data(puechabonsp)
mimage(puechabonsp$map)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 143: adehabitatLT.Rnw:2034-2035
###################################################
ii <- rasterize.ltraj(puechcirc, puechabonsp$map)


###################################################
### code chunk number 144: adehabitatLT.Rnw:2042-2044
###################################################
tr1 <- ii[[1]]
head(tr1)


###################################################
### code chunk number 145: ksdfjk (eval = FALSE)
###################################################
## plot(tr1)
## points(tr1[tr1[[1]]==3,], col="red")


###################################################
### code chunk number 146: adehabitatLT.Rnw:2056-2057 (eval = FALSE)
###################################################
## plot(tr1)
## points(tr1[tr1[[1]]==3,], col="red")


###################################################
### code chunk number 147: adehabitatLT.Rnw:2061-2064
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(tr1)
points(tr1[tr1[[1]]==3,], col="red")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 148: jklmmre (eval = FALSE)
###################################################
## ov <- over(tr1, puechabonsp$map)
## mo <- tapply(ov[[1]], tr1[[1]], mean)
## plot(mo, ty="l")


###################################################
### code chunk number 149: adehabitatLT.Rnw:2080-2081 (eval = FALSE)
###################################################
## ov <- over(tr1, puechabonsp$map)
## mo <- tapply(ov[[1]], tr1[[1]], mean)
## plot(mo, ty="l")


###################################################
### code chunk number 150: adehabitatLT.Rnw:2085-2088
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
ov <- over(tr1, puechabonsp$map)
mo <- tapply(ov[[1]], tr1[[1]], mean)
plot(mo, ty="l")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 151: adehabitatLT.Rnw:2101-2127
###################################################
val <- lapply(1:length(ii), function(i) {

    ## get the rasterized trajectory
    tr <- ii[[i]]

    ## over with the map
    ov <- over(tr, puechabonsp$map)

    ## calculate the mean elevation
    mo <- tapply(ov[[1]], tr[[1]], mean)

    ## prepare the output
    elev <- rep(NA, nrow(puechcirc[[i]]))

    ## place the average values at the right place
    ## names(mo) contains the step number (i.e. relocation
    ## number +1)
    elev[as.numeric(names(mo))+1] <- mo
    df <- data.frame(elevation = elev)

    ## same row.names as in the ltraj
    row.names(df) <- row.names(puechcirc[[i]])

    return(df)
})



###################################################
### code chunk number 152: adehabitatLT.Rnw:2132-2133
###################################################
infolocs(puechcirc) <- val


###################################################
### code chunk number 153: ksdfjeiii (eval = FALSE)
###################################################
## plotltr(puechcirc, "elevation")


###################################################
### code chunk number 154: adehabitatLT.Rnw:2143-2144 (eval = FALSE)
###################################################
## plotltr(puechcirc, "elevation")


###################################################
### code chunk number 155: adehabitatLT.Rnw:2148-2151
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plotltr(puechcirc, "elevation")
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 156: adehabitatLT.Rnw:2192-2194
###################################################
sim <- simm.crw(1:1000, r=0.95)
sim


###################################################
### code chunk number 157: plotcrwsim (eval = FALSE)
###################################################
## plot(sim, addp=FALSE)


###################################################
### code chunk number 158: adehabitatLT.Rnw:2207-2208 (eval = FALSE)
###################################################
## plot(sim, addp=FALSE)


###################################################
### code chunk number 159: adehabitatLT.Rnw:2212-2215
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(sim, addp=FALSE)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 160: figoriboar (eval = FALSE)
###################################################
## data(puechcirc)
## data(puechabonsp)
## boar1 <- puechcirc[1]
## xo <- coordinates(puechabonsp$map)
## ## Note that xo is a matrix containing the coordinates of the
## ## pixels of the elevation map (we will use it later to define
## ## the limits of the study area).
## 
## plot(boar1, spixdf=puechabonsp$map, xlim=range(xo[,1]), ylim=range(xo[,2]))


###################################################
### code chunk number 161: adehabitatLT.Rnw:2275-2276 (eval = FALSE)
###################################################
## data(puechcirc)
## data(puechabonsp)
## boar1 <- puechcirc[1]
## xo <- coordinates(puechabonsp$map)
## ## Note that xo is a matrix containing the coordinates of the
## ## pixels of the elevation map (we will use it later to define
## ## the limits of the study area).
## 
## plot(boar1, spixdf=puechabonsp$map, xlim=range(xo[,1]), ylim=range(xo[,2]))


###################################################
### code chunk number 162: adehabitatLT.Rnw:2280-2283
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
data(puechcirc)
data(puechabonsp)
boar1 <- puechcirc[1]
xo <- coordinates(puechabonsp$map)
## Note that xo is a matrix containing the coordinates of the
## pixels of the elevation map (we will use it later to define
## the limits of the study area).

plot(boar1, spixdf=puechabonsp$map, xlim=range(xo[,1]), ylim=range(xo[,2]))
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 163: adehabitatLT.Rnw:2364-2371
###################################################
plotfun <- function(x, par)
{
    image(par)
    points(x[,1:2], pch=16)
    lines(x[,1:2])
    return(x)
}


###################################################
### code chunk number 164: adehabitatLT.Rnw:2379-2384
###################################################
nmo <- NMs.randomShiftRotation(na.omit(boar1), rshift = TRUE, rrot = TRUE,
                               rx = range(xo[,1]), ry = range(xo[,2]),
                               treatment.func = plotfun,
                               treatment.par = puechabonsp$map[,1], nrep=9)
nmo


###################################################
### code chunk number 165: figtestNM (eval = FALSE)
###################################################
## suppressWarnings(RNGversion("3.5.0"))
## set.seed(90909)
## par(mfrow=c(3,3), mar=c(0,0,0,0))
## resu <- testNM(nmo, count = FALSE)


###################################################
### code chunk number 166: adehabitatLT.Rnw:2406-2407 (eval = FALSE)
###################################################
## suppressWarnings(RNGversion("3.5.0"))
## set.seed(90909)
## par(mfrow=c(3,3), mar=c(0,0,0,0))
## resu <- testNM(nmo, count = FALSE)


###################################################
### code chunk number 167: adehabitatLT.Rnw:2411-2414
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
suppressWarnings(RNGversion("3.5.0"))
set.seed(90909)
par(mfrow=c(3,3), mar=c(0,0,0,0))
resu <- testNM(nmo, count = FALSE)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 168: adehabitatLT.Rnw:2435-2446
###################################################
confun <- function(x, par)
{
    ## Define a SpatialPointsDataFrame from the trajectory
    coordinates(x) <- x[,1:2]
    ## overlap the relocations x to the elevation map par
    jo <- join(x, par)
    ## checks that there are no missing value
    res <- all(!is.na(jo))
    ## return this check
    return(res)
}


###################################################
### code chunk number 169: adehabitatLT.Rnw:2452-2459
###################################################
nmo2 <- NMs.randomShiftRotation(na.omit(boar1), rshift = TRUE, rrot = TRUE,
                                rx = range(xo[,1]), ry = range(xo[,2]),
                                treatment.func = plotfun,
                                treatment.par = puechabonsp$map[,1],
                                constraint.func = confun,
                                constraint.par = puechabonsp$map[,1],
                                nrep=9)


###################################################
### code chunk number 170: figtestnm2 (eval = FALSE)
###################################################
## suppressWarnings(RNGversion("3.5.0"))
## set.seed(90909)
## par(mfrow=c(3,3), mar=c(0,0,0,0))
## resu <- testNM(nmo2, count = FALSE)


###################################################
### code chunk number 171: adehabitatLT.Rnw:2472-2473 (eval = FALSE)
###################################################
## suppressWarnings(RNGversion("3.5.0"))
## set.seed(90909)
## par(mfrow=c(3,3), mar=c(0,0,0,0))
## resu <- testNM(nmo2, count = FALSE)


###################################################
### code chunk number 172: adehabitatLT.Rnw:2477-2480
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
suppressWarnings(RNGversion("3.5.0"))
set.seed(90909)
par(mfrow=c(3,3), mar=c(0,0,0,0))
resu <- testNM(nmo2, count = FALSE)
dev.null <- dev.off()
cat("\\includegraphics[height=10cm,width=10cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 173: adehabitatLT.Rnw:2494-2500
###################################################
varel <- function(x, par)
{
    coordinates(x) <- x[,1:2]
    jo <- join(x, par)
    return(var(jo))
}


###################################################
### code chunk number 174: adehabitatLT.Rnw:2509-2516
###################################################
nmo3 <- NMs.randomShiftRotation(na.omit(boar1), rshift = TRUE, rrot = TRUE,
                                rx = range(xo[,1]), ry = range(xo[,2]),
                                treatment.func = varel,
                                treatment.par = puechabonsp$map[,1],
                                constraint.func = confun,
                                constraint.par = puechabonsp$map[,1],
                                nrep=99)


###################################################
### code chunk number 175: adehabitatLT.Rnw:2523-2524 (eval = FALSE)
###################################################
## sim <- testNM(nmo3, count=FALSE)


###################################################
### code chunk number 176: adehabitatLT.Rnw:2527-2528
###################################################
load("sim.Rdata")


###################################################
### code chunk number 177: adehabitatLT.Rnw:2533-2534
###################################################
(obs <- varel(na.omit(boar1)[[1]], puechabonsp$map[,1]))


###################################################
### code chunk number 178: adehabitatLT.Rnw:2541-2543
###################################################
(ran <- as.randtest(unlist(sim[[1]]), obs))
plot(ran)


###################################################
### code chunk number 179: adehabitatLT.Rnw:2555-2557
###################################################
puechcirc
plot(puechcirc)


###################################################
### code chunk number 180: adehabitatLT.Rnw:2564-2565
###################################################
(boar <- bindltraj(puechcirc))


###################################################
### code chunk number 181: adehabitatLT.Rnw:2575-2582
###################################################
nmo4 <- NMs.randomShiftRotation(na.omit(boar), rshift = TRUE, rrot = TRUE,
                                rx = range(xo[,1]), ry = range(xo[,2]),
                                treatment.func = varel,
                                treatment.par = puechabonsp$map[,1],
                                constraint.func = confun,
                                constraint.par = puechabonsp$map[,1],
                                nrep=99)


###################################################
### code chunk number 182: adehabitatLT.Rnw:2587-2588 (eval = FALSE)
###################################################
## sim2 <- testNM(nmo4, count=FALSE)


###################################################
### code chunk number 183: adehabitatLT.Rnw:2591-2592
###################################################
load("sim2.Rdata")


###################################################
### code chunk number 184: adehabitatLT.Rnw:2601-2604
###################################################
(obs <- lapply(na.omit(boar), function(x) {
  varel(x, puechabonsp$map[,1])
}))


###################################################
### code chunk number 185: adehabitatLT.Rnw:2609-2612
###################################################
lapply(1:2, function(i) {
  as.randtest(unlist(sim2[[i]]), obs[[i]])
})


###################################################
### code chunk number 186: adehabitatLT.Rnw:2659-2668
###################################################
meanvarel <- function(x, par)
{
    livar <- lapply(x, function(y) {
        coordinates(y) <- y[,1:2]
        jo <- join(y, par)
        return(var(jo))
    })
    mean(unlist(livar))
}


###################################################
### code chunk number 187: adehabitatLT.Rnw:2675-2677
###################################################
nmo5 <- NMs2NMm(nmo4, treatment.func = meanvarel,
                treatment.par = puechabonsp$map, nrep = 99)


###################################################
### code chunk number 188: adehabitatLT.Rnw:2686-2687 (eval = FALSE)
###################################################
## sim3 <- testNM(nmo5, count=FALSE)


###################################################
### code chunk number 189: adehabitatLT.Rnw:2690-2691
###################################################
load("sim3.Rdata")


###################################################
### code chunk number 190: adehabitatLT.Rnw:2701-2702
###################################################
(obs <- meanvarel(na.omit(boar), puechabonsp$map))


###################################################
### code chunk number 191: adehabitatLT.Rnw:2708-2710
###################################################
(ran <- as.randtest(unlist(sim3), obs))
plot(ran)


