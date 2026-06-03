### R code from vignette source 'adehabitatHS.Rnw'

###################################################
### code chunk number 1: adehabitatHS.Rnw:29-37
###################################################
owidth <- getOption("width")
options("width"=80)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0
wi <- 480
pt <- 30
reso <- 40


###################################################
### code chunk number 2: afig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
##           .PngNo, ".png", sep="")
## png(file=file, width = wi, height = wi, pointsize = pt, res=reso)


###################################################
### code chunk number 3: afigi (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
##           .PngNo, ".png", sep="")
## png(file=file, width = wi, height = wi, pointsize = pt, res=20)


###################################################
### code chunk number 4: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 5: zfigg (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 6: adehabitatHS.Rnw:112-113
###################################################
library(adehabitatHS)


###################################################
### code chunk number 7: adehabitatHS.Rnw:116-117
###################################################
set.seed(13431)


###################################################
### code chunk number 8: strdesI (eval = FALSE)
###################################################
## par(mar=c(0,0,0,0))
## plot(0,0, xlim=c(0,1), ylim=c(0,1), ty="n")
## polygon(c(0.1, 0.5, 0.5, 0.1, 0.1),
##         c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
## 
## polygon(c(0.6, 0.65, 0.65, 0.6, 0.6),
##         c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
## 
## polygon(c(0.8, 0.85, 0.85, 0.8, 0.8),
##         c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
## text(0.07, 0.2, "N", font=3)
## text(0.5, 0.83, "P", font=3)
## text(0.3, 0.5, "X", font=2, cex=2)
## text(0.625, 0.5, "a", font=2, cex=2)
## text(0.825, 0.5, "u", font=2, cex=2)
## 


###################################################
### code chunk number 9: adehabitatHS.Rnw:212-215
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
par(mar=c(0,0,0,0))
plot(0,0, xlim=c(0,1), ylim=c(0,1), ty="n")
polygon(c(0.1, 0.5, 0.5, 0.1, 0.1),
        c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)

polygon(c(0.6, 0.65, 0.65, 0.6, 0.6),
        c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)

polygon(c(0.8, 0.85, 0.85, 0.8, 0.8),
        c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
text(0.07, 0.2, "N", font=3)
text(0.5, 0.83, "P", font=3)
text(0.3, 0.5, "X", font=2, cex=2)
text(0.625, 0.5, "a", font=2, cex=2)
text(0.825, 0.5, "u", font=2, cex=2)

dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 10: bastrdesII (eval = FALSE)
###################################################
## par(mar=c(0,0,0,0))
## plot(0,0, xlim=c(0,1), ylim=c(0,1), ty="n")
## polygon(c(0.1, 0.5, 0.5, 0.1, 0.1),
##         c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
## polygon(c(0.55, 0.6, 0.6, 0.55, 0.55),
##         c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
## polygon(c(0.65, 1, 1, 0.65, 0.65),
##         c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
## text(0.07, 0.2, "N", font=3)
## text(0.5, 0.83, "P", font=3)
## text(0.3, 0.5, "X", font=2, cex=2)
## text(0.575, 0.5, "a", font=2, cex=2)
## text(0.825, 0.5, "U", font=2, cex=2)
## text(1, 0.83, "K", font=3)


###################################################
### code chunk number 11: adehabitatHS.Rnw:254-257
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
par(mar=c(0,0,0,0))
plot(0,0, xlim=c(0,1), ylim=c(0,1), ty="n")
polygon(c(0.1, 0.5, 0.5, 0.1, 0.1),
        c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
polygon(c(0.55, 0.6, 0.6, 0.55, 0.55),
        c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
polygon(c(0.65, 1, 1, 0.65, 0.65),
        c(0.2, 0.2, 0.8, 0.8, 0.2), col="grey", lwd=2)
text(0.07, 0.2, "N", font=3)
text(0.5, 0.83, "P", font=3)
text(0.3, 0.5, "X", font=2, cex=2)
text(0.575, 0.5, "a", font=2, cex=2)
text(0.825, 0.5, "U", font=2, cex=2)
text(1, 0.83, "K", font=3)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 12: essaidesIII (eval = FALSE)
###################################################
## par(mar=c(0,0,0,0))
## plot(0,0, xlim=c(0,1), ylim=c(0,1), ty="n")
## polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
##         c(0.1, 0.1, 0.3, 0.3, 0.3), col="grey", lwd=2)
## text(0.2, 0.2, expression(X[k]), font=2, cex=2)
## text(0.07, 0.1, expression(N[k]), font=3)
## 
## polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
##         c(0.3, 0.3, 0.5, 0.5, 0.5), col="grey", lwd=2)
## text(0.2, 0.4, "...", font=2, cex=2)
## 
## polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
##         c(0.5, 0.5, 0.7, 0.7, 0.7), col="grey", lwd=2)
## text(0.2, 0.6, expression(X[3]), font=2, cex=2)
## text(0.07, 0.5, expression(N[3]), font=3)
## 
## polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
##         c(0.7, 0.7, 0.85, 0.85, 0.85), col="grey", lwd=2)
## text(0.2, 0.775, expression(X[2]), font=2, cex=2)
## text(0.07, 0.7, expression(N[2]), font=3)
## 
## polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
##         c(0.85, 0.85, 0.95, 0.95, 0.95), col="grey", lwd=2)
## text(0.2, 0.9, expression(X[1]), font=2, cex=2)
## text(0.07, 0.85, expression(N[1]), font=3)
## text(0.3, 0.98, "P", font=3)
## 
## ## Le vecteur
## polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
##         c(0.1, 0.1, 0.3, 0.3, 0.3), col="grey", lwd=2)
## text(0.52, 0.2, expression(a[k]), font=2, cex=1.5)
## 
## polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
##         c(0.3, 0.3, 0.5, 0.5, 0.5), col="grey", lwd=2)
## text(0.52, 0.4, "...", font=2, cex=1.5)
## 
## polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
##         c(0.5, 0.5, 0.7, 0.7, 0.7), col="grey", lwd=2)
## text(0.52, 0.6, expression(a[3]), font=2, cex=1.5)
## 
## polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
##         c(0.7, 0.7, 0.85, 0.85, 0.85), col="grey", lwd=2)
## text(0.52, 0.775, expression(a[2]), font=2, cex=1.5)
## 
## polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
##         c(0.85, 0.85, 0.95, 0.95, 0.95), col="grey", lwd=2)
## text(0.52, 0.9, expression(a[1]), font=2, cex=1.5)
## 
## 
## 
## ## Le vecteur
## polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
##         c(0.1, 0.1, 0.3, 0.3, 0.3), col="grey", lwd=2)
## text(0.72, 0.2, expression(u[k]), font=2, cex=1.5)
## 
## polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
##         c(0.3, 0.3, 0.5, 0.5, 0.5), col="grey", lwd=2)
## text(0.72, 0.4, "...", font=2, cex=1.5)
## 
## polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
##         c(0.5, 0.5, 0.7, 0.7, 0.7), col="grey", lwd=2)
## text(0.72, 0.6, expression(u[3]), font=2, cex=1.5)
## 
## polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
##         c(0.7, 0.7, 0.85, 0.85, 0.85), col="grey", lwd=2)
## text(0.72, 0.775, expression(u[2]), font=2, cex=1.5)
## 
## polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
##         c(0.85, 0.85, 0.95, 0.95, 0.95), col="grey", lwd=2)
## text(0.72, 0.9, expression(u[1]), font=2, cex=1.5)
## 


###################################################
### code chunk number 13: adehabitatHS.Rnw:357-360
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
par(mar=c(0,0,0,0))
plot(0,0, xlim=c(0,1), ylim=c(0,1), ty="n")
polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
        c(0.1, 0.1, 0.3, 0.3, 0.3), col="grey", lwd=2)
text(0.2, 0.2, expression(X[k]), font=2, cex=2)
text(0.07, 0.1, expression(N[k]), font=3)

polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
        c(0.3, 0.3, 0.5, 0.5, 0.5), col="grey", lwd=2)
text(0.2, 0.4, "...", font=2, cex=2)

polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
        c(0.5, 0.5, 0.7, 0.7, 0.7), col="grey", lwd=2)
text(0.2, 0.6, expression(X[3]), font=2, cex=2)
text(0.07, 0.5, expression(N[3]), font=3)

polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
        c(0.7, 0.7, 0.85, 0.85, 0.85), col="grey", lwd=2)
text(0.2, 0.775, expression(X[2]), font=2, cex=2)
text(0.07, 0.7, expression(N[2]), font=3)

polygon(c(0.1, 0.3, 0.3, 0.1, 0.1),
        c(0.85, 0.85, 0.95, 0.95, 0.95), col="grey", lwd=2)
text(0.2, 0.9, expression(X[1]), font=2, cex=2)
text(0.07, 0.85, expression(N[1]), font=3)
text(0.3, 0.98, "P", font=3)

## Le vecteur
polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
        c(0.1, 0.1, 0.3, 0.3, 0.3), col="grey", lwd=2)
text(0.52, 0.2, expression(a[k]), font=2, cex=1.5)

polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
        c(0.3, 0.3, 0.5, 0.5, 0.5), col="grey", lwd=2)
text(0.52, 0.4, "...", font=2, cex=1.5)

polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
        c(0.5, 0.5, 0.7, 0.7, 0.7), col="grey", lwd=2)
text(0.52, 0.6, expression(a[3]), font=2, cex=1.5)

polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
        c(0.7, 0.7, 0.85, 0.85, 0.85), col="grey", lwd=2)
text(0.52, 0.775, expression(a[2]), font=2, cex=1.5)

polygon(c(0.48, 0.55, 0.55, 0.48, 0.48),
        c(0.85, 0.85, 0.95, 0.95, 0.95), col="grey", lwd=2)
text(0.52, 0.9, expression(a[1]), font=2, cex=1.5)



## Le vecteur
polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
        c(0.1, 0.1, 0.3, 0.3, 0.3), col="grey", lwd=2)
text(0.72, 0.2, expression(u[k]), font=2, cex=1.5)

polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
        c(0.3, 0.3, 0.5, 0.5, 0.5), col="grey", lwd=2)
text(0.72, 0.4, "...", font=2, cex=1.5)

polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
        c(0.5, 0.5, 0.7, 0.7, 0.7), col="grey", lwd=2)
text(0.72, 0.6, expression(u[3]), font=2, cex=1.5)

polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
        c(0.7, 0.7, 0.85, 0.85, 0.85), col="grey", lwd=2)
text(0.72, 0.775, expression(u[2]), font=2, cex=1.5)

polygon(c(0.68, 0.75, 0.75, 0.68, 0.68),
        c(0.85, 0.85, 0.95, 0.95, 0.95), col="grey", lwd=2)
text(0.72, 0.9, expression(u[1]), font=2, cex=1.5)

dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 14: ploniche (eval = FALSE)
###################################################
## set.seed(321)
## par(mar=c(0,0,0,0))
## plot(0,0, xlim=c(-1,1), ylim=c(-1,1), ty="n")
## arrows(0,0, 0, 0.9, lwd=2)
## arrows(0,0, -0.9, -0.5, lwd=2)
## arrows(0,0, 0.9, -0.5, lwd=2)
## points(data.frame(rnorm(500,sd=0.25), rnorm(500, sd=0.25)),
##        pch=21, bg="grey")
## points(data.frame(rnorm(50, mean=0.3, sd=0.1),
##                   rnorm(50, mean=0.3, sd=0.1)),
##        pch=21, bg="red", cex=2)
## text(0,0.96, expression(V[3]), cex=1.5)
## text(-0.96, -0.5, expression(V[1]), cex=1.5)
## text(0.96, -0.5, expression(V[2]), cex=1.5)
## 


###################################################
### code chunk number 15: adehabitatHS.Rnw:416-419
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
set.seed(321)
par(mar=c(0,0,0,0))
plot(0,0, xlim=c(-1,1), ylim=c(-1,1), ty="n")
arrows(0,0, 0, 0.9, lwd=2)
arrows(0,0, -0.9, -0.5, lwd=2)
arrows(0,0, 0.9, -0.5, lwd=2)
points(data.frame(rnorm(500,sd=0.25), rnorm(500, sd=0.25)),
       pch=21, bg="grey")
points(data.frame(rnorm(50, mean=0.3, sd=0.1),
                  rnorm(50, mean=0.3, sd=0.1)),
       pch=21, bg="red", cex=2)
text(0,0.96, expression(V[3]), cex=1.5)
text(-0.96, -0.5, expression(V[1]), cex=1.5)
text(0.96, -0.5, expression(V[2]), cex=1.5)

dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 16: adehabitatHS.Rnw:583-585
###################################################
data(bauges)
names(bauges)


###################################################
### code chunk number 17: fig1 (eval = FALSE)
###################################################
## map <- bauges$map
## mimage(map)


###################################################
### code chunk number 18: adehabitatHS.Rnw:596-597 (eval = FALSE)
###################################################
## map <- bauges$map
## mimage(map)


###################################################
### code chunk number 19: adehabitatHS.Rnw:601-604
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
map <- bauges$map
mimage(map)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 20: fig2 (eval = FALSE)
###################################################
## image(map)
## locs <- bauges$locs
## points(locs, pch=3)


###################################################
### code chunk number 21: adehabitatHS.Rnw:621-622 (eval = FALSE)
###################################################
## image(map)
## locs <- bauges$locs
## points(locs, pch=3)


###################################################
### code chunk number 22: adehabitatHS.Rnw:626-629
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
image(map)
locs <- bauges$locs
points(locs, pch=3)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 23: adehabitatHS.Rnw:640-641
###################################################
cp <- count.points(locs, map)


###################################################
### code chunk number 24: adehabitatHS.Rnw:648-650
###################################################
tab <- slot(map, "data")
pr <- slot(cp, "data")[,1]


###################################################
### code chunk number 25: fig3 (eval = FALSE)
###################################################
## histniche(tab, pr)


###################################################
### code chunk number 26: adehabitatHS.Rnw:669-670 (eval = FALSE)
###################################################
## histniche(tab, pr)


###################################################
### code chunk number 27: adehabitatHS.Rnw:674-677
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
histniche(tab, pr)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 28: adehabitatHS.Rnw:693-694
###################################################
args(gnesfa)


###################################################
### code chunk number 29: adehabitatHS.Rnw:750-751
###################################################
pc <- dudi.pca(tab, scannf=FALSE)


###################################################
### code chunk number 30: adehabitatHS.Rnw:779-780
###################################################
gn <- gnesfa(pc, Focus = pr, scan=FALSE, nfFirst=1, nfLast=1)


###################################################
### code chunk number 31: fig5 (eval = FALSE)
###################################################
## barplot(gn$eig)


###################################################
### code chunk number 32: adehabitatHS.Rnw:804-807
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
barplot(gn$eig)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 33: fig6 (eval = FALSE)
###################################################
## barplot(1/gn$eig)


###################################################
### code chunk number 34: adehabitatHS.Rnw:821-824
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
barplot(1/gn$eig)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 35: adehabitatHS.Rnw:833-834
###################################################
gn


###################################################
### code chunk number 36: fig7 (eval = FALSE)
###################################################
## scatterniche(gn$li, pr, pts=TRUE)


###################################################
### code chunk number 37: adehabitatHS.Rnw:845-846 (eval = FALSE)
###################################################
## scatterniche(gn$li, pr, pts=TRUE)


###################################################
### code chunk number 38: adehabitatHS.Rnw:850-853
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
scatterniche(gn$li, pr, pts=TRUE)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 39: adehabitatHS.Rnw:884-885
###################################################
gn2 <- gnesfa(pc, Reference = pr, scan=FALSE, nfFirst=2, nfLast=0)


###################################################
### code chunk number 40: fig8 (eval = FALSE)
###################################################
## barplot(gn2$eig)


###################################################
### code chunk number 41: adehabitatHS.Rnw:893-896
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
barplot(gn2$eig)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 42: fig9 (eval = FALSE)
###################################################
## scatterniche(gn2$li, pr)


###################################################
### code chunk number 43: adehabitatHS.Rnw:909-910 (eval = FALSE)
###################################################
## scatterniche(gn2$li, pr)


###################################################
### code chunk number 44: adehabitatHS.Rnw:914-917
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
scatterniche(gn2$li, pr)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 45: fig10 (eval = FALSE)
###################################################
## s.arrow(gn2$co)


###################################################
### code chunk number 46: adehabitatHS.Rnw:932-933 (eval = FALSE)
###################################################
## s.arrow(gn2$co)


###################################################
### code chunk number 47: adehabitatHS.Rnw:937-940
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
s.arrow(gn2$co)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 48: fig11 (eval = FALSE)
###################################################
## s.arrow(gn2$cor)


###################################################
### code chunk number 49: adehabitatHS.Rnw:958-959 (eval = FALSE)
###################################################
## s.arrow(gn2$cor)


###################################################
### code chunk number 50: adehabitatHS.Rnw:963-966
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
s.arrow(gn2$cor)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 51: adehabitatHS.Rnw:998-999
###################################################
(mad <- madifa(pc, pr, scan=FALSE))


###################################################
### code chunk number 52: adehabitatHS.Rnw:1005-1007
###################################################
mad$eig
gn2$eig


###################################################
### code chunk number 53: adehabitatHS.Rnw:1015-1016
###################################################
pt <- 12


###################################################
### code chunk number 54: fig14 (eval = FALSE)
###################################################
## plot(mad, map)


###################################################
### code chunk number 55: adehabitatHS.Rnw:1023-1024 (eval = FALSE)
###################################################
## plot(mad, map)


###################################################
### code chunk number 56: adehabitatHS.Rnw:1028-1031
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
plot(mad, map)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 57: adehabitatHS.Rnw:1035-1036
###################################################
pt <- 20


###################################################
### code chunk number 58: fig16 (eval = FALSE)
###################################################
## MD <- mahasuhab(map, locs)
## image(MD)


###################################################
### code chunk number 59: adehabitatHS.Rnw:1083-1084 (eval = FALSE)
###################################################
## MD <- mahasuhab(map, locs)
## image(MD)


###################################################
### code chunk number 60: adehabitatHS.Rnw:1088-1091
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=20)
MD <- mahasuhab(map, locs)
image(MD)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 61: figRMD (eval = FALSE)
###################################################
## RMD1 <- data.frame(RMD1=mad$li[,1]^2)
## coordinates(RMD1) <- coordinates(map)
## gridded(RMD1) <- TRUE
## image(RMD1)


###################################################
### code chunk number 62: adehabitatHS.Rnw:1115-1116 (eval = FALSE)
###################################################
## RMD1 <- data.frame(RMD1=mad$li[,1]^2)
## coordinates(RMD1) <- coordinates(map)
## gridded(RMD1) <- TRUE
## image(RMD1)


###################################################
### code chunk number 63: adehabitatHS.Rnw:1120-1123
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=20)
RMD1 <- data.frame(RMD1=mad$li[,1]^2)
coordinates(RMD1) <- coordinates(map)
gridded(RMD1) <- TRUE
image(RMD1)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 64: figRMD2 (eval = FALSE)
###################################################
## RMD2 <- data.frame(RMD2 = apply(mad$li[,1:2], 1, function(x) sum(x^2)))
## coordinates(RMD2) <- coordinates(map)
## gridded(RMD2) <- TRUE
## image(RMD2)


###################################################
### code chunk number 65: adehabitatHS.Rnw:1138-1139 (eval = FALSE)
###################################################
## RMD2 <- data.frame(RMD2 = apply(mad$li[,1:2], 1, function(x) sum(x^2)))
## coordinates(RMD2) <- coordinates(map)
## gridded(RMD2) <- TRUE
## image(RMD2)


###################################################
### code chunk number 66: adehabitatHS.Rnw:1143-1146
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=20)
RMD2 <- data.frame(RMD2 = apply(mad$li[,1:2], 1, function(x) sum(x^2)))
coordinates(RMD2) <- coordinates(map)
gridded(RMD2) <- TRUE
image(RMD2)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 67: adehabitatHS.Rnw:1189-1193
###################################################
en1 <- enfa(pc, pr, scan=FALSE)
gn3 <- gnesfa(pc, Reference=pr, scan=FALSE, centering="twice")
en1$s
gn3$eig


###################################################
### code chunk number 68: baploenfa (eval = FALSE)
###################################################
## barplot(en1$s)


###################################################
### code chunk number 69: adehabitatHS.Rnw:1202-1203 (eval = FALSE)
###################################################
## barplot(en1$s)


###################################################
### code chunk number 70: adehabitatHS.Rnw:1207-1210
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
barplot(en1$s)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 71: scattenfa (eval = FALSE)
###################################################
## scatter(en1)


###################################################
### code chunk number 72: adehabitatHS.Rnw:1233-1234 (eval = FALSE)
###################################################
## scatter(en1)


###################################################
### code chunk number 73: adehabitatHS.Rnw:1238-1241
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
scatter(en1)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 74: plotenmad (eval = FALSE)
###################################################
## plot(mad$li[,1], en1$li[,1], xlab="First axis of the MADIFA",
##      ylab="Marginality axis")


###################################################
### code chunk number 75: adehabitatHS.Rnw:1272-1273 (eval = FALSE)
###################################################
## plot(mad$li[,1], en1$li[,1], xlab="First axis of the MADIFA",
##      ylab="Marginality axis")


###################################################
### code chunk number 76: adehabitatHS.Rnw:1277-1280
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
plot(mad$li[,1], en1$li[,1], xlab="First axis of the MADIFA",
     ylab="Marginality axis")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 77: adehabitatHS.Rnw:1364-1365
###################################################
dun <- dunnfa(pc, pr, scann=FALSE)


###################################################
### code chunk number 78: bpdun (eval = FALSE)
###################################################
## barplot(dun$eig)


###################################################
### code chunk number 79: adehabitatHS.Rnw:1374-1375 (eval = FALSE)
###################################################
## barplot(dun$eig)


###################################################
### code chunk number 80: adehabitatHS.Rnw:1379-1382
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
barplot(dun$eig)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 81: hindun (eval = FALSE)
###################################################
## histniche(data.frame(dun$liA[,1]), pr, main="First Axis")


###################################################
### code chunk number 82: adehabitatHS.Rnw:1394-1395 (eval = FALSE)
###################################################
## histniche(data.frame(dun$liA[,1]), pr, main="First Axis")


###################################################
### code chunk number 83: adehabitatHS.Rnw:1399-1402
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
histniche(data.frame(dun$liA[,1]), pr, main="First Axis")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 84: scordun (eval = FALSE)
###################################################
## s.arrow(dun$cor)


###################################################
### code chunk number 85: adehabitatHS.Rnw:1421-1422 (eval = FALSE)
###################################################
## s.arrow(dun$cor)


###################################################
### code chunk number 86: adehabitatHS.Rnw:1426-1429
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
s.arrow(dun$cor)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 87: figdunmad (eval = FALSE)
###################################################
## plot(dun$liA[,1], mad$li[,1], xlab="DUNNFA 1", ylab="MADIFA 1")


###################################################
### code chunk number 88: adehabitatHS.Rnw:1442-1443 (eval = FALSE)
###################################################
## plot(dun$liA[,1], mad$li[,1], xlab="DUNNFA 1", ylab="MADIFA 1")


###################################################
### code chunk number 89: adehabitatHS.Rnw:1447-1450
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
plot(dun$liA[,1], mad$li[,1], xlab="DUNNFA 1", ylab="MADIFA 1")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 90: adehabitatHS.Rnw:1489-1492
###################################################
elev <- map[,1]
av <- factor(cut(slot(elev, "data")[,1], 4),
             labels=c("Low","Medium","High","Very High"))


###################################################
### code chunk number 91: adehabitatHS.Rnw:1497-1498
###################################################
(tav <- table(av))


###################################################
### code chunk number 92: cartelev (eval = FALSE)
###################################################
## slot(elev, "data")[,1] <- as.numeric(av)
## image(elev)


###################################################
### code chunk number 93: adehabitatHS.Rnw:1508-1509 (eval = FALSE)
###################################################
## slot(elev, "data")[,1] <- as.numeric(av)
## image(elev)


###################################################
### code chunk number 94: adehabitatHS.Rnw:1513-1516
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=20)
slot(elev, "data")[,1] <- as.numeric(av)
image(elev)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 95: adehabitatHS.Rnw:1524-1528
###################################################
us <- join(locs, elev)
tus <- table(us)
names(tus) <- names(tav)
tus


###################################################
### code chunk number 96: adehabitatHS.Rnw:1555-1559
###################################################
class(tus) <- NULL
tav <- tav/sum(tav)
class(tav) <- NULL
(Wi <- widesI(tus, tav))


###################################################
### code chunk number 97: adehabitatHS.Rnw:1568-1569 (eval = FALSE)
###################################################
## plot(Wi)


###################################################
### code chunk number 98: adehabitatHS.Rnw:1589-1592
###################################################
dis <- acm.disjonctif(slot(elev, "data"))
pc2 <- dudi.pca(dis, scan=FALSE)
gnf <- gnesfa(pc2, pr, scan=FALSE)


###################################################
### code chunk number 99: adehabitatHS.Rnw:1603-1605
###################################################
sum(gnf$eig)
sum(Wi$wi)


###################################################
### code chunk number 100: figs (eval = FALSE)
###################################################
## set.seed(321)
## par(mar=c(0,0,0,0))
## plot(0,0, xlim=c(-1,1), ylim=c(-1,1), ty="n", axes=FALSE)
## arrows(0,0, 0, 0.9, lwd=2)
## arrows(0,0, -0.9, -0.5, lwd=2)
## arrows(0,0, 0.9, -0.5, lwd=2)
## points(data.frame(rnorm(500,sd=0.25), rnorm(500, sd=0.25)),
##        pch=16, col="grey")
## points(data.frame(rnorm(50, mean=0.3, sd=0.1),
##                   rnorm(50, mean=0.3, sd=0.1)),
##        pch=21, bg="red", cex=2)
## points(data.frame(rnorm(50, mean=0, sd=0.1),
##                   rnorm(50, mean=0.4, sd=0.1)),
##        pch=21, bg="green", cex=2)
## points(data.frame(rnorm(50, mean=-0.3, sd=0.1),
##                   rnorm(50, mean=0.3, sd=0.1)),
##        pch=21, bg="blue", cex=2)
## text(0,0.96, expression(V[3]), cex=1.5)
## text(-0.96, -0.5, expression(V[1]), cex=1.5)
## text(0.96, -0.5, expression(V[2]), cex=1.5)
## arrows(0,0,0.3,0.3, lwd=4)
## arrows(0,0,0,0.4, lwd=4)
## arrows(0,0,-0.3,0.3, lwd=4)
## arrows(0,0,0.3,0.3, lwd=2, col="pink")
## arrows(0,0,0,0.4, lwd=2, col="lightgreen")
## arrows(0,0,-0.3,0.3, lwd=2, col="lightblue")
## 


###################################################
### code chunk number 101: adehabitatHS.Rnw:1661-1664
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
set.seed(321)
par(mar=c(0,0,0,0))
plot(0,0, xlim=c(-1,1), ylim=c(-1,1), ty="n", axes=FALSE)
arrows(0,0, 0, 0.9, lwd=2)
arrows(0,0, -0.9, -0.5, lwd=2)
arrows(0,0, 0.9, -0.5, lwd=2)
points(data.frame(rnorm(500,sd=0.25), rnorm(500, sd=0.25)),
       pch=16, col="grey")
points(data.frame(rnorm(50, mean=0.3, sd=0.1),
                  rnorm(50, mean=0.3, sd=0.1)),
       pch=21, bg="red", cex=2)
points(data.frame(rnorm(50, mean=0, sd=0.1),
                  rnorm(50, mean=0.4, sd=0.1)),
       pch=21, bg="green", cex=2)
points(data.frame(rnorm(50, mean=-0.3, sd=0.1),
                  rnorm(50, mean=0.3, sd=0.1)),
       pch=21, bg="blue", cex=2)
text(0,0.96, expression(V[3]), cex=1.5)
text(-0.96, -0.5, expression(V[1]), cex=1.5)
text(0.96, -0.5, expression(V[2]), cex=1.5)
arrows(0,0,0.3,0.3, lwd=4)
arrows(0,0,0,0.4, lwd=4)
arrows(0,0,-0.3,0.3, lwd=4)
arrows(0,0,0.3,0.3, lwd=2, col="pink")
arrows(0,0,0,0.4, lwd=2, col="lightgreen")
arrows(0,0,-0.3,0.3, lwd=2, col="lightblue")

dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 102: adehabitatHS.Rnw:1697-1699
###################################################
data(puech)
names(puech)


###################################################
### code chunk number 103: imdes2 (eval = FALSE)
###################################################
## maps <- puech$maps
## mimage(maps)


###################################################
### code chunk number 104: adehabitatHS.Rnw:1710-1711 (eval = FALSE)
###################################################
## maps <- puech$maps
## mimage(maps)


###################################################
### code chunk number 105: adehabitatHS.Rnw:1715-1718
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
maps <- puech$maps
mimage(maps)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 106: ploloc (eval = FALSE)
###################################################
## locs <- puech$relocations
## image(maps)
## points(locs, col=as.numeric(slot(locs, "data")[,1]), pch=16)


###################################################
### code chunk number 107: adehabitatHS.Rnw:1731-1732 (eval = FALSE)
###################################################
## locs <- puech$relocations
## image(maps)
## points(locs, col=as.numeric(slot(locs, "data")[,1]), pch=16)


###################################################
### code chunk number 108: adehabitatHS.Rnw:1737-1740
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
locs <- puech$relocations
image(maps)
points(locs, col=as.numeric(slot(locs, "data")[,1]), pch=16)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 109: cplocs (eval = FALSE)
###################################################
## cp <- count.points(locs, maps)
## mimage(cp)


###################################################
### code chunk number 110: adehabitatHS.Rnw:1753-1754 (eval = FALSE)
###################################################
## cp <- count.points(locs, maps)
## mimage(cp)


###################################################
### code chunk number 111: adehabitatHS.Rnw:1758-1761
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
cp <- count.points(locs, maps)
mimage(cp)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 112: adehabitatHS.Rnw:1768-1770
###################################################
X <- slot(maps, "data")
U <- slot(cp, "data")


###################################################
### code chunk number 113: adehabitatHS.Rnw:1808-1809
###################################################
pc <- dudi.pca(X, scannf=FALSE)


###################################################
### code chunk number 114: adehabitatHS.Rnw:1814-1815
###################################################
(ni <- niche(pc, U, scannf=FALSE))


###################################################
### code chunk number 115: plotnichee (eval = FALSE)
###################################################
## plot(ni)


###################################################
### code chunk number 116: adehabitatHS.Rnw:1824-1825 (eval = FALSE)
###################################################
## plot(ni)


###################################################
### code chunk number 117: adehabitatHS.Rnw:1829-1832
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
plot(ni)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 118: adehabitatHS.Rnw:1843-1844
###################################################
ni$eig[1]/sum(ni$eig)


###################################################
### code chunk number 119: adehabitatHS.Rnw:1852-1853
###################################################
ni$eig[2]/sum(ni$eig)


###################################################
### code chunk number 120: cartaxeniche (eval = FALSE)
###################################################
## ls <- ni$ls
## coordinates(ls) <- coordinates(maps)
## gridded(ls) <- TRUE
## mimage(ls)


###################################################
### code chunk number 121: adehabitatHS.Rnw:1910-1911 (eval = FALSE)
###################################################
## ls <- ni$ls
## coordinates(ls) <- coordinates(maps)
## gridded(ls) <- TRUE
## mimage(ls)


###################################################
### code chunk number 122: adehabitatHS.Rnw:1915-1918
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
ls <- ni$ls
coordinates(ls) <- coordinates(maps)
gridded(ls) <- TRUE
mimage(ls)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 123: figokcanomi (eval = FALSE)
###################################################
## par(mar=c(0.1,0.1,0.1,0.1))
## xx <- cbind(rnorm(1000), rnorm(1000, sd=0.3))
## plot(xx, asp=1, bg="grey", pch=21, axes=FALSE)
## xx1 <- cbind(rnorm(100, mean=-1.5, sd=0.15),
##              rnorm(100, mean=0.2, sd=0.05))
## xx2 <- cbind(rnorm(100, mean=-0.5, sd=0.15),
##              rnorm(100, mean=0.5, sd=0.05))
## xx3 <- cbind(rnorm(100, mean=0.5, sd=0.15),
##              rnorm(100, mean=0.5, sd=0.05))
## xx4 <- cbind(rnorm(100, mean=1.5, sd=0.15),
##              rnorm(100, mean=0.2, sd=0.05))
## points(xx1, pch=21, bg="red", cex=1.5)
## points(xx2, pch=21, bg="yellow", cex=1.5)
## points(xx3, pch=21, bg="green", cex=1.5)
## points(xx4, pch=21, bg="blue", cex=1.5)
## arrows(0,0,-1.5, 0.2, lwd=5)
## arrows(0,0,-1.5, 0.2, lwd=3, col="red")
## arrows(0,0,-0.5, 0.5, lwd=5)
## arrows(0,0,-0.5, 0.5, lwd=3, col="yellow")
## arrows(0,0,0.5, 0.5, lwd=5)
## arrows(0,0,0.5, 0.5, lwd=3, col="green")
## arrows(0,0,1.5, 0.2, lwd=5)
## arrows(0,0,1.5, 0.2, lwd=3, col="blue")
## box()
## abline(v=0, h=0)
## 


###################################################
### code chunk number 124: adehabitatHS.Rnw:2019-2022
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
par(mar=c(0.1,0.1,0.1,0.1))
xx <- cbind(rnorm(1000), rnorm(1000, sd=0.3))
plot(xx, asp=1, bg="grey", pch=21, axes=FALSE)
xx1 <- cbind(rnorm(100, mean=-1.5, sd=0.15),
             rnorm(100, mean=0.2, sd=0.05))
xx2 <- cbind(rnorm(100, mean=-0.5, sd=0.15),
             rnorm(100, mean=0.5, sd=0.05))
xx3 <- cbind(rnorm(100, mean=0.5, sd=0.15),
             rnorm(100, mean=0.5, sd=0.05))
xx4 <- cbind(rnorm(100, mean=1.5, sd=0.15),
             rnorm(100, mean=0.2, sd=0.05))
points(xx1, pch=21, bg="red", cex=1.5)
points(xx2, pch=21, bg="yellow", cex=1.5)
points(xx3, pch=21, bg="green", cex=1.5)
points(xx4, pch=21, bg="blue", cex=1.5)
arrows(0,0,-1.5, 0.2, lwd=5)
arrows(0,0,-1.5, 0.2, lwd=3, col="red")
arrows(0,0,-0.5, 0.5, lwd=5)
arrows(0,0,-0.5, 0.5, lwd=3, col="yellow")
arrows(0,0,0.5, 0.5, lwd=5)
arrows(0,0,0.5, 0.5, lwd=3, col="green")
arrows(0,0,1.5, 0.2, lwd=5)
arrows(0,0,1.5, 0.2, lwd=3, col="blue")
box()
abline(v=0, h=0)

dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 125: canom (eval = FALSE)
###################################################
## par(mar=c(0.1,0.1,0.1,0.1))
## plot(xx[,1]*0.3, xx[,2], asp=1, bg="grey", pch=21, axes=FALSE)
## xx1 <- cbind(rnorm(100, mean=-1.5, sd=0.15),
##              rnorm(100, mean=0.2, sd=0.05))
## xx2 <- cbind(rnorm(100, mean=-0.5, sd=0.15),
##              rnorm(100, mean=0.5, sd=0.05))
## xx3 <- cbind(rnorm(100, mean=0.5, sd=0.15),
##              rnorm(100, mean=0.5, sd=0.05))
## xx4 <- cbind(rnorm(100, mean=1.5, sd=0.15),
##              rnorm(100, mean=0.2, sd=0.05))
## points(xx1[,1]*0.3, xx1[,2], pch=21, bg="red", cex=1.5)
## points(xx2[,1]*0.3, xx2[,2], pch=21, bg="yellow", cex=1.5)
## points(xx3[,1]*0.3, xx3[,2], pch=21, bg="green", cex=1.5)
## points(xx4[,1]*0.3, xx4[,2], pch=21, bg="blue", cex=1.5)
## arrows(0,0,-1.5*0.3, 0.2, lwd=5)
## arrows(0,0,-1.5*0.3, 0.2, lwd=3, col="red")
## arrows(0,0,-0.5*0.3, 0.5, lwd=5)
## arrows(0,0,-0.5*0.3, 0.5, lwd=3, col="yellow")
## arrows(0,0,0.5*0.3, 0.5, lwd=5)
## arrows(0,0,0.5*0.3, 0.5, lwd=3, col="green")
## arrows(0,0,1.5*0.3, 0.2, lwd=5)
## arrows(0,0,1.5*0.3, 0.2, lwd=3, col="blue")
## box()
## abline(v=0, h=0)


###################################################
### code chunk number 126: adehabitatHS.Rnw:2085-2088
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
par(mar=c(0.1,0.1,0.1,0.1))
plot(xx[,1]*0.3, xx[,2], asp=1, bg="grey", pch=21, axes=FALSE)
xx1 <- cbind(rnorm(100, mean=-1.5, sd=0.15),
             rnorm(100, mean=0.2, sd=0.05))
xx2 <- cbind(rnorm(100, mean=-0.5, sd=0.15),
             rnorm(100, mean=0.5, sd=0.05))
xx3 <- cbind(rnorm(100, mean=0.5, sd=0.15),
             rnorm(100, mean=0.5, sd=0.05))
xx4 <- cbind(rnorm(100, mean=1.5, sd=0.15),
             rnorm(100, mean=0.2, sd=0.05))
points(xx1[,1]*0.3, xx1[,2], pch=21, bg="red", cex=1.5)
points(xx2[,1]*0.3, xx2[,2], pch=21, bg="yellow", cex=1.5)
points(xx3[,1]*0.3, xx3[,2], pch=21, bg="green", cex=1.5)
points(xx4[,1]*0.3, xx4[,2], pch=21, bg="blue", cex=1.5)
arrows(0,0,-1.5*0.3, 0.2, lwd=5)
arrows(0,0,-1.5*0.3, 0.2, lwd=3, col="red")
arrows(0,0,-0.5*0.3, 0.5, lwd=5)
arrows(0,0,-0.5*0.3, 0.5, lwd=3, col="yellow")
arrows(0,0,0.5*0.3, 0.5, lwd=5)
arrows(0,0,0.5*0.3, 0.5, lwd=3, col="green")
arrows(0,0,1.5*0.3, 0.2, lwd=5)
arrows(0,0,1.5*0.3, 0.2, lwd=3, col="blue")
box()
abline(v=0, h=0)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 127: adehabitatHS.Rnw:2105-2106
###################################################
com <- canomi(pc, U, scannf=FALSE)


###################################################
### code chunk number 128: plotcan (eval = FALSE)
###################################################
## plot(com)


###################################################
### code chunk number 129: adehabitatHS.Rnw:2115-2116 (eval = FALSE)
###################################################
## plot(com)


###################################################
### code chunk number 130: adehabitatHS.Rnw:2120-2123
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
plot(com)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 131: adehabitatHS.Rnw:2144-2146
###################################################
pt <- 8
wi <- 900


###################################################
### code chunk number 132: plocompcomni (eval = FALSE)
###################################################
## par(mfrow=c(2,1))
## plot(ni$ls[,1], com$ls[,1],
##      xlab="Scores on the first axis of the OMI",
##      ylab="Scores on the first axis of the can. OMI")
## plot(ni$ls[,2], com$ls[,2],
##      xlab="Scores on the first axis of the OMI",
##      ylab="Scores on the first axis of the can. OMI")


###################################################
### code chunk number 133: adehabitatHS.Rnw:2159-2160 (eval = FALSE)
###################################################
## par(mfrow=c(2,1))
## plot(ni$ls[,1], com$ls[,1],
##      xlab="Scores on the first axis of the OMI",
##      ylab="Scores on the first axis of the can. OMI")
## plot(ni$ls[,2], com$ls[,2],
##      xlab="Scores on the first axis of the OMI",
##      ylab="Scores on the first axis of the can. OMI")


###################################################
### code chunk number 134: adehabitatHS.Rnw:2164-2167
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
par(mfrow=c(2,1))
plot(ni$ls[,1], com$ls[,1],
     xlab="Scores on the first axis of the OMI",
     ylab="Scores on the first axis of the can. OMI")
plot(ni$ls[,2], com$ls[,2],
     xlab="Scores on the first axis of the OMI",
     ylab="Scores on the first axis of the can. OMI")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 135: adehabitatHS.Rnw:2171-2173
###################################################
pt <- 20
wi <- 480


###################################################
### code chunk number 136: adehabitatHS.Rnw:2193-2198
###################################################
slope <- maps[,8]
sl <- slot(slope, "data")[,1]
av <- factor(cut(sl, c(-0.1, 2, 5, 12, 50)),
             labels=c("Low","Medium","High","Very High"))



###################################################
### code chunk number 137: adehabitatHS.Rnw:2203-2204
###################################################
(tav <- table(av))


###################################################
### code chunk number 138: carteslope (eval = FALSE)
###################################################
## slot(slope, "data")[,1] <- as.numeric(av)
## image(slope)


###################################################
### code chunk number 139: adehabitatHS.Rnw:2214-2215 (eval = FALSE)
###################################################
## slot(slope, "data")[,1] <- as.numeric(av)
## image(slope)


###################################################
### code chunk number 140: adehabitatHS.Rnw:2219-2222
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=20)
slot(slope, "data")[,1] <- as.numeric(av)
image(slope)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 141: adehabitatHS.Rnw:2230-2236
###################################################
us <- join(locs, slope)
tus <- table(slot(locs,"data")[,1],us)
class(tus) <- NULL
tus <- as.data.frame(tus)
colnames(tus) <- names(tav)
tus


###################################################
### code chunk number 142: adehabitatHS.Rnw:2245-2249
###################################################
tav2 <- matrix(rep(tav, nrow(tus)), nrow=nrow(tus), byrow=TRUE)
colnames(tav2) <- names(tav)
compana(tus, tav2, test = "randomisation",
        rnv = 0.01, nrep = 500, alpha = 0.1)


###################################################
### code chunk number 143: adehabitatHS.Rnw:2265-2268
###################################################
tav <- as.vector(tav)
names(tav) <- names(tus)
(WiII <- widesII(tus, tav))


###################################################
### code chunk number 144: adehabitatHS.Rnw:2309-2310
###################################################
(eis <- eisera(tus,tav2, scannf=FALSE))


###################################################
### code chunk number 145: ploteiser (eval = FALSE)
###################################################
## barplot(eis$eig)


###################################################
### code chunk number 146: adehabitatHS.Rnw:2320-2321 (eval = FALSE)
###################################################
## barplot(eis$eig)


###################################################
### code chunk number 147: adehabitatHS.Rnw:2325-2328
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
barplot(eis$eig)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 148: scateis (eval = FALSE)
###################################################
## scatter(eis)


###################################################
### code chunk number 149: adehabitatHS.Rnw:2338-2339 (eval = FALSE)
###################################################
## scatter(eis)


###################################################
### code chunk number 150: adehabitatHS.Rnw:2343-2346
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
scatter(eis)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 151: adehabitatHS.Rnw:2365-2367
###################################################
avdis <- acm.disjonctif(data.frame(av))
head(avdis)


###################################################
### code chunk number 152: adehabitatHS.Rnw:2372-2374
###################################################
pc <- dudi.pca(avdis, scannf=FALSE)
com2 <- canomi(pc, U, scannf=FALSE)


###################################################
### code chunk number 153: adehabitatHS.Rnw:2380-2381
###################################################
eis$eig/com2$eig


###################################################
### code chunk number 154: adehabitatHS.Rnw:2391-2392
###################################################
eis$li[,1]/com2$li[,1]


###################################################
### code chunk number 155: plotnidesIII (eval = FALSE)
###################################################
## par(mar=c(0.1,0.1,0.1,0.1))
## 
## set.seed(2351)
## plot(1,1, asp=1, bg="grey", pch=21, axes=FALSE, ty="n",
##      xlim=c(-3.5,3.5), ylim=c(-3.5,3.5))
## xx1 <- cbind(rnorm(100, mean=-1.5, sd=0.25),
##              rnorm(100, mean=1.5, sd=0.25))
## points(xx1, bg="grey", col="blue", pch=21)
## set.seed(2)
## xx1b <- cbind(rnorm(20, mean=-1.3, sd=0.1),
##               rnorm(20, mean=1.8, sd=0.1))
## points(xx1b, bg="blue", col="black", pch=21, cex=1.5)
## arrows(-1.5, 1.5, -1.3, 1.8, lwd=5, length=0.1, angle=20)
## arrows(-1.5, 1.5, -1.3, 1.8, lwd=3, col="red", length=0.1, angle=20)
## 
## 
## set.seed(240)
## xx2 <- cbind(rnorm(100, mean=1.5, sd=0.25),
##              rnorm(100, mean=1.5, sd=0.25))
## points(xx2, bg="grey", col="green", pch=21)
## set.seed(2488)
## xx2b <- cbind(rnorm(20, mean=1.7, sd=0.1),
##               rnorm(20, mean=1.8, sd=0.1))
## points(xx2b, bg="green", col="black", pch=21, cex=1.5)
## arrows(1.5, 1.5, 1.7, 1.8, lwd=5, length=0.1, angle=20)
## arrows(1.5, 1.5, 1.7, 1.8, lwd=3, col="red", length=0.1, angle=20)
## 
## 
## set.seed(240987)
## xx2 <- cbind(rnorm(100, mean=0, sd=0.25),
##              rnorm(100, mean=-1.5, sd=0.25))
## points(xx2, bg="grey", col="orange", pch=21)
## set.seed(2488980)
## xx2b <- cbind(rnorm(20, mean=0.3, sd=0.1),
##               rnorm(20, mean=-1.6, sd=0.1))
## points(xx2b, bg="orange", col="black", pch=21, cex=1.5)
## arrows(0, -1.5, 0.3, -1.6, lwd=5, length=0.1, angle=20)
## arrows(0, -1.5, 0.3, -1.6, lwd=3, col="red", length=0.1, angle=20)
## 
## arrows(0,0,0,3)
## arrows(0,0,-sqrt(4.5),-sqrt(4.5))
## arrows(0,0,sqrt(4.5),-sqrt(4.5))
## 
## text(0,3.2, "V3")
## text(-sqrt(4.5)-0.2,-sqrt(4.5)-0.2, "V1")
## text(sqrt(4.5)+0.2,-sqrt(4.5)-0.2, "V2")
## 


###################################################
### code chunk number 156: adehabitatHS.Rnw:2460-2463
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
par(mar=c(0.1,0.1,0.1,0.1))

set.seed(2351)
plot(1,1, asp=1, bg="grey", pch=21, axes=FALSE, ty="n",
     xlim=c(-3.5,3.5), ylim=c(-3.5,3.5))
xx1 <- cbind(rnorm(100, mean=-1.5, sd=0.25),
             rnorm(100, mean=1.5, sd=0.25))
points(xx1, bg="grey", col="blue", pch=21)
set.seed(2)
xx1b <- cbind(rnorm(20, mean=-1.3, sd=0.1),
              rnorm(20, mean=1.8, sd=0.1))
points(xx1b, bg="blue", col="black", pch=21, cex=1.5)
arrows(-1.5, 1.5, -1.3, 1.8, lwd=5, length=0.1, angle=20)
arrows(-1.5, 1.5, -1.3, 1.8, lwd=3, col="red", length=0.1, angle=20)


set.seed(240)
xx2 <- cbind(rnorm(100, mean=1.5, sd=0.25),
             rnorm(100, mean=1.5, sd=0.25))
points(xx2, bg="grey", col="green", pch=21)
set.seed(2488)
xx2b <- cbind(rnorm(20, mean=1.7, sd=0.1),
              rnorm(20, mean=1.8, sd=0.1))
points(xx2b, bg="green", col="black", pch=21, cex=1.5)
arrows(1.5, 1.5, 1.7, 1.8, lwd=5, length=0.1, angle=20)
arrows(1.5, 1.5, 1.7, 1.8, lwd=3, col="red", length=0.1, angle=20)


set.seed(240987)
xx2 <- cbind(rnorm(100, mean=0, sd=0.25),
             rnorm(100, mean=-1.5, sd=0.25))
points(xx2, bg="grey", col="orange", pch=21)
set.seed(2488980)
xx2b <- cbind(rnorm(20, mean=0.3, sd=0.1),
              rnorm(20, mean=-1.6, sd=0.1))
points(xx2b, bg="orange", col="black", pch=21, cex=1.5)
arrows(0, -1.5, 0.3, -1.6, lwd=5, length=0.1, angle=20)
arrows(0, -1.5, 0.3, -1.6, lwd=3, col="red", length=0.1, angle=20)

arrows(0,0,0,3)
arrows(0,0,-sqrt(4.5),-sqrt(4.5))
arrows(0,0,sqrt(4.5),-sqrt(4.5))

text(0,3.2, "V3")
text(-sqrt(4.5)-0.2,-sqrt(4.5)-0.2, "V1")
text(sqrt(4.5)+0.2,-sqrt(4.5)-0.2, "V2")

dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 157: pplomcp (eval = FALSE)
###################################################
## pcc <- mcp(locs)
## image(maps)
## plot(pcc, col=rainbow(6), add=TRUE)


###################################################
### code chunk number 158: adehabitatHS.Rnw:2485-2486 (eval = FALSE)
###################################################
## pcc <- mcp(locs)
## image(maps)
## plot(pcc, col=rainbow(6), add=TRUE)


###################################################
### code chunk number 159: adehabitatHS.Rnw:2490-2493
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=20)
pcc <- mcp(locs)
image(maps)
plot(pcc, col=rainbow(6), add=TRUE)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 160: plohr (eval = FALSE)
###################################################
## hr <- do.call("data.frame", lapply(1:nrow(pcc), function(i) {
##     over(maps,geometry(pcc[i,]))
## }))
## names(hr) <- slot(pcc, "data")[,1]
## coordinates(hr) <- coordinates(maps)
## gridded(hr) <- TRUE
## mimage(hr)


###################################################
### code chunk number 161: adehabitatHS.Rnw:2509-2510 (eval = FALSE)
###################################################
## hr <- do.call("data.frame", lapply(1:nrow(pcc), function(i) {
##     over(maps,geometry(pcc[i,]))
## }))
## names(hr) <- slot(pcc, "data")[,1]
## coordinates(hr) <- coordinates(maps)
## gridded(hr) <- TRUE
## mimage(hr)


###################################################
### code chunk number 162: adehabitatHS.Rnw:2514-2517
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
hr <- do.call("data.frame", lapply(1:nrow(pcc), function(i) {
    over(maps,geometry(pcc[i,]))
}))
names(hr) <- slot(pcc, "data")[,1]
coordinates(hr) <- coordinates(maps)
gridded(hr) <- TRUE
mimage(hr)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 163: adehabitatHS.Rnw:2548-2550
###################################################
pks <- prepksel(maps, hr, cp)
names(pks)


###################################################
### code chunk number 164: adehabitatHS.Rnw:2564-2565
###################################################
pc <- dudi.pca(pks$tab, scannf=FALSE)


###################################################
### code chunk number 165: plotksel (eval = FALSE)
###################################################
## ksel <- kselect(pc, pks$factor, pks$weight, scann=FALSE)
## plot(ksel)


###################################################
### code chunk number 166: adehabitatHS.Rnw:2575-2576 (eval = FALSE)
###################################################
## ksel <- kselect(pc, pks$factor, pks$weight, scann=FALSE)
## plot(ksel)


###################################################
### code chunk number 167: adehabitatHS.Rnw:2580-2583
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt, res=reso)
ksel <- kselect(pc, pks$factor, pks$weight, scann=FALSE)
plot(ksel)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


