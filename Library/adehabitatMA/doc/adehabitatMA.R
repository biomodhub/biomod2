### R code from vignette source 'adehabitatMA.Rnw'

###################################################
### code chunk number 1: adehabitatMA.Rnw:28-32
###################################################
oldopt <- options(width=80, warn=-1)
.PngNo <- 0
wi <- 480
pt <- 20


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
### code chunk number 5: adehabitatMA.Rnw:103-104
###################################################
library(adehabitatMA)


###################################################
### code chunk number 6: adehabitatMA.Rnw:107-110
###################################################
suppressWarnings(RNGversion("3.5.0"))
set.seed(13431)
adeoptions(shortprint=TRUE)


###################################################
### code chunk number 7: adehabitatMA.Rnw:130-131
###################################################
data(lynxjura)


###################################################
### code chunk number 8: adehabitatMA.Rnw:143-144
###################################################
head(lynxjura$locs)


###################################################
### code chunk number 9: adehabitatMA.Rnw:151-152
###################################################
lynxjura$map


###################################################
### code chunk number 10: figu1 (eval = FALSE)
###################################################
## mimage(lynxjura$map)


###################################################
### code chunk number 11: adehabitatMA.Rnw:162-163 (eval = FALSE)
###################################################
## mimage(lynxjura$map)


###################################################
### code chunk number 12: adehabitatMA.Rnw:167-170
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
mimage(lynxjura$map)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 13: adehabitatMA.Rnw:187-188 (eval = FALSE)
###################################################
## explore(lynxjura$map)


###################################################
### code chunk number 14: adehabitatMA.Rnw:199-200 (eval = FALSE)
###################################################
## explore(lynxjura$map, panel.last=function() points(lynxjura$locs, pch=3))


###################################################
### code chunk number 15: adehabitatMA.Rnw:210-211
###################################################
map <- lynxjura$map


###################################################
### code chunk number 16: sodssss (eval = FALSE)
###################################################
## hist(map)


###################################################
### code chunk number 17: adehabitatMA.Rnw:220-221 (eval = FALSE)
###################################################
## hist(map)


###################################################
### code chunk number 18: adehabitatMA.Rnw:225-228
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
hist(map)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 19: qsfddfdd (eval = FALSE)
###################################################
## forest <- map[,1]
## forest[[1]][forest[[1]]<95] <- NA
## image(forest, col="green")


###################################################
### code chunk number 20: adehabitatMA.Rnw:243-244 (eval = FALSE)
###################################################
## forest <- map[,1]
## forest[[1]][forest[[1]]<95] <- NA
## image(forest, col="green")


###################################################
### code chunk number 21: adehabitatMA.Rnw:249-252
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
forest <- map[,1]
forest[[1]][forest[[1]]<95] <- NA
image(forest, col="green")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 22: ksdkksss (eval = FALSE)
###################################################
## lab <- labcon(forest)
## image(lab)


###################################################
### code chunk number 23: adehabitatMA.Rnw:266-267 (eval = FALSE)
###################################################
## lab <- labcon(forest)
## image(lab)


###################################################
### code chunk number 24: adehabitatMA.Rnw:271-274
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
lab <- labcon(forest)
image(lab)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 25: adehabitatMA.Rnw:283-285
###################################################
lab
max(lab[[1]])


###################################################
### code chunk number 26: adehabitatMA.Rnw:291-292
###################################################
gridparameters(lab)


###################################################
### code chunk number 27: adehabitatMA.Rnw:297-298
###################################################
table(lab[[1]])*500*500


###################################################
### code chunk number 28: adehabitatMA.Rnw:309-311
###################################################
fullgrid(lab) <- TRUE
fullgrid(map) <- TRUE


###################################################
### code chunk number 29: adehabitatMA.Rnw:316-317
###################################################
mean(map[[2]][lab[[1]]==1], na.rm=TRUE)


###################################################
### code chunk number 30: sldsss (eval = FALSE)
###################################################
## comp1 <- map[2]
## comp1[[1]][map[[1]]<95] <- NA
## comp1[[1]][lab[[1]]!=1] <- NA
## image(comp1)


###################################################
### code chunk number 31: adehabitatMA.Rnw:330-331 (eval = FALSE)
###################################################
## comp1 <- map[2]
## comp1[[1]][map[[1]]<95] <- NA
## comp1[[1]][lab[[1]]!=1] <- NA
## image(comp1)


###################################################
### code chunk number 32: adehabitatMA.Rnw:335-338
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
comp1 <- map[2]
comp1[[1]][map[[1]]<95] <- NA
comp1[[1]][lab[[1]]!=1] <- NA
image(comp1)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 33: skks (eval = FALSE)
###################################################
## image(forest, col="red")


###################################################
### code chunk number 34: adehabitatMA.Rnw:352-353 (eval = FALSE)
###################################################
## image(forest, col="red")


###################################################
### code chunk number 35: adehabitatMA.Rnw:357-360
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(forest, col="red")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 36: adehabitatMA.Rnw:367-368
###################################################
con <- getcontour(forest)


###################################################
### code chunk number 37: ssskkkq (eval = FALSE)
###################################################
## plot(con, col="green")


###################################################
### code chunk number 38: adehabitatMA.Rnw:395-396 (eval = FALSE)
###################################################
## plot(con, col="green")


###################################################
### code chunk number 39: adehabitatMA.Rnw:400-403
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(con, col="green")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 40: adehabitatMA.Rnw:426-428
###################################################
for1 <- morphology(forest, "dilate", nt=1)
for1


###################################################
### code chunk number 41: sdkskss (eval = FALSE)
###################################################
## image(for1, col="blue")
## image(forest, col="yellow", add=TRUE)


###################################################
### code chunk number 42: adehabitatMA.Rnw:440-441 (eval = FALSE)
###################################################
## image(for1, col="blue")
## image(forest, col="yellow", add=TRUE)


###################################################
### code chunk number 43: adehabitatMA.Rnw:445-448
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(for1, col="blue")
image(forest, col="yellow", add=TRUE)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 44: skkssss (eval = FALSE)
###################################################
## plot(getcontour(for1), col="green")


###################################################
### code chunk number 45: adehabitatMA.Rnw:460-461 (eval = FALSE)
###################################################
## plot(getcontour(for1), col="green")


###################################################
### code chunk number 46: adehabitatMA.Rnw:465-468
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(getcontour(for1), col="green")
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 47: skdskdk (eval = FALSE)
###################################################
## map <- lynxjura$map
## mimage(map)


###################################################
### code chunk number 48: adehabitatMA.Rnw:492-493 (eval = FALSE)
###################################################
## map <- lynxjura$map
## mimage(map)


###################################################
### code chunk number 49: adehabitatMA.Rnw:497-500
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
map <- lynxjura$map
mimage(map)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 50: adehabitatMA.Rnw:506-507
###################################################
gridparameters(map)


###################################################
### code chunk number 51: sssshhh (eval = FALSE)
###################################################
## mimage(lowres(map, 10))


###################################################
### code chunk number 52: adehabitatMA.Rnw:520-521 (eval = FALSE)
###################################################
## mimage(lowres(map, 10))


###################################################
### code chunk number 53: adehabitatMA.Rnw:525-528
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
mimage(lowres(map, 10))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 54: sqkdksssq (eval = FALSE)
###################################################
## map[[1]] <- as.numeric(cut(map[[1]],3))
## image(map, 1)


###################################################
### code chunk number 55: adehabitatMA.Rnw:543-544 (eval = FALSE)
###################################################
## map[[1]] <- as.numeric(cut(map[[1]],3))
## image(map, 1)


###################################################
### code chunk number 56: adehabitatMA.Rnw:548-551
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
map[[1]] <- as.numeric(cut(map[[1]],3))
image(map, 1)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 57: sdskkkk (eval = FALSE)
###################################################
## image(lowres(map, 10, which.fac=1))


###################################################
### code chunk number 58: adehabitatMA.Rnw:567-568 (eval = FALSE)
###################################################
## image(lowres(map, 10, which.fac=1))


###################################################
### code chunk number 59: adehabitatMA.Rnw:572-575
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(lowres(map, 10, which.fac=1))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 60: skssjdx (eval = FALSE)
###################################################
## image(forest, col="green")
## box()


###################################################
### code chunk number 61: adehabitatMA.Rnw:590-591 (eval = FALSE)
###################################################
## image(forest, col="green")
## box()


###################################################
### code chunk number 62: adehabitatMA.Rnw:595-598
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(forest, col="green")
box()
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 63: adehabitatMA.Rnw:605-606 (eval = FALSE)
###################################################
## for2 <- subsetmap(forest)


###################################################
### code chunk number 64: adehabitatMA.Rnw:612-614
###################################################
for2 <- subsetmap(forest, xlim=c(850254.2, 878990.2),
                  ylim=c(2128744, 2172175))


###################################################
### code chunk number 65: skksks (eval = FALSE)
###################################################
## image(for2, col="green")
## box()


###################################################
### code chunk number 66: adehabitatMA.Rnw:622-623 (eval = FALSE)
###################################################
## image(for2, col="green")
## box()


###################################################
### code chunk number 67: adehabitatMA.Rnw:627-630
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(for2, col="green")
box()
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 68: adehabitatMA.Rnw:666-671
###################################################
data(lynxjura)
map <- lynxjura$map
class(map)
locs <- lynxjura$loc
class(locs)


###################################################
### code chunk number 69: flkqskqfkjc (eval = FALSE)
###################################################
## image(map, 1)
## points(locs, pch=3)


###################################################
### code chunk number 70: adehabitatMA.Rnw:681-682 (eval = FALSE)
###################################################
## image(map, 1)
## points(locs, pch=3)


###################################################
### code chunk number 71: adehabitatMA.Rnw:686-689
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(map, 1)
points(locs, pch=3)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 72: adehabitatMA.Rnw:695-696
###################################################
cp <- count.points(locs, map)


###################################################
### code chunk number 73: ssckcc (eval = FALSE)
###################################################
## image(cp)


###################################################
### code chunk number 74: adehabitatMA.Rnw:717-718 (eval = FALSE)
###################################################
## image(cp)


###################################################
### code chunk number 75: adehabitatMA.Rnw:722-725
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(cp)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 76: adehabitatMA.Rnw:738-739
###################################################
head(locs)


###################################################
### code chunk number 77: adehabitatMA.Rnw:749-751
###################################################
cpr <- count.points(locs[,"Type"], map)
cpr


###################################################
### code chunk number 78: sdkkckck (eval = FALSE)
###################################################
## mimage(cpr)


###################################################
### code chunk number 79: adehabitatMA.Rnw:760-761 (eval = FALSE)
###################################################
## mimage(cpr)


###################################################
### code chunk number 80: adehabitatMA.Rnw:765-768
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
mimage(cpr)
dev.null <- dev.off()
cat("\\includegraphics[height=14cm,width=14cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 81: adehabitatMA.Rnw:782-783
###################################################
df <- join(locs, map)


###################################################
### code chunk number 82: adehabitatMA.Rnw:789-790
###################################################
head(df)


###################################################
### code chunk number 83: adehabitatMA.Rnw:809-811
###################################################
asc <- ascgen(locs, cellsize=5000)
asc


###################################################
### code chunk number 84: llldcvdv (eval = FALSE)
###################################################
## image(asc)


###################################################
### code chunk number 85: adehabitatMA.Rnw:822-823 (eval = FALSE)
###################################################
## image(asc)


###################################################
### code chunk number 86: adehabitatMA.Rnw:827-830
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(asc)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 87: adehabitatMA.Rnw:843-844
###################################################
po <- locs[locs[["Type"]]=="O",]


###################################################
### code chunk number 88: sfjkfc (eval = FALSE)
###################################################
## image(buffer(po, map, 3000))


###################################################
### code chunk number 89: adehabitatMA.Rnw:854-855 (eval = FALSE)
###################################################
## image(buffer(po, map, 3000))


###################################################
### code chunk number 90: adehabitatMA.Rnw:859-862
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(buffer(po, map, 3000))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 91: sdskwkl (eval = FALSE)
###################################################
## plot(con)


###################################################
### code chunk number 92: adehabitatMA.Rnw:873-874 (eval = FALSE)
###################################################
## plot(con)


###################################################
### code chunk number 93: adehabitatMA.Rnw:878-881
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
plot(con)
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 94: dfkwkfckj (eval = FALSE)
###################################################
## image(buffer(con, map, 3000))


###################################################
### code chunk number 95: adehabitatMA.Rnw:892-893 (eval = FALSE)
###################################################
## image(buffer(con, map, 3000))


###################################################
### code chunk number 96: adehabitatMA.Rnw:897-900
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
image(buffer(con, map, 3000))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 97: qksdk (eval = FALSE)
###################################################
## sl <- as(con, "SpatialLines")
## image(buffer(sl, map, 500))


###################################################
### code chunk number 98: adehabitatMA.Rnw:914-915 (eval = FALSE)
###################################################
## sl <- as(con, "SpatialLines")
## image(buffer(sl, map, 500))


###################################################
### code chunk number 99: adehabitatMA.Rnw:919-922
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-",
          .PngNo, ".png", sep="")
png(file=file, width = wi, height = wi, pointsize = pt)
sl <- as(con, "SpatialLines")
image(buffer(sl, map, 500))
dev.null <- dev.off()
cat("\\includegraphics[height=7cm,width=7cm]{", file, "}\n\n", sep="")


###################################################
### code chunk number 100: adehabitatMA.Rnw:959-960
###################################################
options(oldopt)


