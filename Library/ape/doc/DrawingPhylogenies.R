### R code from vignette source 'DrawingPhylogenies.Rnw'

###################################################
### code chunk number 1: DrawingPhylogenies.Rnw:16-17
###################################################
options(width = 80, prompt = "> ")


###################################################
### code chunk number 2: DrawingPhylogenies.Rnw:61-63
###################################################
library(ape)
mytr <- read.tree(text = "((Pan:5,Homo:5):2,Gorilla:7);")


###################################################
### code chunk number 3: DrawingPhylogenies.Rnw:69-75
###################################################
foo <- function() {
    col <- "green"
    for (i in 1:2)
        axis(i, col = col, col.ticks = col, col.axis = col, las = 1)
    box(lty = "19")
}


###################################################
### code chunk number 4: DrawingPhylogenies.Rnw:82-89
###################################################
layout(matrix(1:4, 2, 2, byrow = TRUE))
plot(mytr); foo()
plot(mytr, "c", FALSE); foo()
plot(mytr, "u"); foo()
par(xpd = TRUE)
plot(mytr, "f"); foo()
box("outer")


###################################################
### code chunk number 5: DrawingPhylogenies.Rnw:135-137
###################################################
tr <- compute.brlen(stree(8, "l"), 0.1)
tr$tip.label[] <- ""


###################################################
### code chunk number 6: DrawingPhylogenies.Rnw:144-150
###################################################
foo <- function() {
    col <- "green"
    axis(1, col = col, col.ticks = col, col.axis = col)
    axis(2, col = col, col.ticks = col, col.axis = col, at = 1:Ntip(tr), las = 1)
    box(lty = "19")
}


###################################################
### code chunk number 7: DrawingPhylogenies.Rnw:152-152
###################################################



###################################################
### code chunk number 8: DrawingPhylogenies.Rnw:155-165
###################################################
layout(matrix(1:12, 6, 2))
par(mar = c(2, 2, 0.3, 0))
for (type in c("p", "c")) {
   plot(tr, type); foo()
   plot(tr, type, node.pos = 2); foo()
   plot(tr, type, FALSE); foo()
   plot(tr, type, FALSE, node.pos = 1, node.depth = 2); foo()
   plot(tr, type, FALSE, node.pos = 2); foo()
   plot(tr, type, FALSE, node.pos = 2, node.depth = 2); foo()
}


###################################################
### code chunk number 9: DrawingPhylogenies.Rnw:173-178
###################################################
foo <- function() {
    col <- "green"
    for (i in 1:2) axis(i, col = col, col.ticks = col, col.axis = col, las = 1)
    box(lty = "19")
}


###################################################
### code chunk number 10: DrawingPhylogenies.Rnw:181-187
###################################################
layout(matrix(1:4, 2, 2))
par(las = 1)
plot(tr, "u"); foo()
plot(tr, "u", FALSE); foo()
plot(tr, "f"); foo()
plot(tr, "f", FALSE); foo()


###################################################
### code chunk number 11: DrawingPhylogenies.Rnw:257-259
###################################################
args(node.depth.edgelength)
args(node.depth)


###################################################
### code chunk number 12: DrawingPhylogenies.Rnw:271-272
###################################################
args(node.height)


###################################################
### code chunk number 13: DrawingPhylogenies.Rnw:288-289
###################################################
args(unrooted.xy)


###################################################
### code chunk number 14: DrawingPhylogenies.Rnw:296-299
###################################################
args(phylogram.plot)
args(cladogram.plot)
args(circular.plot)


###################################################
### code chunk number 15: DrawingPhylogenies.Rnw:315-316
###################################################
args(plot.phylo)


###################################################
### code chunk number 16: DrawingPhylogenies.Rnw:369-370
###################################################
geo <- factor(c("Africa", "World", "Africa"))


###################################################
### code chunk number 17: DrawingPhylogenies.Rnw:378-380
###################################################
(mycol <- c("blue", "red")[geo])
plot(mytr, tip.color = mycol)


###################################################
### code chunk number 18: DrawingPhylogenies.Rnw:392-394
###################################################
par(xpd = TRUE)
plot(mytr, tip.color = mycol, cex = c(1, 1, 1.5))


###################################################
### code chunk number 19: DrawingPhylogenies.Rnw:401-402
###################################################
args(def)


###################################################
### code chunk number 20: DrawingPhylogenies.Rnw:413-415
###################################################
mycol2 <- def(mytr$tip.label, Homo = "red", default = "blue")
identical(mycol, mycol2)


###################################################
### code chunk number 21: DrawingPhylogenies.Rnw:436-442
###################################################
(i <- which.edge(mytr, c("Homo", "Pan")))
co <- rep("black", Nedge(mytr))
co[i] <- "blue"
par(xpd = TRUE)
plot(mytr, edge.col = co)
edgelabels()


###################################################
### code chunk number 22: DrawingPhylogenies.Rnw:480-481
###################################################
args(axisPhylo)


###################################################
### code chunk number 23: DrawingPhylogenies.Rnw:493-494
###################################################
args(add.scale.bar)


###################################################
### code chunk number 24: DrawingPhylogenies.Rnw:506-510
###################################################
plot(mytr)
add.scale.bar()
axisPhylo()
axis(3)


###################################################
### code chunk number 25: DrawingPhylogenies.Rnw:520-523
###################################################
args(nodelabels)
args(tiplabels)
args(edgelabels)


###################################################
### code chunk number 26: DrawingPhylogenies.Rnw:535-540
###################################################
par(xpd = TRUE)
plot(mytr)
nodelabels()
tiplabels()
edgelabels()


###################################################
### code chunk number 27: DrawingPhylogenies.Rnw:563-566
###################################################
plot(rep(1:5, 5), rep(1:5, each = 5), pch = 1:25, xlim = c(1, 5.2),
     col = "blue", bg = "yellow", cex = 2)
text(rep(1:5, 5) + 0.2, rep(1:5, each = 5), 1:25)


###################################################
### code chunk number 28: DrawingPhylogenies.Rnw:583-593
###################################################
v <- c(0, 0.5, 1)
layout(matrix(1:9, 3, 3))
par(mar = c(3, 3, 3, 0), las = 1)
for (i in v) {
    for (j in v) {
        plot(0, 0, "n", main = paste0("adj = c(", i, ", ", j, ")"))
        abline(v = 0, h = 0, lty = 3)
        text(0, 0, "Gorilla", adj = c(i, j))
    }
}


###################################################
### code chunk number 29: DrawingPhylogenies.Rnw:601-605
###################################################
plot(rep(0, 3), 0:2, "n", las = 1)
abline(v = 0, h = 0:2, lty = 3)
text(0, 0:2, mytr$tip.label, adj = -1)
text(0, 0:2, mytr$tip.label, adj = 2, font = 2)


###################################################
### code chunk number 30: DrawingPhylogenies.Rnw:623-625
###################################################
args(phydataplot)
args(ring)


###################################################
### code chunk number 31: DrawingPhylogenies.Rnw:677-679
###################################################
x <- matrix(1:6, 3, 10)
dimnames(x) <- list(c("Homo", "Gorilla", "Pan"), LETTERS[1:ncol(x)])


###################################################
### code chunk number 32: DrawingPhylogenies.Rnw:685-687
###################################################
x
mytr$tip.label


###################################################
### code chunk number 33: DrawingPhylogenies.Rnw:699-705
###################################################
par(mar = c(10, 2, 5, 5))
plot(mytr, x.lim = 30)
phydataplot(x, mytr, "m", border = "white", offset = 2, width = 1)
phydataplot(x, mytr, "m", border = "white", offset = 15,
            width = 1, continuous = 2)



###################################################
### code chunk number 34: DrawingPhylogenies.Rnw:718-728
###################################################
n <- 100
p <- 3
set.seed(3)
tr <- rcoal(n)
x <- matrix(runif(p * n), n, p)
rownames(x) <- tr$tip.label
COL <- c("red", "green", "blue")
par(xpd = TRUE)
plot(tr, "f", x.lim = c(-5, 5), show.tip.label = FALSE, no.margin = TRUE)
for (j in 1:p) ring(x[, j], tr, offset = j - 1 + 0.1, col = COL[j])


###################################################
### code chunk number 35: DrawingPhylogenies.Rnw:754-759
###################################################
tree <- rtree(10)
tree$edge.length[1] <- 1000
layout(matrix(1:2, 1))
plot(tree); axisPhylo()
plotBreakLongEdges(tree); axisPhylo()


###################################################
### code chunk number 36: DrawingPhylogenies.Rnw:767-769
###################################################
el <- tree$edge.length
sum(el > mean(el) + 2 * sd(el))


###################################################
### code chunk number 37: DrawingPhylogenies.Rnw:774-780
###################################################
tree$edge.length[8] <- 1000
el <- tree$edge.length
sum(el > mean(el) + 2 * sd(el))
layout(matrix(1:2, 1))
plotBreakLongEdges(tree)
plotBreakLongEdges(tree, n = 2)


###################################################
### code chunk number 38: DrawingPhylogenies.Rnw:804-806
###################################################
TR <- replicate(10, rcoal(sample(11:20, size = 1)), simplify = FALSE)
kronoviz(TR, type = "c", show.tip.label = FALSE)


###################################################
### code chunk number 39: DrawingPhylogenies.Rnw:817-820
###################################################
dates <- as.Date(.leap.seconds)
tr <- rcoal(length(dates))
plotTreeTime(tr, dates)


###################################################
### code chunk number 40: DrawingPhylogenies.Rnw:833-840
###################################################
par(oma = rep(2, 4))
plot(mytr)
box(lty = 3)
box("figure", col = "blue")
for (i in 1:4)
    mtext("Outer margin", i, 0.5, outer = TRUE, font = 2, col = "blue")
box("outer")


###################################################
### code chunk number 41: DrawingPhylogenies.Rnw:856-867
###################################################
layout(matrix(1:4, 2, 2))
par(oma = rep(2, 4))
for (i in 1:4) {
    par(mar = rep(i, 4))
    plot(mytr)
    box(lty = 3)
    box("figure", col = "blue")
}
for (i in 1:4)
    mtext("Outer margin", i, 0.5, outer = TRUE, font = 2, col = "blue")
box("outer")


###################################################
### code chunk number 42: DrawingPhylogenies.Rnw:885-887
###################################################
layout(matrix(c(1, 2, 1, 3), 2, 2))
for (i in 1:3) plot(rtree(5))


###################################################
### code chunk number 43: DrawingPhylogenies.Rnw:906-907
###################################################
(pp <- get("last_plot.phylo", envir = .PlotPhyloEnv))


###################################################
### code chunk number 44: DrawingPhylogenies.Rnw:972-1029
###################################################
## split the device into 2:
layout(matrix(1:2, 2))
## plot the 1st tree
plot(rtree(4))
## get the parameters of this first plot:
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
## keep the coordinates of the first tip:
x <- pp$xx[1]
y <- pp$yy[1]
## extremes of the coordinates on both axes:
(pu <- par("usr"))
## fraction of the plot region:
(x - pu[1]) / (pu[2] - pu[1])
(y - pu[3]) / (pu[4] - pu[3])
## get the dimensions of plotting region and margins in inches:
pi <- par("pin")
mi <- par("mai")
## convert the coordinates into inches:
xi1 <- (x - pu[1]) / (pu[2] - pu[1]) * pi[1] + mi[2]
yi1 <- (y - pu[3]) / (pu[4] - pu[3]) * pi[2] + mi[1]

## xi1 and yi1 are the final coordinates of this tip in inches

## plot the 2nd tree:
plot(rtree(4))
## same as above:
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
## ... except we take the coordinates of the root:
x <- pp$xx[5]
y <- pp$yy[5]
pu <- par("usr")
pi <- par("pin")
mi <- par("mai")
xi2 <- (x - pu[1]) / (pu[2] - pu[1]) * pi[1] + mi[2]
yi2 <- (y - pu[3]) / (pu[4] - pu[3]) * pi[2] + mi[1]

## xi2 and yi2 are the final coordinates of this root in inches

## we add the height of this second plot to the 'y' coordinate
## of the first tip of the first tree which is above the second
## one according to layout()
yi1 <- yi1 + par("fin")[2]
## => this operation depends on the specific layout of plots

## get the dimension of the device in inches:
di <- par("din")

## reset the layout
layout(1)
## set new = TRUE and the margins to zero:
par(new = TRUE, mai = rep(0, 4))
## set the scales to be [0,1] on both axes (in user coordinates):
plot(NA, type = "n", ann = FALSE, axes = FALSE, xlim = 0:1,
     ylim = 0:1, xaxs = "i", yaxs = "i")
## graphical elements can now be added:
fancyarrows(xi1/di[1], yi1/di[2], xi2/di[1], yi2/di[2], 1,
            lwd = 10, col = rgb(1, .5, 0, .5), type = "h")


###################################################
### code chunk number 45: DrawingPhylogenies.Rnw:1039-1074
###################################################
foo <- function(phy1, phy2, from)
{
    layout(matrix(1:2, 2))
    plot(phy1, font = 1)
    pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    from <- which(phy1$tip.label == from)
    x <- pp$xx[from]
    y <- pp$yy[from]
    ## fraction of the plot region:
    pu <- par("usr")
    ## convert into inches:
    pi <- par("pin")
    mi <- par("mai")
    xi1 <- (x - pu[1]) / (pu[2] - pu[1]) * pi[1] + mi[2]
    yi1 <- (y - pu[3]) / (pu[4] - pu[3]) * pi[2] + mi[1]
    plot(phy2)
    pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    to <- Ntip(phy2) + 1
    x <- pp$xx[to]
    y <- pp$yy[to]
    ## same as above:
    pu <- par("usr")
    pi <- par("pin")
    mi <- par("mai")
    xi2 <- (x - pu[1]) / (pu[2] - pu[1]) * pi[1] + mi[2]
    yi2 <- (y - pu[3]) / (pu[4] - pu[3]) * pi[2] + mi[1]
    yi1 <- yi1 + par("fin")[2]
    di <- par("din")
    layout(1)
    par(new = TRUE, mai = rep(0, 4))
    plot(NA, type = "n", ann = FALSE, axes = FALSE, xlim = 0:1,
         ylim = 0:1, xaxs = "i", yaxs = "i")
    fancyarrows(xi1/di[1], yi1/di[2], xi2/di[1], yi2/di[2], 1,
                lwd = 10, col = rgb(1, .5, 0, .5), type = "h")
}


###################################################
### code chunk number 46: DrawingPhylogenies.Rnw:1080-1083
###################################################
trb <- read.tree(text = "((Primates,Carnivora),Monotremata);")
par(xpd = TRUE)
foo(trb, mytr, "Primates")


###################################################
### code chunk number 47: DrawingPhylogenies.Rnw:1100-1106
###################################################
plot(mytr, plot = FALSE)
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
rect(pp$xx[5] - 0.1, pp$yy[1] - 0.1, pp$xx[1] + 2, pp$yy[2] + 0.1,
     col = "yellow", border = NA)
par(new = TRUE)
plot(mytr)


