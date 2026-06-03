### R code from vignette source 'intro-vegan.Rnw'

###################################################
### code chunk number 1: intro-vegan.Rnw:18-23
###################################################
par(mfrow=c(1,1))
options(width=72)
figset <- function() par(mar=c(4,4,1,1)+.1)
options(SweaveHooks = list(fig = figset))
options("prompt" = "> ", "continue" = "  ")


###################################################
### code chunk number 2: intro-vegan.Rnw:74-77
###################################################
library(vegan)
data(dune)
ord <- decorana(dune)


###################################################
### code chunk number 3: intro-vegan.Rnw:80-81
###################################################
ord


###################################################
### code chunk number 4: intro-vegan.Rnw:102-104
###################################################
ord <- metaMDS(dune, trace = FALSE)
ord


###################################################
### code chunk number 5: a
###################################################
plot(ord)


###################################################
### code chunk number 6: intro-vegan.Rnw:119-120
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(ord)


###################################################
### code chunk number 7: a
###################################################
plot(ord, type = "n") |>
    points("sites", cex = 0.8, pch=21, col="red", bg="yellow") |>
    text("species", cex=0.7, col="blue")


###################################################
### code chunk number 8: intro-vegan.Rnw:143-144
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(ord, type = "n") |>
    points("sites", cex = 0.8, pch=21, col="red", bg="yellow") |>
    text("species", cex=0.7, col="blue")


###################################################
### code chunk number 9: intro-vegan.Rnw:233-235
###################################################
data(dune.env)
attach(dune.env)


###################################################
### code chunk number 10: a
###################################################
plot(ord, disp="sites", type="n")
ordihull(ord, Management, col=1:4, lwd=3)
ordiellipse(ord, Management, col=1:4, kind = "ehull", lwd=3)
ordiellipse(ord, Management, col=1:4, draw="polygon")
ordispider(ord, Management, col=1:4, label = TRUE)
points(ord, disp="sites", pch=21, col="red", bg="yellow", cex=1.3)


###################################################
### code chunk number 11: intro-vegan.Rnw:246-247
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(ord, disp="sites", type="n")
ordihull(ord, Management, col=1:4, lwd=3)
ordiellipse(ord, Management, col=1:4, kind = "ehull", lwd=3)
ordiellipse(ord, Management, col=1:4, draw="polygon")
ordispider(ord, Management, col=1:4, label = TRUE)
points(ord, disp="sites", pch=21, col="red", bg="yellow", cex=1.3)


###################################################
### code chunk number 12: intro-vegan.Rnw:277-279
###################################################
ord.fit <- envfit(ord ~ A1 + Management, data=dune.env, perm=999)
ord.fit


###################################################
### code chunk number 13: a
###################################################
plot(ord, dis="site")
plot(ord.fit, bg = "yellow")


###################################################
### code chunk number 14: b
###################################################
ordisurf(ord, A1, add=TRUE)


###################################################
### code chunk number 15: intro-vegan.Rnw:295-297
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(ord, dis="site")
plot(ord.fit, bg = "yellow")
ordisurf(ord, A1, add=TRUE)


###################################################
### code chunk number 16: intro-vegan.Rnw:317-319
###################################################
ord <- cca(dune ~ A1 + Management, data=dune.env)
ord


###################################################
### code chunk number 17: a
###################################################
plot(ord, spe.par = list(optimize = TRUE), sit.par = list(type="p"))


###################################################
### code chunk number 18: intro-vegan.Rnw:326-327
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(ord, spe.par = list(optimize = TRUE), sit.par = list(type="p"))


###################################################
### code chunk number 19: intro-vegan.Rnw:347-348
###################################################
cca(dune ~ ., data=dune.env)


###################################################
### code chunk number 20: intro-vegan.Rnw:367-368
###################################################
anova(ord)


###################################################
### code chunk number 21: intro-vegan.Rnw:376-377
###################################################
anova(ord, by="term")


###################################################
### code chunk number 22: intro-vegan.Rnw:382-383
###################################################
anova(ord, by="margin")


###################################################
### code chunk number 23: intro-vegan.Rnw:390-392
###################################################
ord <- cca(dune ~ A1 + Management + Condition(Moisture), data=dune.env)
ord


###################################################
### code chunk number 24: intro-vegan.Rnw:397-398
###################################################
anova(ord, by="term")


###################################################
### code chunk number 25: intro-vegan.Rnw:406-408
###################################################
how <- how(nperm=999, plots = Plots(strata=dune.env$Moisture))
anova(ord, by="term", permutations = how)


###################################################
### code chunk number 26: intro-vegan.Rnw:412-413
###################################################
detach(dune.env)


