### R code from vignette source 'MoranI.Rnw'

###################################################
### code chunk number 1: MoranI.Rnw:19-20
###################################################
options(width = 80, prompt = "> ")


###################################################
### code chunk number 2: MoranI.Rnw:119-123
###################################################
body <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
longevity <- c(4.74493, 3.3322, 3.3673, 2.89037, 2.30259)
names(body) <- names(longevity) <- c("Homo", "Pongo", "Macaca",
                                     "Ateles", "Galago")


###################################################
### code chunk number 3: MoranI.Rnw:129-135
###################################################
library(ape)
trnwk <- "((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles:0.62)"
trnwk[2] <- ":0.38,Galago:1.00);"
tr <- read.tree(text = trnwk)
plot(tr)
axisPhylo()


###################################################
### code chunk number 4: MoranI.Rnw:141-143
###################################################
w <- 1/cophenetic(tr)
w


###################################################
### code chunk number 5: MoranI.Rnw:147-148
###################################################
diag(w) <- 0


###################################################
### code chunk number 6: MoranI.Rnw:152-153
###################################################
Moran.I(body, w)


###################################################
### code chunk number 7: MoranI.Rnw:173-174
###################################################
Moran.I(body, w, alt = "greater")


###################################################
### code chunk number 8: MoranI.Rnw:179-180
###################################################
Moran.I(longevity, w)


###################################################
### code chunk number 9: MoranI.Rnw:241-244
###################################################
data(carnivora)
carnivora$log10SW <- log10(carnivora$SW)
carnivora$log10FW <- log10(carnivora$FW)


###################################################
### code chunk number 10: MoranI.Rnw:248-251
###################################################
fm1.carn <- log10SW ~ Order/SuperFamily/Family/Genus
co1 <- correlogram.formula(fm1.carn, data = carnivora)
plot(co1)


###################################################
### code chunk number 11: MoranI.Rnw:266-269
###################################################
fm2.carn <- log10SW + log10FW ~ Order/SuperFamily/Family/Genus
co2 <- correlogram.formula(fm2.carn, data = carnivora)
print(plot(co2))


###################################################
### code chunk number 12: MoranI.Rnw:277-278
###################################################
plot(co2, FALSE)


