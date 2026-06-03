## ----setup, echo=FALSE--------------------------------------------------------
# set global chunk options: images will be bigger
knitr::opts_chunk$set(fig.width=6, fig.height=4)
#, global.par=TRUE
options(digits = 4)

## ----generate data------------------------------------------------------------
library(phangorn)
data <- matrix(c("r","a","y","g","g","a","c","-","c","t","c","g",
    "a","a","t","g","g","a","t","-","c","t","c","a",
    "a","a","t","-","g","a","c","c","c","t","?","g"),
    dimnames = list(c("t1", "t2", "t3"),NULL), nrow=3, byrow=TRUE)
data

## ----dna----------------------------------------------------------------------
gapsdata1 <- phyDat(data)
gapsdata1

## ----5state-------------------------------------------------------------------
gapsdata2 <- phyDat(data, type="USER", levels=c("a","c","g","t","-"),
    ambiguity = c("?", "n"))
gapsdata2

## ----contrast-----------------------------------------------------------------
contrast <- matrix(data = c(1,0,0,0,0,
    0,1,0,0,0,
    0,0,1,0,0,
    0,0,0,1,0,
    1,0,1,0,0,
    0,1,0,1,0,
    0,0,0,0,1,
    1,1,1,1,0,
    1,1,1,1,1),
    ncol = 5, byrow = TRUE)
dimnames(contrast) <- list(c("a","c","g","t","r","y","-","n","?"),
    c("a", "c", "g", "t", "-"))
contrast
gapsdata3 <- phyDat(data, type="USER", contrast=contrast)
gapsdata3

## ----optim.pml subs-----------------------------------------------------------
library(ape)
tree <- unroot(rtree(3))
fit <- pml(tree, gapsdata3)
fit <- optim.pml(fit, optQ=TRUE, subs=c(1,0,1,2,1,0,2,1,2,2),
    control=pml.control(trace=0))
fit

## ----read codon data----------------------------------------------------------
fdir <- system.file("extdata/trees", package = "phangorn")
hiv_2_nef <- read.phyDat(file.path(fdir, "seqfile.txt"), format="sequential")
tree <- read.tree(file.path(fdir, "tree.txt"))

## ----echo=FALSE---------------------------------------------------------------
load("AF.RData")

## ----codonTest, eval=FALSE----------------------------------------------------
#  cdn <- codonTest(tree, hiv_2_nef)
#  cdn

## ----codonTest_cheat, echo=FALSE----------------------------------------------
cdn

## ----plot_codon---------------------------------------------------------------
plot(cdn, "M1a")
plot(cdn, "M2a")

## ----M0-----------------------------------------------------------------------
treeM0 <- cdn$estimates[["M0"]]$tree # tree with edge lengths
M0 <- pml(treeM0, dna2codon(hiv_2_nef), bf="F3x4")
M0 <- optim.pml(M0, model="codon1", control=pml.control(trace=0))
M0

## ----M0+F3x4, eval=FALSE------------------------------------------------------
#  M0_opt <- optim.pml(M0, model="codon1", optBf=TRUE, control=pml.control(trace=0))
#  M0_opt

## ----M0+F3x4_cheat, echo=FALSE------------------------------------------------
M0_opt

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

