## ----setup, echo=FALSE--------------------------------------------------------
# set global chunk options: images will be bigger
knitr::opts_chunk$set(fig.width=6, fig.height=4)
#, global.par=TRUE
options(digits = 4)
suppressPackageStartupMessages(library(phangorn))

## ----parsimony----------------------------------------------------------------
library(phangorn)
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")
tree <- pratchet(primates, trace=0) |> acctran(primates) |> makeNodeLabel()
parsimony(tree, primates)

## ----ancestral_reconstruction-------------------------------------------------
anc.pars <- anc_pars(tree, primates)

## ----seqLogo, fig.cap="Fig 1. Ancestral reconstruction for a node.", eval=TRUE----
plotSeqLogo(anc.pars, node=getRoot(tree), 1, 20)

## ----MPR, fig.cap="Fig 2. Ancestral reconstruction using MPR."----------------
plotAnc(anc.pars, 17)
title("MPR")

## ----fit_ML-------------------------------------------------------------------
fit <- pml(tree, primates)
fit <- optim.pml(fit, model="F81", control = pml.control(trace=0))

## ----ML_reconstruction--------------------------------------------------------
anc.ml <- anc_pml(fit)

## ----plotML, fig.cap="Fig 4. Ancestral reconstruction the using the maximum likelihood."----
plotAnc(anc.ml, 17)
title("ML")

## ----plotB, fig.cap="Fig 5. Ancestral reconstruction using (empirical) Bayes."----
#plotAnc(anc.bayes, 17)
#title("Bayes")

## ----read_geospiza_data-------------------------------------------------------
data("bird.orders")
x <- c(rep(0, 5), rep(1, 18))
x[c(20,22,23)] <- 2
x <- factor(x)
names(x) <- bird.orders$tip.label
dat <- phyDat(x, "USER", levels=c(0,1,2))

## ----ER_model-----------------------------------------------------------------
fit <- pml(bird.orders, dat)
fit_ER <- optim.pml(fit, optEdge = FALSE, optRate=TRUE, 
                    control = pml.control(trace=0))
fit_ER

## ----SYM_model----------------------------------------------------------------
fit_SYM <- optim.pml(fit, optEdge = FALSE, optRate=TRUE, model="SYM", 
                    control = pml.control(trace=0))
fit_SYM

## ----ace----------------------------------------------------------------------
fit_ace <- ace(x, bird.orders, model="SYM", type = "d")

## ----comparison---------------------------------------------------------------
fit_SYM$logLik
fit_ace$loglik+log(1/3)
all.equal(fit_SYM$logLik, fit_ace$loglik+log(1/3))

## ----SYM_reconstruction-------------------------------------------------------
anc_SYM <- anc_pml(fit_SYM)
plotAnc(anc_SYM)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

