## ----setup, echo=FALSE--------------------------------------------------------
# set global chunk options: images will be bigger
knitr::opts_chunk$set(fig.width=8, fig.height=6)
options(digits = 4)

## ----load packages------------------------------------------------------------
library(ape)
library(phangorn)
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")

## ----nj tree------------------------------------------------------------------
dm <- dist.ml(primates)
treeNJ  <- NJ(dm)

## ----pml----------------------------------------------------------------------
fit <- pml(treeNJ, data=primates)
fit

## ----methods pml--------------------------------------------------------------
methods(class="pml")

## ----optim.pml----------------------------------------------------------------
fitJC  <- optim.pml(fit, rearrangement="NNI")
logLik(fitJC)

## ----F81+G+I, cache=TRUE------------------------------------------------------
fitF81 <- update(fitJC, k=4, inv=0.2, bf=baseFreq(primates))
fitF81

## ----GTR+G+I, cache=TRUE------------------------------------------------------
fitGTR <- optim.pml(fitF81, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR

## ----stochastic, cache=TRUE---------------------------------------------------
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR

## ----anova--------------------------------------------------------------------
anova(fitJC, fitGTR)

## ----SH_test------------------------------------------------------------------
SH.test(fitGTR, fitJC)

## ----AIC----------------------------------------------------------------------
AIC(fitJC)
AIC(fitGTR)
AICc(fitGTR)
BIC(fitGTR)

## ----echo=TRUE, eval=TRUE, cache=TRUE-----------------------------------------
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE,
    control = pml.control(trace = 0))

## ----plotBS, fig.cap="Tree with bootstrap support. Unrooted tree (midpoint rooted) with bootstrap support values.", echo=TRUE----
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

## ----ConsensusNet, fig.cap="ConsensusNet from the bootstrap sample.", echo=TRUE, eval=TRUE----
cnet <- consensusNet(bs, p=0.2)
plot(cnet, show.edge.label=TRUE)

## ----allTrees-----------------------------------------------------------------
trees <- allTrees(5)
par(mfrow=c(3,5), mar=rep(0,4))
for(i in 1:15)plot(trees[[i]], cex=1, type="u")

## ----nni----------------------------------------------------------------------
nni(trees[[1]])

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

