## ----eval=TRUE----------------------------------------------------------------
library(phangorn)
data(Laurasiatherian)
data(yeast)

## ----eval=TRUE----------------------------------------------------------------
set.seed(1)
bs <- bootstrap.phyDat(yeast, FUN = function(x)nj(dist.hamming(x)), 
    bs=100)
tree <- nj(dist.hamming(yeast))
par("mar" = rep(1, 4))
tree <- plotBS(tree, bs, "phylogram")
cnet <- consensusNet(bs, .3)
plot(cnet, show.edge.label=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  plot(cnet, "3D")
#  # rotate 3d plot
#  play3d(spin3d(axis=c(0,1,0), rpm=6), duration=10)
#  # create animated gif file
#  movie3d(spin3d(axis=c(0,1,0), rpm=6), duration=10)

## ----eval=TRUE----------------------------------------------------------------
dm <- dist.hamming(yeast)
nnet <- neighborNet(dm)
par("mar" = rep(1, 4))
plot(nnet)

## ----eval=TRUE----------------------------------------------------------------
nnet <- addConfidences(nnet, tree)
par("mar" = rep(1, 4))
plot(nnet, show.edge.label=TRUE)

## ----eval=TRUE----------------------------------------------------------------
tree2 <- rNNI(tree, 2)
tree2 <- addConfidences(tree2, tree)
# several support values are missing
par("mar" = rep(1, 4))
plot(tree2, show.node.label=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  cnet <- nnls.networx(cnet, dm)
#  par("mar" = rep(1, 4))
#  plot(cnet, show.edge.label=TRUE)

## ----eval=TRUE----------------------------------------------------------------
sessionInfo()

