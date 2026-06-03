## ----setup, echo=FALSE--------------------------------------------------------
# set global chunk options: images will be bigger
knitr::opts_chunk$set(fig.width=6, fig.height=4)
options(digits = 2)

## ----load packages------------------------------------------------------------
library(phangorn)
set.seed(9)

## ----load data----------------------------------------------------------------
fdir <- system.file("extdata", package = "phangorn")
mm <- read.csv(file.path(fdir, "mites.csv"), row.names = 1)
mm_pd <- phyDat(as.matrix(mm), type = "USER", levels = 0:7)

## ----write_nexus, eval=FALSE--------------------------------------------------
#  write.phyDat(mm_pd, file.path(fdir, "mites.nex"), format = "nexus")

## ----read_nexus---------------------------------------------------------------
mm_pd <- read.phyDat(file.path(fdir, "mites.nex"), format = "nexus", type = "STANDARD")

## ----contrast_matrix----------------------------------------------------------
contrast <- matrix(data = c(1,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,
    0,0,1,0,0,0,0,0,0,
    0,0,0,1,0,0,0,0,0,
    0,0,0,0,1,0,0,0,0,
    0,0,0,0,0,1,0,0,0,
    0,0,0,0,0,0,1,0,0,
    0,0,0,0,0,0,0,1,0,
    0,0,0,0,0,0,0,0,1,
    1,1,1,1,1,1,1,1,1),
    ncol = 9, byrow = TRUE)
dimnames(contrast) <- list(c(0:7,"-","?"),
    c(0:7, "-"))
contrast
mm_pd <- phyDat(mm_pd, type="USER", contrast=contrast)

## ----random addition----------------------------------------------------------
mm_start <- random.addition(mm_pd)

## ----pratchet, cache=TRUE-----------------------------------------------------
mm_tree <- pratchet(mm_pd, start = mm_start, minit = 1000, maxit = 10000,
                    all = TRUE, trace = 0)
mm_tree

## ----edge lengths-------------------------------------------------------------
mm_tree <- acctran(mm_tree, mm_pd)

## ----bab----------------------------------------------------------------------
mm_bab <- bab(mm_pd, trace = 0)
mm_bab

## ----root trees, message=FALSE------------------------------------------------
mm_tree_rooted <- root(mm_tree, outgroup = "C._cymba", resolve.root = TRUE,
                       edgelabel = TRUE)

## ----plot_trees, eval=FALSE---------------------------------------------------
#  # subsetting for tree nr. 9
#  plotBS(mm_tree_rooted[[9]], digits = 2)
#  
#  # save plot as pdf
#  pdf(file = "mm_rooted.pdf")
#  plotBS(mm_tree_rooted[[9]], digits = 2)
#  dev.off()

## ----consensus tree-----------------------------------------------------------
# unrooted pratchet tree
mm_cons <- consensus(mm_tree)

# rooted pratchet tree
mm_cons_root <- consensus(mm_tree_rooted, rooted = TRUE)

# branch and bound, we root the consensus tree in the same step
mm_bab_cons <- root(consensus(mm_bab), outgroup = "C._cymba",
                    resolve.root = TRUE, edgelabel = TRUE)

## ----plot_cons_tree, fig.cap="Unrooted and rooted consensus trees of the mites dataset with MP.", fig.show="hold", out.width="33%"----
plot(mm_cons, main="Unrooted pratchet consensus tree")
plot(mm_cons_root, main="Rooted pratchet consensus tree")
plot(mm_bab_cons, main="Rooted bab consensus tree")

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

