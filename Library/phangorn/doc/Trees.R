## ----setup, echo=FALSE--------------------------------------------------------
# set global chunk options: images will be bigger
knitr::opts_chunk$set(fig.width=8, fig.height=6)

## ----load_packages------------------------------------------------------------
library(ape)
library(phangorn)
fdir <- system.file("extdata/trees", package = "phangorn")
primates <- read.phyDat(file.path(fdir, "primates.dna"),
                        format = "interleaved")

## ----distance_calculation-----------------------------------------------------
dm  <- dist.ml(primates)
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)

## ----plot1, fig.cap="Rooted UPGMA tree.", echo=TRUE---------------------------
plot(treeUPGMA, main="UPGMA")

## ----plot2, fig.cap="Unrooted NJ tree.", echo=TRUE----------------------------
plot(treeNJ, "unrooted", main="NJ")

## ----bootstrap_dist, echo=TRUE------------------------------------------------
fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(primates,  fun)

## ----bootstrap_dist_new, echo=TRUE, eval=FALSE, cache=TRUE--------------------
#  bs_upgma <- bootstrap.phyDat(primates,  \(x){dist.ml(x) |> upgma})

## ----plot_bs, fig.cap="Rooted UPGMA tree.", echo=TRUE-------------------------
plotBS(treeUPGMA, bs_upgma, main="UPGMA")

## ----pars_calc----------------------------------------------------------------
parsimony(treeUPGMA, primates)
parsimony(treeNJ, primates)

## ----pratchet-----------------------------------------------------------------
treeRatchet  <- pratchet(primates, trace = 0, minit=100)
parsimony(treeRatchet, primates)

## ----acctran------------------------------------------------------------------
treeRatchet  <- acctran(treeRatchet, primates)

## ----di2multi-----------------------------------------------------------------
treeRatchet  <- di2multi(treeRatchet)

## ----unique_trees-------------------------------------------------------------
if(inherits(treeRatchet, "multiPhylo")){
  treeRatchet <- unique(treeRatchet)
}

## ----midpoint-----------------------------------------------------------------
plotBS(midpoint(treeRatchet), type="phylogram")
add.scale.bar()

## ----random_addition----------------------------------------------------------
treeRA <- random.addition(primates)
treeSPR  <- optim.parsimony(treeRA, primates)
parsimony(c(treeRA, treeSPR), primates)

## ----bab----------------------------------------------------------------------
(trees <- bab(primates[1:10,]))

## ----mt, echo=TRUE, eval=FALSE------------------------------------------------
#  mt <- modelTest(primates)

## ----echo=FALSE---------------------------------------------------------------
load("Trees.RData")

## ----mt_selected, echo=TRUE, eval=FALSE---------------------------------------
#  mt <- modelTest(primates, model=c("JC", "F81", "K80", "HKY", "SYM", "GTR"),
#                  control = pml.control(trace = 0))

## ----echo=FALSE---------------------------------------------------------------
library(knitr)
kable(mt, digits=2)

## ----as.pml, echo=TRUE--------------------------------------------------------
fit <- as.pml(mt, "HKY+G(4)+I")
fit <- as.pml(mt, "BIC")

## ----pml_bb_modelTest, cache=TRUE---------------------------------------------
fit_mt <- pml_bb(mt, control = pml.control(trace = 0))
fit_mt

## ----pml_bb GTR, eval=FALSE---------------------------------------------------
#  fitGTR <- pml_bb(primates, model="GTR+G(4)+I")

## ----bs, echo=TRUE, eval=FALSE------------------------------------------------
#  bs <- bootstrap.pml(fit_mt, bs=100, optNni=TRUE,
#      control = pml.control(trace = 0))

## ----plotBS_ultrafast_bs, fig.cap="Unrooted tree (midpoint rooted) with ultrafast, standard and transfer bootstrap support values.", echo=TRUE, fig.show="hold", out.width="33%"----
plotBS(midpoint(fit_mt$tree), p = .5, type="p", digits=2, main="Ultrafast bootstrap")

plotBS(midpoint(fit_mt$tree), bs, p = 50, type="p", main="Standard bootstrap")

plotBS(midpoint(fit_mt$tree), bs, p = 50, type="p", digits=0, method = "TBE", main="Transfer bootstrap")

## ----assign_bs_values, eval=FALSE---------------------------------------------
#  # assigning standard bootstrap values to our tree; this is the default method
#  tree_stdbs <- plotBS(fit_mt$tree, bs, type = "n")
#  
#  # assigning transfer bootstrap values to our tree
#  tree_tfbs <- plotBS(fit_mt$tree, bs, type = "n", method = "TBE")

## ----ConsensusNet, fig.cap="ConsensusNet from the standard bootstrap sample.", echo=TRUE----
cnet <- consensusNet(bs, p=0.2)
plot(cnet, show.edge.label=TRUE)

## ----write_tree, eval=FALSE---------------------------------------------------
#  # tree with ultrafast bootstrap values
#  write.tree(fit_mt$tree, "primates.tree")
#  
#  # tree with standard bootstrap values
#  write.tree(tree_stdbs, "primates.tree")
#  
#  # tree with transfer bootstrap values
#  write.tree(tree_tfbs, "primates.tree")

## ----strict_primates, echo=TRUE, cache=TRUE-----------------------------------
fit_strict <- pml_bb(primates, model="HKY+G(4)", method="ultrametric",
                     rearrangement="NNI", control = pml.control(trace = 0))

## ----plot_strict_primates-----------------------------------------------------
plot(fit_strict)

## ----tipdated_data------------------------------------------------------------
fdir <- system.file("extdata/trees", package = "phangorn")
tmp <- read.csv(file.path(fdir,"H3N2_NA_20.csv"))
H3N2 <- read.phyDat(file.path(fdir,"H3N2_NA_20.fasta"), format="fasta")

## ----tipdated_processing------------------------------------------------------
dates <- setNames(tmp$numdate_given, tmp$name)
head(dates)

## ----tipdated_fit-------------------------------------------------------------
fit_td <- pml_bb(H3N2, model="HKY+I", method="tipdated", tip.dates=dates, 
               rearrangement="NNI", control = pml.control(trace = 0))
fit_td

## ----tipdated_plot------------------------------------------------------------
plot(fit_td, align.tip.label=TRUE)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

