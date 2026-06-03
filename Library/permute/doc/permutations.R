## ----setup, echo = FALSE, results = "hide", message = FALSE, cache = FALSE----
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
library("permute")

## ----load_jackal--------------------------------------------------------------
library("permute")
data(jackal)
jackal

## ----ttest_jackal-------------------------------------------------------------
jack.t <- t.test(Length ~ Sex, data = jackal, var.equal = TRUE,
                 alternative = "greater")
jack.t

## ----meanFun------------------------------------------------------------------
meanDif <- function(x, grp) {
  mean(x[grp == "Male"]) - mean(x[grp == "Female"])
}

## ----randJackal---------------------------------------------------------------
Djackal <- numeric(length = 5000)
N <- nrow(jackal)
set.seed(42)
for(i in seq_len(length(Djackal) - 1)) {
    perm <- shuffle(N)
    Djackal[i] <- with(jackal, meanDif(Length, Sex[perm]))
}
Djackal[5000] <- with(jackal, meanDif(Length, Sex))

## ----hist_jackal, fig=FALSE, echo=TRUE, eval=FALSE----------------------------
# hist(Djackal, main = "",
#      xlab = expression("Mean difference (Male - Female) in mm"))
# rug(Djackal[5000], col = "red", lwd = 2)

## -----------------------------------------------------------------------------
(Dbig <- sum(Djackal >= Djackal[5000]))

## -----------------------------------------------------------------------------
Dbig / length(Djackal)

## ----draw_hist_jackal, fig=TRUE, echo=FALSE, fig.cap="Distribution of the difference of mean mandible length in random allocations, ten to each sex."----
hist(Djackal, main = "",
     xlab = expression("Mean difference (Male - Female) in mm"))
rug(Djackal[5000], col = "red", lwd = 2)

## -----------------------------------------------------------------------------
choose(20, 10)

## ----show_args----------------------------------------------------------------
args(shuffle)

## ----show_str-----------------------------------------------------------------
str(how())

## ----compare_shuffle_sample---------------------------------------------------
set.seed(2)
(r1 <- shuffle(10))
set.seed(2)
(r2 <- sample(1:10, 10, replace = FALSE))
all.equal(r1, r2)

## ----series1------------------------------------------------------------------
set.seed(4)
x <- 1:10
CTRL <- how(within = Within(type = "series"))
perm <- shuffle(10, control = CTRL)
perm
x[perm] ## equivalent

## ----grid1--------------------------------------------------------------------
set.seed(4)
plt <- gl(3, 9)
CTRL <- how(within = Within(type = "grid", ncol = 3, nrow = 3),
            plots = Plots(strata = plt))
perm <- shuffle(length(plt), control = CTRL)
perm

## ----vis_grid1, keep.source=TRUE----------------------------------------------
## Original
lapply(split(seq_along(plt), plt), matrix, ncol = 3)
## Shuffled
lapply(split(perm, plt), matrix, ncol = 3)

## ----grid_2, keep.source=TRUE-------------------------------------------------
set.seed(4)
CTRL <- how(within = Within(type = "grid", ncol = 3, nrow = 3,
                            constant = TRUE),
            plots = Plots(strata = plt))
perm2 <- shuffle(length(plt), control = CTRL)
lapply(split(perm2, plt), matrix, ncol = 3)

## ----series_2, results="hide"-------------------------------------------------
how(nperm = 10, within = Within(type = "series"))

## ----shuffleSet_1-------------------------------------------------------------
set.seed(4)
CTRL <- how(within = Within(type = "series"))
pset <- shuffleSet(10, nset = 5, control = CTRL)
pset

## ----results="hide"-----------------------------------------------------------
how(nperm = 999)

## ----withinArgs, echo=FALSE---------------------------------------------------
args(Within)

## ----ptest-fun----------------------------------------------------------------
pt.test <- function(x, group, nperm = 199) {
    ## mean difference function
    meanDif <- function(i, x, grp) {
        grp <- grp[i]
        mean(x[grp == "Male"]) - mean(x[grp == "Female"])
    }
    ## check x and group are of same length
    stopifnot(all.equal(length(x), length(group)))
    ## number of observations
    N <- nobs(x)
    ## generate the required set of permutations
    pset <- shuffleSet(N, nset = nperm)
    ## iterate over the set of permutations applying meanDif
    D <- apply(pset, 1, meanDif, x = x, grp = group)
    ## add on the observed mean difference
    D <- c(meanDif(seq_len(N), x, group), D)
    ## compute & return the p-value
    Ds <- sum(D >= D[1]) # how many >= to the observed diff?
    Ds / (nperm + 1)     # what proportion of perms is this (the pval)?
}

## ----run-ptest----------------------------------------------------------------
set.seed(42) ## same seed as earlier
pval <- with(jackal, pt.test(Length, Sex, nperm = 4999))
pval

## ----parallel-ptest-fun-------------------------------------------------------
ppt.test <- function(x, group, nperm = 199, cores = 2) {
    ## mean difference function
    meanDif <- function(i, .x, .grp) {
        .grp <- .grp[i]
        mean(.x[.grp == "Male"]) - mean(.x[.grp == "Female"])
    }
    ## check x and group are of same length
    stopifnot(all.equal(length(x), length(group)))
    ## number of observations
    N <- nobs(x)
    ## generate the required set of permutations
    pset <- shuffleSet(N, nset = nperm)
    if (cores > 1) {
        ## initiate a cluster
        cl <- makeCluster(cores)
        on.exit(stopCluster(cl = cl))
        ## iterate over the set of permutations applying meanDif
        D <- parRapply(cl, pset, meanDif, .x = x, .grp = group)
    } else {
        D <- apply(pset, 1, meanDif, .x = x, .grp = group)
    }
    ## add on the observed mean difference
    D <- c(meanDif(seq_len(N), x, group), D)
    ## compute & return the p-value
    Ds <- sum(D >= D[1]) # how many >= to the observed diff?
    Ds / (nperm + 1)     # what proportion of perms is this (the pval)?
}

## ----run-pptest---------------------------------------------------------------
require("parallel")
set.seed(42)
system.time(ppval <- ppt.test(jackal$Length, jackal$Sex, nperm = 9999,
                              cores = 2))
ppval

## ----run-pptest2--------------------------------------------------------------
set.seed(42)
system.time(ppval2 <- ppt.test(jackal$Length, jackal$Sex, nperm = 9999,
                               cores = 1))
ppval2

## ----get-set-eg0--------------------------------------------------------------
hh <- how()

## ----get-set-eg1--------------------------------------------------------------
getNperm(hh)

## ----<get-set-eg2-------------------------------------------------------------
getCall(hh)
setNperm(hh) <- 999
getNperm(hh)
getCall(hh)

## ----get-set-eg3--------------------------------------------------------------
hh <- how(within = Within(type = "series"),
          plots = Plots(type = "series", strata = gl(10, 5)),
          blocks = gl(5, 10))

## ----get-set-eg4--------------------------------------------------------------
pl <- getPlots(hh)
setType(pl) <- "free"
setPlots(hh) <- pl

## ----get-set-eg5--------------------------------------------------------------
getType(hh, which = "plots")

## ----get-set-eg6--------------------------------------------------------------
getCall(getPlots(hh))

## ----get-set-eg7--------------------------------------------------------------
hh <- update(hh, plots = update(getPlots(hh), type = "series"))
getType(hh, which = "plots")

## ----seesionInfo, echo=FALSE--------------------------------------------------
sessioninfo::session_info()

