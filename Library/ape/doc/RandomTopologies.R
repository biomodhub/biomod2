### R code from vignette source 'RandomTopologies.Rnw'

###################################################
### code chunk number 1: RandomTopologies.Rnw:15-16
###################################################
options(width = 80, prompt = "> ")


###################################################
### code chunk number 2: RandomTopologies.Rnw:49-52
###################################################
library(ape)
n <- 1:10
data.frame(N.topo = sapply(n, howmanytrees))


###################################################
### code chunk number 3: RandomTopologies.Rnw:60-61
###################################################
howmanytrees(25)


###################################################
### code chunk number 4: RandomTopologies.Rnw:68-69
###################################################
howmanytrees(25) * 0.21


###################################################
### code chunk number 5: RandomTopologies.Rnw:80-81
###################################################
data.frame(N.topo = sapply(n, howmanytrees, labeled = FALSE))


###################################################
### code chunk number 6: RandomTopologies.Rnw:83-86
###################################################
n2 <- c(50, 100, 200)
data.frame(N.topo = sapply(n2, howmanytrees, labeled = FALSE),
           row.names = n2)


###################################################
### code chunk number 7: RandomTopologies.Rnw:104-105
###################################################
1 / howmanytrees(50)


###################################################
### code chunk number 8: RandomTopologies.Rnw:110-111
###################################################
1 / sapply(n, howmanytrees)


###################################################
### code chunk number 9: RandomTopologies.Rnw:116-117
###################################################
1 / sapply(n, howmanytrees, labeled = FALSE)


###################################################
### code chunk number 10: RandomTopologies.Rnw:129-132
###################################################
library(phangorn)
TR <- allTrees(4, rooted = TRUE)
TR


###################################################
### code chunk number 11: RandomTopologies.Rnw:135-138
###################################################
layout(matrix(1:15, 3, 5))
par(mar = rep(1.5, 4), xpd = TRUE)
for (i in 1:15) plot(TR[[i]], "c", cex = 1.2, font = 1)


###################################################
### code chunk number 12: RandomTopologies.Rnw:160-165
###################################################
K <- 100
n <- 50
N <- howmanytrees(n, labeled = FALSE)
p <- 1/N
1 - (dbinom(0, K, p) + dbinom(1, K, p))


###################################################
### code chunk number 13: RandomTopologies.Rnw:172-173
###################################################
1 - dbinom(0, K, p) - dbinom(1, K, p)


###################################################
### code chunk number 14: RandomTopologies.Rnw:179-189
###################################################
f <- function(n, K) {
    N <- howmanytrees(n, labeled = FALSE)
    p <- 1/N
    1 - (dbinom(0, K, p) + dbinom(1, K, p))
}
DF <- expand.grid(n = 4:20, K = c(2, 5, 10, 100, 1e3, 1e4))
DF$P <- mapply(f, n = DF$n, K = DF$K)
library(lattice)
xyplot(P ~ n, DF, groups = K, panel = panel.superpose,
     type = "b", auto.key = list(x = .8, y = .9, title = "K"))


###################################################
### code chunk number 15: RandomTopologies.Rnw:195-196
###################################################
subset(DF, n == 20)


###################################################
### code chunk number 16: RandomTopologies.Rnw:198-199
###################################################
f(20, 1e5)


###################################################
### code chunk number 17: RandomTopologies.Rnw:201-202
###################################################
f(30, 1e6)


###################################################
### code chunk number 18: RandomTopologies.Rnw:268-272
###################################################
txt <- c("((,),((,),));", "(((,),(,)),);", "((((,),),),);")
TR5 <- read.tree(text = txt)
layout(matrix(1:3, 1))
for (i in 1:3) plot(TR5[[i]], "c", cex = 1.2, main = LETTERS[i])


###################################################
### code chunk number 19: RandomTopologies.Rnw:304-305
###################################################
howmanytrees(1, labeled = FALSE) * howmanytrees(4, labeled = FALSE)


###################################################
### code chunk number 20: RandomTopologies.Rnw:307-308
###################################################
howmanytrees(2, labeled = FALSE) * howmanytrees(3, labeled = FALSE)


###################################################
### code chunk number 21: RandomTopologies.Rnw:314-315
###################################################
howmanytrees(400, labeled = FALSE) * howmanytrees(403, labeled = FALSE)


###################################################
### code chunk number 22: RandomTopologies.Rnw:371-374
###################################################
TR <- list(rmtopology(1000, 4, TRUE, br = NULL),
           rmtree(1000, 4, TRUE, br = NULL),
           rmtree(1000, 4, TRUE, br = NULL, equiprob = TRUE))


###################################################
### code chunk number 23: RandomTopologies.Rnw:382-394
###################################################
foo.test <- function(x, labeled = TRUE) {
    uTR <- unique(x, use.tip.label = labeled)
    f <- integer(length(uTR))
    for (j in seq_along(uTR)) {
        for (i in seq_along(x)) {
            if (all.equal(x[[i]], uTR[[j]], use.tip.label = labeled))
                f[j] <- f[j] + 1L
        }
    }
    print(f)
    chisq.test(f)
}


###################################################
### code chunk number 24: RandomTopologies.Rnw:399-400
###################################################
lapply(TR, foo.test)


###################################################
### code chunk number 25: RandomTopologies.Rnw:402-403
###################################################
lapply(TR, foo.test, labeled = FALSE)


###################################################
### code chunk number 26: RandomTopologies.Rnw:435-436
###################################################
.Machine$double.xmax


###################################################
### code chunk number 27: RandomTopologies.Rnw:441-442
###################################################
howmanytrees(151)


###################################################
### code chunk number 28: RandomTopologies.Rnw:444-445
###################################################
howmanytrees(792, labeled = FALSE)


###################################################
### code chunk number 29: RandomTopologies.Rnw:488-489
###################################################
(N <- howmanytrees(2e6))


###################################################
### code chunk number 30: RandomTopologies.Rnw:494-495
###################################################
N[2] + 1


###################################################
### code chunk number 31: RandomTopologies.Rnw:500-501
###################################################
ceiling((N[2] + 1) / 3000)


###################################################
### code chunk number 32: RandomTopologies.Rnw:511-515
###################################################
n <- 1:792
N <- howmanytrees(792, labeled = FALSE, detail = TRUE)
lN <- log10(N)
summary(lm(lN ~ n))


###################################################
### code chunk number 33: RandomTopologies.Rnw:520-522
###################################################
s <- n > 700
summary(lm(lN[s] ~ n[s]))


###################################################
### code chunk number 34: RandomTopologies.Rnw:532-535
###################################################
n <- 788:792
log10(sapply(n, howmanytrees, labeled = FALSE))
0.3941 * n - 4.153


###################################################
### code chunk number 35: RandomTopologies.Rnw:594-600
###################################################
xi <- 2.477993
lxi <- log10(xi)
n <- 1000
(n - 1) * lxi # k = 1
for (k in 2:11)
    print((n - k) * lxi + log10(sum(xi^(0:(k - 1)))))


