### R code from vignette source 'expm.Rnw'

###################################################
### code chunk number 1: expm.Rnw:49-55
###################################################
library(expm)
m <- matrix(c(4, 1, 1, 2, 4, 1, 0, 1, 4), 3, 3)
expm(m)
dimnames(m) <- list(letters[1:3], LETTERS[1:3])
m
expm(m)


