### R code from vignette source 'MVT_Rnews.Rnw'

###################################################
### code chunk number 1: prelim
###################################################
set.seed(290875)
### options(width=60, prompt="R> ")


###################################################
### code chunk number 2: smallex
###################################################
library("mvtnorm")
m <- 3
sigma <- diag(3)
sigma[2,1] <- 3/5
sigma[3,1] <- 1/3
sigma[3,2] <- 11/15
(prb <- pmvnorm(lower = rep(-Inf, m), upper = c(1, 4, 2),
                mean = rep(0, m), 
                corr = sigma) ### only lower triangular 
                              ### part used
)


###################################################
### code chunk number 3: cats
###################################################
n <- c(26, 24, 20, 33, 32)
V <- diag(1/n)
df <- sum(n) - length(n)
C <- matrix(c( 1, 1, 1, 0, 0,
              -1, 0, 0, 1, 0,
               0,-1, 0, 0, 1,
               0, 0, 0,-1,-1,
               0, 0,-1, 0 ,0), 
            ncol = length(n))
### covariance matrix 
cv <- tcrossprod(C %*% V, C)
### correlation matrix
cr <- cov2cor(cv)
delta <- rep(0, 5)
(qnt <- qmvt(0.95, df = df, delta = delta, corr = cr, 
             abseps = 0.0001, maxpts = 100000, tail = "both"))
