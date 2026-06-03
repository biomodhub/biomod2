#### "demo":  Make these (function) definitions easily available to useRs
####  ----   --> we use them in tests in ../tests/exact-ex.R
##                                       ~~~~~~~~~~~~~~~~~~~

### For nilpotent matrices A, exp(A) is polynomial in A
###  Mathworld gives the example of  the general  3 x 3  upper triangle
nilA3 <- function(x,y,z) {
    ## Purpose: simple nilpotent matrix 3x3  A (with A^n = 0 for n >= 3)
    ##          / 0 x z \
    ##     A = [  0 0 y  ]
    ##          \ 0 0 0 /
    ## and its exact matrix exponential
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Jan 2008
    stopifnot((n <- length(x)) == 1, length(y) == 1, length(z) == 1,
              is.numeric(x), is.numeric(y), is.numeric(z))
    list(A = cbind(0, rbind(matrix(c(x,0,z,y), 2,2), 0)),
         expA = cbind(c(1,0,0), c(x,1,0), c(z + x*y/2, y, 1)))
}


## The relative error typically returned by all.equal -- simplified here
relErr <- function(target, current)
    mean(abs(target - current)) / mean(abs(target))

facMat <- function(n, R_FUN, ev = R_FUN(n), M = rMat(n, R_FUN = R_FUN))
{
    ## Purpose: Construct random matrix x of which we "know" expm(x)
    ## because we set  x :=  M %*% diag(ev)  %*% solve(M)
    ## ----------------------------------------------------------------------
    ## Arguments: n:     dimension of matrices
    ##            R_FUN: random number generator function (n)
    ##            ev:    numeric length-n vector of eigenvalues
    ##            M:     n x n matrix. Note that the default,
    ##                   rMat() will give matrices ``not close to singular''
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: Feb 2008 / Aug. 2024 for <complex>
    R_FUN <- match.fun(R_FUN)
    num <- is.numeric(ev)
    stopifnot(n > 0, num || is.complex(ev), length(ev) == n,
              dim(M) == c(n,n), is.numeric(M) || is.complex(M))

    iM <- solve(M)
    ## D <- diag(ev); A = M %*% D %*% iM
    list(A    = M %*% (ev      * iM),
         expA = M %*% (exp(ev) * iM))
}

### --- The 2x2 example with bad condition , see A3 in ./ex2.R
m2ex3 <- function(eps = 0) {
    stopifnot(is.numeric(eps), length(eps) == 1)
    A <- rbind(c(-1,     1),
               c(eps^2, -1))
    I.e <- 1 - eps^2 / 2
    V <- I.e* rbind(    c(-1, 1),
                    eps*c( 1, 1))
    D <- c(-1-eps, -1+eps)
    iV <- ## solve(V)
        rbind(c(-1, 1/eps),
              c( 1, 1/eps)) / (2 * I.e)
    ## NOTE:  kappa(V) = condition_number(V) == 1/eps exactly
    useTol <- 2e-16 / eps
    stopifnot(all.equal(diag(2), V %*% iV,	 tolerance=useTol),
	      all.equal(A, V %*% diag(D) %*% iV, tolerance=useTol) )
    ch.e <- cosh(eps)
    sh.e <- sinh(eps)
    list(A = A,
         expA = exp(-1) *
         rbind(c( ch.e,  sh.e/eps),
               c(sh.e*eps, ch.e  )))
}

###---

rnilMat <- function(n, R_FUN = function(n) rpois(n, lambda=5))
{
    ## random upper triangular (zero-diagonal) nilpotent  n x n matrix
    m <- matrix(0, n,n)
    ut <- upper.tri(m)
    R_FUN <- match.fun(R_FUN)
    m[ut] <- R_FUN(sum(ut))
    m
}

set.seed(17)
m <- rnilMat(10)
if(FALSE)
    Matrix(m)
## 10 x 10 sparse Matrix of class "dtCMatrix"
##
##  . 3 10 7 3 4 9 5 9 6
##  . .  5 4 3 . 5 6 3 6
##  . .  . 5 7 7 3 7 5 6
##  . .  . . 3 7 6 8 2 7
##  . .  . . . 9 5 2 7 6
##  . .  . . . . 8 5 4 6
##  . .  . . . . . 5 5 3
##  . .  . . . . . . 3 5
##  . .  . . . . . . . 3
##  . .  . . . . . . . .

## An interesting example, rounded from above {see ../tests/exact-ex.R} :
dN <- 9*7*320 # 20160
EmN <- matrix(c(dN, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                3*dN, dN, 0, 0, 0, 0, 0, 0, 0, 0,
                352800, 5*dN, dN, 0, 0, 0, 0, 0, 0, 0,
                1018080, 332640, 5*dN, dN, 0, 0, 0, 0, 0, 0,
                2235240, 786240, 292320, 3*dN, dN, 0, 0, 0, 0, 0,
                9368520, 3483480, 1582560, 413280, 181440, dN, 0, 0, 0, 0,
                24676176, 9598680, 5073600, 1562400, 826560, 161280, dN, 0,0,0,
                43730160, 17451000, 10051440, 3430560, 1955520, 504000,
                5*dN, dN, 0, 0,
                68438436, 27747480, 16853760, 6036240, 3638880, 1038240,
                252000, 3*dN, dN, 0,
                119725855, 49165892, 31046760, 11652480, 7198800, 2264640,
                614880, 191520, 3*dN, dN),
              10, 10)
xct10 <- list(m = m, expm = EmN / dN, expmNum = EmN, expmDen = dN)
