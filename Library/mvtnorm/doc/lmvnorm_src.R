### R code from vignette source 'lmvnorm_src.Rnw'

###################################################
### code chunk number 1: mvtnorm-citation
###################################################
year <- substr(packageDescription("mvtnorm")$Date, 1, 4)
version <- packageDescription("mvtnorm")$Version


###################################################
### code chunk number 2: digits
###################################################
options(digits = 4)


###################################################
### code chunk number 3: chk
###################################################
chk <- function(...) stopifnot(isTRUE(all.equal(...)))


###################################################
### code chunk number 4: example
###################################################
library("mvtnorm")
set.seed(290875)
N <- 4L
J <- 5L
rn <- paste0("C_", 1:N)
nm <- LETTERS[1:J]
Jn <- J * (J - 1) / 2
## data
xn <- matrix(runif(N * Jn), ncol = N)
colnames(xn) <- rn
xd <- matrix(runif(N * (Jn + J)), ncol = N)
colnames(xd) <- rn

(lxn <- ltMatrices(xn, byrow = TRUE, names = nm))
dim(lxn)
dimnames(lxn)
lxd <- ltMatrices(xd, byrow = TRUE, diag = TRUE, names = nm)
dim(lxd)
dimnames(lxd)

lxn <- as.syMatrices(lxn)
lxn


###################################################
### code chunk number 5: ex-reorder
###################################################
## constructor + .reorder + as.array
a <- as.array(ltMatrices(xn, byrow = TRUE))
b <- as.array(ltMatrices(ltMatrices(xn, byrow = TRUE), 
                         byrow = FALSE))
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = FALSE))
b <- as.array(ltMatrices(ltMatrices(xn, byrow = FALSE), 
                         byrow = TRUE))
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE))
b <- as.array(ltMatrices(ltMatrices(xd, byrow = TRUE, diag = TRUE), 
                         byrow = FALSE))
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE))
b <- as.array(ltMatrices(ltMatrices(xd, byrow = FALSE, diag = TRUE), 
                         byrow = TRUE))
chk(a, b)


###################################################
### code chunk number 6: ex-subset (eval = FALSE)
###################################################
## ## subset
## a <- as.array(ltMatrices(xn, byrow = FALSE, names = nm)[i, j])
## b <- as.array(ltMatrices(xn, byrow = FALSE, names = nm))[j, j, i]
## chk(a, b)
## 
## a <- as.array(ltMatrices(xn, byrow = TRUE, names = nm)[i, j])
## b <- as.array(ltMatrices(xn, byrow = TRUE, names = nm))[j, j, i]
## chk(a, b)
## 
## a <- as.array(ltMatrices(xd, byrow = FALSE, 
##                          diag = TRUE, names = nm)[i, j])
## b <- as.array(ltMatrices(xd, byrow = FALSE, 
##                          diag = TRUE, names = nm))[j, j, i]
## chk(a, b)
## 
## a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
##                          names = nm)[i, j])
## b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
##                          names = nm))[j, j, i]
## chk(a, b)


###################################################
### code chunk number 7: ex-subset-1
###################################################
i <- colnames(xn)[1:2]
j <- 2:4
## subset
a <- as.array(ltMatrices(xn, byrow = FALSE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = FALSE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm))[j, j, i]
chk(a, b)


###################################################
### code chunk number 8: ex-subset-2
###################################################
i <- 1:2
j <- nm[2:4]
## subset
a <- as.array(ltMatrices(xn, byrow = FALSE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = FALSE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm))[j, j, i]
chk(a, b)


###################################################
### code chunk number 9: ex-subset-3
###################################################
j <- c(1, 3, 5)
## subset
a <- as.array(ltMatrices(xn, byrow = FALSE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = FALSE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm))[j, j, i]
chk(a, b)


###################################################
### code chunk number 10: ex-subset-4
###################################################
j <- nm[c(1, 3, 5)]
## subset
a <- as.array(ltMatrices(xn, byrow = FALSE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = FALSE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm))[j, j, i]
chk(a, b)


###################################################
### code chunk number 11: ex-subset-5
###################################################
j <- -c(1, 3, 5)
## subset
a <- as.array(ltMatrices(xn, byrow = FALSE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = FALSE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xn, byrow = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = FALSE, 
                         diag = TRUE, names = nm))[j, j, i]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm)[i, j])
b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE, 
                         names = nm))[j, j, i]
chk(a, b)


###################################################
### code chunk number 12: ex-subset-6
###################################################
## subset
j <- nm[sample(1:J)]
ltM <- ltMatrices(xn, byrow = FALSE, names = nm)
try(ltM[i, j])
ltM <- as.syMatrices(ltM)
a <- as.array(ltM[i, j])
b <- as.array(ltM)[j, j, i]
chk(a, b)


###################################################
### code chunk number 13: ex-Lower_tri
###################################################
## J <- 4
M <- ltMatrices(matrix(1:10, nrow = 10, ncol = 2), diag = TRUE)
Lower_tri(M, diag = FALSE)
Lower_tri(M, diag = TRUE)
M <- ltMatrices(matrix(1:6, nrow = 6, ncol = 2), diag = FALSE)
Lower_tri(M, diag = FALSE)
Lower_tri(M, diag = TRUE)
## multiple symmetric matrices
Lower_tri(invchol2cor(M))


###################################################
### code chunk number 14: ex-diag
###################################################
all(diagonals(ltMatrices(xn, byrow = TRUE)) == 1L)


###################################################
### code chunk number 15: ex-addiag
###################################################
lxd2 <- lxn
diagonals(lxd2) <- 1
chk(as.array(lxd2), as.array(lxn))


###################################################
### code chunk number 16: ex-diagJ
###################################################
(I5 <- diagonals(5L))
diagonals(I5) <- 1:5
I5


###################################################
### code chunk number 17: ex-mult
###################################################
lxn <- ltMatrices(xn, byrow = TRUE)
lxd <- ltMatrices(xd, byrow = TRUE, diag = TRUE)
y <- matrix(runif(N * J), nrow = J)
a <- Mult(lxn, y)
A <- as.array(lxn)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(A[,,i] %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

a <- Mult(lxd, y)
A <- as.array(lxd)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(A[,,i] %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

### recycle C
chk(Mult(lxn[rep(1, N),], y), Mult(lxn[1,], y), check.attributes = FALSE)

### recycle y
chk(Mult(lxn, y[,1]), Mult(lxn, y[,rep(1, N)]))

### tcrossprod as multiplication
i <- sample(1:N)[1]
M <- t(as.array(lxn)[,,i])
a <- sapply(1:J, function(j) Mult(lxn[i,], M[,j,drop = FALSE]))
rownames(a) <- colnames(a) <- dimnames(lxn)[[2L]]
b <- as.array(Tcrossprod(lxn[i,]))[,,1]
chk(a, b, check.attributes = FALSE)


###################################################
### code chunk number 18: ex-tmult
###################################################
a <- Mult(lxn, y, transpose = TRUE)
A <- as.array(lxn)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(t(A[,,i]) %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

a <- Mult(lxd, y, transpose = TRUE)
A <- as.array(lxd)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(t(A[,,i]) %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

### recycle C
chk(Mult(lxn[rep(1, N),], y, transpose = TRUE), 
    Mult(lxn[1,], y, transpose = TRUE), check.attributes = FALSE)

### recycle y
chk(Mult(lxn, y[,1], transpose = TRUE), 
    Mult(lxn, y[,rep(1, N)], transpose = TRUE))


###################################################
### code chunk number 19: ex-symult
###################################################
J <- 5
N1 <- 10
ex <- expression({
  C <- syMatrices(matrix(runif(N2 * J * (J + c(-1, 1)[DIAG + 1L] ) / 2),
                          ncol = N2), 
                  diag = DIAG)
  x <- matrix(runif(N1 * J), nrow = J)
  Ca <- as.array(C)
  p1 <- do.call("cbind", lapply(1:N1, function(i) 
      Ca[,,c(1,i)[(N2 > 1) + 1]] %*% x[,i]))
  p2 <- Mult(C, x)
  chk(p1, p2)
})

N2 <- N1
DIAG <- TRUE
eval(ex)
N2 <- 1
DIAG <- TRUE
eval(ex)
N2 <- 1
DIAG <- FALSE
eval(ex)
N2 <- N1
DIAG <- FALSE
eval(ex)


###################################################
### code chunk number 20: ex-solve
###################################################
## solve
A <- as.array(lxn)
a <- solve(lxn)
a <- as.array(a)
b <- array(apply(A, 3L, function(x) solve(x), simplify = TRUE), 
           dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

A <- as.array(lxd)
a <- as.array(solve(lxd))
b <- array(apply(A, 3L, function(x) solve(x), simplify = TRUE), 
           dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

chk(solve(lxn, y), Mult(solve(lxn), y))
chk(solve(lxd, y), Mult(solve(lxd), y))

### recycle C
chk(solve(lxn[1,], y), as.array(solve(lxn[1,]))[,,1] %*% y)
chk(solve(lxn[rep(1, N),], y), solve(lxn[1,], y), check.attributes = FALSE)

### recycle y
chk(solve(lxn, y[,1]), solve(lxn, y[,rep(1, N)]))


###################################################
### code chunk number 21: ex-tsolve
###################################################
chk(solve(lxn[1,], y, transpose = TRUE), 
    t(as.array(solve(lxn[1,]))[,,1]) %*% y)


###################################################
### code chunk number 22: ex-logdet
###################################################
chk(logdet(lxn), colSums(log(diagonals(lxn))))
chk(logdet(lxd[1,]), colSums(log(diagonals(lxd[1,]))))
chk(logdet(lxd), colSums(log(diagonals(lxd))))
lxd2 <- ltMatrices(lxd, byrow = !attr(lxd, "byrow"))
chk(logdet(lxd2), colSums(log(diagonals(lxd2))))


###################################################
### code chunk number 23: ex-tcrossprod
###################################################
## Tcrossprod
a <- as.array(Tcrossprod(lxn))
b <- array(apply(as.array(lxn), 3L, function(x) tcrossprod(x), simplify = TRUE), 
           dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Tcrossprod(lxn, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Tcrossprod(lxn)))

a <- as.array(Tcrossprod(lxd))
b <- array(apply(as.array(lxd), 3L, function(x) tcrossprod(x), simplify = TRUE), 
           dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Tcrossprod(lxd, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Tcrossprod(lxd)))


###################################################
### code chunk number 24: ex-crossprod
###################################################
## Crossprod
a <- as.array(Crossprod(lxn))
b <- array(apply(as.array(lxn), 3L, function(x) crossprod(x), simplify = TRUE), 
           dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Crossprod(lxn, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Crossprod(lxn)))

a <- as.array(Crossprod(lxd))
b <- array(apply(as.array(lxd), 3L, function(x) crossprod(x), simplify = TRUE), 
           dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Crossprod(lxd, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Crossprod(lxd)))


###################################################
### code chunk number 25: ex-tcrossprod-2
###################################################
### tcrossprod
a <- as.array(tcrossprod(lxn))
b <- array(apply(as.array(lxn), 3L, function(x) tcrossprod(x), simplify = TRUE), 
           dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

a <- as.array(tcrossprod(lxd))
b <- array(apply(as.array(lxd), 3L, function(x) tcrossprod(x), simplify = TRUE), 
           dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

## crossprod
a <- as.array(crossprod(lxn))
b <- array(apply(as.array(lxn), 3L, function(x) crossprod(x), simplify = TRUE), 
           dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

a <- as.array(crossprod(lxd))
b <- array(apply(as.array(lxd), 3L, function(x) crossprod(x), simplify = TRUE), 
           dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)


###################################################
### code chunk number 26: ex-Mult-2
###################################################
a <- lxn %*% y
A <- as.array(lxn)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(A[,,i] %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

a <- lxd %*% y
A <- as.array(lxd)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(A[,,i] %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

### recycle C
chk(lxn[rep(1, N),] %*% y, lxn[1,] %*% y, check.attributes = FALSE)

### recycle y
chk(lxn %*% y[,1], lxn %*% y[,rep(1, N)])

### tcrossprod as multiplication
i <- sample(1:N)[1]
M <- t(as.array(lxn)[,,i])
a <- sapply(1:J, function(j) lxn[i,] %*% M[,j,drop = FALSE])
rownames(a) <- colnames(a) <- dimnames(lxn)[[2L]]
b <- as.array(tcrossprod(lxn[i,]))[,,1]
chk(a, b, check.attributes = FALSE)

a <- crossprod(lxn, y)
A <- as.array(lxn)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(t(A[,,i]) %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

a <- crossprod(lxd, y)
A <- as.array(lxd)
b <- do.call("rbind", lapply(1:ncol(y), 
    function(i) t(t(A[,,i]) %*% y[,i,drop = FALSE])))
chk(a, t(b), check.attributes = FALSE)

### recycle C
chk(crossprod(lxn[rep(1, N),], y), 
    crossprod(lxn[1,], y), check.attributes = FALSE)

### recycle y
chk(crossprod(lxn, y[,1]), 
    crossprod(lxn, y[,rep(1, N)]))


###################################################
### code chunk number 27: chol
###################################################
Sigma <- tcrossprod(lxd)
chk(chol(Sigma), lxd)
Sigma <- tcrossprod(lxn)
## Sigma and chol(Sigma) always have diagonal, lxn doesn't
chk(as.array(chol(Sigma)), as.array(lxn))


###################################################
### code chunk number 28: kronecker
###################################################
J <- 10

d <- TRUE
L <- diag(J)
L[lower.tri(L, diag = d)] <- prm <- runif(J * (J + c(-1, 1)[d + 1]) / 2)

C <- solve(L)

D <- -kronecker(t(C), C)

S <- diag(J)
S[lower.tri(S, diag = TRUE)] <- x <- runif(J * (J + 1) / 2)

SD0 <- matrix(c(S) %*% D, ncol = J)

SD1 <- -crossprod(C, tcrossprod(S, C))

a <- ltMatrices(C[lower.tri(C, diag = TRUE)], diag = TRUE, byrow = FALSE)
b <- ltMatrices(x, diag = TRUE, byrow = FALSE)

SD2 <- -vectrick(a, b, a)
SD2a <- -vectrick(a, b)
chk(SD2, SD2a)

chk(SD0[lower.tri(SD0, diag = d)], 
    SD1[lower.tri(SD1, diag = d)])
chk(SD0[lower.tri(SD0, diag = d)],
    c(unclass(SD2)))

### same; but SD2 is vec(SD0)
S <- t(matrix(as.array(b), byrow = FALSE, nrow = 1))
SD2 <- -vectrick(a, S, a)
SD2a <- -vectrick(a, S)
chk(SD2, SD2a)

chk(c(SD0), c(SD2))

### N > 1
N <- 4L
prm <- runif(J * (J - 1) / 2)
C <- ltMatrices(prm)
S <- matrix(runif(J^2 * N), ncol = N)
A <- vectrick(C, S, C)
Cx <- as.array(C)[,,1]
B <- apply(S, 2, function(x) t(Cx) %*% matrix(x, ncol = J) %*% t(Cx))
chk(A, B)

A <- vectrick(C, S, C, transpose = c(FALSE, FALSE))
Cx <- as.array(C)[,,1]
B <- apply(S, 2, function(x) Cx %*% matrix(x, ncol = J) %*% Cx)
chk(A, B)


###################################################
### code chunk number 29: conv-ex-1
###################################################
prec2pc <- function(x) {
    ret <- -cov2cor(x)
    diag(ret) <- 0
    ret
}
L <- lxn
Sigma <- apply(as.array(L), 3, 
               function(x) tcrossprod(solve(x)), simplify = FALSE)
Prec <- lapply(Sigma, solve)
Corr <- lapply(Sigma, cov2cor)
CP <- lapply(Corr, solve)
PC <- lapply(Prec, function(x) prec2pc(x))
chk(unlist(Sigma), c(as.array(invchol2cov(L))), 
    check.attributes = FALSE)
chk(unlist(Prec), c(as.array(invchol2pre(L))), 
    check.attributes = FALSE)
chk(unlist(Corr), c(as.array(invchol2cor(L))), 
    check.attributes = FALSE)
chk(unlist(CP), c(as.array(crossprod(invcholD(L)))), 
    check.attributes = FALSE)
chk(unlist(PC), c(as.array(invchol2pc(L))), 
    check.attributes = FALSE)


###################################################
### code chunk number 30: conv-ex-2
###################################################
C <- lxn
Sigma <- apply(as.array(C), 3, 
               function(x) tcrossprod(x), simplify = FALSE)
Prec <- lapply(Sigma, solve)
Corr <- lapply(Sigma, cov2cor)
CP <- lapply(Corr, solve)
PC <- lapply(Prec, function(x) prec2pc(x))
chk(unlist(Sigma), c(as.array(chol2cov(C))), 
    check.attributes = FALSE)
chk(unlist(Prec), c(as.array(chol2pre(C))), 
    check.attributes = FALSE)
chk(unlist(Corr), c(as.array(chol2cor(C))), 
    check.attributes = FALSE)
chk(unlist(CP), c(as.array(crossprod(solve(Dchol(C))))), 
    check.attributes = FALSE)
chk(unlist(PC), c(as.array(chol2pc(C))), 
    check.attributes = FALSE)


###################################################
### code chunk number 31: conv-ex-3
###################################################
L <- lxd
Sigma <- apply(as.array(L), 3, 
               function(x) tcrossprod(solve(x)), simplify = FALSE)
Prec <- lapply(Sigma, solve)
Corr <- lapply(Sigma, cov2cor)
CP <- lapply(Corr, solve)
PC <- lapply(Prec, function(x) prec2pc(x))
chk(unlist(Sigma), c(as.array(invchol2cov(L))), 
    check.attributes = FALSE)
chk(unlist(Prec), c(as.array(invchol2pre(L))), 
    check.attributes = FALSE)
chk(unlist(Corr), c(as.array(invchol2cor(L))), 
    check.attributes = FALSE)
chk(unlist(CP), c(as.array(crossprod(invcholD(L)))), 
    check.attributes = FALSE)
chk(unlist(PC), c(as.array(invchol2pc(L))), 
    check.attributes = FALSE)


###################################################
### code chunk number 32: conv-ex-4
###################################################
C <- lxd
Sigma <- apply(as.array(C), 3, 
               function(x) tcrossprod(x), simplify = FALSE)
Prec <- lapply(Sigma, solve)
Corr <- lapply(Sigma, cov2cor)
CP <- lapply(Corr, solve)
PC <- lapply(Prec, function(x) prec2pc(x))
chk(unlist(Sigma), c(as.array(chol2cov(C))), 
    check.attributes = FALSE)
chk(unlist(Prec), c(as.array(chol2pre(C))), 
    check.attributes = FALSE)
chk(unlist(Corr), c(as.array(chol2cor(C))), 
    check.attributes = FALSE)
chk(unlist(CP), c(as.array(crossprod(solve(Dchol(C))))), 
    check.attributes = FALSE)
chk(unlist(PC), c(as.array(chol2pc(C))), 
    check.attributes = FALSE)


###################################################
### code chunk number 33: aperm-tests
###################################################
L <- as.invchol(lxn)
J <- dim(L)[2L]
Lp <- aperm(a = L, perm = p <- sample(1:J))
chk(invchol2cov(L)[,p], invchol2cov(Lp))

C <- as.chol(lxn)
J <- dim(C)[2L]
Cp <- aperm(a = C, perm = p <- sample(1:J))
chk(chol2cov(C)[,p], chol2cov(Cp))


###################################################
### code chunk number 34: marg
###################################################
Sigma <- tcrossprod(lxd)
j <- 1:3
chk(Sigma[,j], tcrossprod(marg_mvnorm(chol = lxd, which = j)$chol))
j <- 2:4
chk(Sigma[,j], tcrossprod(marg_mvnorm(chol = lxd, which = j)$chol))

Sigma <- tcrossprod(solve(lxd))
j <- 1:3
chk(Sigma[,j], tcrossprod(solve(marg_mvnorm(invchol = lxd, which = j)$invchol)))
j <- 2:4
chk(Sigma[,j], tcrossprod(solve(marg_mvnorm(invchol = lxd, which = j)$invchol)))


###################################################
### code chunk number 35: cond-general
###################################################
Sigma <- as.array(tcrossprod(lxd))[,,1]
j <- 2:4
y <- matrix(c(-1, 2, 1), nrow = 3)

cm <- Sigma[-j, j,drop = FALSE] %*% solve(Sigma[j,j]) %*%  y
colnames(cm) <- dimnames(lxd)[[1L]][1L]
cS <- Sigma[-j, -j] - Sigma[-j,j,drop = FALSE] %*% 
      solve(Sigma[j,j]) %*% Sigma[j,-j,drop = FALSE]

cmv <- cond_mvnorm(chol = lxd[1,], which_given = j, given = y)

chk(cm, cmv$mean)
chk(cS, as.array(tcrossprod(cmv$chol))[,,1])

Sigma <- as.array(tcrossprod(solve(lxd)))[,,1]
j <- 2:4
y <- matrix(c(-1, 2, 1), nrow = 3)

cm <- Sigma[-j, j,drop = FALSE] %*% solve(Sigma[j,j]) %*%  y
colnames(cm) <- dimnames(lxd)[[1L]][1L]
cS <- Sigma[-j, -j] - Sigma[-j,j,drop = FALSE] %*% 
      solve(Sigma[j,j]) %*% Sigma[j,-j,drop = FALSE]

cmv <- cond_mvnorm(invchol = lxd[1,], which_given = j, given = y)

chk(cm, cmv$mean)
chk(cS, as.array(tcrossprod(solve(cmv$invchol)))[,,1])


###################################################
### code chunk number 36: cond-simple
###################################################
Sigma <- as.array(tcrossprod(lxd))[,,1]
j <- 1:3
y <- matrix(c(-1, 2, 1), nrow = 3)

cm <- Sigma[-j, j,drop = FALSE] %*% solve(Sigma[j,j]) %*%  y
colnames(cm) <- dimnames(lxd)[[1L]][1L]
cS <- Sigma[-j, -j] - Sigma[-j,j,drop = FALSE] %*% 
      solve(Sigma[j,j]) %*% Sigma[j,-j,drop = FALSE]

cmv <- cond_mvnorm(chol = lxd[1,], which_given = j, given = y)

chk(c(cm), c(cmv$mean))
chk(cS, as.array(tcrossprod(cmv$chol))[,,1])

Sigma <- as.array(tcrossprod(solve(lxd)))[,,1]
j <- 1:3
y <- matrix(c(-1, 2, 1), nrow = 3)

cm <- Sigma[-j, j,drop = FALSE] %*% solve(Sigma[j,j]) %*%  y
colnames(cm) <- dimnames(lxd)[[1L]][1L]
cS <- Sigma[-j, -j] - Sigma[-j,j,drop = FALSE] %*% 
      solve(Sigma[j,j]) %*% Sigma[j,-j,drop = FALSE]

cmv <- cond_mvnorm(invchol = lxd[1,], which_given = j, given = y)

chk(c(cm), c(cmv$mean))
chk(cS, as.array(tcrossprod(solve(cmv$invchol)))[,,1])


###################################################
### code chunk number 37: ex-MV
###################################################
N <- 1000L
J <- 50L
lt <- ltMatrices(matrix(runif(N * J * (J + 1) / 2) + 1, ncol = N), 
                 diag = TRUE, byrow = FALSE)
Z <- matrix(rnorm(N * J), ncol = N)
Y <- solve(lt, Z)
ll1 <- sum(dnorm(lt %*% Y, log = TRUE)) + sum(log(diagonals(lt)))

S <- as.array(tcrossprod(solve(lt)))
ll2 <- sum(sapply(1:N, function(i) 
                           dmvnorm(x = Y[,i], sigma = S[,,i], log = TRUE)))
chk(ll1, ll2)


###################################################
### code chunk number 38: ex-MV-d
###################################################
ll3 <- ldmvnorm(obs = Y, invchol = lt)
chk(ll1, ll3)


###################################################
### code chunk number 39: ex-MV-mc
###################################################
## marginal of and conditional on these
(j <- 1:5 * 10)
md <- marg_mvnorm(invchol = lt, which = j)
cd <- cond_mvnorm(invchol = lt, which_given = j, given = Y[j,])

ll3 <- sum(dnorm(md$invchol %*% Y[j,], log = TRUE)) + 
       sum(log(diagonals(md$invchol))) +
       sum(dnorm(cd$invchol %*% (Y[-j,] - cd$mean), log = TRUE)) + 
       sum(log(diagonals(cd$invchol)))
chk(ll1, ll3)


###################################################
### code chunk number 40: chapterseed
###################################################
set.seed(270312)


###################################################
### code chunk number 41: fct-lpmvnormR
###################################################

lpmvnormR <- function(lower, upper, mean = 0, invcholmean, center = NULL, chol, logLik = TRUE, ...) {

    ### only needed for input checks
    stopifnot(missing(invcholmean))

    ### get access to internal mvtnorm functions
    .check_obs_mean <- mvtnorm:::.check_obs_mean

    
    if (!is.matrix(lower)) lower <- matrix(lower, ncol = 1)
    if (!is.matrix(upper)) upper <- matrix(upper, ncol = 1)
    stopifnot(isTRUE(all.equal(dim(lower), dim(upper))))

    stopifnot(is.ltMatrices(chol))		### NOTE: replace with is.chol
    byrow_orig <- attr(chol, "byrow")
    chol <- ltMatrices(chol, byrow = TRUE)
    d <- dim(chol)
    ### allow single matrix C
    N <- ifelse(d[1L] == 1, ncol(lower), d[1L])
    J <- d[2L]

    stopifnot(nrow(lower) == J && ncol(lower) == N)
    stopifnot(nrow(upper) == J && ncol(upper) == N)

    if (!missing(mean)) {
        lower <- .check_obs_mean(lower, mean, J = J, N = N)
        upper <- .check_obs_mean(upper, mean, J = J, N = N)
    }

    if (!missing(invcholmean)) {
        stopifnot(.check_obs_invcholmean(lower, invcholmean, J = J, N = N))
        center <- - invcholmean
    }

    if (!is.null(center)) {
        if (!is.matrix(center)) center <- matrix(center, ncol = 1)
        stopifnot(nrow(center) == J && ncol(center == N))
    }
    

    sigma <- Tcrossprod(chol)
    S <- as.array(sigma)
    idx <- 1

    ret <- error <- numeric(N)
    for (i in 1:N) {
        if (dim(sigma)[[1L]] > 1) idx <- i
        tmp <- pmvnorm(lower = lower[,i], upper = upper[,i], sigma = S[,,idx], ...)
        ret[i] <- tmp
        error[i] <- attr(tmp, "error")
    }
    attr(ret, "error") <- error

    if (logLik)
        return(sum(log(pmax(ret, .Machine$double.eps))))

    ret
}



###################################################
### code chunk number 42: ex-lpmvnorm_R
###################################################
J <- 5L
N <- 10L

x <- matrix(runif(N * J * (J + 1) / 2), ncol = N)
lx <- ltMatrices(x, byrow = TRUE, diag = TRUE)

a <- matrix(runif(N * J), nrow = J) - 2
a[sample(J * N)[1:2]] <- -Inf
b <- a + 2 + matrix(runif(N * J), nrow = J)
b[sample(J * N)[1:2]] <- Inf

(phat <- c(lpmvnormR(a, b, chol = lx, logLik = FALSE)))


###################################################
### code chunk number 43: ex-again
###################################################
phat
exp(lpmvnorm(a, b, chol = lx, M = 25000, logLik = FALSE, fast = TRUE))
exp(lpmvnorm(a, b, chol = lx, M = 25000, logLik = FALSE, fast = FALSE))


###################################################
### code chunk number 44: ex-lpmvnorm
###################################################
M <- 10000L
if (require("qrng", quietly = TRUE)) {
    ### quasi-Monte-Carlo
    W <- t(ghalton(M, d = J - 1))
} else {
    ### Monte-Carlo
    ### byrow = TRUE because adding / removing dimensions
    ### keeps the MC points for the remaining dimensions constant
    W <- matrix(runif(M * (J - 1)), nrow = J - 1, byrow = TRUE)
}

### Genz & Bretz, 2002, without early stopping (really?)
pGB <- lpmvnormR(a, b, chol = lx, logLik = FALSE, 
                algorithm = GenzBretz(maxpts = M, abseps = 0, releps = 0))
### Genz 1992 with quasi-Monte-Carlo, fast pnorm
pGqf <- exp(lpmvnorm(a, b, chol = lx, w = W, M = M, logLik = FALSE, 
                     fast = TRUE))
### Genz 1992, original Monte-Carlo, fast pnorm
pGf <- exp(lpmvnorm(a, b, chol = lx, w = NULL, M = M, logLik = FALSE, 
                    fast = TRUE))
### Genz 1992 with quasi-Monte-Carlo, R::pnorm
pGqs <- exp(lpmvnorm(a, b, chol = lx, w = W, M = M, logLik = FALSE, 
                     fast = FALSE))
### Genz 1992, original Monte-Carlo, R::pnorm
pGs <- exp(lpmvnorm(a, b, chol = lx, w = NULL, M = M, logLik = FALSE, 
                    fast = FALSE))

cbind(pGB, pGqf, pGf, pGqs, pGs)


###################################################
### code chunk number 45: ex-uni
###################################################
### test univariate problem
### call pmvnorm
pGB <- lpmvnormR(a[1,,drop = FALSE], b[1,,drop = FALSE], chol = lx[,1], 
                logLik = FALSE, 
                algorithm = GenzBretz(maxpts = M, abseps = 0, releps = 0))
### call lpmvnorm
pGq <- exp(lpmvnorm(a[1,,drop = FALSE], b[1,,drop = FALSE], chol = lx[,1], 
                   logLik = FALSE))
### ground truth
ptr <- pnorm(b[1,] / c(unclass(lx[,1]))) - pnorm(a[1,] / c(unclass(lx[,1])))

cbind(c(ptr), pGB, pGq)


###################################################
### code chunk number 46: ex-score
###################################################
J <- 5L
N <- 4L

S <- crossprod(matrix(runif(J^2), nrow = J))
prm <- t(chol(S))[lower.tri(S, diag = TRUE)]

### define C
mC <- ltMatrices(matrix(prm, ncol = 1), diag = TRUE)

a <- matrix(runif(N * J), nrow = J) - 2
b <- a + 4
a[2,] <- -Inf
b[3,] <- Inf

M <- 10000L
W <- matrix(runif(M * (J - 1)), ncol = M)

lli <- c(lpmvnorm(a, b, chol = mC, w = W, M = M, logLik = FALSE))

fC <- function(prm) {
    C <- ltMatrices(matrix(prm, ncol = 1), diag = TRUE)
    lpmvnorm(a, b, chol = C, w = W, M = M)
}

sC <- slpmvnorm(a, b, chol = mC, w = W, M = M)

chk(lli, sC$logLik)

if (require("numDeriv", quietly = TRUE))
    chk(grad(fC, unclass(mC)), rowSums(unclass(sC$chol)), 
        check.attributes = FALSE)


###################################################
### code chunk number 47: ex-Lscore
###################################################
mL <- solve(mC)

lliL <- c(lpmvnorm(a, b, invchol = mL, w = W, M = M, logLik = FALSE))

chk(lli, lliL)

fL <- function(prm) {
    L <- ltMatrices(matrix(prm, ncol = 1), diag = TRUE)
    lpmvnorm(a, b, invchol = L, w = W, M = M)
}

sL <- slpmvnorm(a, b, invchol = mL, w = W, M = M)

chk(lliL, sL$logLik)

if (require("numDeriv", quietly = TRUE))
    chk(grad(fL, unclass(mL)), rowSums(unclass(sL$invchol)),
        check.attributes = FALSE)


###################################################
### code chunk number 48: ex-uni-score
###################################################
ptr <- pnorm(b[1,] / c(unclass(mC[,1]))) - pnorm(a[1,] / c(unclass(mC[,1])))
log(ptr)
lpmvnorm(a[1,,drop = FALSE], b[1,,drop = FALSE], chol = mC[,1], logLik = FALSE)
lapply(slpmvnorm(a[1,,drop = FALSE], b[1,,drop = FALSE], chol = mC[,1], 
                 logLik = TRUE), unclass)
sd1 <- c(unclass(mC[,1]))
(dnorm(b[1,] / sd1) * b[1,] - dnorm(a[1,] / sd1) * a[1,]) * (-1) / sd1^2 / ptr


###################################################
### code chunk number 49: chapterseed
###################################################
set.seed(110515)


###################################################
### code chunk number 50: ex-ML-dgp
###################################################
J <- 4
R <- diag(J)
R[1,2] <- R[2,1] <- .25
R[1,3] <- R[3,1] <- .5
R[2,4] <- R[4,2] <- .75
Sigma <- diag(sqrt(1:J / 2)) %*% R %*% diag(sqrt(1:J / 2))
C <- t(chol(Sigma))


###################################################
### code chunk number 51: ex-ML-C
###################################################
prm <- C[lower.tri(C, diag = TRUE)]
lt <- ltMatrices(matrix(prm, ncol = 1L), 
                 diag = TRUE,    ### has diagonal elements
                 byrow = FALSE)  ### prm is column-major
BYROW <- FALSE   ### later checks
lt <- ltMatrices(lt, 
                 byrow = BYROW)   ### convert to row-major
chk(C, as.array(lt)[,,1], check.attributes = FALSE)
chk(Sigma, as.array(tcrossprod(lt))[,,1], check.attributes = FALSE)


###################################################
### code chunk number 52: ex-ML-data
###################################################
N <- 100L
Z <- matrix(rnorm(N * J), nrow = J)
Y <- lt %*% Z + (mn <- 1:J)


###################################################
### code chunk number 53: ex-ML-mu-vcov
###################################################
rowMeans(Y)
(Shat <- var(t(Y)) * (N - 1) / N)


###################################################
### code chunk number 54: ex-ML-clogLik
###################################################
Yc <- Y - rowMeans(Y)

ll <- function(parm) {
    C <- ltMatrices(parm, diag = TRUE, byrow = BYROW)
    -ldmvnorm(obs = Yc, chol = C)
}

sc <- function(parm) {
    C <- ltMatrices(parm, diag = TRUE, byrow = BYROW)
    -rowSums(unclass(sldmvnorm(obs = Yc, chol = C)$chol))
}


###################################################
### code chunk number 55: ex-ML-const
###################################################
llim <- rep(-Inf, J * (J + 1) / 2)
llim[which(rownames(unclass(lt)) %in% paste(1:J, 1:J, sep = "."))] <- 1e-4


###################################################
### code chunk number 56: ex-ML-c
###################################################
if (BYROW) {
  cML <- chol(Shat)[upper.tri(Shat, diag = TRUE)]
} else {
  cML <- t(chol(Shat))[lower.tri(Shat, diag = TRUE)]
}
ll(cML)
start <- runif(length(cML))
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll, start), sc(start), check.attributes = FALSE)


###################################################
### code chunk number 57: ex-ML-coptim
###################################################
op <- optim(start, fn = ll, gr = sc, method = "L-BFGS-B", 
            lower = llim, control = list(trace = FALSE))
## ML numerically
ltMatrices(op$par, diag = TRUE, byrow = BYROW)
ll(op$par)
## ML analytically
t(chol(Shat))
ll(cML)
## true C matrix
lt


###################################################
### code chunk number 58: ex-ML-cens
###################################################
prb <- 1:9 / 10
sds <- sqrt(diag(Sigma))
ct <- sapply(1:J, function(j) qnorm(prb, mean = mn[j], sd = sds[j])) 
lwr <- upr <- Y
for (j in 1:J) {
    f <- cut(Y[j,], breaks = c(-Inf, ct[,j], Inf))
    lwr[j,] <- c(-Inf, ct[,j])[f]
    upr[j,] <- c(ct[,j], Inf)[f]
}


###################################################
### code chunk number 59: ex-ML-chk (eval = FALSE)
###################################################
## M <- floor(exp(0:25/10) * 1000)
## lGB <- sapply(M, function(m) {
##     st <- system.time(ret <- 
##         lpmvnormR(lwr, upr, mean = mn, chol = lt, algorithm = 
##                   GenzBretz(maxpts = m, abseps = 0, releps = 0)))
##     return(c(st["user.self"], ll = ret))
## })
## lH <- sapply(M, function(m) {
##     W <- NULL
##     if (require("qrng", quietly = TRUE))
##         W <- t(ghalton(m, d = J - 1))
##     st <- system.time(ret <- lpmvnorm(lwr, upr, mean = mn, 
##                                       chol = lt, w = W, M = m))
##     return(c(st["user.self"], ll = ret))
## })
## lHf <- sapply(M, function(m) {
##     W <- NULL
##     if (require("qrng", quietly = TRUE))
##         W <- t(ghalton(m, d = J - 1))
##     st <- system.time(ret <- lpmvnorm(lwr, upr, mean = mn, chol = lt, 
##                                       w = W, M = m, fast = TRUE))
##     return(c(st["user.self"], ll = ret))
## })


###################################################
### code chunk number 60: ex-ML-fig-data
###################################################
### use pre-computed data, otherwise CRAN complains.
M <-
c(1000, 1105, 1221, 1349, 1491, 1648, 1822, 2013, 2225, 2459, 
2718, 3004, 3320, 3669, 4055, 4481, 4953, 5473, 6049, 6685, 7389, 
8166, 9025, 9974, 11023, 12182)
lGB <- matrix(c(0.054, -880.492612, 0.054, -880.492426, 0.054, -880.492996, 0.054,
-880.492629, 0.054, -880.490231, 0.055, -880.492784, 0.054, -880.492632,
0.055, -880.489297, 0.054, -880.492516, 0.054, -880.491339, 0.054,
-880.492091, 0.11, -880.491601, 0.114, -880.493553, 0.111, -880.49125,
0.108, -880.492151, 0.108, -880.492275, 0.109, -880.491879, 0.109,
-880.492008, 0.192, -880.492132, 0.195, -880.491839, 0.194, -880.492139,
0.194, -880.491042, 0.198, -880.492198, 0.328, -880.4916, 0.323,
-880.491941, 0.323, -880.491698), nrow = 2)
rownames(lGB) <- c("user.self", "ll")
lH <- matrix(c(0.023, -880.480296, 0.027, -880.496166, 0.029, -880.488683,
0.032, -880.496171, 0.035, -880.485597, 0.039, -880.491333, 0.043,
-880.494557, 0.048, -880.495429, 0.053, -880.494391, 0.06, -880.485546,
0.064, -880.491455, 0.071, -880.494138, 0.079, -880.491619, 0.087,
-880.493393, 0.095, -880.492541, 0.106, -880.491649, 0.118, -880.492508,
0.129, -880.492558, 0.141, -880.492509, 0.157, -880.490448, 0.173,
-880.491686, 0.193, -880.491178, 0.211, -880.492286, 0.233, -880.491511,
0.258, -880.49153, 0.287, -880.491929), nrow = 2) 
rownames(lH) <- c("user.self", "ll")
lHf <- matrix(c(0.018, -880.487067, 0.019, -880.488639, 0.022, -880.488569,
0.024, -880.49393, 0.026, -880.486029, 0.029, -880.491563, 0.033,
-880.499415, 0.035, -880.494457, 0.038, -880.493954, 0.043, -880.493648,
0.047, -880.492955, 0.052, -880.494667, 0.059, -880.493745, 0.065,
-880.494195, 0.07, -880.49333, 0.078, -880.491451, 0.086, -880.492379,
0.094, -880.490392, 0.106, -880.491061, 0.115, -880.491577, 0.129,
-880.492523, 0.142, -880.491027, 0.158, -880.492086, 0.171, -880.492069,
0.189, -880.492251, 0.208, -880.492347), nrow = 2) 
rownames(lHf) <- c("user.self", "ll")


###################################################
### code chunk number 61: ex-ML-fig
###################################################
layout(matrix(1:2, nrow = 1))
plot(M, lGB["ll",], ylim = range(c(lGB["ll",], lH["ll",], lHf["ll",])), ylab = "Log-likelihood")
points(M, lH["ll",], pch = 4)
points(M, lHf["ll",], pch = 5)
plot(M, lGB["user.self",], ylim = c(0, max(lGB["user.self",])), ylab = "Time (in sec)")
points(M, lH["user.self",], pch = 4)
points(M, lHf["user.self",], pch = 5)
legend("bottomright", legend = c("pmvnorm", "lpmvnorm", "lpmvnorm(fast)"), pch = c(1, 4, 5), bty = "n")


###################################################
### code chunk number 62: ex-ML-ll
###################################################
M <- 500 
if (require("qrng", quietly = TRUE)) {
    ### quasi-Monte-Carlo
    W <- t(ghalton(M, d = J - 1))
} else {
    ### Monte-Carlo
    W <- matrix(runif(M * (J - 1)), nrow = J - 1, byrow = TRUE)
}
ll <- function(parm, J) {
     m <- parm[1:J]             ### mean parameters
     parm <- parm[-(1:J)]       ### chol parameters
     C <- matrix(c(parm), ncol = 1L)
     C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
     -lpmvnorm(lower = lwr, upper = upr, mean = m, chol = C, 
               w = W, M = M, logLik = TRUE)
}


###################################################
### code chunk number 63: ex-ML-check
###################################################
prm <- c(mn, unclass(lt))
ll(prm, J = J)
### ATLAS gives -880.4908, M1mac gives -880.4911
round(lpmvnormR(lwr, upr, mean = mn, chol = lt, 
                algorithm = GenzBretz(maxpts = M, abseps = 0, releps = 0)), 3)
(llprm <- lpmvnorm(lwr, upr, mean = mn, chol = lt, w = W, M = M))
chk(llprm, sum(lpmvnorm(lwr, upr, mean = mn, chol = lt, w = W, 
                        M = M, logLik = FALSE)))


###################################################
### code chunk number 64: ex-ML-sc
###################################################
sc <- function(parm, J) {
    m <- parm[1:J]             ### mean parameters
    parm <- parm[-(1:J)]       ### chol parameters
    C <- matrix(c(parm), ncol = 1L)
    C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
    ret <- slpmvnorm(lower = lwr, upper = upr, mean = m, chol = C, 
                     w = W, M = M, logLik = TRUE)
    return(-c(rowSums(ret$mean), rowSums(unclass(ret$chol))))
}


###################################################
### code chunk number 65: ex-ML-sc-chk
###################################################
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll, prm, J = J), sc(prm, J = J), check.attributes = FALSE)


###################################################
### code chunk number 66: ex-ML
###################################################
llim <- rep(-Inf, J + J * (J + 1) / 2)
llim[J + which(rownames(unclass(lt)) %in% paste(1:J, 1:J, sep = "."))] <- 1e-4

if (BYROW) {
  start <- c(rowMeans(Y), chol(Shat)[upper.tri(Shat, diag = TRUE)])
} else {
  start <- c(rowMeans(Y), t(chol(Shat))[lower.tri(Shat, diag = TRUE)])
}

ll(start, J = J)

op <- optim(start, fn = ll, gr = sc, J = J, method = "L-BFGS-B", 
            lower = llim, control = list(trace = FALSE))

op$value ## compare with 
ll(prm, J = J)


###################################################
### code chunk number 67: ex-ML-C
###################################################
(C <- ltMatrices(matrix(op$par[-(1:J)], ncol = 1), 
                 diag = TRUE, byrow = BYROW))
lt


###################################################
### code chunk number 68: ex-ML-mu
###################################################
op$par[1:J]
mn


###################################################
### code chunk number 69: ex-ML-Shat
###################################################
### ATLAS print issues
round(tcrossprod(lt), 4)  ### true Sigma
tcrossprod(C)             ### interval-censored obs
Shat                      ### "exact" obs


###################################################
### code chunk number 70: regressions
###################################################
c(cond_mvnorm(chol = C, which_given = 2:J, given = diag(J - 1))$mean)


###################################################
### code chunk number 71: regressionsC
###################################################
c(cond_mvnorm(chol = aperm(as.chol(C), perm = c(2:J, 1)),
              which_given = 1:(J - 1), given = diag(J - 1))$mean)


###################################################
### code chunk number 72: regressionsP
###################################################
x <- as.array(chol2pre(aperm(as.chol(C), perm = c(2:J, 1))))[J,,1]
c(-x[-J] / x[J])


###################################################
### code chunk number 73: lm-ex
###################################################
dY <- as.data.frame(t(Y))
colnames(dY) <- paste0("Y", 1:J)
coef(m1 <- lm(Y1 ~ ., data = dY))[-1L]


###################################################
### code chunk number 74: hessian
###################################################
H <- optim(op$par, fn = ll, gr = sc, J = J, method = "L-BFGS-B", 
           lower = llim, hessian = TRUE, 
           control = list(trace = FALSE))$hessian


###################################################
### code chunk number 75: ML-sample
###################################################
L <- try(t(chol(H)))
### some check on r-oldrel-macos-arm64
if (inherits(L, "try-error"))
    L <- t(chol(H + 1e-4 * diag(nrow(H))))
L <- ltMatrices(L[lower.tri(L, diag = TRUE)], diag = TRUE)
Nsim <- 50000
Z <- matrix(rnorm(Nsim * nrow(H)), ncol = Nsim)
rC <- solve(L, Z)[-(1:J),] + op$par[-(1:J)] ### remove mean parameters


###################################################
### code chunk number 76: ML-check
###################################################
c(sqrt(rowMeans((rC - rowMeans(rC))^2)))
c(sqrt(diagonals(crossprod(solve(L)))))


###################################################
### code chunk number 77: rC
###################################################
rC <- ltMatrices(rC, diag = TRUE)


###################################################
### code chunk number 78: ML-beta
###################################################
rbeta <- cond_mvnorm(chol = rC, which_given = 2:J, given = diag(J - 1))$mean
sqrt(rowMeans((rbeta - rowMeans(rbeta))^2))


###################################################
### code chunk number 79: se-ex
###################################################
sqrt(diag(vcov(m1)))[-1L]


###################################################
### code chunk number 80: ex-ML-cd
###################################################
ic <- 1:2 	### position of continuous variables
ll_cd <- function(parm, J) {
     m <- parm[1:J]             ### mean parameters
     parm <- parm[-(1:J)]       ### chol parameters
     C <- matrix(c(parm), ncol = 1L)
     C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
     -ldpmvnorm(obs = Y[ic,], lower = lwr[-ic,], 
                upper = upr[-ic,], mean = m, chol = C, 
                w = W[-ic,,drop = FALSE], M = M)
}
sc_cd <- function(parm, J) {
    m <- parm[1:J]             ### mean parameters
    parm <- parm[-(1:J)]       ### chol parameters
    C <- matrix(c(parm), ncol = 1L)
    C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
    ret <- sldpmvnorm(obs = Y[ic,], lower = lwr[-ic,],
                      upper = upr[-ic,], mean = m, chol = C, 
                      w = W[-ic,,drop = FALSE], M = M)
    return(-c(rowSums(ret$mean),
              rowSums(Lower_tri(ret$chol, diag = TRUE))))
}


###################################################
### code chunk number 81: ex-ML-cd-score
###################################################
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll_cd, start, J = J), sc_cd(start, J = J), 
        check.attributes = FALSE, tolerance = 1e-6)


###################################################
### code chunk number 82: ex-ML-cd-optim
###################################################
op <- optim(start, fn = ll_cd, gr = sc_cd, J = J, 
            method = "L-BFGS-B", lower = llim, 
            control = list(trace = FALSE))
## estimated C
ltMatrices(matrix(op$par[-(1:J)], ncol = 1), 
           diag = TRUE, byrow = BYROW)
## compare with true C
lt
## estimated means
op$par[1:J]
## compare with true means
mn


###################################################
### code chunk number 83: ex-ML-ap
###################################################
### discrete variables first
perm <- c((1:J)[-ic], ic)
ll_ap <- function(parm, J) {
     m <- parm[1:J]             ### mean parameters; NOT permuted
     parm <- parm[-(1:J)]       ### chol parameters
     C <- matrix(c(parm), ncol = 1L)
     C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
     Ct <- aperm(as.chol(C), perm = perm)
     -ldpmvnorm(obs = Y[ic,], lower = lwr[-ic,], 
                upper = upr[-ic,], mean = m, chol = Ct, 
                w = W[-ic,,drop = FALSE], M = M)
}


###################################################
### code chunk number 84: ex-ML-ap-score
###################################################
sc_ap <- function(parm, J) {
    m <- parm[1:J]               ### mean parameters; NOT permuted
    parm <- parm[-(1:J)]         ### chol parameters
    C <- matrix(c(parm), ncol = 1L)
    C <- ltMatrices(C, diag = TRUE, byrow = BYROW)
    ### permutation
    Ct <- aperm(as.chol(C), perm = perm)
    ret <- sldpmvnorm(obs = Y[ic,], lower = lwr[-ic,],
                      upper = upr[-ic,], mean = m, chol = Ct, 
                      w = W[-ic,,drop = FALSE], M = M)
    ### undo permutation for chol
    retC <- deperma(chol = C, permuted_chol = Ct, 
                   perm = perm, score_schol = ret$chol)
    return(-c(rowSums(ret$mean),
              rowSums(Lower_tri(retC, diag = TRUE))))
}


###################################################
### code chunk number 85: ex-ML-ap-grad
###################################################
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll_ap, start, J = J), sc_ap(start, J = J), 
        check.attributes = FALSE, tolerance = 1e-6)


###################################################
### code chunk number 86: ex-ML-ap-optim-
###################################################
op <- optim(start, fn = ll_ap, gr = sc_ap, J = J, 
            method = "L-BFGS-B", lower = llim, 
            control = list(trace = FALSE))
## estimated C for (X, Y)
ltMatrices(matrix(op$par[-(1:J)], ncol = 1), 
           diag = TRUE, byrow = BYROW)
## compare with true _permuted_ C for (X, Y)
round(as.array(aperm(as.chol(lt), perm = perm)), 4)


###################################################
### code chunk number 87: ex-stand
###################################################
C <- ltMatrices(runif(10))
chk(as.array(chol2cov(standardize(chol = C))),
    as.array(chol2cor(standardize(chol = C))))
L <- solve(C)
chk(as.array(invchol2cov(standardize(invchol = L))),
    as.array(invchol2cor(standardize(invchol = L))))


###################################################
### code chunk number 88: gc-classical
###################################################
data("iris", package = "datasets")
J <- 4
Z <- t(qnorm(do.call("cbind", lapply(iris[1:J], rank, ties.method = "max")) / 
       (nrow(iris) + 1)))
(CR <- cor(t(Z)))
ll <- function(parm) {
    C <- ltMatrices(parm)
    Cs <- standardize(chol = C)
    -ldmvnorm(obs = Z, chol = Cs)
}
sc <- function(parm) {
    C <- ltMatrices(parm)
    Cs <- standardize(chol = C)
    -rowSums(Lower_tri(destandardize(chol = C, 
        score_schol = sldmvnorm(obs = Z, chol = Cs)$chol)))
}
start <- t(chol(CR))
start <- start[lower.tri(start)]
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll, start), sc(start), check.attributes = FALSE)
op <- optim(start, fn = ll, gr = sc, method = "BFGS", 
            control = list(trace = FALSE), hessian = TRUE)
op$value
S_ML <- chol2cov(standardize(chol = ltMatrices(op$par)))


###################################################
### code chunk number 89: gc-NPML
###################################################
lwr <- do.call("cbind", lapply(iris[1:J], rank, ties.method = "min")) - 1L
upr <- do.call("cbind", lapply(iris[1:J], rank, ties.method = "max"))
lwr <- t(qnorm(lwr / nrow(iris)))
upr <- t(qnorm(upr / nrow(iris)))

M <- 500 
if (require("qrng", quietly = TRUE)) {
    ### quasi-Monte-Carlo
    W <- t(ghalton(M, d = J - 1))
} else {
    ### Monte-Carlo
    W <- matrix(runif(M * (J - 1)), nrow = J - 1, byrow = TRUE)
}

ll <- function(parm) {
    C <- ltMatrices(parm)
    Cs <- standardize(chol = C)
    -lpmvnorm(lower = lwr, upper = upr, chol = Cs, M = M, w = W)
}
sc <- function(parm) {
    C <- ltMatrices(parm)
    Cs <- standardize(chol = C)
    -rowSums(Lower_tri(destandardize(chol = C, 
        score_schol = slpmvnorm(lower = lwr, upper = upr, chol = Cs, 
                               M = M, w = W)$chol)))
}
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll, start), sc(start), check.attributes = FALSE)
op2 <- optim(start, fn = ll, gr = sc, method = "BFGS", 
             control = list(trace = FALSE), hessian = TRUE)
S_NPML <- chol2cov(standardize(chol = ltMatrices(op2$par)))


###################################################
### code chunk number 90: gc
###################################################
S_ML
S_NPML


###################################################
### code chunk number 91: gc-se
###################################################
sd_ML <- ltMatrices(sqrt(diag(solve(op$hessian))))
diagonals(sd_ML) <- 0
sd_NPML <- try(ltMatrices(sqrt(diag(solve(op2$hessian)))))
if (!inherits(sd_NPML, "try-error")) {
    diagonals(sd_NPML) <- 0
    print(sd_ML)
    print(sd_NPML)
}


###################################################
### code chunk number 92: concave-L
###################################################
J <- 5
N <- 100
### mean
m <- rnorm(J)
L <- ltMatrices(prm <- runif(J * (J + 1) / 2), diag = TRUE)
Z <- matrix(rnorm(N * J), nrow = J)
Y <- solve(L, Z) + m
### scaled mean; FIXME: update example using invcholmean args
d <- L %*% m

nll <- function(parm) {
    d <- parm[seq_len(J)]
    L <- ltMatrices(parm[-seq_len(J)], diag = TRUE)
    -ldmvnorm(obs = Y, mean = solve(L, d), invchol = L)
}

start <- c(d, prm)

nll(start)
### identical
-ldmvnorm(obs = Y, mean = m, invchol = L)


###################################################
### code chunk number 93: scores-L
###################################################
nsc <- function(parm) {
    d <- parm[seq_len(J)]
    L <- ltMatrices(parm[-seq_len(J)], diag = TRUE)
    ret <- sldmvnorm(obs = Y, mean = solve(L, d), invchol = L)
    C <- solve(L)

    J <- dim(L)[2L]
    M <- matrix(seq_len(J^2), nrow = J, byrow = FALSE)
    idx <- M[lower.tri(M, diag = TRUE)]
   
    X <- -ret$obs
    Y <- matrix(d, nrow = nrow(X), ncol = ncol(X))
    A <- X[rep(1:nrow(X), times = nrow(X)),,drop = FALSE] * 
         Y[rep(1:nrow(Y), each = nrow(X)),,drop = FALSE]

    scL <- - vectrick(C, A)
    scL <- scL[idx,,drop = FALSE]
    scL <- rowSums(unclass(scL) + unclass(ret$invchol))
    - c(rowSums(solve(L, -ret$obs, transpose = TRUE)), 
        scL)
}
chk(unname(nsc(start)), grad(nll, start))


###################################################
### code chunk number 94: concave-C
###################################################
C <- ltMatrices(prm <- runif(J * (J + 1) / 2), diag = TRUE)
Z <- matrix(rnorm(N * J), nrow = J)
Y <- C %*% Z + m
### scaled mean
d <- solve(C, m)

nll <- function(parm) {
    d <- parm[seq_len(J)]
    C <- ltMatrices(parm[-seq_len(J)], diag = TRUE)
    -ldmvnorm(obs = Y, mean = C %*% d, chol = C)
}

start <- c(d, prm)

nll(start)
### identical
-ldmvnorm(obs = Y, mean = m, chol = C)


###################################################
### code chunk number 95: scores-C
###################################################
nsc <- function(parm) {
    d <- parm[seq_len(J)]
    C <- ltMatrices(parm[-seq_len(J)], diag = TRUE)
    ret <- sldmvnorm(obs = Y, mean = C %*% d, chol = C)

    J <- dim(C)[2L]
    M <- matrix(seq_len(J^2), nrow = J, byrow = FALSE)
    idx <- M[lower.tri(M, diag = TRUE)]

    X <- -ret$obs
    Y <- matrix(d, nrow = nrow(X), ncol = ncol(X))
    A <- X[rep(1:nrow(X), times = nrow(X)),,drop = FALSE] * 
         Y[rep(1:nrow(Y), each = nrow(X)),,drop = FALSE]

    scC <- A[idx,,drop = FALSE]
    scC <- rowSums(unclass(scC) + unclass(ret$chol))
    - c(rowSums(crossprod(C, -ret$obs)), 
        scC)
}
chk(unname(nsc(start)), grad(nll, start))


###################################################
### code chunk number 96: iris-model
###################################################
data("iris", package = "datasets")
vars <- names(iris)[-5L]
N <- nrow(iris)
m <- colMeans(iris[,vars])
V <- var(iris[,vars]) * (N - 1) / N
iris_mvn <- mvnorm(mean = m, chol = t(chol(V)))
iris_var <- simulate(iris_mvn, nsim = nrow(iris))


###################################################
### code chunk number 97: iris-mc
###################################################
j <- 3:4
margDist(iris_mvn, which = vars[j])
gm <- t(iris[,vars[-(j)]])
iris_cmvn <- condDist(iris_mvn, which_given = vars[j], given = gm)


###################################################
### code chunk number 98: iris-ll
###################################################
logLik(object = iris_cmvn, obs = t(iris[,vars[-j]]))


###################################################
### code chunk number 99: iris-ll-perm
###################################################
logLik(object = iris_cmvn, obs = t(iris[,rev(vars[-j])]))


###################################################
### code chunk number 100: iris-lLgrad
###################################################
J <- length(vars)
obs <- t(iris[, vars])
ll <- function(parm) {
    C <- ltMatrices(parm[-(1:J)], diag = TRUE, names = vars)
    x <- mvnorm(mean = parm[1:J], chol = C)
    -logLik(object = x, obs = obs)
}
sc <- function(parm) {
    C <- ltMatrices(parm[-(1:J)], diag = TRUE, names = vars)
    x <- mvnorm(mean = parm[1:J], chol = C)
    ret <- lLgrad(object = x, obs = obs)
    -c(rowSums(ret$mean), rowSums(Lower_tri(ret$scale, diag = TRUE)))
}


###################################################
### code chunk number 101: iris-ML
###################################################
### don't start at the solution
start <- round(c(c(iris_mvn$mean), 
                 Lower_tri(iris_mvn$scale, diag = TRUE)), 2)
llim <- rep(-Inf, J + J * (J + 1) / 2)
llim[J + c(diagonals(ltMatrices(seq_len(J * (J + 1) / 2), diag = TRUE)))] <- 1e-4
op <- optim(start, fn = ll, gr = sc, method = "L-BFGS-B", 
            lower = llim, 
            control = list(trace = FALSE, 
                           factr = 1e-6)) ### noLD machines
Chat <- ltMatrices(op$par[-(1:J)], diag = TRUE, names = vars)
ML <- mvnorm(mean = op$par[1:J], chol = Chat)


###################################################
### code chunk number 102: iris-ML-hat
###################################################
### covariance (noLD brings up small differences)
round(vcov(ML), 3)
V
### mean
mean(ML)[,,drop = TRUE]
m


###################################################
### code chunk number 103: iris-lLgrad-nu
###################################################
ll <- function(parm, logLik = TRUE) {
    L <- ltMatrices(parm[-(1:J)], diag = TRUE, names = vars)
    x <- mvnorm(invcholmean = parm[1:J], invchol = L)
    if (!logLik) return(x)
    -logLik(object = x, obs = obs)
}
sc <- function(parm) {
    x <- ll(parm, logLik = FALSE)
    ret <- lLgrad(object = x, obs = obs)
    -c(rowSums(ret$invcholmean), rowSums(Lower_tri(ret$scale, diag = TRUE)))
}
### note: This is a convex problem now, so (here incorrect) 
### starting values shouldn't matter
opL <- optim(start, fn = ll, gr = sc, method = "L-BFGS-B", 
            lower = llim, control = list(trace = FALSE, factr = 1e-6))
MLL <- ll(opL$par, logLik = FALSE)


###################################################
### code chunk number 104: iris-ML-hat-nu
###################################################
### log-likelihood
op$value
opL$value
### covariance
round(vcov(MLL), 3)
V
### mean
mean(MLL)[,,drop = TRUE]
m


###################################################
### code chunk number 105: iris-interval
###################################################
v1 <- vars[1]
q1 <- quantile(iris[[v1]], probs = 1:4 / 5)
head(f1 <- cut(iris[[v1]], breaks = c(-Inf, q1, Inf)))
lower <- matrix(c(-Inf, q1)[f1], nrow = 1)
upper <- matrix(c(q1, Inf)[f1], nrow = 1)
rownames(lower) <- rownames(upper) <- v1
obs <- obs[!rownames(obs) %in% v1,,drop = FALSE]


###################################################
### code chunk number 106: iris-MLi
###################################################
ll <- function(parm, logLik = TRUE) {
    L <- ltMatrices(parm[-(1:J)], diag = TRUE, names = vars)
    x <- mvnorm(invcholmean = parm[1:J], invchol = L)
    if (!logLik) return(x)
    -logLik(object = x, obs = obs, lower = lower, upper = upper, 
            tol = 1e-6) ### probs < tol are considered 0
}
sc <- function(parm) {
    x <- ll(parm, logLik = FALSE)
    ret <- lLgrad(object = x, obs = obs, lower = lower, upper = upper, 
                  tol = 1e-6) ### probs < tol are considered 0
    -c(rowSums(ret$invcholmean, na.rm = TRUE), 
       rowSums(Lower_tri(ret$scale, diag = TRUE), na.rm = TRUE))
}


###################################################
### code chunk number 107: iris-MLi-opt
###################################################
start <- round(opL$par, 2)
if (require("numDeriv", quietly = TRUE))
    chk(grad(ll, start), sc(start), check.attributes = FALSE)
opi <- optim(start, fn = ll, gr = sc, method = "L-BFGS-B", 
             lower = llim, control = list(trace = FALSE, factr = 1e-6))
MLi <- ll(opi$par, logLik = FALSE)


###################################################
### code chunk number 108: iris-MLi-hat
###################################################
### covariance
round(vcov(MLi), 3)
round(vcov(MLL), 3)
round(vcov(ML), 3)
### mean
mean(MLi)[,,drop = TRUE]
mean(MLL)[,,drop = TRUE]
mean(ML)[,,drop = TRUE]


###################################################
### code chunk number 109: iris-lm
###################################################
### least-squares coefficients
coef(irislm <- lm(Petal.Width ~ Sepal.Length + Sepal.Width + Petal.Length, 
                  data = iris))
### residual standard deviation
summary(irislm)$sigma
### compare with 
round(coef(ML, which = "Petal.Width"), 3)


###################################################
### code chunk number 110: lmvnorm_src.Rnw:9120-9122
###################################################
# slightly different results on noLD machines
cat("> ## IGNORE_RDIFF_BEGIN\n")


###################################################
### code chunk number 111: iris-lm-iL
###################################################
### nu, L for exact observations
round(coef(MLL, which = "Petal.Width"), 3)
### nu, L for censored observations
round(coef(MLi, which = "Petal.Width"), 3)


###################################################
### code chunk number 112: lmvnorm_src.Rnw:9130-9132
###################################################
## slightly different results on noLD machines
cat("> ## IGNORE_RDIFF_END\n")


###################################################
### code chunk number 113: marginB
###################################################
N <- 3
J <- 4
L <- ltMatrices(runif(J * (J + 1) / 2), diag = TRUE, names = LETTERS[1:J])
Z <- matrix(rnorm(J * N), nrow = J)
Y <- solve(L, Z)

lwrA <- matrix(-1, nrow = 1, ncol = N)
uprA <- matrix(1, nrow = 1, ncol = N)
rownames(lwrA) <- rownames(uprA) <- "A"

lwrB <- matrix(-Inf, nrow = 1, ncol = N)
uprB <- matrix(Inf, nrow = 1, ncol = N)
rownames(lwrB) <- rownames(uprB) <- "B"

lwr <- rbind(lwrA, lwrB)
upr <- rbind(uprA, uprB)
obs <- Y[rev(LETTERS[3:J]),]    ### change order of dimensions


###################################################
### code chunk number 114: marginBllsc
###################################################
w <- matrix(runif(1000), nrow = 1, byrow = TRUE)
lABCD <- logLik(mvnorm(invchol = L), obs = obs, lower = lwr, upper = upr, w = w)
sABCD <- lLgrad(mvnorm(invchol = L), obs = obs, lower = lwr, upper = upr, w = w)


###################################################
### code chunk number 115: marginllsc
###################################################
lACD <- logLik(mvnorm(invchol = L), obs = obs, lower = lwrA, upper = uprA)
sACD <- lLgrad(mvnorm(invchol = L), obs = obs, lower = lwrA, upper = uprA)


###################################################
### code chunk number 116: marginchk
###################################################
chk(lABCD, lACD)
nm <- names(sABCD)
nm <- nm[!nm %in% c("lower", "upper")]
chk(sABCD[nm], sACD[nm])


###################################################
### code chunk number 117: marginsc
###################################################
chk(sABCD$lower["A",,drop = FALSE], sACD$lower)
chk(sABCD$upper["A",,drop = FALSE], sACD$upper)
sABCD$lower["B",]	### zero
sABCD$upper["B",]	### zero


###################################################
### code chunk number 118: RR-ll
###################################################
J <- 6
K <- 3
B <- matrix(rnorm(J * K), nrow = J)
D <- runif(J)
S <- tcrossprod(B) + diag(D)
Linv <- t(chol(S))
Linv <- ltMatrices(Linv[lower.tri(Linv, diag = TRUE)], diag = TRUE)
a <- -(2 + runif(J))
b <- 2 + runif(J)
M <- 1e6
dim(w <- matrix(runif((J - 1) * M), nrow = J - 1, byrow = TRUE))
lpmvnorm(lower = a, upper = b, chol = Linv, w = w)
dim(Z <- matrix(rnorm(K * M), nrow = K))
lpRR(lower = a, upper = b, B = B, D = D, Z = Z)


###################################################
### code chunk number 119: RR-sc
###################################################
smv <- slpmvnorm(lower = a, upper = b, chol = Linv, w = w)
sRR <- slpRR(lower = a, upper = b, B = B, D = D, Z = Z)
chk(c(smv$lower), sRR$lower, tolerance = 1e-2)
chk(c(smv$upper), sRR$upper, tolerance = 1e-2)
chk(c(smv$mean), sRR$mean, tolerance = 1e-2 * 2)


###################################################
### code chunk number 120: RR-sc-BD
###################################################
Z <- matrix(rnorm(K * 1000), nrow = K)
lB <- function(B) lpRR(lower = a, upper = b, B = B, D = D, Z = Z)
gB <- grad(lB, B)
sRR <- slpRR(lower = a, upper = b, B = B, D = D, Z = Z)
chk(gB, c(sRR$B), tolerance = 1e-3)
lD <- function(D) lpRR(lower = a, upper = b, B = B, D = D, Z = Z)
gD <- grad(lD, D)
chk(gD, c(sRR$D), tolerance = 1e-3)
### while we are at it, check lower and again
llwr <- function(a) lpRR(lower = a, upper = b, B = B, D = D, Z = Z)
glwr <- grad(llwr, a)
chk(glwr, c(sRR$lower))
lupr <- function(b) lpRR(lower = a, upper = b, B = B, D = D, Z = Z)
gupr <- grad(lupr, b)
chk(gupr, c(sRR$upper))
