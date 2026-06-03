### R code from vignette source 'vinny.Rnw'

###################################################
### code chunk number 1: try1a
###################################################
library(rcdd)
d <- 3
# unit simplex in H-representation
qux <- makeH(- diag(d), rep(0, d), rep(1, d), 1)
print(qux)


###################################################
### code chunk number 2: try1b
###################################################
# unit simplex in V-representation
out <- scdd(qux)
print(out)


###################################################
### code chunk number 3: try1c
###################################################
# unit simplex in H-representation
# note: different from original, but equivalent
out <- scdd(out$output)
print(out)


###################################################
### code chunk number 4: try2a
###################################################
# add equality constraint
quux <- addHeq(1:d, 2.2, qux)
print(quux)
out <- scdd(quux)
print(out)


###################################################
### code chunk number 5: try3a
###################################################
quuxq <- d2q(quux)
print(quuxq)


###################################################
### code chunk number 6: try3b
###################################################
bar <- as.numeric(unlist(strsplit(quuxq[5,2], "/")))
print(bar)
bar[1] / bar[2]


###################################################
### code chunk number 7: try3c
###################################################
q2d(quuxq)


###################################################
### code chunk number 8: try3d
###################################################
outq <- scdd(quuxq)
print(outq)


###################################################
### code chunk number 9: try3d
###################################################
print(q2d(outq$output))


###################################################
### code chunk number 10: try3e
###################################################
quuxq <- z2q(round(quux * 10), rep(10, length(quux)))
print(quuxq)
outq <- scdd(quuxq)
print(outq)


###################################################
### code chunk number 11: try3f
###################################################
qmq(outq$output, out$output)


###################################################
### code chunk number 12: try4a
###################################################
d <- 4
n <- 100
set.seed(42)
x <- matrix(rnorm(d * n), nrow = n)
foo <- makeV(d2q(x))
out <- scdd(foo)
l <- out$output[ , 1]
b <- out$output[ , 2]
v <- out$output[ , - c(1, 2)]
a <- qneg(v)


###################################################
### code chunk number 13: try4b
###################################################
axb <- qmatmult(a, t(x))
axb <- sweep(axb, 1, b, FUN = qmq)
fred <- apply(axb, 2, function(foo) max(qsign(foo)))

all(fred <= 0)
sum(fred < 0)
sum(fred == 0)


###################################################
### code chunk number 14: try4c
###################################################
y <- matrix(rnorm(2 * n * d), nrow = 2 * n)
ayb <- qmatmult(a, t(d2q(y)))
ayb <- sweep(ayb, 1, b, FUN = qmq)
sally <- apply(ayb, 2, function(foo) max(qsign(foo)))

sum(sally < 0)
sum(sally == 0)
sum(sally > 0)


###################################################
### code chunk number 15: try5a
###################################################
hrep <- rbind(c("0", "0", "1", "1", "0", "0"),
              c("0", "0", "0", "2", "0", "0"),
              c("1", "3", "0", "-1", "0", "0"),
              c("1", "9/2", "0", "0", "-1", "-1"))
print(hrep)
a <- c("2", "3/5", "0", "0")
out <- lpcdd(hrep, a)
print(out)


###################################################
### code chunk number 16: try5a-chk1
###################################################
qsum(qxq(a, out$primal.solution))


###################################################
### code chunk number 17: try5a-chk2
###################################################
xbar <- out$primal.solution
foo <- qmatmult(hrep[ , - c(1, 2)], cbind(xbar))
foo <- qpq(hrep[ , 2], foo)
print(foo)


###################################################
### code chunk number 18: try5a-chk3
###################################################
qxq(foo, out$dual.solution)


###################################################
### code chunk number 19: try5a-chk4
###################################################
qpq(a, qmatmult(rbind(out$dual.solution), hrep[ , -c(1, 2)]))


###################################################
### code chunk number 20: try5b
###################################################
hrep <- rbind(c("0", "0", "1", "0"),
              c("0", "0", "0", "1"),
              c("0", "-2", "-1", "-1"))
print(hrep)
a <- c("1", "1")
out <- lpcdd(hrep, a)
print(out)


###################################################
### code chunk number 21: try5c
###################################################
hrep <- rbind(c("0", "0", "1", "0"),
              c("0", "0", "0", "1"))
print(hrep)
a <- c("1", "1")
out <- lpcdd(hrep, a, minimize = FALSE)
print(out)


###################################################
### code chunk number 22: try5c-chk1
###################################################
qmatmult(hrep[ , - c(1, 2)], cbind(out$primal.direction))


###################################################
### code chunk number 23: try5c-chk2
###################################################
qsum(qxq(a, out$primal.direction))


###################################################
### code chunk number 24: interior-lp
###################################################
xin <- x[fred < 0, , drop = FALSE]
qin <- xin[sample(nrow(xin), 1), ]
qin
hrep <- cbind(0, 0, 1, - x)
hrep <- rbind(hrep, c(0, 1, 1, - qin))

out <- lpcdd(d2q(hrep), d2q(c(-1, qin)), minimize = FALSE)
out$optimal.value


###################################################
### code chunk number 25: interior-lp
###################################################
yout <- y[sally > 0, , drop = FALSE]
qout <- yout[sample(nrow(yout), 1), ]
qout
hrep <- cbind(0, 0, 1, - x)
hrep <- rbind(hrep, c(0, 1, 1, - qout))

out <- lpcdd(d2q(hrep), d2q(c(-1, qout)), minimize = FALSE)
out$optimal.value


###################################################
### code chunk number 26: toy
###################################################
hrep <- rbind(c(0, 0,  1,  1,  0),
              c(0, 0, -1,  0,  0),
              c(0, 0,  0, -1,  0),
              c(0, 0,  0,  0, -1),
              c(0, 0, -1, -1, -1))
print(hrep)
redundant(hrep, representation = "H")


###################################################
### code chunk number 27: vertex
###################################################
foo <- makeV(points = d2q(x))
out <- redundant(foo)
nrow(out$output)
all((out$new.position == 0) == (fred < 0))


###################################################
### code chunk number 28: faces
###################################################
vrep <- rbind(c(0, 1,  1,  1, 0),
              c(0, 1,  1, -1, 0),
              c(0, 1, -1,  1, 0),
              c(0, 1, -1, -1, 0),
              c(0, 1,  0,  0, 1))
print(vrep)
hrep <- scdd(vrep, rep = "V")$output
print(hrep)


###################################################
### code chunk number 29: faces-numbers
###################################################
out <- allfaces(hrep)
d <- unlist(out$dimension)
nd <- tabulate(d + 1)
names(nd) <- seq(0, 3)
print(nd)


###################################################
### code chunk number 30: faces-numbers
###################################################
asl <- sapply(out$active.set, paste, collapse = " ")
names(asl) <- d
asl <- asl[order(d)]
print(asl)


###################################################
### code chunk number 31: agresti-one
###################################################
x <- seq(10, 90, 10)
x <- x[x != 50]
x
y <- as.numeric(x > 50)
y


###################################################
### code chunk number 32: agresti-one-generators
###################################################
yy <- matrix(0:1, nrow = 2, ncol = length(x))
colnames(yy) <- paste0("y", x)
yy <- expand.grid(as.data.frame(yy))
head(yy)
nrow(yy)


###################################################
### code chunk number 33: agresti-one-v-rep
###################################################
yy <- as.matrix(yy) # was data frame
yy.trans = yy %*% cbind(1, x)
dim(yy.trans)


###################################################
### code chunk number 34: agresti-one-convex-hull
###################################################
foo <- makeV(points = d2q(yy.trans))
out <- redundant(foo)
nrow(out$output)
yy.trans <- out$output[ , - c(1, 2)]
dim(yy.trans)


###################################################
### code chunk number 35: fig1plot
###################################################
plot(yy.trans[ , 1], yy.trans[ , 2], xlab = "sum(y)", ylab = "sum(x * y)")


###################################################
### code chunk number 36: fig1
###################################################
plot(yy.trans[ , 1], yy.trans[ , 2], xlab = "sum(y)", ylab = "sum(x * y)")


###################################################
### code chunk number 37: agresti-one-v-rep
###################################################
out <- scdd(out$output)
nrow(out$output)


