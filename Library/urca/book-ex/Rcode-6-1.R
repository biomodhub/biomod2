set.seed(123456)
e <- rnorm(500)
## trend
trd <- 1:500
S <- c(rep(0, 249), rep(1, 251))
## random walk with drift
y1 <- 0.1*trd + cumsum(e)
## random walk with drift and shift
y2 <- 0.1*trd + 10*S + cumsum(e)
