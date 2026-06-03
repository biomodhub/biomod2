vecm <- ca.jo(Canada[, c("rw", "prod", "e", "U")],
              type = "trace", ecdet = "trend",
              K = 3, spec = "transitory")
vecm.r1 <- cajorls(vecm, r = 1)
alpha <- coef(vecm.r1$rlm)[1, ]
beta <- vecm.r1$beta
resids <- resid(vecm.r1$rlm)
N <- nrow(resids)
sigma <- crossprod(resids) / N
## t-stats for alpha 
alpha.se <- sqrt(solve(crossprod(
                 cbind(vecm@ZK %*% beta, vecm@Z1)))
                 [1, 1]* diag(sigma))
alpha.t <- alpha / alpha.se
## t-stats for beta
beta.se <- sqrt(diag(kronecker(solve(
                crossprod(vecm@RK[, -1])),
                solve(t(alpha) %*% solve(sigma)
                %*% alpha))))
beta.t <- c(NA, beta[-1] / beta.se)
