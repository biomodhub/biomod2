library(dse1)
library(vars)
## A-model
Apoly   <- array(c(1.0, -0.5, 0.3, 0.8,
                   0.2, 0.1, -0.7, -0.2,
                   0.7, 1, 0.5, -0.3) ,
                 c(3, 2, 2))
## Setting covariance to identity-matrix
B <- diag(2)
## Generating the VAR(2) model 
svarA  <- ARMA(A = Apoly, B = B)
## Simulating 500 observations
svarsim <- simulate(svarA, sampleT = 500,
                    rng = list(seed = c(123456)))
## Obtaining the generated series
svardat <- matrix(svarsim$output, nrow = 500, ncol = 2)
colnames(svardat) <- c("y1", "y2")
## Estimating the VAR
varest <- VAR(svardat, p = 2, type = "none")
## Setting up matrices for A-model
Amat <- diag(2)
Amat[2, 1] <- NA
Amat[1, 2] <- NA
## Estimating the SVAR A-type by direct maximisation
## of the log-likelihood
args(SVAR)
svar.A <- SVAR(varest, estmethod = "direct",
               Amat = Amat, hessian = TRUE)
