library(dse1)
library(vars)
## B-model
Apoly   <- array(c(1.0, -0.5, 0.3, 0,
                   0.2, 0.1, 0, -0.2,
                   0.7, 1, 0.5, -0.3) ,
                 c(3, 2, 2))
## Setting covariance to identity-matrix
B <- diag(2)
B[2, 1] <- -0.8
## Generating the VAR(2) model 
svarB  <- ARMA(A = Apoly, B = B)
## Simulating 500 observations
svarsim <- simulate(svarB, sampleT = 500,
                    rng = list(seed = c(123456)))
svardat <- matrix(svarsim$output, nrow = 500, ncol = 2)
colnames(svardat) <- c("y1", "y2")
varest <- VAR(svardat, p = 2, type = "none")
## Estimating the SVAR B-type by scoring algorithm
## Setting up the restriction matrix and vector
## for B-model
Bmat <- diag(2)
Bmat[2, 1] <- NA
svar.B <- SVAR(varest, estmethod = "scoring",
               Bmat = Bmat, max.iter = 200)
