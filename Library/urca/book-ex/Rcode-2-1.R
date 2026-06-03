## Simulate VAR(2)-data
library(dse1)
library(vars)
## Setting the lag-polynomial A(L) 
Apoly   <- array(c(1.0, -0.5, 0.3, 0,
                   0.2, 0.1, 0, -0.2,
                   0.7, 1, 0.5, -0.3) ,
                 c(3, 2, 2))
## Setting Covariance to identity-matrix
B <- diag(2)
## Setting constant term to 5 and 10
TRD <- c(5, 10)
## Generating the VAR(2) model 
var2  <- ARMA(A = Apoly, B = B, TREND = TRD)
## Simulating 500 observations
varsim <- simulate(var2, sampleT = 500,
                   noise = list(w = matrix(rnorm(1000),
nrow = 500, ncol = 2)), rng = list(seed = c(123456))) 
## Obtaining the generated series
vardat <- matrix(varsim$output, nrow = 500, ncol = 2)
colnames(vardat) <- c("y1", "y2")
## Plotting the series
plot.ts(vardat, main = "", xlab = "")
## Determining an appropriate lag-order
infocrit <- VARselect(vardat, lag.max = 3,
                      type = "const")
## Estimating the model
varsimest <- VAR(vardat, p = 2, type = "const",
                 season = NULL, exogen = NULL)
## Alternatively, selection according to AIC
varsimest <- VAR(vardat, type = "const",
                 lag.max = 3, ic = "SC")
## Checking the roots
roots <- roots(varsimest)
