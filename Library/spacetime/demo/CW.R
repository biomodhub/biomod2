# chapter 3, page 90:
# example of Yuler-Walker equations -- find the
# coefficients of an AR(3) model, using correlations:
library(gstat)
data(wind)
# take station in column 5:
x = wind[[5]]
plot(x, type='l')
acf4 = acf(x, plot = FALSE)$acf[1:4]
acf4
C = matrix(NA,3,3)
diag(C) = 1
C[,1] = acf4[1:3]
C[1,] = acf4[1:3]
C[2,3]=C[3,2]=acf4[2]
C
# coeficients:
solve(C,acf4[-1])
# check:
arima(x, c(3,0,0))

# Numerical experiment, ch3, page 105
Y = c(-.0831,.7187,.2182,.3068,-.9568,-.4170)
plot(Y, type='l')
t = 1:6
phi0 = rep(1/sqrt(6),6)
phi1 = sqrt(2/6)*cos(2*pi*t/6)
phi2 = sqrt(2/6)*sin(2*pi*t/6)
phi3 = sqrt(2/6)*cos(4*pi*t/6)
phi4 = sqrt(2/6)*sin(4*pi*t/6)
phi5 = sqrt(1/6)*cos(pi*t)
Phi = cbind(phi0,phi1,phi2,phi3,phi4,phi5)
alpha = t(Phi) %*% Y
A = c(alpha[1], sqrt(sum(alpha[2:3]^2)), sqrt(sum(alpha[4:5]^2)), alpha[6])
# answer: A1, A2, A3, A4:
round(A, 2)
# create reduced series, taking out A0 and A2:
alphaR = alpha; alphaR[c(1,4,5)] = 0
YR = Phi %*% alphaR
round(YR, 4)
plot(Y, type = 'l', col = 'blue')
lines(YR, type = 'l', col = 'green', lty = 2)

#### CHAPTER 4: the difference between "kriging" Z or Y. -- p. 137.
require(gstat)

# kriging on a grid (no points coincide with grid cell centers):
loadMeuse()
m1 <- vgm(.59, "Sph", 874, .04)
z <- krige(log(zinc)~1, meuse, meuse.grid, model = m1)
m2 <- vgm(.59, "Sph", 874, add.to = vgm(.04, "Err", 0))
y <- krige(log(zinc)~1, meuse, meuse.grid, model = m2)
# the difference in terms of predictions:
summary(z$var1.pred - y$var1.pred) # no difference
# the difference in terms of prediction variances:
summary(z$var1.var - y$var1.var) # difference 0.04

# kriging on data points:
z <- krige(log(zinc)~1, meuse, meuse[1:10,], model = m1)
y <- krige(log(zinc)~1, meuse, meuse[1:10,], model = m2)
# z equals the data, y does not:
summary(log(meuse$zinc)[1:10] - z[[1]])
summary(log(meuse$zinc)[1:10] - y[[1]])
# z has zero variance, y does not:
summary(z$var1.var)
summary(y$var1.var)
# this causes a difference in terms of predictions:
summary(z$var1.pred - y$var1.pred) # small, but varying, as
# and difference in terms of prediction variances:
summary(z$var1.var - y$var1.var)
