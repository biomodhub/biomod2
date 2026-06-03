library(urca)
data(npext)
npext
y <- ts(na.omit(npext$unemploy), start=1890, end=1988,
        frequency=1)
op <- par(no.readonly=TRUE)
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow=TRUE))
plot(y, ylab="unemployment rate (logarithm)")
acf(y, main='Autocorrelations', ylab='', ylim=c(-1, 1))
pacf(y, main='Partial Autocorrelations', ylab='',
     ylim=c(-1, 1))
par(op)
## tentative ARMA(2,0)
arma20 <- arima(y, order=c(2, 0, 0))
ll20 <- logLik(arma20)
aic20 <- arma20$aic
res20 <- residuals(arma20)
Box.test(res20, lag = 20, type =  "Ljung-Box")
shapiro.test(res20)
## alternative specifications
## ARMA(3,0)
arma30 <- arima(y, order=c(3, 0, 0))
ll30 <- logLik(arma30)
aic30 <- arma30$aic
lrtest <- as.numeric(2*(ll30 - ll20))
chi.pval <- pchisq(lrtest, df = 1, lower.tail = FALSE)
## ARMA(1,1)
arma11 <- arima(y, order = c(1, 0, 1))
ll11 <- logLik(arma11)
aic11 <- arma11$aic
tsdiag(arma11)
res11 <- residuals(arma11)
Box.test(res11, lag = 20, type =  "Ljung-Box")
shapiro.test(res11)
tsdiag(arma11)
## Using auto.arima()
library(forecast)
auto.arima(y, max.p = 3, max.q = 3, start.p = 1,
           start.q = 1, ic = "aic")
