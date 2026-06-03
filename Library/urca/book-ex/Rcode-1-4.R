## Forecasts
arma11.pred <- predict(arma11, n.ahead = 10)
predict <- ts(c(rep(NA, length(y) - 1), y[length(y)],
                arma11.pred$pred), start = 1909,
              frequency = 1)
upper <- ts(c(rep(NA, length(y) - 1), y[length(y)],
              arma11.pred$pred + 2 * arma11.pred$se),
            start = 1909, frequency = 1)
lower <- ts(c(rep(NA, length(y) - 1), y[length(y)],
              arma11.pred$pred - 2 * arma11.pred$se),
            start = 1909, frequency = 1)
observed <- ts(c(y, rep(NA, 10)), start=1909,
               frequency = 1)
## Plot of actual and forecasted values
plot(observed, type = "l",
     ylab = "Actual and predicted values", xlab = "")
lines(predict, col = "blue", lty = 2)
lines(lower, col = "red", lty = 5)
lines(upper, col = "red", lty = 5)
abline(v = 1988, col = "gray", lty = 3)
