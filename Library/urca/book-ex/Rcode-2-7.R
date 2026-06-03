## Forecast error variance decomposition
fevd.var2 <- fevd(varsimest, n.ahead = 10)
args(vars:::plot.varfevd)
plot(fevd.var2, addbars = 2)
