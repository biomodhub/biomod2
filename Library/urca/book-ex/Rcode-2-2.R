## testing serial correlation
args(serial.test)
## Portmanteau-Test
var2c.serial <- serial.test(varsimest, lags.pt = 16,
                            type = "PT.asymptotic")
var2c.serial
plot(var2c.serial, names = "y1")
plot(var2c.serial, names = "y2")
## testing heteroscedasticity
args(arch.test)
var2c.arch <- arch.test(varsimest, lags.multi = 5,
                        multivariate.only = TRUE)
var2c.arch
## testing for normality
args(normality.test)
var2c.norm <- normality.test(varsimest,
                             multivariate.only = TRUE)
var2c.norm
## class and methods for diganostic tests
class(var2c.serial)
class(var2c.arch)
class(var2c.norm)
methods(class = "varcheck")
## Plot of objects "varcheck"
args(vars:::plot.varcheck)
plot(var2c.serial, names = "y1")
