### R code from vignette source 'sto.Rnw'

###################################################
### code chunk number 1: sto.Rnw:89-90 (eval = FALSE)
###################################################
## over(x, geometry(y))


###################################################
### code chunk number 2: sto.Rnw:97-98 (eval = FALSE)
###################################################
## over(x, y)


###################################################
### code chunk number 3: sto.Rnw:106-107 (eval = FALSE)
###################################################
## over(x, y, returnList = TRUE)


###################################################
### code chunk number 4: sto.Rnw:118-121
###################################################
library(sp)
library(spacetime)
showMethods(over)


###################################################
### code chunk number 5: sto.Rnw:133-134 (eval = FALSE)
###################################################
## over(x, y)


###################################################
### code chunk number 6: sto.Rnw:182-183 (eval = FALSE)
###################################################
## aggregate(x, by, FUN, ..., simplify = TRUE)


###################################################
### code chunk number 7: sto.Rnw:199-203
###################################################
data(air)
rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
class(rural)
class(DE_NUTS1)


###################################################
### code chunk number 8: sto.Rnw:216-218
###################################################
x = as(rural[,"2008"], "xts")
apply(x, 1, mean, na.rm=TRUE)[1:5]


###################################################
### code chunk number 9: sto.Rnw:224-229
###################################################
dim(rural[,"2008"])
x = aggregate(rural[,"2008"], DE_NUTS1, mean, na.rm=TRUE)
dim(x)
summary(x)
stplot(x, mode = "tp")


###################################################
### code chunk number 10: sto.Rnw:233-234 (eval = FALSE)
###################################################
## stplot(x, mode = "tp", par.strip.text = list(cex=.5))


###################################################
### code chunk number 11: sto.Rnw:238-239
###################################################
print(stplot(x, mode = "tp", par.strip.text = list(cex=.5)))


###################################################
### code chunk number 12: sto.Rnw:249-252
###################################################
x = aggregate(rural[,"2008"], DE, mean, na.rm=TRUE)
class(x)
plot(x[,"PM10"])


###################################################
### code chunk number 13: sto.Rnw:257-258
###################################################
plot(x[,"PM10"])


###################################################
### code chunk number 14: sto.Rnw:269-271
###################################################
x = as(rural[,"2008"], "xts")
apply(x, 2, mean, na.rm=TRUE)[1:5]


###################################################
### code chunk number 15: sto.Rnw:278-281
###################################################
sel = which(!apply(as(rural[,"2008"], "xts"), 2, function(x) all(is.na(x))))
x = aggregate(rural[sel,"2008"], "month", mean, na.rm=TRUE)
stplot(x, mode = "tp")


###################################################
### code chunk number 16: sto.Rnw:286-287
###################################################
print(stplot(x, mode = "tp", par.strip.text = list(cex=.5)))


###################################################
### code chunk number 17: sto.Rnw:303-306
###################################################
library(zoo)
x = aggregate(rural[sel,"2005::2011"], as.yearqtr, median, na.rm=TRUE)
stplot(x, mode = "tp")


###################################################
### code chunk number 18: sto.Rnw:311-312
###################################################
as.year <- function(x) as.numeric(floor(as.yearmon(x)))


###################################################
### code chunk number 19: sto.Rnw:318-319
###################################################
print(stplot(x, mode = "tp", par.strip.text = list(cex=.5)))


###################################################
### code chunk number 20: sto.Rnw:329-331
###################################################
DE.years = STF(DE, as.Date(c("2008-01-01", "2009-01-01")))
aggregate(rural[,"2008::2009"], DE.years, mean, na.rm=TRUE)


