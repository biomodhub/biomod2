## ----load_forecast, echo=FALSE, message=FALSE---------------------------------
library('forecast')

## ----load_expsmooth, echo=FALSE, message=FALSE, eval=FALSE--------------------
# library('expsmooth')

## ----expsmooth_datsets, echo=FALSE, message=FALSE-----------------------------
bonds <-
structure(c(5.83, 6.06, 6.58, 7.09, 7.31, 7.23, 7.43, 7.37, 7.6,
7.89, 8.12, 7.96, 7.93, 7.61, 7.33, 7.18, 6.74, 6.27, 6.38, 6.6,
6.3, 6.13, 6.02, 5.79, 5.73, 5.89, 6.37, 6.62, 6.85, 7.03, 6.99,
6.75, 6.95, 6.64, 6.3, 6.4, 6.69, 6.52, 6.8, 7.01, 6.82, 6.6,
6.32, 6.4, 6.11, 5.82, 5.87, 5.89, 5.63, 5.65, 5.73, 5.72, 5.73,
5.58, 5.53, 5.41, 4.87, 4.58, 4.89, 4.69, 4.78, 4.99, 5.23, 5.18,
5.54, 5.9, 5.8, 5.94, 5.91, 6.1, 6.03, 6.26, 6.66, 6.52, 6.26,
6, 6.42, 6.1, 6.04, 5.83, 5.8, 5.74, 5.72, 5.23, 5.14, 5.1, 4.89,
5.13, 5.37, 5.26, 5.23, 4.97, 4.76, 4.55, 4.61, 5.07, 5, 4.9,
5.28, 5.21, 5.15, 4.9, 4.62, 4.24, 3.88, 3.91, 4.04, 4.03, 4.02,
3.9, 3.79, 3.94, 3.56, 3.32, 3.93, 4.44, 4.29, 4.27, 4.29, 4.26,
4.13, 4.06, 3.81, 4.32, 4.7), .Tsp = c(1994, 2004.33333333333,
12), class = "ts")
usnetelec <-
structure(c(296.1, 334.1, 375.3, 403.8, 447, 476.3, 550.3, 603.9,
634.6, 648.5, 713.4, 759.2, 797.1, 857.9, 920, 987.2, 1058.4,
1147.5, 1217.8, 1332.8, 1445.5, 1535.1, 1615.9, 1753, 1864.1,
1870.3, 1920.8, 2040.9, 2127.4, 2209.4, 2250.7, 2289.6, 2298,
2244.4, 2313.4, 2419.5, 2473, 2490.5, 2575.3, 2707.4, 2967.3,
3038, 3073.8, 3083.9, 3197.2, 3247.5, 3353.5, 3444.2, 3492.2,
3620.3, 3694.8, 3802.1, 3736.6, 3858.5, 3848), .Tsp = c(1949,
2003, 1), class = "ts")
ukcars <-
structure(c(330.371, 371.051, 270.67, 343.88, 358.491, 362.822,
261.281, 240.355, 325.382, 316.7, 171.153, 257.217, 298.127,
251.464, 181.555, 192.598, 245.652, 245.526, 225.261, 238.211,
257.385, 228.461, 175.371, 226.462, 266.15, 287.251, 225.883,
265.313, 272.759, 234.134, 196.462, 205.551, 291.283, 284.422,
221.571, 250.697, 253.757, 267.016, 220.388, 277.801, 283.233,
302.072, 259.72, 297.658, 306.129, 322.106, 256.723, 341.877,
356.004, 361.54, 270.433, 311.105, 326.688, 327.059, 274.257,
367.606, 346.163, 348.211, 250.008, 292.518, 343.318, 343.429,
275.386, 329.747, 364.521, 378.448, 300.798, 331.757, 362.536,
389.133, 323.322, 391.832, 421.646, 416.823, 311.713, 381.902,
422.982, 427.722, 376.85, 458.58, 436.225, 441.487, 369.566,
450.723, 462.442, 468.232, 403.636, 413.948, 460.496, 448.932,
407.787, 469.408, 494.311, 433.24, 335.106, 378.795, 387.1, 372.395,
335.79, 397.08, 449.755, 402.252, 391.847, 385.89, 424.325, 433.28,
391.213, 408.74, 445.458, 428.202, 379.048, 394.042, 432.796), .Tsp = c(1977,
2005, 4), class = "ts")
visitors <-
structure(c(75.7, 75.4, 83.1, 82.9, 77.3, 105.7, 121.9, 150,
98, 118, 129.5, 110.6, 91.7, 94.8, 109.5, 105.1, 95, 130.3, 156.7,
190.1, 139.7, 147.8, 145.2, 132.7, 120.7, 116.5, 142, 140.4,
128, 165.7, 183.1, 222.8, 161.3, 180.4, 185.2, 160.5, 157.1,
163.8, 203.3, 196.9, 179.6, 207.3, 208, 245.8, 168.9, 191.1,
180, 160.1, 136.6, 142.7, 175.4, 161.4, 149.9, 174.1, 192.7,
247.4, 176.2, 192.8, 189.1, 181.1, 149.9, 157.3, 185.3, 178.2,
162.7, 190.6, 198.6, 253.1, 177.4, 190.6, 189.2, 168, 161.4,
172.2, 208.3, 199.3, 197.4, 216, 223.9, 266.8, 196.1, 238.2,
217.8, 203.8, 175.2, 176.9, 219.3, 199.1, 190, 229.3, 255, 302.4,
242.8, 245.5, 257.9, 226.3, 213.4, 204.6, 244.6, 239.9, 224,
267.2, 285.9, 344, 250.5, 304.3, 307.4, 255.1, 214.9, 230.9,
282.5, 265.4, 254, 301.6, 311, 384, 303.8, 319.1, 313.5, 294.2,
244.8, 261.4, 329.7, 304.9, 268.6, 320.7, 342.9, 422.3, 317.2,
392.7, 365.6, 333.2, 261.5, 306.9, 358.2, 329.2, 309.2, 350.4,
375.6, 465.2, 342.9, 408, 390.9, 325.9, 289.1, 308.2, 397.4,
330.4, 330.9, 366.5, 379.5, 448.3, 346.2, 353.6, 338.6, 341.1,
283.4, 304.2, 372.3, 323.7, 323.9, 354.8, 367.9, 457.6, 351,
398.6, 389, 334.1, 298.1, 317.1, 388.5, 355.6, 353.1, 397, 416.7,
460.8, 360.8, 434.6, 411.9, 405.6, 319.3, 347.9, 429, 372.9,
403, 426.5, 459.9, 559.9, 416.6, 429.2, 428.7, 405.4, 330.2,
370, 446.9, 384.6, 366.3, 378.5, 376.2, 523.2, 379.3, 437.2,
446.5, 360.3, 329.9, 339.4, 418.2, 371.9, 358.6, 428.9, 437,
534, 396.6, 427.5, 392.5, 321.5, 260.9, 308.3, 415.5, 362.2,
385.6, 435.3, 473.3, 566.6, 420.2, 454.8, 432.3, 402.8, 341.3,
367.3, 472, 405.8, 395.6, 449.9, 479.9, 593.1, 462.4, 501.6,
504.7, 409.5), .Tsp = c(1985.33333333333, 2005.25, 12), class = "ts")

## ----etsexamples, fig.height=7, fig.width=9, echo=FALSE, fig.cap="Four time series showing point forecasts and 80\\% \\& 95\\% prediction intervals obtained using exponential smoothing state space models."----
par(mfrow = c(2,2))
mod1 <- ets(bonds)
mod2 <- ets(usnetelec)
mod3 <- ets(ukcars)
mod4 <- ets(visitors)

plot(forecast(mod1), main="(a) US 10-year bonds yield", xlab="Year", ylab="Percentage per annum")
plot(forecast(mod2), main="(b) US net electricity generation", xlab="Year", ylab="Billion kwh")
plot(forecast(mod3), main="(c) UK passenger motor vehicle production", xlab="Year", ylab="Thousands of cars")
plot(forecast(mod4), main="(d) Overseas visitors to Australia", xlab="Year", ylab="Thousands of people")

## ----etsnames, echo=FALSE-----------------------------------------------------
etsnames <- c(mod1$method, mod2$method, mod3$method, mod4$method)
etsnames <- gsub("Ad","A\\\\damped",etsnames)

## ----ets-usnetelec, echo=TRUE-------------------------------------------------
etsfit <- ets(usnetelec)

## ----ets-usnetelec-print,echo=TRUE--------------------------------------------
etsfit

## ----ets-usnetelec-accuracy,eval=TRUE,echo=TRUE-------------------------------
accuracy(etsfit)

## ----ets-usnetelec-fcast, fig.height=5, fig.width=8, message=FALSE, warning=FALSE, include=FALSE, output=FALSE----
fcast <- forecast(etsfit)
plot(fcast)

## ----ets-usnetelec-fcast-print,eval=TRUE,echo=TRUE----------------------------
fcast

## ----ets-usnetelec-newdata,eval=FALSE,echo=TRUE-------------------------------
# fit <- ets(usnetelec[1:45])
# test <- ets(usnetelec[46:55], model = fit)
# accuracy(test)

## ----ets-usnetelec-fcast-accuracy,eval=FALSE,echo=TRUE------------------------
# accuracy(forecast(fit,10), usnetelec[46:55])

## ----arimaexamples, fig.height=7, fig.width=9, echo=FALSE, fig.cap="Four time series showing point forecasts and 80\\% \\& 95\\% prediction intervals obtained using ARIMA models."----
mod1 <- auto.arima(bonds, seasonal=FALSE, approximation=FALSE)
mod2 <- auto.arima(usnetelec)
mod3 <- auto.arima(ukcars)
mod4 <- auto.arima(visitors)
par(mfrow = c(2,2))
plot(forecast(mod1), main="(a) US 10-year bonds yield", xlab="Year", ylab="Percentage per annum")
plot(forecast(mod2), main="(b) US net electricity generation", xlab="Year", ylab="Billion kwh")
plot(forecast(mod3), main="(c) UK passenger motor vehicle production", xlab="Year", ylab="Thousands of cars")
plot(forecast(mod4), main="(d) Overseas visitors to Australia", xlab="Year", ylab="Thousands of people")

## ----arima-auto-fcast,eval=TRUE,echo=TRUE,fig.show="hide"---------------------
arimafit <- auto.arima(usnetelec)
fcast <- forecast(arimafit)
plot(fcast)

## ----arimanames, echo=FALSE---------------------------------------------------
# Convert character strings to latex
arimanames <- c(as.character(mod1),
  as.character(mod2),
  as.character(mod3),
  as.character(mod4))
arimanames <-
    gsub("\\[([0-9]*)\\]", "$_{\\1}$", arimanames)

## ----arimafcastsummary, echo=TRUE, message=FALSE, warning=FALSE, as.is=TRUE----
summary(fcast)

