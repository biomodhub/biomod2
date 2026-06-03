Canada <- Canada[, c("prod", "e", "U", "rw")]
p1ct <- VAR(Canada, p = 1, type = "both")
p2ct <- VAR(Canada, p = 2, type = "both")
p3ct <- VAR(Canada, p = 3, type = "both")
## Serial
serial.test(p3ct, lags.pt = 16,
            type = "PT.asymptotic")
serial.test(p2ct, lags.pt = 16,
            type = "PT.asymptotic")
serial.test(p1ct, lags.pt = 16,
            type = "PT.asymptotic")
serial.test(p3ct, lags.pt = 16,
            type = "PT.adjusted")
serial.test(p2ct, lags.pt = 16,
            type = "PT.adjusted")
serial.test(p1ct, lags.pt = 16,
            type = "PT.adjusted")
## JB
normality.test(p3ct)
normality.test(p2ct)
normality.test(p1ct)
## ARCH
arch.test(p3ct, lags.multi = 5)
arch.test(p2ct, lags.multi = 5)
arch.test(p1ct, lags.multi = 5)
## Stability (Recursive CUSUM)
plot(stability(p3ct), nc = 2)
plot(stability(p2ct), nc = 2)
plot(stability(p1ct), nc = 2)
