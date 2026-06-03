library(urca)
library(uroot)
data(UKconinc)
incl <- ts(UKconinc$incl, start = c(1955,1),
           end = c(1984,4), frequency = 4)
HEGY000 <- HEGY.test(wts = incl, itsd = c(0, 0, c(0)),
                     selectlags = list(mode = c(1,4,5)))
HEGY100 <- HEGY.test(wts = incl, itsd = c(1, 0, c(0)),
                     selectlags = list(mode = c(1,4,5)))
HEGY110 <- HEGY.test(wts = incl, itsd = c(1, 1, c(0)),
                     selectlags = list(mode = c(1,4,5)))
HEGY101 <- HEGY.test(wts = incl,
                     itsd = c(1, 0, c(1, 2, 3)),
                     selectlags = list(mode = c(1,4,5)))
HEGY111 <- HEGY.test(wts = incl,
                     itsd = c(1, 1, c(1, 2, 3)),
                     selectlags = list(mode = c(1,4,5)))
