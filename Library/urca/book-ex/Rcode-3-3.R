library(fracdiff)
set.seed(123456)
# ARFIMA(0.0,0.3,0.0)
y <- fracdiff.sim(n=1000, ar=0.0, ma=0.0, d=0.3)
# Get the data series, demean this if necessary
y.dm <- y$series 
max.y <- max(cumsum(y.dm))
min.y <- min(cumsum(y.dm))
sd.y <- sd(y$series)
RS <- (max.y - min.y)/sd.y 
H <- log(RS)/log(1000)
d <- H - 0.5
