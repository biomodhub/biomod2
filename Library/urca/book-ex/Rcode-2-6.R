## Impulse response analysis
irf.y1 <- irf(varsimest, impulse = "y1",
              response = "y2", n.ahead = 10,
              ortho = FALSE, cumulative = FALSE,
              boot = FALSE, seed = 12345)
args(vars:::plot.varirf)
plot(irf.y1)
irf.y2 <- irf(varsimest, impulse = "y2",
              response = "y1", n.ahead = 10,
              ortho = TRUE, cumulative = TRUE,
              boot = FALSE, seed = 12345)
plot(irf.y2)
