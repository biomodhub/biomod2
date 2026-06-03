## Impulse response analysis of SVAR A-type model
args(vars:::irf.svarest)
irf.svara <- irf(svar.A, impulse = "y1",
                 response = "y2", boot = FALSE)
args(vars:::plot.varirf)
plot(irf.svara)
