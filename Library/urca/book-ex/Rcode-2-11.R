## FEVD analysis of SVAR B-type model
args(vars:::fevd.svarest)
fevd.svarb <- fevd(svar.B, n.ahead = 5)
class(fevd.svarb)
methods(class = "varfevd")
plot(fevd.svarb)
