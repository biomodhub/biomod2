svec.irf <- irf(svec, response = "U",
                n.ahead = 48, boot = TRUE)
svec.irf
plot(svec.irf)


