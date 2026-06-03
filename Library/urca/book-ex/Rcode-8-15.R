LR[3, 3] <- 0
LR
svec.oi <- SVEC(vecm, LR = LR, SR = SR, r = 1,
                lrtest = TRUE, boot = FALSE)
svec.oi <- update(svec, LR = LR, lrtest = TRUE,
                  boot = FALSE)
svec.oi$LRover
