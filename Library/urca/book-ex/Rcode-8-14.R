vecm <- ca.jo(Canada[, c("prod", "e", "U", "rw")],
              type = "trace", ecdet = "trend",
              K = 3, spec = "transitory")
SR <- matrix(NA, nrow = 4, ncol = 4)
SR[4, 2] <- 0
SR
LR <- matrix(NA, nrow = 4, ncol = 4)
LR[1, 2:4] <- 0
LR[2:4, 4] <- 0
LR
svec <- SVEC(vecm, LR = LR, SR = SR, r = 1,
             lrtest = FALSE,
             boot = TRUE, runs = 100)
svec
svec$SR / svec$SRse
svec$LR / svec$LRse
