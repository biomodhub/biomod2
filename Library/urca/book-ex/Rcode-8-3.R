A1 <- matrix(c(1,0,0,0,0, 0,0,1,0,0,
               0,0,0,1,0, 0,0,0,0,1),
             nrow=5, ncol=4)
A2 <- matrix(c(1,0,0,0,0, 0,1,0,0,0,
               0,0,1,0,0, 0,0,0,1,0),
             nrow=5, ncol=4)
H41 <- summary(alrtest(z = H1, A = A1, r = 2))
H42 <- summary(alrtest(z = H1, A = A2, r = 2))
