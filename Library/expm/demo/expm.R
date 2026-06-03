require("expm")


# diagonalisable matrix
T <-  rbind(c(-2,  2, 0),
            c(-3, -2, 2),
            c( 2,  1,-2))

expm(T)


# numerically singular matrix
T <- rbind(c(-2, 2, 0),
           c( 0,-2, 2),
           c( 0, 0,-2))

expm(T)

## logm(), the inverse of expm() :
T - logm(expm(T)) # small ( ~ 1e-13 )
stopifnot(all.equal(T, logm(expm(T))))

# solve shows T is numerically singular
try(solve(eigen(T)$vectors))

# singular matrix
T <- rbind(c(0, 2, 1),
           c(0, 0, 2),
           c(0, 0, 0))

expm(T)
stopifnot(all.equal(logm(expm(T)), T))
## and show how close it is
all.equal(logm(expm(T)), T,  tolerance=0)# 2.39e-15 {64b ubuntu 12-04}
