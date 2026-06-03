### R code from vignette source 'prs.Rnw'

###################################################
### code chunk number 1: prs.Rnw:48-54
###################################################
library(gstat)
cluster = read.table(system.file("external/cluster.txt", package="gstat"),
	header = TRUE)
summary(cluster)
library(sp)
coordinates(cluster) = ~X+Y


###################################################
### code chunk number 2: prs.Rnw:59-61
###################################################
bnd = c(0,2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5)
variogram(Primary~1, cluster, boundaries = bnd)


###################################################
### code chunk number 3: prs.Rnw:66-67
###################################################
variogram(Primary~1, cluster, boundaries=bnd, PR = TRUE)


###################################################
### code chunk number 4: prs.Rnw:71-75
###################################################
pl1 = plot(variogram(Primary~1, cluster, boundaries=bnd, PR = FALSE))
pl2 = plot(variogram(Primary~1, cluster, boundaries=bnd, PR = TRUE))
print(pl1, split = c(1,1,2,1), more = TRUE)
print(pl2, split = c(2,1,2,1), more = FALSE)


###################################################
### code chunk number 5: prs.Rnw:86-96
###################################################
z = cluster$Primary
d = spDists(cluster)
zd = outer(z, z, "-")
zs = outer(z, z, "+")
pr = (2 * zd / zs )^2
prv = as.vector(pr)
dv = as.vector(d)
mean(prv[dv > 0 & dv < 2.5])/2
mean(prv[dv > 2.5 & dv < 7.5])/2
mean(prv[dv > 7.5 & dv < 12.5])/2


