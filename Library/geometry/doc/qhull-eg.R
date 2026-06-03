### R code from vignette source 'qhull-eg.Rnw'

###################################################
### code chunk number 1: qhull-eg.Rnw:27-31
###################################################
library(geometry)
ps <-matrix(rnorm(30), , 2)
ch <- convhulln(ps)
head(ch)


###################################################
### code chunk number 2: qhull-eg.Rnw:43-47
###################################################
ps <-matrix(rnorm(30), , 2)
ch <- convhulln(ps, options="FA")
print(ch$area)
print(ch$vol)


###################################################
### code chunk number 3: qhull-eg.Rnw:51-52
###################################################
plot(ch)


###################################################
### code chunk number 4: qhull-eg.Rnw:56-58
###################################################
ch <- convhulln(ps, options="n")
head(ch$normals)


###################################################
### code chunk number 5: qhull-eg.Rnw:70-75
###################################################
tp <- rbox(n=200, D=2, B=4)
in_ch <- inhulln(ch, tp)
plot(tp[!in_ch,], col="gray")
points(tp[in_ch,], col="red")
plot(ch, add=TRUE)


###################################################
### code chunk number 6: qhull-eg.Rnw:86-91
###################################################
ps <- rbox(n=10, D=2)
dt <- delaunayn(ps)
head(dt)
trimesh(dt, ps)
points(ps)


###################################################
### code chunk number 7: qhull-eg.Rnw:102-106
###################################################
dt2 <- delaunayn(ps, options="Fa")
print(dt2$areas)
dt2 <- delaunayn(ps, options="Fn")
print(dt2$neighbours)


