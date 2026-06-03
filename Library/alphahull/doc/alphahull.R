### R code from vignette source 'alphahull.rnw'

###################################################
### code chunk number 1: alphahull.rnw:77-78
###################################################
library(alphahull)


###################################################
### code chunk number 2: alphahull.rnw:118-125
###################################################
# Uniform sample of size n=300 in the disc B(c,0.5)\B(c,0.25) 
# with c=(0.5,0.5).
n<-150
theta<-runif(n,0,2*pi)
r<-sqrt(runif(n,0.25^2,0.5^2))
x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
alpha<-0.1


###################################################
### code chunk number 3: alphahull.rnw:131-134
###################################################
par(mfrow=c(1,1))
print(plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F))
print(plot(ahull(x,alpha=alpha),col=c(6,1,1),xlab="",ylab="",add=T))


###################################################
### code chunk number 4: alphahull.rnw:149-152
###################################################
par(mfrow=c(1,1))
print(plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F))
print(plot(ashape(x,alpha=alpha),col=c(4,1),xlab="",ylab="",add=T))


###################################################
### code chunk number 5: alphahull.rnw:161-167
###################################################
# Uniform sample of size n=300 in the disc B(c,0.5)\B(c,0.25) 
# with c=(0.5,0.5).
n<-300
theta<-runif(n,0,2*pi)
r<-sqrt(runif(n,0.25^2,0.5^2))
x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))


###################################################
### code chunk number 6: alphahull.rnw:172-180
###################################################
par(mfrow=c(1,3))
alpha1=0.02
alpha2=0.25
alpha3=1
plot(ashape(x,alpha=alpha1),col=c(4,1),xlab="",ylab="",main=expression(paste(alpha," = 0.02 ")))
plot(ashape(x,alpha=alpha2),col=c(4,1),xlab="",ylab="",main=expression(paste(alpha," = 0.25")))
plot(ashape(x,alpha=alpha3),col=c(4,1),xlab="",ylab="",main=expression(paste(alpha," = 1")))
par(mfrow=c(1,1))


###################################################
### code chunk number 7: alphahull.rnw:195-199
###################################################
#x<-matrix(runif(20),nc=2)
x1<-c(0.5915,0.6230,0.9689,0.8248,0.9392,0.8156,0.2050,0.9757,0.0957,0.4139)
y1<-c(0.472,0.619,0.304,0.197,0.716,0.575,0.507,0.574,0.996,0.893)
x<-cbind(x1,y1)


###################################################
### code chunk number 8: alphahull.rnw:205-208
###################################################
par(mfrow=c(1,1))
print(plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F))
print(plot(delvor(x),col=1:3,xlab="",ylab="",add=T))


###################################################
### code chunk number 9: alphahull.rnw:223-227
###################################################
x <- c(0.905, 0.606, 0.458, 0.988, 0.744)
y <- c(0.763, 0.937, 0.095, 0.259, 0.731) 
dv <- delvor(x,y)
dv


###################################################
### code chunk number 10: alphahull.rnw:233-234
###################################################
plot(dv, main = "Delaunay triangulation and Voronoi diagram",col = 1:3, xlab = "x-coordinate", ylab = "y-coordinate",xlim = c(-0.5, 1.5), ylim = c(-0.5, 1.5), number = TRUE)


###################################################
### code chunk number 11: alphahull.rnw:241-243
###################################################
par(mfrow=c(1,1))
print(plot(dv, main = "Delaunay triangulation and Voronoi diagram",col = 1:3, xlab = "x-coordinate", ylab = "y-coordinate",xlim = c(-0.5, 1.5), ylim = c(-0.5, 1.5), number = TRUE))


###################################################
### code chunk number 12: alphahull.rnw:275-278
###################################################
x<-matrix(runif(100),ncol=2)
alpha <- 0.2
alphashape <- ashape(x, alpha = alpha)


###################################################
### code chunk number 13: alphahull.rnw:281-282
###################################################
names(alphashape)


###################################################
### code chunk number 14: alphahull.rnw:287-288
###################################################
plot(alphashape, col = c(4, 1), xlab = "x-coordinate", ylab = "y-coordinate", main = expression(paste(alpha, "-shape")))


###################################################
### code chunk number 15: alphahull.rnw:294-296
###################################################
par(mfrow=c(1,1))
print(plot(alphashape, col = c(4, 1), xlab = "x-coordinate", ylab = "y-coordinate", main = expression(paste(alpha, "-shape"))))


###################################################
### code chunk number 16: alphahull.rnw:304-305
###################################################
plot(alphashape, wlines = "del", col = c(4, 1, 2), xlab = "x-coordinate", ylab = "y-coordinate")


###################################################
### code chunk number 17: alphahull.rnw:311-313
###################################################
par(mfrow=c(1,1))
print(plot(alphashape, wlines = "del", col = c(4, 1, 2), xlab = "x-coordinate",ylab = "y-coordinate",main = expression(paste(alpha, "-shape and Delaunay triangulation"))))


###################################################
### code chunk number 18: alphahull.rnw:367-393
###################################################
par(mfrow=c(1,1))
plot(0,type="n",axes=FALSE,xlim=c(0,0.5),ylim=c(0,0.5),xlab="",ylab="")
r<-0.5
t<-0
segments(0,0,r*cos(t),r*sin(t),col=4,lty=2)
points(r*cos(t),r*sin(t),pch=19,col=4)

t<-pi/6
arrows(0,0,0.3*cos(t),0.3*sin(t))
v<-c(0.3*cos(t),0.3*sin(t))
v<-v/sqrt(sum(v^2))
arc(c(0,0),0.5,v,pi/6,col=4,lwd=2)
t<-pi/3
segments(0,0,r*cos(t),r*sin(t),col=4,lty=2)
points(r*cos(t),r*sin(t),pch=19,col=4)

t<-pi/12
v<-c(0.3*cos(t),0.3*sin(t))
v<-v/sqrt(sum(v^2))
arc(c(0,0),0.3,v,pi/12,col=1)


text(0.15,0.11,expression(italic(v)),cex=1.5)
text(0.31,0.07,expression(italic(theta)),cex=1.5)
text(-0.015,0,expression(italic(c)),cex=1.5)
text(0.15,0.3,expression(italic(r)),cex=1.5)


###################################################
### code chunk number 19: alphahull.rnw:402-408
###################################################
n <- 200
theta<-runif(n,0,2*pi)
r<-sqrt(runif(n,0.25^2,0.5^2))
x<-cbind(0.5+r*cos(theta),0.5+r*sin(theta))
alpha <- 0.15
alphahull <- ahull(x, alpha = alpha)


###################################################
### code chunk number 20: alphahull.rnw:411-415
###################################################
names(alphahull)
alphahull$complement[1:5,1:3]
alphahull$arcs[1:5,]
alphahull$length


###################################################
### code chunk number 21: alphahull.rnw:421-422
###################################################
plot(alphahull, col = c(6, rep(1, 5)), xlab = "x-coordinate",ylab = "y-coordinate", main = expression(paste(alpha, "-hull")))


###################################################
### code chunk number 22: alphahull.rnw:428-430
###################################################
par(mfrow=c(1,1))
print(plot(alphahull,col=c(6,rep(1,5)), xlab = "x-coordinate",ylab = "y-coordinate", main = expression(paste(alpha, "-hull"))))


###################################################
### code chunk number 23: alphahull.rnw:438-443
###################################################
plot(alphahull, col = c(6, rep(1, 5)), xlab = "x-coordinate", ylab = "y-coordinate", main = expression(paste(alpha, "-hull")))
warcs<- which(alphahull$arcs[,3]>0)
for (i in warcs) {
arc(alphahull$arcs[i, 1:2], alphahull$arcs[i,3], c(0,1), pi, col = "gray", lty = 2)
}


###################################################
### code chunk number 24: alphahull.rnw:448-449
###################################################
plot(alphahull, do.shape = TRUE, col = c(6, 4, rep(1, 4)), xlab = "x-coordinate", ylab = "y-coordinate", main = expression(paste(alpha, "-hull and ",alpha, "-shape")))


###################################################
### code chunk number 25: alphahull.rnw:455-463
###################################################
par(mfrow=c(1,2))
print(plot(alphahull, col = c(6, rep(1, 5)), xlab = "x-coordinate", ylab = "y-coordinate", main = expression(paste(alpha, "-hull"))))
warcs<- which(alphahull$arcs[,3]>0)
for (i in warcs) {
arc(alphahull$arcs[i, 1:2], alphahull$arcs[i,3], c(0,1), 2 * pi, col = "gray", lty = 2)
}
print(plot(alphahull, do.shape = TRUE, col = c(6, 4, rep(1, 4)),xlab = "x-coordinate", ylab = "y-coordinate",main = expression(paste(alpha, "-hull and ",alpha, "-shape")))
)


