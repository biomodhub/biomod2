## ----echo=FALSE, message=FALSE------------------------------------------------
knitr::opts_chunk$set(global.par=TRUE, collapse=TRUE, comment="#>", fig.width=5, fig.height=5, fig.align="center", dpi=96)
options(tibble.print_min=3L, tibble.print_max=3L)

## ----message=FALSE------------------------------------------------------------
library(ks)
mus <- rbind(c(-2,2), c(0,0), c(2,-2))
Sigmas <- rbind(diag(2), matrix(c(0.8, -0.72, -0.72, 0.8), nrow=2), diag(2))
cwt <- 3/11
props <- c((1-cwt)/2, cwt, (1-cwt)/2)
plotmixt(mus=mus, Sigmas=Sigmas, props=props, display="filled.contour", xlim=c(-4,4), ylim=c(-4,4))
set.seed(88192)
samp <- 2000
x <- rmvnorm.mixt(n=samp, mus=mus, Sigmas=Sigmas, props=props)
colnames(x) <- c("x","y")
plot(x, pch=16, cex=0.5, xlim=c(-4,4), ylim=c(-4,4))

## -----------------------------------------------------------------------------
Hpi1 <- Hpi(x=x)
Hpi2 <- Hpi.diag(x=x)

## -----------------------------------------------------------------------------
fhat.pi1 <- kde(x=x)         ## same as Hpi(x=x, H=Hpi1)
fhat.pi2 <- kde(x=x, H=Hpi2)

## -----------------------------------------------------------------------------
plot(fhat.pi1, main="Plug-in", display="filled.contour", xlim=c(-4,4), ylim=c(-4,4))
plot(fhat.pi2, main="Plug-in diagonal", display="filled.contour", xlim=c(-4,4), ylim=c(-4,4)) 

## -----------------------------------------------------------------------------
Hscv1 <- Hscv(x=x)
Hscv2 <- Hscv.diag(x=x)
fhat.cv1 <- kde(x=x, H=Hscv1)
fhat.cv2 <- kde(x=x, H=Hscv2)
plot(fhat.cv1, main="SCV", display="filled.contour", xlim=c(-4,4), ylim=c(-4,4))
plot(fhat.cv2, main="SCV diagonal", display="filled.contour", xlim=c(-4,4), ylim=c(-4,4))

