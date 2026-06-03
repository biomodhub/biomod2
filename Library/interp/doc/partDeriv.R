### R code from vignette source 'partDeriv.Rnw'

###################################################
### code chunk number 1: init
###################################################
set.seed(42)
options(width=80)
options(continue=" ")
options(SweaveHooks=list(fig=function()
    par(mar=c(5.1, 4.1, 1.1, 2.1))))
library(interp)
library(Deriv)
library(Ryacas)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)


###################################################
### code chunk number 2: partDeriv.Rnw:415-416
###################################################
ng <- 11


###################################################
### code chunk number 3: partDeriv.Rnw:421-422
###################################################
knl <- "gaussian"


###################################################
### code chunk number 4: partDeriv.Rnw:429-431
###################################################
bwg <- 0.33  
bwl <- 0.11  


###################################################
### code chunk number 5: partDeriv.Rnw:444-445
###################################################
dg=3


###################################################
### code chunk number 6: partDeriv.Rnw:448-449
###################################################
f <- function(x,y) (x-0.5)*(x-0.2)*(y-0.6)*y*(x-1)


###################################################
### code chunk number 7: helperR2Yacas
###################################################
# helper functions for translation between R and Yacas
fn_y  <- function(f){
    b <- toString(as.expression(body(f)))
    b <- stringr::str_replace_all(b,"cos","Cos")
    b <- stringr::str_replace_all(b,"sin","Sin")
    b <- stringr::str_replace_all(b,"exp","Exp")
    b <- stringr::str_replace_all(b,"log","Log")
    b <- stringr::str_replace_all(b,"sqrt","Sqrt")
    b
}


###################################################
### code chunk number 8: helperYacas2R
###################################################
ys_fn  <- function(f){
    f <- stringr::str_replace_all(f,"Cos","cos")
    f <- stringr::str_replace_all(f,"Sin","sin")
    f <- stringr::str_replace_all(f,"Exp","exp")
    f <- stringr::str_replace_all(f,"Log","log")
    f <- stringr::str_replace_all(f,"Sqrt","sqrt")
    f
}


###################################################
### code chunk number 9: helperDerivs
###################################################
derivs <- function(f,dg){
    ret<-list(f=f,
              f_str=ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),""),")"))))

    if(dg>0){

        ret$fx <- function(x,y){
            myfx <- Deriv(f,"x");
            tmp <- myfx(x,y);
            if(length(tmp)==1)
                return(rep(tmp,length(x)))
            else
                return(tmp)
        }
        ret$fx_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)"),")")))


        ret$fy <- function(x,y){
            myfy <- Deriv(f,"y");
            tmp <- myfy(x,y);
            if(length(tmp)==1)
                return(rep(tmp,length(x)))
            else
                return(tmp)
        }
        ret$fy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(y)"),")")))


        if(dg>1){
            ret$fxy <- function(x,y){
                myfxy <- Deriv(Deriv(f,"y"),"x");
                tmp <- myfxy(x,y);
                if(length(tmp)==1)
                    return(rep(tmp,length(x)))
                else
                    return(tmp)
            }
            ret$fxy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(y)"),")")))

            ret$fxx <- function(x,y){
                myfxx <- Deriv(Deriv(f,"x"),"x");
                tmp <- myfxx(x,y);
                if(length(tmp)==1)
                    return(rep(tmp,length(x)))
                else
                    return(tmp)
            }
            ret$fxx_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(x)"),")")))

            ret$fyy <- function(x,y){
                myfyy <- Deriv(Deriv(f,"y"),"y");
                tmp <- myfyy(x,y);
                if(length(tmp)==1)
                    return(rep(tmp,length(x)))
                else
                    return(tmp)
            }
            ret$fyy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(y)D(y)"),")")))

            if(dg>2){
                ret$fxxy <- function(x,y){
                    myfxxy <- Deriv(Deriv(Deriv(f,"y"),"x"),"x");
                    tmp <- myfxxy(x,y);
                    if(length(tmp)==1)
                        return(rep(tmp,length(x)))
                    else
                        return(tmp)
                }
                ret$fxxy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(x)D(y)"),")")))

                ret$fxyy <- function(x,y){
                    myfxyy <- Deriv(Deriv(Deriv(f,"y"),"y"),"x");
                    tmp <- myfxyy(x,y);
                    if(length(tmp)==1)
                        return(rep(tmp,length(x)))
                    else
                        return(tmp)
                }
                ret$fxyy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(y)D(y)"),")")))

                ret$fxxx <- function(x,y){
                    myfxxx <- Deriv(Deriv(Deriv(f,"x"),"x"),"x");
                    tmp <- myfxxx(x,y);
                    if(length(tmp)==1)
                        return(rep(tmp,length(x)))
                    else
                        return(tmp)
                }
                ret$fxxx_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(x)D(x)"),")")))

                ret$fyyy <- function(x,y){
                    myfyyy <- Deriv(Deriv(Deriv(f,"y"),"y"),"y");
                    tmp <- myfyyy(x,y);
                    if(length(tmp)==1)
                        return(rep(tmp,length(x)))
                    else
                        return(tmp)
                }
                ret$fyyy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(y)D(y)D(y)"),")")))
            }
        }
    }
    ret
}


###################################################
### code chunk number 10: partDeriv.Rnw:583-584
###################################################
df <- derivs(f,dg)


###################################################
### code chunk number 11: partDeriv.Rnw:587-590
###################################################
xg <- seq(0,1,length=ng)
yg <- seq(0,1,length=ng)
xyg <- expand.grid(xg,yg)


###################################################
### code chunk number 12: partDeriv.Rnw:592-593
###################################################
af=4


###################################################
### code chunk number 13: partDeriv.Rnw:597-601
###################################################
af <- 4
xfg <- seq(0,1,length=af*ng)
yfg <- seq(0,1,length=af*ng)
xyfg <- expand.grid(xfg,yfg)


###################################################
### code chunk number 14: partDeriv.Rnw:604-608
###################################################
nx <- length(xg)
ny <- length(yg)
xx <- t(matrix(rep(xg,ny),nx,ny))
yy <- matrix(rep(yg,nx),ny,nx)


###################################################
### code chunk number 15: helperGrid
###################################################
# for plots of exact values
fgrid <- function(f,xg,yg,dg){
  ret <- list(f=outer(xg,yg,f))
  df <- derivs(f,dg)
  if(dg>0){
    ret$fx  <- outer(xg,yg,df$fx)
    ret$fy  <- outer(xg,yg,df$fy)
    if(dg>1){
      ret$fxy <- outer(xg,yg,df$fxy)
      ret$fxx <- outer(xg,yg,df$fxx)
      ret$fyy <- outer(xg,yg,df$fyy)
      if(dg>2){
        ret$fxxy <- outer(xg,yg,df$fxxy)
        ret$fxyy <- outer(xg,yg,df$fxyy)
        ret$fxxx <- outer(xg,yg,df$fxxx)
        ret$fyyy <- outer(xg,yg,df$fyyy)
      }
    }
  }
  ret
}


###################################################
### code chunk number 16: partDeriv.Rnw:636-640
###################################################
## data for local regression
fg   <- outer(xg,yg,f)
## data for exact plots on fine grid
ffg <- fgrid(f,xfg,yfg,dg)


###################################################
### code chunk number 17: partDeriv.Rnw:644-648
###################################################
## global bandwidth:
pdg <- interp::locpoly(xg,yg,fg, input="grid", pd="all", h=c(bwg,bwg), solver="QR", degree=dg,kernel=knl,nx=af*ng,ny=af*ng)
## local bandwidth:
pdl <- interp::locpoly(xg,yg,fg, input="grid", pd="all", h=bwl, solver="QR", degree=dg,kernel=knl,nx=af*ng,ny=af*ng)


###################################################
### code chunk number 18: helperSplit
###################################################
split_str <- function(txt,l){
  start <- seq(1, nchar(txt), l)
  stop <- seq(l, nchar(txt)+l, l)[1:length(start)]
  substring(txt, start, stop)
}


###################################################
### code chunk number 19: helperImage
###################################################
grid2df <- function(x,y,z)
    subset(data.frame(x = rep(x, nrow(z)),
                      y = rep(y, each = ncol(z)),
                      z = as.numeric(z)),
           !is.na(z))

gg1image2contours <- function(x,y,z1,z2,z3,xyg,ttl=""){
    breaks <- pretty(seq(min(z1,na.rm=T),max(z1,na.rm=T),length=11))
    griddf1 <- grid2df(x,y,z1)
    griddf2 <- grid2df(x,y,z2)
    griddf3 <- grid2df(x,y,z3)
    griddf  <- data.frame(x=griddf1$x,y=griddf1$y,z1=griddf1$z,z2=griddf2$z,z3=griddf3$z)
    ggplot(griddf, aes(x=x, y=y, z = z1)) +
        ggtitle(ttl) +
        theme(plot.title = element_text(size = 6, face = "bold"),
              axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +
        geom_contour_filled(breaks=breaks) +
        scale_fill_brewer(palette = "YlOrRd") +
        geom_contour(aes(z=z2),breaks=breaks,color="green",lty="dashed",lwd=0.5) +
        geom_contour(aes(z=z3),breaks=breaks,color="blue",lty="dotted",lwd=0.5) +
        theme(legend.position="none") +
        geom_point(data=xyg, aes(x=Var1,y=Var2), inherit.aes = FALSE,size=1,pch="+")
}


###################################################
### code chunk number 20: helperPrint
###################################################
print_deriv <- function(txt,l,at=42){
    ret<-""
    for(t in txt){
        if(stringi::stri_length(t)<at)
            btxt <- t
        else
            btxt <- split_str(t,at)
        ftxt <- rep(paste(rep(" ",stringi::stri_length(l)),sep="",collapse=""),length(btxt))
        ftxt[1] <- l
        ret <- paste(ret,paste(ftxt,btxt,sep="",collapse = "\n"),sep="",collapse = "\n")
    }
    ret
}
print_f <- function(f,df,dg,offset=0.8){
  lns <- c(print_deriv(df$f_str,"f(x,y) ="))
  if(dg>=1)
    lns <- c(lns,
    print_deriv(df$fx_str,"f_x(x,y) ="),
    print_deriv(df$fy_str,"f_y(x,y) ="))
  if(dg>=2)
    lns <- c(lns,
    print_deriv(df$fxx_str,"f_xx(x,y) ="),
    print_deriv(df$fyy_str,"f_yy(x,y) ="),
    print_deriv(df$fxy_str,"f_xy(x,y) ="))
  if(dg>=3)
    lns <- c(lns,
    print_deriv(df$fxxx_str,"f_xxx(x,y) ="),
    print_deriv(df$fyyy_str,"f_yyy(x,y) ="),
    print_deriv(df$fxxy_str,"f_xxy(x,y) ="),
    print_deriv(df$fxyy_str,"f_xyy(x,y) ="))
  txt <- grid.text(paste(lns,
    collapse="\n"),gp=gpar(fontsize=8),
    x=0,y=offset,draw=FALSE,
    just = c("left","top"))
  txt
}


###################################################
### code chunk number 21: prepare1
###################################################
t1 <- grid.text(paste(c(paste("regular data grid",nx,"x",ny),
  "colors = exaxt values",
  "dashed green = global bw",
  "dotted blue = local bw",
  "crosses: data points"),collapse="\n"),
  gp=gpar(fontsize=8),
  x=0,y=0.8,draw=FALSE,
  just = c("left","top"))

t3 <- grid.text(paste(c(paste("kernel:",knl),
                        paste("global bandwidth",bwg*100,"%"),
                        paste("local bandwidth",bwl*100,"%")),
                      collapse="\n"),
                gp=gpar(fontsize=8),x=0,y=0.8,draw=FALSE,
                just = c("left","top"))



###################################################
### code chunk number 22: partDeriv.Rnw:746-759
###################################################
pf <- gg1image2contours(xfg,yfg,ffg$f,pdg$z,pdl$z,xyg,"f")
pfx <- gg1image2contours(xfg,yfg,ffg$fx,pdg$zx,pdl$zx,xyg,"f_x")
pfy <- gg1image2contours(xfg,yfg,ffg$fy,pdg$zy,pdl$zy,xyg,"f_x")
pfxx <- gg1image2contours(xfg,yfg,ffg$fxx,pdg$zxx,pdl$zxx,xyg,"f_xx")
pfxy <- gg1image2contours(xfg,yfg,ffg$fxy,pdg$zxy,pdl$zxy,xyg,"f_xy")
pfyy <- gg1image2contours(xfg,yfg,ffg$fyy,pdg$zyy,pdl$zyy,xyg,"f_yy")
pfxxx <- gg1image2contours(xfg,yfg,ffg$fxxx,pdg$zxxx,pdl$zxxx,xyg,"f_xxx")
pfxxy <- gg1image2contours(xfg,yfg,ffg$fxxy,pdg$zxxy,pdl$zxxy,xyg,"f_xxy")
pfxyy <- gg1image2contours(xfg,yfg,ffg$fxyy,pdg$zxyy,pdl$zxyy,xyg,"f_xyy")
pfyyy <- gg1image2contours(xfg,yfg,ffg$fyyy,pdg$zyyy,pdl$zyyy,xyg,"f_yyy")
## t1 and t3 contain pure texts generated hidden in this Sweave file.
## t2 contains aas much of the symbolic computation output as possible:
t2 <- print_f(f,df,3)


###################################################
### code chunk number 23: plotbicubic
###################################################
lay<-rbind(c( 1, 2, 3, 3),
           c( 4, 5, 3, 3),
           c( 6, 7, 8, 9),
           c(10,11,12,13))
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 24: partDeriv.Rnw:782-783
###################################################
getOption("SweaveHooks")[["fig"]]()
lay<-rbind(c( 1, 2, 3, 3),
           c( 4, 5, 3, 3),
           c( 6, 7, 8, 9),
           c(10,11,12,13))
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 25: partDeriv.Rnw:791-795
###################################################
f <- function(x,y) 0.75*exp(-((9*x-2)^2+(9*y-2)^2)/4)+0.75*exp(-((9*x+1)^2)/49-(9*y+1)/10)+0.5*exp(-((9*x-7)^2+(9*y-3)^2)/4)-0.2*exp(-(9*x-4)^2-(9*y-7)^2)
fg  <- outer(xg,yg,f)
ffg <- fgrid(f,xfg,yfg,dg)
df  <- derivs(f,dg)


###################################################
### code chunk number 26: partDeriv.Rnw:798-802
###################################################
## global bw,
pdg <- interp::locpoly(xg,yg,fg, input="grid", pd="all", h=c(bwg,bwg), solver="QR", degree=dg,kernel=knl,nx=af*ng,ny=af*ng)
## local bw:
pdl <- interp::locpoly(xg,yg,fg, input="grid", pd="all", h=bwl, solver="QR", degree=dg,kernel=knl,nx=af*ng,ny=af*ng)


###################################################
### code chunk number 27: plotfranke1
###################################################
pf <- gg1image2contours(xfg,yfg,ffg$f,pdg$z,pdl$z,xyg,"f")
pfx <- gg1image2contours(xfg,yfg,ffg$fx,pdg$zx,pdl$zx,xyg,"f_x")
pfy <- gg1image2contours(xfg,yfg,ffg$fy,pdg$zy,pdl$zy,xyg,"f_x")
pfxx <- gg1image2contours(xfg,yfg,ffg$fxx,pdg$zxx,pdl$zxx,xyg,"f_xx")
pfxy <- gg1image2contours(xfg,yfg,ffg$fxy,pdg$zxy,pdl$zxy,xyg,"f_xy")
pfyy <- gg1image2contours(xfg,yfg,ffg$fyy,pdg$zyy,pdl$zyy,xyg,"f_yy")
pfxxx <- gg1image2contours(xfg,yfg,ffg$fxxx,pdg$zxxx,pdl$zxxx,xyg,"f_xxx")
pfxxy <- gg1image2contours(xfg,yfg,ffg$fxxy,pdg$zxxy,pdl$zxxy,xyg,"f_xxy")
pfxyy <- gg1image2contours(xfg,yfg,ffg$fxyy,pdg$zxyy,pdl$zxyy,xyg,"f_xyy")
pfyyy <- gg1image2contours(xfg,yfg,ffg$fyyy,pdg$zyyy,pdl$zyyy,xyg,"f_yyy")

t2 <- print_f(f,df,1,0.9)

gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix=lay)


###################################################
### code chunk number 28: partDeriv.Rnw:825-826
###################################################
getOption("SweaveHooks")[["fig"]]()
pf <- gg1image2contours(xfg,yfg,ffg$f,pdg$z,pdl$z,xyg,"f")
pfx <- gg1image2contours(xfg,yfg,ffg$fx,pdg$zx,pdl$zx,xyg,"f_x")
pfy <- gg1image2contours(xfg,yfg,ffg$fy,pdg$zy,pdl$zy,xyg,"f_x")
pfxx <- gg1image2contours(xfg,yfg,ffg$fxx,pdg$zxx,pdl$zxx,xyg,"f_xx")
pfxy <- gg1image2contours(xfg,yfg,ffg$fxy,pdg$zxy,pdl$zxy,xyg,"f_xy")
pfyy <- gg1image2contours(xfg,yfg,ffg$fyy,pdg$zyy,pdl$zyy,xyg,"f_yy")
pfxxx <- gg1image2contours(xfg,yfg,ffg$fxxx,pdg$zxxx,pdl$zxxx,xyg,"f_xxx")
pfxxy <- gg1image2contours(xfg,yfg,ffg$fxxy,pdg$zxxy,pdl$zxxy,xyg,"f_xxy")
pfxyy <- gg1image2contours(xfg,yfg,ffg$fxyy,pdg$zxyy,pdl$zxyy,xyg,"f_xyy")
pfyyy <- gg1image2contours(xfg,yfg,ffg$fyyy,pdg$zyyy,pdl$zyyy,xyg,"f_yyy")

t2 <- print_f(f,df,1,0.9)

gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix=lay)


###################################################
### code chunk number 29: partDeriv.Rnw:838-839
###################################################
n <- ng*ng


###################################################
### code chunk number 30: partDeriv.Rnw:842-843
###################################################
f <- function(x,y) (x-0.5)*(x-0.2)*(y-0.6)*y*(x-1)


###################################################
### code chunk number 31: partDeriv.Rnw:846-851
###################################################
## random irregular data
x<-runif(n)
y<-runif(n)
xy<-data.frame(Var1=x,Var2=y)
z <- f(x,y)


###################################################
### code chunk number 32: partDeriv.Rnw:854-856
###################################################
ffg <- fgrid(f,xfg,yfg,dg)
df <- derivs(f,dg)


###################################################
### code chunk number 33: partDeriv.Rnw:859-863
###################################################
## global bandwidth
pdg <- interp::locpoly(x,y,z, xfg,yfg, pd="all", h=c(bwg,bwg), solver="QR", degree=dg,kernel=knl)
## local bandwidth:
pdl <- interp::locpoly(x,y,z, xfg,yfg, pd="all", h=bwl, solver="QR", degree=dg,kernel=knl)


###################################################
### code chunk number 34: partDeriv.Rnw:868-889
###################################################
pf <- gg1image2contours(xfg,yfg,ffg$f,pdg$z,pdl$z,xy,"f")
pfx <- gg1image2contours(xfg,yfg,ffg$fx,pdg$zx,pdl$zx,xy,"f_x")
pfy <- gg1image2contours(xfg,yfg,ffg$fy,pdg$zy,pdl$zy,xy,"f_x")
pfxx <- gg1image2contours(xfg,yfg,ffg$fxx,pdg$zxx,pdl$zxx,xy,"f_xx")
pfxy <- gg1image2contours(xfg,yfg,ffg$fxy,pdg$zxy,pdl$zxy,xy,"f_xy")
pfyy <- gg1image2contours(xfg,yfg,ffg$fyy,pdg$zyy,pdl$zyy,xy,"f_yy")
pfxxx <- gg1image2contours(xfg,yfg,ffg$fxxx,pdg$zxxx,pdl$zxxx,xy,"f_xxx")
pfxxy <- gg1image2contours(xfg,yfg,ffg$fxxy,pdg$zxxy,pdl$zxxy,xy,"f_xxy")
pfxyy <- gg1image2contours(xfg,yfg,ffg$fxyy,pdg$zxyy,pdl$zxyy,xy,"f_xyy")
pfyyy <- gg1image2contours(xfg,yfg,ffg$fyyy,pdg$zyyy,pdl$zyyy,xy,"f_yyy")

t1 <- grid.text(paste(c(paste("irregular data grid",n,"pts"),
  "colors = exaxt values",
  "dashed green = global bw",
  "dotted blue = local bw",
  "crosses: data points"),collapse="\n"),
  gp=gpar(fontsize=8),
  x=0,y=0.8,draw=FALSE,
  just = c("left","top"))

t2 <- print_f(f,df,3)


###################################################
### code chunk number 35: plotbicubic2
###################################################
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 36: partDeriv.Rnw:896-897
###################################################
getOption("SweaveHooks")[["fig"]]()
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 37: partDeriv.Rnw:903-904
###################################################
f <- function(x,y) 0.75*exp(-((9*x-2)^2+(9*y-2)^2)/4)+0.75*exp(-((9*x+1)^2)/49-(9*y+1)/10)+0.5*exp(-((9*x-7)^2+(9*y-3)^2)/4)-0.2*exp(-(9*x-4)^2-(9*y-7)^2)


###################################################
### code chunk number 38: partDeriv.Rnw:906-910
###################################################
z <- f(x,y)
fg  <- outer(xg,yg,f)
ffg <- fgrid(f,xfg,yfg,dg)
df <- derivs(f,dg)


###################################################
### code chunk number 39: partDeriv.Rnw:912-917
###################################################
## global bandwidth:
ttg <- system.time(pdg <- interp::locpoly(x,y,z, xfg,yfg, pd="all", h=c(bwg,bwg), solver="QR", degree=dg,kernel=knl))

## local bandwidth:
ttl <- system.time(pdl <- interp::locpoly(x,y,z, xfg,yfg, pd="all", h=bwl, solver="QR", degree=dg,kernel=knl))


###################################################
### code chunk number 40: partDeriv.Rnw:919-931
###################################################
pf <- gg1image2contours(xfg,yfg,ffg$f,pdg$z,pdl$z,xy,"f")
pfx <- gg1image2contours(xfg,yfg,ffg$fx,pdg$zx,pdl$zx,xy,"f_x")
pfy <- gg1image2contours(xfg,yfg,ffg$fy,pdg$zy,pdl$zy,xy,"f_x")
pfxx <- gg1image2contours(xfg,yfg,ffg$fxx,pdg$zxx,pdl$zxx,xy,"f_xx")
pfxy <- gg1image2contours(xfg,yfg,ffg$fxy,pdg$zxy,pdl$zxy,xy,"f_xy")
pfyy <- gg1image2contours(xfg,yfg,ffg$fyy,pdg$zyy,pdl$zyy,xy,"f_yy")
pfxxx <- gg1image2contours(xfg,yfg,ffg$fxxx,pdg$zxxx,pdl$zxxx,xy,"f_xxx")
pfxxy <- gg1image2contours(xfg,yfg,ffg$fxxy,pdg$zxxy,pdl$zxxy,xy,"f_xxy")
pfxyy <- gg1image2contours(xfg,yfg,ffg$fxyy,pdg$zxyy,pdl$zxyy,xy,"f_xyy")
pfyyy <- gg1image2contours(xfg,yfg,ffg$fyyy,pdg$zyyy,pdl$zyyy,xy,"f_yyy")

t2 <- print_f(f,df,1,0.9)


###################################################
### code chunk number 41: plotfranke12
###################################################
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 42: partDeriv.Rnw:938-939
###################################################
getOption("SweaveHooks")[["fig"]]()
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 43: partDeriv.Rnw:950-954
###################################################
## global bandwidth:
pdg <- interp::locpoly(x,y,z, xfg,yfg, pd="all", h=c(bwg,bwg), solver="QR", degree=dg,kernel="uniform")
## local bandwidth:
pdl <- interp::locpoly(x,y,z, xfg,yfg, pd="all", h=bwl, solver="QR", degree=dg,kernel="uniform")


###################################################
### code chunk number 44: partDeriv.Rnw:956-975
###################################################
pf <- gg1image2contours(xfg,yfg,ffg$f,pdg$z,pdl$z,xy,"f")
pfx <- gg1image2contours(xfg,yfg,ffg$fx,pdg$zx,pdl$zx,xy,"f_x")
pfy <- gg1image2contours(xfg,yfg,ffg$fy,pdg$zy,pdl$zy,xy,"f_x")
pfxx <- gg1image2contours(xfg,yfg,ffg$fxx,pdg$zxx,pdl$zxx,xy,"f_xx")
pfxy <- gg1image2contours(xfg,yfg,ffg$fxy,pdg$zxy,pdl$zxy,xy,"f_xy")
pfyy <- gg1image2contours(xfg,yfg,ffg$fyy,pdg$zyy,pdl$zyy,xy,"f_yy")
pfxxx <- gg1image2contours(xfg,yfg,ffg$fxxx,pdg$zxxx,pdl$zxxx,xy,"f_xxx")
pfxxy <- gg1image2contours(xfg,yfg,ffg$fxxy,pdg$zxxy,pdl$zxxy,xy,"f_xxy")
pfxyy <- gg1image2contours(xfg,yfg,ffg$fxyy,pdg$zxyy,pdl$zxyy,xy,"f_xyy")
pfyyy <- gg1image2contours(xfg,yfg,ffg$fyyy,pdg$zyyy,pdl$zyyy,xy,"f_yyy")

t2 <- print_f(f,df,1,0.9)
t3 <- grid.text(paste(c(paste("kernel:","uniform"),
                        paste("global bandwidth",bwg*100,"%"),
                        paste("local bandwidth",bwl*100,"%")),
                      collapse="\n"),
                gp=gpar(fontsize=8),x=0,y=0.8,draw=FALSE,
                just = c("left","top"))



###################################################
### code chunk number 45: plotfranke12unif
###################################################
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 46: partDeriv.Rnw:982-983
###################################################
getOption("SweaveHooks")[["fig"]]()
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 47: partDeriv.Rnw:988-992
###################################################
## global bandwidth:
pdg <- interp::locpoly(x,y,z, xfg,yfg, pd="all", h=c(bwg,bwg), solver="QR", degree=dg,kernel="epanechnikov")
## local bandwidth:
pdl <- interp::locpoly(x,y,z, xfg,yfg, pd="all", h=bwl, solver="QR", degree=dg,kernel="epanechnikov")


###################################################
### code chunk number 48: partDeriv.Rnw:994-1012
###################################################
pf <- gg1image2contours(xfg,yfg,ffg$f,pdg$z,pdl$z,xy,"f")
pfx <- gg1image2contours(xfg,yfg,ffg$fx,pdg$zx,pdl$zx,xy,"f_x")
pfy <- gg1image2contours(xfg,yfg,ffg$fy,pdg$zy,pdl$zy,xy,"f_x")
pfxx <- gg1image2contours(xfg,yfg,ffg$fxx,pdg$zxx,pdl$zxx,xy,"f_xx")
pfxy <- gg1image2contours(xfg,yfg,ffg$fxy,pdg$zxy,pdl$zxy,xy,"f_xy")
pfyy <- gg1image2contours(xfg,yfg,ffg$fyy,pdg$zyy,pdl$zyy,xy,"f_yy")
pfxxx <- gg1image2contours(xfg,yfg,ffg$fxxx,pdg$zxxx,pdl$zxxx,xy,"f_xxx")
pfxxy <- gg1image2contours(xfg,yfg,ffg$fxxy,pdg$zxxy,pdl$zxxy,xy,"f_xxy")
pfxyy <- gg1image2contours(xfg,yfg,ffg$fxyy,pdg$zxyy,pdl$zxyy,xy,"f_xyy")
pfyyy <- gg1image2contours(xfg,yfg,ffg$fyyy,pdg$zyyy,pdl$zyyy,xy,"f_yyy")

t2 <- print_f(f,df,1,0.9)
t3 <- grid.text(paste(c(paste("kernel:","epanechnikov"),
                        paste("global bandwidth",bwg*100,"%"),
                        paste("local bandwidth",bwl*100,"%")),
                      collapse="\n"),
                gp=gpar(fontsize=8),x=0,y=0.8,draw=FALSE,
                just = c("left","top"))


###################################################
### code chunk number 49: plotfranke12epa
###################################################
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 50: partDeriv.Rnw:1019-1020
###################################################
getOption("SweaveHooks")[["fig"]]()
gg <- grid.arrange(grobs=gList(ggplotGrob(pf),t1,t2,ggplotGrob(pfx),ggplotGrob(pfy),ggplotGrob(pfxx),ggplotGrob(pfxy),ggplotGrob(pfyy),t3,ggplotGrob(pfxxx),ggplotGrob(pfxxy),ggplotGrob(pfxyy),ggplotGrob(pfyyy)),layout_matrix = lay)


###################################################
### code chunk number 51: partDeriv.Rnw:1035-1037
###################################################
# helper functions for translation between R and Yacas
fn_y  <- function(f){
    b <- toString(as.expression(body(f)))
    b <- stringr::str_replace_all(b,"cos","Cos")
    b <- stringr::str_replace_all(b,"sin","Sin")
    b <- stringr::str_replace_all(b,"exp","Exp")
    b <- stringr::str_replace_all(b,"log","Log")
    b <- stringr::str_replace_all(b,"sqrt","Sqrt")
    b
}
ys_fn  <- function(f){
    f <- stringr::str_replace_all(f,"Cos","cos")
    f <- stringr::str_replace_all(f,"Sin","sin")
    f <- stringr::str_replace_all(f,"Exp","exp")
    f <- stringr::str_replace_all(f,"Log","log")
    f <- stringr::str_replace_all(f,"Sqrt","sqrt")
    f
}


###################################################
### code chunk number 52: partDeriv.Rnw:1044-1045
###################################################
derivs <- function(f,dg){
    ret<-list(f=f,
              f_str=ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),""),")"))))

    if(dg>0){

        ret$fx <- function(x,y){
            myfx <- Deriv(f,"x");
            tmp <- myfx(x,y);
            if(length(tmp)==1)
                return(rep(tmp,length(x)))
            else
                return(tmp)
        }
        ret$fx_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)"),")")))


        ret$fy <- function(x,y){
            myfy <- Deriv(f,"y");
            tmp <- myfy(x,y);
            if(length(tmp)==1)
                return(rep(tmp,length(x)))
            else
                return(tmp)
        }
        ret$fy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(y)"),")")))


        if(dg>1){
            ret$fxy <- function(x,y){
                myfxy <- Deriv(Deriv(f,"y"),"x");
                tmp <- myfxy(x,y);
                if(length(tmp)==1)
                    return(rep(tmp,length(x)))
                else
                    return(tmp)
            }
            ret$fxy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(y)"),")")))

            ret$fxx <- function(x,y){
                myfxx <- Deriv(Deriv(f,"x"),"x");
                tmp <- myfxx(x,y);
                if(length(tmp)==1)
                    return(rep(tmp,length(x)))
                else
                    return(tmp)
            }
            ret$fxx_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(x)"),")")))

            ret$fyy <- function(x,y){
                myfyy <- Deriv(Deriv(f,"y"),"y");
                tmp <- myfyy(x,y);
                if(length(tmp)==1)
                    return(rep(tmp,length(x)))
                else
                    return(tmp)
            }
            ret$fyy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(y)D(y)"),")")))

            if(dg>2){
                ret$fxxy <- function(x,y){
                    myfxxy <- Deriv(Deriv(Deriv(f,"y"),"x"),"x");
                    tmp <- myfxxy(x,y);
                    if(length(tmp)==1)
                        return(rep(tmp,length(x)))
                    else
                        return(tmp)
                }
                ret$fxxy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(x)D(y)"),")")))

                ret$fxyy <- function(x,y){
                    myfxyy <- Deriv(Deriv(Deriv(f,"y"),"y"),"x");
                    tmp <- myfxyy(x,y);
                    if(length(tmp)==1)
                        return(rep(tmp,length(x)))
                    else
                        return(tmp)
                }
                ret$fxyy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(y)D(y)"),")")))

                ret$fxxx <- function(x,y){
                    myfxxx <- Deriv(Deriv(Deriv(f,"x"),"x"),"x");
                    tmp <- myfxxx(x,y);
                    if(length(tmp)==1)
                        return(rep(tmp,length(x)))
                    else
                        return(tmp)
                }
                ret$fxxx_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(x)D(x)D(x)"),")")))

                ret$fyyy <- function(x,y){
                    myfyyy <- Deriv(Deriv(Deriv(f,"y"),"y"),"y");
                    tmp <- myfyyy(x,y);
                    if(length(tmp)==1)
                        return(rep(tmp,length(x)))
                    else
                        return(tmp)
                }
                ret$fyyy_str  <- ys_fn(yac(paste("Simplify(",y_fn(fn_y(f),"D(y)D(y)D(y)"),")")))
            }
        }
    }
    ret
}


###################################################
### code chunk number 53: partDeriv.Rnw:1051-1052
###################################################
# for plots of exact values
fgrid <- function(f,xg,yg,dg){
  ret <- list(f=outer(xg,yg,f))
  df <- derivs(f,dg)
  if(dg>0){
    ret$fx  <- outer(xg,yg,df$fx)
    ret$fy  <- outer(xg,yg,df$fy)
    if(dg>1){
      ret$fxy <- outer(xg,yg,df$fxy)
      ret$fxx <- outer(xg,yg,df$fxx)
      ret$fyy <- outer(xg,yg,df$fyy)
      if(dg>2){
        ret$fxxy <- outer(xg,yg,df$fxxy)
        ret$fxyy <- outer(xg,yg,df$fxyy)
        ret$fxxx <- outer(xg,yg,df$fxxx)
        ret$fyyy <- outer(xg,yg,df$fyyy)
      }
    }
  }
  ret
}


###################################################
### code chunk number 54: partDeriv.Rnw:1056-1057
###################################################
split_str <- function(txt,l){
  start <- seq(1, nchar(txt), l)
  stop <- seq(l, nchar(txt)+l, l)[1:length(start)]
  substring(txt, start, stop)
}


###################################################
### code chunk number 55: partDeriv.Rnw:1061-1062
###################################################
grid2df <- function(x,y,z)
    subset(data.frame(x = rep(x, nrow(z)),
                      y = rep(y, each = ncol(z)),
                      z = as.numeric(z)),
           !is.na(z))

gg1image2contours <- function(x,y,z1,z2,z3,xyg,ttl=""){
    breaks <- pretty(seq(min(z1,na.rm=T),max(z1,na.rm=T),length=11))
    griddf1 <- grid2df(x,y,z1)
    griddf2 <- grid2df(x,y,z2)
    griddf3 <- grid2df(x,y,z3)
    griddf  <- data.frame(x=griddf1$x,y=griddf1$y,z1=griddf1$z,z2=griddf2$z,z3=griddf3$z)
    ggplot(griddf, aes(x=x, y=y, z = z1)) +
        ggtitle(ttl) +
        theme(plot.title = element_text(size = 6, face = "bold"),
              axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +
        geom_contour_filled(breaks=breaks) +
        scale_fill_brewer(palette = "YlOrRd") +
        geom_contour(aes(z=z2),breaks=breaks,color="green",lty="dashed",lwd=0.5) +
        geom_contour(aes(z=z3),breaks=breaks,color="blue",lty="dotted",lwd=0.5) +
        theme(legend.position="none") +
        geom_point(data=xyg, aes(x=Var1,y=Var2), inherit.aes = FALSE,size=1,pch="+")
}


###################################################
### code chunk number 56: partDeriv.Rnw:1066-1067
###################################################
print_deriv <- function(txt,l,at=42){
    ret<-""
    for(t in txt){
        if(stringi::stri_length(t)<at)
            btxt <- t
        else
            btxt <- split_str(t,at)
        ftxt <- rep(paste(rep(" ",stringi::stri_length(l)),sep="",collapse=""),length(btxt))
        ftxt[1] <- l
        ret <- paste(ret,paste(ftxt,btxt,sep="",collapse = "\n"),sep="",collapse = "\n")
    }
    ret
}
print_f <- function(f,df,dg,offset=0.8){
  lns <- c(print_deriv(df$f_str,"f(x,y) ="))
  if(dg>=1)
    lns <- c(lns,
    print_deriv(df$fx_str,"f_x(x,y) ="),
    print_deriv(df$fy_str,"f_y(x,y) ="))
  if(dg>=2)
    lns <- c(lns,
    print_deriv(df$fxx_str,"f_xx(x,y) ="),
    print_deriv(df$fyy_str,"f_yy(x,y) ="),
    print_deriv(df$fxy_str,"f_xy(x,y) ="))
  if(dg>=3)
    lns <- c(lns,
    print_deriv(df$fxxx_str,"f_xxx(x,y) ="),
    print_deriv(df$fyyy_str,"f_yyy(x,y) ="),
    print_deriv(df$fxxy_str,"f_xxy(x,y) ="),
    print_deriv(df$fxyy_str,"f_xyy(x,y) ="))
  txt <- grid.text(paste(lns,
    collapse="\n"),gp=gpar(fontsize=8),
    x=0,y=offset,draw=FALSE,
    just = c("left","top"))
  txt
}


