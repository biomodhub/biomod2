collinChk <- function(i,j,k,x,y,eps=NULL) {
xt <- x[c(i,j,k)]
yt <- y[c(i,j,k)]
xx <- .Fortran(
          "cross",
          xt=as.double(xt),
          yt=as.double(yt),
          ijk=as.integer(0),
          cprd=double(1),
          PACKAGE='deldir'
       )
cprd <- xx$cprd
if(is.null(eps)) eps <- sqrt(.Machine$double.eps)
rslt <- abs(cprd) < eps
attr(rslt,"crossproduct") <- cprd
rslt
}
