### R code from vignette source 'magic.Rnw'

###################################################
### code chunk number 1: magic.Rnw:101-101
###################################################



###################################################
### code chunk number 2: magic.Rnw:102-103
###################################################
require(magic)


###################################################
### code chunk number 3: magic.Rnw:109-110
###################################################
 magic(3) 


###################################################
### code chunk number 4: magic.Rnw:123-124
###################################################
magicplot(magic.2np1(3))


###################################################
### code chunk number 5: magic.Rnw:175-197
###################################################
shadedsquare <- function(m=2){
  n <- 4*m
  jj.1 <- kronecker(diag(2), matrix(1, 2, 2))
  jj <- kronecker(matrix(1, m + 1, m + 1), jj.1)[2:(n + 1), 2:(n + 1)]
  par(xaxt="n",yaxt="n")
  image(1:n,1:n,jj,xlab="",ylab="",asp=1,frame=FALSE,col=c(gray(0.9),gray(0.4)))
  abline(v=0.5+(0:n))
  segments(x0=rep(0.5,n),y0=0.5+(0:n),x1=rep(n+0.5,n),y1=0.5+(0:n))
  return(invisible(jj))
}

jj <- shadedsquare()
#a <- magic(8)
#text(row(a),col(a),as.character(a),col="white")

for(i in 1:8){
  for(j in 1:8){
    if(jj[i,j]==1){
      text(i,j,magic(8)[i,9-j],col="white")
    }
  }
}


###################################################
### code chunk number 6: magic.Rnw:206-216
###################################################
shadedsquare()
for(i in 1:8){
  for(j in 1:8){
    if(jj[i,j]==1){
      text(i,j,magic(8)[i,9-j],col="white")
         } else {
                 text(i,j,magic(8)[i,9-j],col="black")
               }
  }
}


