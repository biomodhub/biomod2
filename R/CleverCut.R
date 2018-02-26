.CleverCut <- function(x){
#### old version
#   switch(EXPR=x,
#          '1' = return(c(1,1)),
#          '2' = return(c(1,2)),
#          '3' = return(c(2,2)),
#          '4' = return(c(2,2)),
#          '5' = return(c(2,3)),
#          '6' = return(c(2,3)),
#          '7' = return(c(3,3)),
#          '8' = return(c(3,3)),
#          return(c(3,3)))
  
  nb_col = ceiling(sqrt(x))
  nb_row = ceiling(x/nb_col)
  return(c(nb_row,nb_col))
}

.bmCat <- function(x=NULL,...){
  if(is.null(x)){
    cat("\n")
    cat(paste(rep("-=", round(.Options$width/2) ), collapse=""))
    cat("\n")
  } else{
    x.length = nchar(x) + 2
    y.length = (.Options$width - x.length) / 2
    cat("\n")
    cat(paste(rep("-=", round(y.length/2) ), collapse=""), x, paste(rep("-=", round(y.length/2) ), collapse=""))
    cat("\n")
  }
}