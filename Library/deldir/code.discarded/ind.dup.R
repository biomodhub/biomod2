ind.dup <- function(x,y,rw=NULL,frac=0.0001) {
#
# Function ind.dup to calculate the indices of data pairs 
# which duplicate earlier ones.  (Returns a logical vector;
# true for such indices, false for the rest.)
#

if(is.null(rw)) rw <- c(0,1,0,1)
n <- length(x)
rslt <- .Fortran(
		'inddup',
		x=as.double(x),
		y=as.double(y),
		n=as.integer(n),
		rw=as.double(rw),
		frac=as.double(frac),
		dup=logical(n),
		PACKAGE='deldir'
	)

rslt$dup
}
