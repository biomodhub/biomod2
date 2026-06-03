#### Will be sourced by several R scripts in ../tests/

source(system.file("test-tools-1.R", package="Matrix"), keep.source=FALSE)

expm.t.identity <- function(x, method,
                            tol = .Machine$double.eps^0.5,
                            check.attributes = FALSE,
                            ...)
{
  ## Purpose: Test the identity   expm(A') = (expm(A))'
  ## ----------------------------------------------------------------------
  ## Arguments: method, ... :          arguments to  expm()
  ##            tol, check.attributes: arguments to  all.equal()
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 23 Feb 2008, 17:26
    ex <- expm::expm(x   , method=method, ...)
    et <- expm::expm(t(x), method=method, ...)
    all.equal(t(ex), et, tolerance = tol, check.attributes = check.attributes)
}


### This is similar to Matrix'  example(spMatrix) :
##' @title random sparse matrix
##' @param nrow,ncol dimension
##' @param ncol
##' @param nnz number of non-zero entries
##' @param density
##' @param rand.x random number generator for 'x' slot
##' @return an  nrow x ncol  matrix
##' @author Martin Maechler, 14.-16. May 2007
rSpMatrix <- function(nrow, ncol = nrow, density, nnz = density*nrow*ncol,
                      sparse = FALSE,
		      rand.x = function(n) round(100 * rnorm(n)))
{
    stopifnot((nnz <- as.integer(nnz)) >= 0,
	      nrow >= 0, ncol >= 0,
	      nnz <= nrow * ncol)
    xx <- rand.x(nnz)
    ## unfortunately, the two resulting matrices might *not* be identical:
    ## because the x's of repeated  (i,j)'s will be *added* for sparse, but not dense:
    ## set.seed(11); m <- rSpMatrix(12, density = 1/10)
    ## set.seed(11); M <- rSpMatrix(12, density = 1/10, sparse=TRUE)
    if(sparse)
	spMatrix(nrow, ncol,
		 i = sample(nrow, nnz, replace = TRUE),
		 j = sample(ncol, nnz, replace = TRUE), x = xx)
    else {
	m <- matrix(0, nrow, ncol)
	m[cbind(i = sample(nrow, nnz, replace = TRUE),
		j = sample(ncol, nnz, replace = TRUE))] <- xx
	m
    }
}

zeroTrace <- function(m) {
    ## Make the {average} trace to 0  -- as it is inside expm(. "Ward77")
    ## This version also works for 'Matrices'
    stopifnot(length(dim(m)) == 2,
              is.numeric(dd <- diag(m)))
    diag(m) <- dd - mean(dd)
    m
}

uniqEntries <- function(m, diagS = FALSE)
{
    ## Purpose: make the non-zero entries of matrix 'm' ``unique''
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 26 Feb 2008, 14:40
    m[m > 0] <-  seq_len(sum(m > 0))
    m[m < 0] <- -seq_len(sum(m < 0))
    if(diagS)
        diag(m) <- 10 * sign(diag(m))
    m
}


## This needs "Matrix" package
rMat <- function(n, R_FUN = rnorm,
                 rcondMin = 1.4 * n ^ -1.6226,
                 iterMax = 100)
{
    ## Purpose: random square matrix "not close to singular"
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## NOTE: needs  Matrix::rcond() -- 2023-11: WHY? {it has more norm = "<c>", but..}
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Jan 2008
    ##
    ##--> /u/maechler/R/MM/Pkg-ex/Matrix/rcondition-numb.R  researches rcond( <random mat>)
    ## Result :
    ##   -log[rcond] = log(Kappa) = 1.051 + 1.6226 * log(n)
    ##   ==================================================
    ##   1/rcond = Kappa = exp(1.051 + 1.6226 * log(n))
    ##                   = 2.8605 * n ^ 1.6226
    ##   ==================================================

    ## since we *search* a bit, take a factor ~ 4  higher rcond:
    ##  4 / 2.8605 ~ 1.4 --> default of rcondMin  above

    stopifnot(require("Matrix")) # needs also as(*, ..) etc
    it <- 1
    rcOpt <- 0
    repeat {
        M <- matrix(R_FUN(n^2), n,n)
        if((rc <- Matrix::rcond(M)) >= rcondMin) break
        if(rc > rcOpt) {
            rcOpt <- rc
            M.Opt <- M
        }
        if((it <- it+1) > iterMax) {
            warning("No Matrix found with rcond() >= ",format(rcondMin),
                    "\n Achieved rcond() = ", format(rcOpt),"\n")
            M <- M.Opt
            break
        }
    }
    M
}

##' call  expm(A, <meth>)  for (all possible) methods <meth> and do catch errors
expmAll <- function(A, meths = eval(formals(expm)$method), errFUN = conditionMessage) {
    if(!missing(meths)) stopifnot(meths %in% eval(formals(expm)$method))
    sapply(meths, simplify = FALSE, function(mtd)
        tryCatch(expm(A, method = mtd), error = errFUN))
}

##' Are they "equal" -- typically applied to result of expmAll()
allEq <- function(Lst, iBest = 1L, check.attributes=FALSE, tol = 1e-10, ...) {
    stopifnot(!is.na(iB <- as.integer(iBest)), length(iB) == 1L, 1L <= iB, iB <= length(Lst),
              is.list(Lst))
    sapply(Lst[-iB], simplify=FALSE, # not vapply() : result TRUE or "..."
           function(R) all.equal(Lst[[iB]], R,
                                 check.attributes=check.attributes, tolerance=tol, ...))
}




doExtras <- interactive() || nzchar(Sys.getenv("R_EXPM_CHECK_EXTRA")) ||
    identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
