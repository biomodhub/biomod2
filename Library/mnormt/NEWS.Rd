\name{NEWS}
\title{News (Changelog) for Package \pkg{mnormt}}
\encoding{UTF-8}

\section{Changes in mnormt version 2.1.2 (2026-01-25)}{
 \itemize{    
   \item Some declarations of variables, notably \code{DOUBLE PRECISION TVTMFN} 
     in \code{tvpack.f}, have been inserted in Fortran source files, as pointed
     out by the CRAN group. 
 \item Fixing some syntax glitches in the Rd documentation files, pointed
    out by the CRAN group.     
  }
}

\section{Changes in mnormt version 2.1.1 (2022-09-26)}{
  \itemize{
    \item When \code{sadmvn, pmnorm, sadmvt, pmt} are called with dimension 
       larger than 20, the returned value is now \code{NA},
       instead of 0, with \code{attr(,"status")="oversize"}.
    \item The NEWS file is now in Rd format.  
  }
}
\section{Changes in mnormt version 2.1.0 (2022-06-07)}{
  \itemize{
    \item Function \code{dmtruncnorm} is now based on Fortran coding of a Gibbs
        sampling, removing an incorrect dependence on package \pkg{tmvnsim}. 
    \item Function \code{plot_fxy} is introduced. 
    \item Now \code{\{d,p,r\}mtruncnorm} can be called with parameters 
        which actually represent a univariate distribution. 
    \item In \code{recintab}, \code{biv.nt.prob} is used instead of \code{sadmvn} 
       when \code{d=2}. 
  }
}
\section{Changes in mnormt version 2.0.2 (2020-08-31)}{
  \itemize{
    \item R version requirement is reduced from 4.0.0 to 2.2.0.
  }
}
\section{Changes in mnormt version 2.0.1 (2020-06-24)}{
  \itemize{
    \item Fixed a bug which effectively prevented using \code{dmtruncnorm} and 
       \code{dmtrunct} with an argument \code{x} of vector type.
    \item Improvements in the documentation of \code{recintab} and \code{mom2cum}.
  }
}
\section{Changes in mnormt version 2.0.0 (2020-06-02)}{
  \itemize{
    \item Support for the truncated versions of the multivariate normal and the 
       \emph{t} distributions is introduced. 
    \item For the truncated normal distribution, functions \code{recintab} and 
      \code{mom2cum} compute moments and cumulants.
    \item Function \code{mom2cum} can also be used for obtaining the 
       cumulants of other distributions of which the moments are known.
    \item Evaluation of the normal and \emph{t} distribution function in the 
       trivariate case now uses specific Fortran code.     
  }
}
\section{Changes in mnormt version 1.5-7 (2020-04-29)}{
  \itemize{
    \item Improved R coding of \code{sadmvn} and \code{sadmvt}, with more 
      extended use of \code{biv.nt.prob}.
    \item Fixed a bug of \code{pmnorm} affecting a peculiar input set.
    \item New function \code{sample_Mardia_measures} is introduced.
  }
}
\section{Changes in mnormt version 1.5-6 (2020-02-02)}{
  \itemize{
    \item Fixed problem connected to use of 'unintialized' variable in Fortan 77 
      code, as claimed by valgrind, with thanks to Professor Brian Ripley for 
      sorting out this issue;
    \item Replaced \code{is.matrix(.)} with code compliant with \R 4.0.0 
  }
}
\section{Changes in mnormt version 1.5-5 (2016-10-15)}{
  \itemize{
    \item Fix a bug in rmt; this drops the earlier property that repeated calls with 
      a larger \code{n} and reset set.seed(.) started with the shorter sequence.
  }
}
\section{Changes in mnormt version 1.5-4 (2016-03-08)}{
  \itemize{
    \item Improved NAMESPACE and documentation.
    \item More flexible use of \code{mean} in \code{\{d,p,r\}mnorm} and in 
      \code{\{d,p,r\}mt}.
    \item \code{pd.solve} transfers \code{dimnames(x)} to output.
  }
}
\section{Changes in mnormt version 1.5-3 (2015-05-25)}{
  \itemize{
    \item \code{rmt} and \code{rmnorm} output with \code{d=1} have now the same 
       structure.
    \item improved checks on arguments supplied to \code{dmnorm} and \code{dmt}.                  
  }
}
\section{Changes in mnormt version 1.5-2 (2015-04-02)}{
  \itemize{
    \item Introduced new argument \code{sqrt} of \code{dmnorm} and \code{dmt}.
    \item Limited check of the arguments is introduced (skipped for speed 
       in earlier versions).
    \item Calls to \code{rmt}  with a smaller/larger of \code{n} generate a 
      sub/superset of values.
  }
}
\section{Changes in mnormt version 1.5-1 (2014-06-30)}{
  \itemize{
    \item fixes a bug of parameter checking in \code{dmnorm} and \code{dmt}.
  }
}
\section{Changes in mnormt version 1.5-0 (2014-06-25)}{
  \itemize{
    \item Improved coding of \code{dmnorm} and \code{dmt} functions, 
       avoiding use of \code{apply}.
    \item Allow calling \code{pmnorm} and \code{pmt} with a matrix argument \code{x}.
  }
}
\section{Changes in mnormt version 1.4-7 (2013-12-04)}{
  \itemize{
    \item \code{.First.lib()} removed on request from CRAN
    \item other minor fixes.
  }
}
\section{Changes in mnormt version 1.4-5 (2012-01-06)}{
  \itemize{
    \item NAMESPACE introduced on request by the R Core Team.
  }
}
\section{Changes in mnormt version 1.4-3 (2011-05-03)}{
  \itemize{
    \item fixed bug similar to the one of 1.4-2 in \code{sadmvn} and \code{sadmvt}.
  }
}
\section{Changes in mnormt version 1.4-2 (2011-04-28)}{
  \itemize{
    \item fixed bug of \code{biv.not.prob} if some components had
        \code{(lower=-Inf & upper=Inf)}.
  }
}
\section{Changes in mnormt version 1.4-1 (2011-04-06)}{
  \itemize{
    \item Fixed inadequate documentation 
    \item coding of \code{dmnorm} and \code{dmt} is straightened;
    \item \code{log.det} in \code{pd.solve} is returned only on request.
  }
}
\section{Changes in mnormt version 1.4-0 (2011-01-16)}{
  \itemize{
    \item Function \code{pd.solve} is introduced
    \item argument \code{mean} can now be a matrix for \code{\{d,p\}mnorm} 
       and \code{\{d,p\}mt}.
  }
}
\section{Changes in mnormt version 1.3-4 (2009-11-25)}{
  \itemize{
    \item fixed a bug affecting the case \code{d>2 & df=Inf}.
  }
}
\section{Changes in mnormt version 1.3-3 (2009-03-21)}{
  \itemize{
    \item removed check on the R version on request of the R Core team.
  }
}
\section{Changes in mnormt version 1.3-2 (2009-01-26)}{
  \itemize{
    \item minor modification to Fortran 77 code to avoid GOTO jumps out of IF blocks.
  }
}
\section{Changes in mnormt version 1.3-0 (2008-06-12)}{
  \itemize{
    \item  evaluation of normal and \emph{t} distribution functions in the 
      bivariate case now uses specific Fortran code 
    \item improved R coding.
  }
}
\section{Changes in mnormt version 1.2-1 (2007-03-16)}{
  \itemize{
    \item some amendments in documentation.
  }
}
\section{Changes in mnormt version 1.2-0 (2006-10-26)}{
  \itemize{
    \item fixed a bug in function \code{rmt}.
  }
}
\section{Changes in mnormt version 1.1-0 (2006-10-09)}{
  \itemize{
    \item amedend Fortran 77 code of Alan Genz (subtroutine TRESTR).
  }
}
\section{Changes in mnormt version 1.0-1 (2006-04-26)}{
  \itemize{
    \item fixed a bug of non-compliance to R programming standards.
  }
}
\section{Changes in mnormt version 1.0-0 (2006-01-23)}{
  \itemize{
    \item first version uploaded on CRAN. 
  }
}
\section{Initial mnormt version 0.0-1 (2005-11-26)}{
  \itemize{
    \item built private version.
  }
}
