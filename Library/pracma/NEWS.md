# pracma NEWS


## pracma 2.4.6 (2025-10-20)

* Corrected erfz() with erfz(0) = 0.
* simpson() DeepSeek code for the 1/3 Simpson rule.

## pracma 2.4.5 (2024-05-23)

* Changed midpoint() to run without Richardson approximation.
* In fmincon(), added '...' variables for the numerical gradient.

## pracma 2.4.4 (2023-11-08)

* gcd(), Lcm() now require integers as input (G. Rücker reported).
* Corrections (missing brackets) in several Rd files.

## pracma 2.4.3 (2023-08-15)

* Corrected a bug in the help file of nelder_mead().
* Added stereographic projection and its inverse as functions
  stereographic() and stereographic_inv().

## pracma 2.4.2 (2022-09-21)

* Corrected a bug in logit() - 'length_1 error' again.
* Replaced 'class() =' with 'inherits(...)' in function ppval().

## pracma 2.4.1 (2022-08-03)

* Added an example of using inequality constraints with anms().
* Added more vectorized 'distmat2' to the help page of 'distmat'.

## pracma 2.4.0 (2022-04-01)

* hessvec() approximates the multiplication of Hessian and vector
  (thanks to Ravi Varadhan for pointing out Pearlmutter's method).
* hessdiag() computes the diagonal of Hessian (finite differences).

## pracma 2.3.9 (2022-03-22)

* Spherical coordinates better explained in help files.
* Deprecated option 'fast' was removed from circlefit().

## pracma 2.3.8 (2022-02-19)

* 'length_1 eror' for logical operators on integral() help page.
* Bug in pinv() for complex input, reported by Stephane Laurent.

## pracma 2.3.7 (2022-01-17)

* Added lu_crout(), Crout's algorithm for LU matrix decomposition.
* Corrected a small bug in randi(), thanks to Luke Cherveny.

## pracma 2.3.6 (2021-12-03)

* Corrected randortho() with return value 'q %*% diag(ph) # %*% q',
  as researched and reported by Daniel Kessler; thanks a lot.
* Function deeve() requires the x-coordinates to be sorted.
* Removed the Nile data (without warning).

## pracma 2.3.5 (2021-07-10)

* Corrected a bug in subspace() - thanks to David Fleischer.
* Corrected a small but annoying bug in steep_descent().
* Allows for complex matrices in pinv() (and mldivide()).

## pracma 2.3.4 (2021-03-18)

* Small correction in 'movavg()': default type now is 's'.

## pracma 2.3.3 (2021-01-22)

* Forgot to correct URL addresses in the Readme.md file.

## pracma 2.3.2 (2021-01-16)

* Removed and corrected URL addresses from 'http' to 'https'.

## pracma 2.3.1 (2021-01-13)

* 'ellipke' help page: compute circumference of an ellipse.
* 'Mode()' now handling all types of NAs (thx. Michael Henry).

## pracma 2.3.0 (2020-04-09)

* circlefit(): option 'fast' is deprecated and will not be used.
* gammainc(0, a) returns 0, thanks to Mark Chappell for reporting.
* ndims() now returns 1 for vectors and 0 for empty objects.

## pracma 2.2.9 (2019-12-15)

* Changed URL reference of Abramowitz and Stegun (link missing).
* Fixed warning with the 'try' construct in several functions
  using 'if(inherits(e, "try-error"))' (help from Bert Gunter).
* Link to R Base HTML help page gave a warning (for Windows).

## pracma 2.2.8 (2019-07-09)

* erfi() returns real values when the input values are all real.
* hypot() now allows for scalar plus numeric vector as inputs.

## pracma 2.2.7 (2019-05-21)

* Polynomial division with polydiv(); for two plynomials
  polygcf() finds the greatest common factor; and rootmult()
  returns the multiplicity of a polynomial root (or 0).
* polyroots() refines the result of roots() in case of roots
  with multiplicities (where roots() is quite inaccurate).

## pracma 2.2.6 (2019-05-02)

* All polynomial functions now accept complex coefficients, 
  esp. roots() finds roots for complex polynomials.
* Fixed a bug in laguerre() for zeros of complex polynomials.

## pracma 2.2.5 (2019-04-08)

* fsolve() and broyden() are no longer applicable to univariate 
  functions (Morrison-Sherman formula not working in this case).
* Alias cgmin() and option 'dfree=F' in fminsearch() are removed,
  both have been deprecated since half a year.

## pracma 2.2.4 (2018-12-12)

* qpsolve() minimizes quadratic forms such as 0.5 t(x) x - d x
  with linear quality and inequality constraints.
* fmincon() now has an 'augmented Lagrangian' option with a
  'variable metric' approach as inner solver. 

## pracma 2.2.3 (2018-12-10)

* linearproj() linear projection onto a linear subspace, and
  affineproj() linear projection onto an affine subspace of R^n.

## pracma 2.2.2 (2018-11-30)

* Corrected "length > 1 in coercion to logical" in expm().

## pracma 2.2.1 (2018-11-30)

* fminunc() unconstrained minimization of nonlinear objective
  function, based on stripped-down 'Rvmmin' code by John Nash.
* 'fmincon()' minimization of nonlinear objective function with
  constraints; wraps suggested package NlcOptim with SQP method.

## pracma 2.2.0 (2018-11-27)

* Reintroduced 'nelder_mead()' and 'hooke_jeeves()'.
* fminsearch() now calls 'Nelder-Mead' or 'Hooke-Jeeves',
  i.e., derivative-free methods only; 'dfree=F' gets deprecated.

## pracma 2.1.9 (2018-11-22)

* Renamed 'cgmin()' to its original name 'fletcher_powell',
  alias 'cgmin' is deprecated since this version.
* Removed alias 'normest2' that was anyway non-existing.

## pracma 2.1.8 (2018-10-16)

* Corrected a bug in hessenberg() reported by Ben Ubah.

## pracma 2.1.7 (2018-09-24)

* Removed the deprecated 'rortho' function, use randortho() instead.

## pracma 2.1.6 (2018-09-09)

* Si(), Ci() sine and cosine integral functions added.
* Added dot notation for brent(), bisect(), newton(), halley(), and
  ridders() -- on request of John Nash for the histRalg project.

## pracma 2.1.4 (2018-01-29)

* shubert() implements one-dimensional Shubert-Piyavskii method.
* fminsearch() and anms() stop for one-dimensional minimization.

## pracma 2.1.3 (2018-01-23)

* bsxfun() now uses sweep() for matrices in search of higher speed.
* direct1d() removed because slow and not effective.

## pracma 2.1.2 (2018-01-21)

* poisson2disk() approximate Poisson disk distribution
* Corrected small bug in findpeaks(), reported by Mike Badescu.

## pracma 2.1.1 (2017-11-21)

* Added a field "Authors@R" in the DESCRIPTION, deleted others.
* Added README.md and NEWS.md (for a future Github repository).
* Needed a new version for resubmitting (because of 'survivalsvm')

## pracma 2.1.0 (2017-11-20)

* Package 'quadprog' is now suggested, not imported; the functions 
  quadprog() and lsqlincon() work only when 'quadprog' is installed.

## pracma 2.0.9 (2017-09-21)

* Package byte-compiled on loading (Requires R version >= 3.1.0).

## pracma 2.0.8 (2017-09-20)

* findpeaks() function not checking for NAs (reported by Wesley Burr).
* fplot() extra parameters were not handed over to plotting routine.

## pracma 2.0.7 (2017-06-17)

* bernstein() generates the Bernstein polynomial B_,_().
* legendre(n,_) corrected for n=0, thanks to Peter W. Marcy.
* cgmin() alias for fletcher_powell(), a constraint gradient method.

## pracma 2.0.6 (2017-06-06)

* polyvalm() evaluates a polynomial in the matrix sense.
* arnoldi() Arnoldi iteration (incl. Hessenberg matrix).

## pracma 2.0.5 (2017-04-30)

* integral() redesigned, less methods, several starting intervals
  with regular or random intermediate nodes (similar to MATLAB).
* quadgr() corrected as functions vectorized with Vectorize() did
  not behave as expected with apply(); still needs vectorization.
* Help page of quadgk() did not mention the need for vectorization.

## pracma 2.0.4 (2017-04-01)

* hessenberg() computes the Hessenberg form of a matrix through
  Householder transformations (this is named hess() in MATLAB).

## pracma 2.0.3 (2017-03-23)

* Corrected functions with conditions in control statements with
  conditions of length greater than one: rem().

## pracma 2.0.2 (2017-02-23)

* isposdef() test for positive definiteness of a (real) matrix.
* hooke_jeeves() removed; similar implementations are available
  in packages 'dfoptim::hjk[b]' and 'adagio::hookejeeves'.

## pracma 2.0.1 (2017-02-06)

* nelder_mead() replaced by an adaptive Nelder-Mead implementation,
  anms(), following F. Gao and L. Han.
* fminsearch() now calls this new version of Nelder-Mead.

## pracma 2.0.0 (2017-01-26)

* incgam(x,a) computes the incomplete upper gamma function using
  the R function pgamma for higher precision than gammainc().
* Corrected a small oversight in hurstexp(), thnx George Ostrouchov.

## pracma 1.9.9 (2017-01-10)

* Slightly changed the description lines on request of CRAN.

## pracma 1.9.8 (2017-01-10)

* whittaker() finally implemented avoiding the sparse matrix package.
* nelder_mead() now applies adaptive parameters for the simplicial
  search, depending on the dimension of the problem space.
* psinc(x,n), the so-called periodic sinc function.

## pracma 1.9.7 (2016-12-14)

* shooting() implements the shooting method for boundary value problems
  of second order differential equations.
* interp2() corrected the help page with size(z) = length(y) length(x).
* Corrected a small oversight on the help page of Gauss-Laguerre.

## pracma 1.9.6 (2016-09-11)

* haversine() Haversine formula for geographical distances on earth.
* trigonometric functions accepting degrees instead of inputs in radians:
  sind   cosd   tand   cotd   secd   cscd
  asind  acosd  atand  acotd  asecd  acscd  atan2d

## pracma 1.9.5 (2016-09-06)

* fprintf() mimicks MATLAB's function of the same name.
* Added ezsurf(), an easy surface plot following MATLAB.
* fplot() is almost an alias for ezplot(); please note that in future
  versions ez...() will be renamed to f...() according MATLAB 2016/17.

## pracma 1.9.4 (2016-07-27)

* rortho() renamed to randortho(), the underlying code was buggy
  (not truely random) and has been replaced, thanks to Jan Tuitman.
* an error in the final step of calculating approx_energy() was
  corrected, thanks to Daniel Krefl.

## pracma 1.9.3 (2016-05-28)

* bvp() now solves boundary value problems for linear 2nd order ODEs
  using a 'finite differences' approach and a tridiagonal solver.
* polyfit2() has been removed, use polyfix() instead.

## pracma 1.9.2 (2016-03-04)

* romberg() corrected an error estimation that diminished the accuracy.
* trapzfun() realizes trapezoidal integration with iterated calculations.

## pracma 1.9.1 (2016-02-15)

* fractalcurve() generates some fractal curves of order n, i.e. the
  Hilbert, Sierpinski, Snowflake, Dragon, and Molecule curves.
* ode23(), ode23s() changed the size of the returned components,
  now it is similar to what is returned by ode45() and ode78().
* arclength() corrected a boundary condition ('on the left'), added
  an example how to generate an arc-length parametrization of a curve.

## pracma 1.9.0 (2015-12-17)

* quadprog() solves quadratic programming problems (QP) with linear
  equality and inequality constraints, based on package 'quadprog'.
* lsqlincon() solves linear least-squares problems with linear equality 
  and inequality constraints (as well as bound constraints).
* pracma now imports package 'quadprog'.

## pracma 1.8.9 (2015-12-05)

* polyfix() fits a polynomial that exactly passes through given
  fixed points. polyfit2() will be deprecated in future versions.
* Important bug fix for polyApprox() (thanks to Max Marchi).

## pracma 1.8.8 (2015-10-28)

* Option 'minpeakdistance' for function findpeaks() added
  (thanks to Razvan Chereji for providing a workable approach).

## pracma 1.8.7 (2015-07-20)

* Removed invperm().
* 'linear' is now the default method for interp1().
* Cases n = 0, 1 for legendre() corrected (thanks to Nuzhdin Yury).

## pracma 1.8.6 (2015-07-11)

* Removed two non-existing links pointing to Gander's pages at the ETHZ.
* Removed a link explaining approximate entropy.

## pracma 1.8.5 (2015-07-07)

* Added 'Imports' field in description and 'import' in namespace,
  as requested for the new R development version.
* strrep() renamed to strRep(), because of a new function in R Base.

## pracma 1.8.4 (2015-06-25)

* bernoulli() calculates the Bernoulli numbers and polynomials.
* factorial2() the product of all even resp. odd integers below n.

## pracma 1.8.3 (2015-02-08)

* Deleted some URLs that were not working properly anymore.

## pracma 1.8.2 (2015-02-07)

* Special functions gathered under topics 'specfun' resp. 'specmat'.

## pracma 1.8.1 (2015-02-06)

* sumalt() accelerating (infinite) alternating sums.
* Option 'fast=FALSE' in circlefit() to avoid optim().
* Added Gauss' AGM-based computation of pi to agmean(). 

## pracma 1.8.0 (2015-01-26)

* hurstexp() amended for vectors of uneven length.

## pracma 1.7.9 (2014-11-15)

* qpspecial() special quadratic programming solver.
* Reintroduces the 'tol' keyword in fminbnd() for compatibility.

## pracma 1.7.8 (2014-11-10)

* bulirsch_stoer() Bulirsch-Stoer method for solving
  ordinary differential equations with high accuracy.
* midpoint() implements the midpoint rule for solving ODEs
  combined with Richardson extrapolation for high accuracy.

## pracma 1.7.7 (2014-11-01)

* lufact() LU factorization with partial pivoting;
  lusys() solves linear systems through Gaussian elimination.

## pracma 1.7.6 (2014-10-30)

* ode23s() for stiff ordinary differential equations refining
  Rosenbrock's method (supply Jacobian if available).
* euler_heun() Euler-Heun ODE solver has been corrected.

## pracma 1.7.5 (2014-10-20)

* fminbnd() much improved implementation of Brent's method;
  added challenging example by Trefethen to the help page.
* lambertWn() for the second (real) branch of Lambert W.
* Function name alias cintegral() removed.

## pracma 1.7.4 (2014-10-13)

* hooke_jeeves() replaced by a much more efficient implementation
  and equipped with a special approach to bound constraints.
* nelder_mead() replaced by a much more efficient implementation
  and utilizing a transformation to handle bound constraints;
  functions nelmin() and nelminb() are not exported anymore.

## pracma 1.7.3 (2014-10-11)

* quadinf() now uses the double exponential method with the
  tanh-sinh quadrature scheme for (semi-)infinite intervals.
* Removed the not-exported and too slow .quadcc() function.
* brent() alias for brentDekker(), newton() for newtonRaphson().

## pracma 1.7.2 (2014-09-08)

* pchipfun() function wrapper around pchip();
  missing error handling in pchip() was added.
* hurst() removed, functionality merged with hurstexp().
* Nile overflow data set 1871--1984 added as time series.

## pracma 1.7.1 (2014-08-12)

* bits() binary representation of a number as string.
* agmean() returns AGM, no of iterations, and estimated precision.
* trapz() tiny improvement on error handling.

## pracma 1.7.0 (2014-06-30)

* ode45() ODE solver using Dormand-Prince (4,5) coefficients.
* ode78() ODE solver using Fehlberg (7,8) coefficients.
* cintegral() renamed to line_integral().

## pracma 1.6.9 (2014-06-14)

* Version 1.6.8 "Failed to build" on R-Forge.
  [Maybe it's time to move pracma to a github repository.]

## pracma 1.6.8 (2014-06-07)

* nelmin() a more efficient and accurate version of Nelder-Mead.
* nelminb() Nelder-Mead in bounded regions (applies a transformation).

## pracma 1.6.7 (2014-05-23)

* trisolve() stopping for singular tridiagonal matrices.
* romberg() slightly improved accuracy and speed.

## pracma 1.6.6 (2014-04-12)

* Corrected rref() (as pointed out by Peter Audano).

## pracma 1.6.5 (2014-02-24)

* lsqnonneg() changed to an active-set approach.
* bisect() trimmed bisection to return almost exact results.

## pracma 1.6.4 (2014-02-05)

* halley() Halley's variant of the Newton-Raphson method.
* numderiv() corrected Richardson's method by breaking the loop.

## pracma 1.6.3 (2014-01-25)

* lambertWp() improved inner accuracy from 1e-12 to 1e-15.
* complexstepJ() renamed to jacobian_csd(); introduced grad_csd().
* hessian_csd() applies Richardson's method as the second step,
  and the same for laplacian_csd().

## pracma 1.6.2 (2014-01-19)

* Removed zeroin(); for fzero() a variation of Brent-Dekker is used,
  that applies cubic instead of quadratic interpolation.
* Corrected an oversight in newtonRaphson().
* brentDekker() returns a list now.

## pracma 1.6.1 (2014-01-14)

* samp_entropy() complements approx_entropy() for short time series.
* Removed NEWS.Rd and NEWS.pdf in favour of NEWS.

## pracma 1.6.0 (2013-12-06)

* integral3() now handles functions as inner interval limits.
* poly_crossings() calculates crossing points of two polygons.
* erfz() complex error function vectorized (thanks to Michael Lachmann).

## pracma 1.5.9 (2013-11-30)

* muller() implements Muller's root-finding method [Mueller, 1956],
  especially suited for polynomials and complex functions.
* Inserted a safeguard for the distmat() function to prevent different
  results on Mac OS X, (Ubuntu) Linux, and Windows operating systems.
* Removed pltcross() and kmeanspp().

## pracma 1.5.8 (2013-11-28)

* interp1() with option method ``spline'' now computes Moler's spline
  functions, for compatibility with MATLAB (hint by Boudewijn Klijn).

## pracma 1.5.7 (2013-10-11)

* Corrected parameter 'waypoints' in cintegral().

## pracma 1.5.6 (2013-09-22)

* odregress() orthogonal distance (or: total least-squares) regression.
* Changed maintainer name to its long form (CRAN request).

## pracma 1.5.5 (2013-09-11)

* L1linreg() L1 (a.k.a. LAD or median) linear regression.
* geo_median() geometric median (minimizes sum of distances).

## pracma 1.5.4 (2013-08-31)

* rectint() rectangular intersection areas (MATLAB style).
* cumtrapz() cumulative trapezoidal integration (MATLAB style).
* Some corrections to help pages and function names.

## pracma 1.5.3 (2013-08-25)

* arclength() length of a parametrized curve in n-dimensional space,
  w/ improved convergence by applying Richardson's extrapolation method.
* legendre() associated Legendre functions (MATLAB style).

## pracma 1.5.2 (2013-08-23)

* poly_center() calculates the center coordinates of a polygon.
* poly_length() calculates the (euclidean) length of a polygon.
* polyarea() corrected, returns the true, not the absolute value.

## pracma 1.5.1 (2013-08-19)

* fsolve() will use broyden() if m = n; fzsolve() the same;
  additionally, improved broyden() and gaussNewton().
* ezplot() can draw markers on the line, with equal distances
  measured along the curve length.

## pracma 1.5.0 (2013-08-08)

* gmres() generalized minimum residual method.
* nearest_spd() finds nearest symmetric positive-definite matrix.
* eps() floating point relative accuracy.

## pracma 1.4.9 (2013-07-16)

* lapacian() now works in n dimensions, not only for n = 2.
* mldivide(), mrdivide() corrected a severe typo.
* numderiv(), numdiff() start with h = 1/2 instead of h = 1.
* figure() platform-independent by using dev.new().

## pracma 1.4.8 (2013-06-17)

* findzeros() now finds 'quadratic' roots, too.
* pdist2() added as an alias for distmat(),
  while pdist(X) now is distmat(X, X) (MATLAB style).

## pracma 1.4.7 (2013-05-20)

* histcc() histogram with optimized number of bins.
* Example of correction term for the trapz() integration.

## pracma 1.4.6 (2013-03-31)

* psi() Psi polygamma function (MATLAB style).
* rosenbrock() and rastrigin() functions removed.

## pracma 1.4.5 (2013-03-21)

* quadcc() new, iterative Clenshaw-Curtis quadrature.
* squareform() formats distance matrix (MATLAB style).

## pracma 1.4.4 (2013-03-12)

* integral2() implements the two-dimensional numerical integration
  approach `TwoD', i.e. Gauss-Kronrod (3, 7)-points on rectangles.
* integral3() three-dimensional integration based on integral2().
* triplequad() 3-dim. integration based on dblquad() (MATLAB style).

## pracma 1.4.3 (2013-03-10)

* integral() combines adaptive numerical integration procedures.
* cintegral() complex line integrals (rectangles and curves).

## pracma 1.4.2 (2013-03-03)

* linprog() linear programming solver for linear equality and
  inequality constraints.

## pracma 1.4.1 (2013-02-20)

* romberg() Romberg integration completely rewritten.
* idivide() integer division with different roundings.

## pracma 1.4.0 (2013-02-10)

* fderiv(), taylor() expanded to higher orders.
* itersolve() iteration methods for solving linear systems.
* lu() LU decomposition with different schemes (w/o pivoting).

## pracma 1.3.9 (2013-01-26)

* pdist() as an alias for distmat() (MATLAB style).
* fftshift(), ifftshift() shifting Fourier frequencies.
* Improved grad(), jacobian(), hessian(), and laplacian().

## pracma 1.3.8 (2013-01-10)

* Smaller corrections, e.g., removed deprecated 'is.real';
  no startup messages anymore.
* geomean(), harmmean(), trimmean() geometric, harmonic, and
  trimmed arithmetic mean (MATLAB style).
* agmean() algebraic-geometric mean.

## pracma 1.3.7 (2013-01-07)

* mexpfit() multi-exponentiell fitting, separating linear and
  nonlinear parts of the problem.

## pracma 1.3.6 (2013-01-06)

* lsqsep() separable least-squares fitting.
* lsqcurvefit() nonlinear least-squares curve fitting.

## pracma 1.3.5 (2013-01-05)

* cd(), pwd() directory functions (MATLAB style).
* rand(), randn() changed to accept size() as input.
* whos(), what() corrected for empty lists resp. directories.

## pracma 1.3.4 (2012-12-19)

* what(), who(), whos(), ver() (MATLAB style).
* semilogx(), semilogy(), loglog() logarithmic plots (MATLAB style)

## pracma 1.3.3 (2012-12-12)

* quadv() vectorized integration.
* ezpolar() easy access to the polar() function.
* sortrows() sorting rows of matrices (MATLAB style).
* null() alias for nullspace function (MATLAB style).
* eigjacobi() Jacobi's method for eigenvalues and eigenvectors.

## pracma 1.3.2 (2012-12-08)

* ellipke(), ellipj() elliptic and Jacobi elliptic integrals.
* expint() implements E1 and Ei, the exponential integrals, with
  aliases expint_E1() and expint_Ei().
* li() the logarithmic integral (w/o offset).

## pracma 1.3.1 (2012-12-06)

* Explicitely listing about 200 MATLAB-emulating function( name)s.
* Dismissed matlab(), using it now for infos only, not assigning any
  MATLAB function names to the environment (because of CRAN policies).

## pracma 1.3.0 (2012-12-05)

* cot(), csc(), sec() cotangens, cosecans, and secans functions.
* acot(), acsc(), asec() inverse cotangens, cosecans, secans functions.
* coth(), csch(), sech() hyperbolic cotangens, cosecans, secans functions.
* acoth(), acsch(), asech() inverse hyperbolic cotangens, cosecans, and
  secans functions.

## pracma 1.2.9 (2012-12-02)

* bvp() changed to solve second order boundary value problems.
* trisolve() solves tridiagonal linear equation systems.
* curvefit() fits points in the plane with a polynomial curve.

## pracma 1.2.8 (2012-11-30)

* lsqlin() least-squares solver with linear equality constraints.
* pinv() now works like MASS::ginv() for singular matrices.
* Added the end-';' feature to str2num().
* toc() added invisible return value.

## pracma 1.2.7 (2012-11-22)

* procrustes() solving the Procrustes problem,
  and kabsch() implements the Kabsch algorithm.
* kriging() ordinary and simple Kriging interpolation.
* Corrected some stupid errors in str2num().

## pracma 1.2.6 (2012-11-11)

* akimaInterp() univariate Akima interpolation.
* Moved transfinite() to package 'adagio'.

## pracma 1.2.5 (2012-10-28)

* histc() Histogram-like counting (MATLAB style).
* Added warning to complexstep() if imaginary part is zero.

## pracma 1.2.4 (2012-10-25)

* Added option 'pinv' to mldivide() to return same results as MATLAB.
* str2num(), num2str() conversion functions (MATLAB style).
* Removed some 'author' entries on help pages.

## pracma 1.2.3 (2012-10-17)

* Renamed mrank() to Rank().
* Corrected nullspace() [thanks to Stephane Laurent], which now agrees
  with Octave's null() function (MASS:Null appears buggy, too).
* Corrected gaussNewton() and fsolve() [thanks to Etienne Chamayou].

## pracma 1.2.2 (2012-10-10)

* bsxfun() apply binary function elementwise (MATLAB style).
* added the analytic solution for the example in bvp().

## pracma 1.2.1 (2012-09-28)

* rosenbrock() added, moved testfunctions to 'adagio' package.
* euler_heun() improved Euler method for solving ODEs.
* logit() function added to sigmoid().
* Keyword 'ode' introduced.

## pracma 1.2.0 (2012-09-27)

* matlab() can reinstall MATLAB function names.

## pracma 1.1.9 (2012-09-25)

* gcd(), lcm() greatest common divisor, least common multiple
  now working on a vector of integers.
* Removed number-theoretic functions: eulersPhi(),
  moebiusFun(), mertensFun(), sigma(), tau(), omega(), Omega(),
  primes2(), twinPrimes(), nextPrime(), previousPrime(),
  modpower(), modorder(), modinv(), modlin(),
  primroot(), contfrac(), coprime(), GCD(), LCM(), extGCD(),
  (these functions are now available in the 'numbers' package).

## pracma 1.1.8 (2012-09-19)

* ezcontour(), ezmesh() wrappers for contour(), image(), persp().
* erfi() imaginary error function.

## pracma 1.1.7 (2012-08-06)

* moler() Moler matrix

## pracma 1.1.6 (2012-07-20)

* Removed '.Rapphistory' from the tests directory (again)
  [and use "--as-cran" for the checks].
* disp() display text or array (MATLAB Style), i.e. cat() with newline.

## pracma 1.1.5 (2012-07-18)

* Renamed functions with capital first letter to avoid name clashes:
  mtrace -> Trace, mdiag -> Diag, strtrim -> strTrim, vnorm -> Norm,
  reshape -> Reshape, find -> finds, fix -> Fix, poly ->Poly,
  mode -> Mode, real -> Real, imag -> Imag, toeplitz -> Toeplitz.

## pracma 1.1.4 (2012-06-26)

* gammainc() (lower and upper) incomplete gamma function, also the
  regularized gamma function, all allowing negative x values.
* polylog() the polylogarithm functions for |z| < 1 and n >= -4 .

## pracma 1.1.3 (2012-06-17)

* fminsearch() now implements Nelder-Mead (similar to optim), and
  Fletcher-Powell when "dfree=FALSE" is chosen.
* Test functions rosenbrock() and rastrigin().

## pracma 1.1.2 (2012-06-13)

* nelder_mead() implements Nelder-Mead for nonlinear optimization.
* hooke-jeeves() Hooke-Jeeves algorithm for direct search.
* fletcher_powell() Davidon-Fletcher-Powell method for function
  minimization (alternative to BFGS approach).
* steep_descent() minimization of functions using steepest descent.

## pracma 1.1.1 (2012-06-10)

* fminbnd() now implements Brent's function minimization algorithm with
  golden section search and parabolic interpolation (same as optimize).
* transfinite() transformation function between bounded and unbounded
  (box constraint) regions.

## pracma 1.1.0 (2012-06-06)

* hurst(), hurstexp() calculate the Hurst exponent of a time series.
* Updated the NEWS.Rd file.

## pracma 1.0.9 (2012-06-03)

* lsqnonneg() solves nonnegative least-squares problems by using the
  trick "x --> exp(x)" and applying lsqnonlin();
  example function lsqcurvefit() for nonlinear curve fitting.
* Renamed ridder() to ridders(), thanks to Robert Monfera for pointing
  it out (he also suggested a multi-dimensional variant).

## pracma 1.0.8 (2012-05-22)

* movavg() moving average of types "simple", "weighted", "modified",
  "exponential" (EMA), or "triangular".
* modlin() solves modular linear equations.

## pracma 1.0.7 (2012-05-11)

* lsqnonlin() solves nonlinear least-squares problems using the
  Levenberg-Marquardt approach.
* renamed froots() to findzeros(), and fmins() to findmins().

## pracma 1.0.6 (2012-04-21)

* fornberg() finite difference (i.e., polynomial) approximation of
  derivatives for unevenly spaced grid points -- Fornberg's method.

## pracma 1.0.5 (2012-04-15)

* randsample() randomly sampling, alias for sample (MATLAB style).
* rands() generates uniform random points on an N-sphere.
* Added tic(), toc() measuring elapsed time (MATLAB style).
* previousPrime() finds the next prime below a number.

## pracma 1.0.4 (2012-04-01)

* invlap() computes the inverse Lapacian numerically.
* ppfit() piecewise polynomial fitting procedure.

## pracma 1.0.3 (2012-03-21)

* cubicspline() interpolating cubic spline (w/ endpoint conditions).
* mkpp() and ppval() for piecewise polynomial structures.

## pracma 1.0.2 (2012-03-17)

* accumarray() resembles the related MATLAB function more closely.
* invperm() returns the inverse of a permutation.
* randperm() changed to make it more MATLAB-like.

## pracma 1.0.1 (2012-03-09)

* plotyy() corrected right ordinate, prettying the labels.
* peaks() peaks function (MATLAB style).

## pracma 1.0.0 (2012-03-01)

* Updated the NEWS.Rd file.

## pracma 0.9.9 (2012-02-29)

* qrSolve solves overdetermined system of linear equations.
* DSCsearch() removed, now in package 'pracopt'.
* randp() found a better, non-selective approach.

## pracma 0.9.8 (2012-02-23)

* gramSchmidt() modified Gram-Schmidt process.
* householder() Householder reflections and QR decomposition.
* givens() Givens rotation and QR decomposition.
* corrected a small error in ridder() (thanks to Roger Harbord); new
  example of how to use ridder() with Rmpfr for multiple precision.

## pracma 0.9.7 (2012-02-17)

* erf() corrected, erfc() and erfcx() as new functions,
  including their inverses erfinv() and erfcinv().
* hypot() now numerically more stable (thanks to Jerry Lewis).

## pracma 0.9.6 (2012-01-25)

* Changed third example for dblquad() [new Windows toolchain problem].
* Deactivated the test for gammaz() because of problems on Solaris.

## pracma 0.9.5 (2012-01-16)

* kmeanspp() kmeans++ clustering algorithm.
* hampel() with new option, fuelled by a blog entry of Ron Pearson.

## pracma 0.9.4 (2012-01-08)

* DSCsearch() Davies-Swann-Campey search in one dimension.
* Improved modpower() through modular exponentiation.
  added lehmann_test() Lehmann's primality test as example.
* Corrected polar() and andrewsplot().

## pracma 0.9.3 (2011-12-27)

* direct1d() one-dimensional version of the DIRECT algorithm for
  global function minimization.

## pracma 0.9.2 (2011-12-26)

* approx_entropy() approximate entropy of a time series.
* circshift() circularly shifting arrays (MATLAB Style).

## pracma 0.9.1 (2011-12-12)

* plotyy() plots curves with y-axes on both left and right side.
* fplot() plots components of a multivariate function.

## pracma 0.9.0 (2011-12-11)

* errorbar() routine for plotting error bars in both directions.
* whittaker() Whittaker-Henderson smoothing** Not yet running** .
* rref() reduced row echelon form.

## pracma 0.8.9 (2011-12-08)

* cutpoints() automatically finds cutting points based on gaps.
* hausdorff_dist calculates the Hausdorff distance / Hausdorff dimension.
* nnz() number of non-zeros elements (MATLAB style).

## pracma 0.8.8 (2011-12-06)

* polar() for polar plots (MATLAB style), see the example plots.
* andrewsplot() plots Andrews curves in polar coordinates.
* Vectorized: cart2sph(), sph2cart(), cart2pol(), pol2cart().

## pracma 0.8.7 (2011-11-30)

* deg2rad(), rad2deg().
* figure() MATLAB style and pltcross() plotting crosses.

## pracma 0.8.6 (2011-11-21)

* ridder() Ridder's method for zero finding of univariate functions.

## pracma 0.8.5 (2011-11-19)

* sqrtm() matrix square root, based on Denman-Beavers iteration,
  rootm() matrix p-th root, computing a complex contour integral,
  signm() matrix sign function.
	o fzero() now uses the new zeroin() function,
	  i.e., a Brent-Dekker approach instead of refering to uniroot().
* twinPrimes() twin primes in a given interval, and nextPrime()
  will find the next higher prime.

## pracma 0.8.4 (2011-11-14)

* Transformations between cartesian, spherical, polar and cylindrical
  coordinate systems: cart2sph(), sph2cart(), cart2pol(), pol2cart().
* randp() uniformly random points in the unit circle.

## pracma 0.8.3 (2011-11-11)

* accumarray() grouping elements and applying a function to each group.
* uniq() MATLAB-style 'unique' function, allsums() in the examples.
* small correction to fsolve(), mentioned on the 'check summary' page.

## pracma 0.8.2 (2011-11-04)

* newmark() Newmark's method for solving second order differential
  equations of the form  y''(t) = f(t, y(t), y'(t)) on [t1, t2].
* cranknic() Crank-Nicolson 'ivp' solver, combining the forward and
  backward Euler methods for ordinary differential equations.

## pracma 0.8.1 (2011-10-30)

* Corrected pinv() for (nearly) singular matrices.
* Renamed ifactor() to factors().

## pracma 0.8.0 (2011-10-27)

* Minor corrections and improvements to the 'pracma.pdf' manual,
  incl. numdiff(), refindall(), trigApprox(), and subspace().

## pracma 0.7.9 (2011-10-22)

* spinterp() monotonic (and later on shape-preserving) interpolation
  following the approach of Delbourgo and Gregory.

## pracma 0.7.8 (2011-10-17)

* bvp() solves boundary value problems of the following kind:
  -u''(x) + c1 u'(x) + c2 u(x) = f(x) for x in [a, b].

## pracma 0.7.7 (2011-10-14)

* primes2(n1, n2) will return all prime numbers betweeen n1 and n2
  (without storing the numbers from sqrt(n2) up to n2).

## pracma 0.7.6 (2011-08-05)

* gaussNewton() for function minimization and solving systems of
  nonlinear equations. fsolve() as a wrapper for it.
* fzsolve() for root finding of complex functions.
* softline() Fletcher's inexact linesearch algorithm.

## pracma 0.7.5 (2011-07-26)

* Put NEWS.Rd in the /inst subdirectory (and NEWS.pdf in /doc),
  thanks to Kurt Hornik; slightly changed the version numbering.

## pracma 0.7-4 (2011-07-22)

* rortho() generate random orthogonal matrix of size n.
* Titanium data set for testing fitting procedures.

## pracma 0.7-3 (2011-07-15)

* erf() and erfc() error and complementary error functions
  (MATLAB style) as (almost) aliases for pnorm().
* erfz() complex error function.

## pracma 0.7-2 (2011-07-11)

* broyden() quasi-Newton root finding method for systems of nonlinear
  equations.

## pracma 0.7-1 (2011-07-09)

* cross() has been vectorized (remark on R-help).

## pracma 0.7-0 (2011-07-07)

* Sigmoid and Einstein functions.

## pracma 0.6-9 (2011-07-06)

* Runge-Kutta-Fehlberg method of order (5,4).

## pracma 0.6-8 (2011-07-05)

* triquad() Gaussian quadrature over triangles.
* cotes() Newton-Cotes integration formulae for 2 to 8 nodes.

## pracma 0.6-7 (2011-07-04)

* lagrangeInterp(), newtonInterp() Lagrange and Newton polynomial
  interpolation, neville() Neville's methods.
* tril(), triu() extracting triangular matrices (MATLAB style).

## pracma 0.6-6 (2011-07-02)

* charpoly() computes the characteristic polynomial, the determinant,
  and the inverse for matrices that are relativly small, applying the
  Faddejew-Leverrier method.
* froots() to find *all* roots (also of second or higher order) of
  a univariate function in a given interval. The same with fmins()
  to find *all* minima.

## pracma 0.6-5 (2011-07-01)

* Adams-Bashford and Adams-Moulton (i.e., multi-step) methods
  for ordinary differential equations in function abm3pc().

## pracma 0.6-4 (2011-06-30)

* Changed the description to be more precise about the package.

## pracma 0.6-3 (2011-06-28)

* rationalfit() rational function approximation
* ratinterp() rational interpolation a la Burlisch-Stoer.

## pracma 0.6-2 (2011-06-26)

* pade() Pade approximation.

## pracma 0.6-1 (2011-06-25)

* quadgk() adaptive Gauss-Kronrod quadrature.

## pracma 0.6-0 (2011-06-24)

* Added differential equation example to expm()'s help page.
* Changed NEWS file to become simpler (no subsections).

## pracma 0.5-9 (2011-06-23)

* quadl() recursive adaptive Gauss-Lobatto quadrature.
* simpadpt() another recursively adaptive Simpson's rule.
* Added testing procedures for all integration routines;
  corrected, refined some of these procedures.

## pracma 0.5-8 (2011-06-20)

* quadgr() Gaussian Quadrature with Richardson extrapolation, can
  handle singularities at endpoints and (half-)infinite intervals.

## pracma 0.5-7 (2011-06-18)

* expm() for matrix exponentials.
* clenshaw_curtis() the Clenshaw-Curtis quadrature formula.

## pracma 0.5-6 (2011-06-17)

* simpson2d() as non-adaptive 2-dimensional Simpson integration.
* dblquad() twofold application of internal function integrate().

## pracma 0.5-5 (2011-06-15)

* gaussHermite() and gaussLaguerre() for infinite intervals.
* Fresnel integrals fresnelS() and frenelC().

## pracma 0.5-4 (2011-06-12)

* gaussLegendre() computes coefficients for Gauss Quadrature,
  and quad2d() uses these weights for 2-dimensional integration.
* quadinf() wrapper for integrate() on infinite intervals.

## pracma 0.5-3 (2011-06-06)

* ode23() solving first order (systems of) differential equations.
* barylag2d() 2-dimensional barycentric Lagrange interpolation.

## pracma 0.5-2 (2011-06-04)

* interp2() for two-dimensional interpolation.
* gradient() now works in two dimensions too.

## pracma 0.5-1 (2011-06-01)

* fzero(), fminbnd(), fminsearch(), fsolve() as aliases for
  uniroot(), optimize(), optim() with Nelder-Mead, newtonsys().

## pracma 0.5-0 (2011-05-31)

* Corrections to help pages.

## pracma 0.4-9 (2011-05-30)

* romberg() and gauss_kronrod() for numerical integration.
* Richardson's extrapolation in numderiv(), numdiff().
* Discrete numerical derivatives (one dimension): gradient().

## pracma 0.4-8 (2011-05-28)

* Numerical function derivatives: fderiv(), grad().
* Specialized operators: hessian(), laplacian().
* Application: taylor().

## pracma 0.4-7 (2011-05-27)

* plot vector fields: quiver() and vectorfield().
* findintervals().
* Corrections in deval(), deeve(), using findintervals().

## pracma 0.4-6 (2011-05-26)

* Laguerre's method laguerre().
* rk4() and rk4sys() classical fourth order Runge-Kutta.
* deval(), deeve() evaluate ODE solutions.

## pracma 0.4-5 (2011-05-24)

* Lebesgue coefficient: lebesgue().
* poly2str() for string representation of a polynomial.

## pracma 0.4-4 (2011-05-23)

* Dirichlet's eta() and Riemann's zeta() function.
* rmserr() different accuracy measures; std_err() standard error.

## pracma 0.4-3 (2011-05-22)

* polypow() and polytrans() for polynomials.
* polyApprox() polynomial approximation using Chebyshev.
* trigPoly(), trigApprox() for trigonometric regression.

## pracma 0.4-2 (2011-05-17)

* segm_intersect() and segm_distance() segment distances.
* inpolygon().

## pracma 0.4-1 (2011-05-13)

* polyadd() polynomial addition.
* conv() and deconv() time series (de)convolution.
* detrend() removes (piecewise) linear trends.
* ifft() for normalized inverse Fast Fourier Transform.

## pracma 0.4-0 (2011-05-10)

* Added tests for functions since version 0.3-7.

## pracma 0.3-9 (2011-05-09)

* and() and or().

## pracma 0.3-8 (2011-05-06)

* pchip() and option `cubic' for interp1() interpolation.
* The complex gamma functions gammaz().
* hadamard() and toeplitz() matrices.

## pracma 0.3-7 (2011-05-04)

* Rank of a matrix, mrank(), and nullspace() for the kernel.
* orth(), orthogonal basis of the image space, and subspace()
  determines the angle between two subspaces.
* normest() for estimating the (Frobenius) norm of a matrix, and
  cond() determines the condition number of a matrix.

## pracma 0.3-6 (2011-04-30)

* fact(), more accurate than the R internal function `factorial'.
* ezplot() as an alias for curve(), but with option ``fill = TRUE''.
* aitken() for accelerating iterations.
* Renamed polycnv() to polymul().
* Renamed outlierMAD() to hampel().

## pracma 0.3-5 (2011-04-23)

* Lambert W function lambertWp() for the real principal branch.
* ``Complex Step'' derivation with complexstep() and complexstepJ().

## pracma 0.3-4 (2011-04-21)

* Barycentric Lagrange interpolation through barylag().
* polyfit2() fits a polynomial that exactly meets one additional point.
* Added more references to the help entry `pracma-package.Rd'.

## pracma 0.3-3 (2011-04-19)

* hornerdefl() for also returning the deflated polynomial.
* newtonHorner() combining Newton's method and the Horner scheme
  for root finding for polynomials.
* jacobian() computes the Jacobian of a function R^n --> R^m as simple
  numerical derivative.
* newtonsys() applies Newton's method to functions R^n --> R^n with
  special application to root finding of complex functions.
* newton() renamed to newtonRaphson().

## pracma 0.3-2 (2011-04-17)

* Sorting functions: bubbleSort(), insertionSort(), selectionSort(),
  shellSort(), heapSort(), mergeSort(), mergeOrdered(), quickSort(),
  quickSortx(), is.sorted(), and testSort().
* Functions from number theory: eulersPhi(), moebiusFun() and the
  mertensFun(), sigma(), tau(), omega(), and Omega().

## pracma 0.3-1 (2011-04-16)

* Chebyshev polynomials of the first kind: chebPoly(), chebCoeff(),
  and chebApprox().

## pracma 0.3-0 (2011-04-09)

* New version of NEWS.Rd, NEWS.pdf.
* More test functions for root finding and quadrature.

## pracma 0.2-9

* fnorm() and the Runge function runge().
* contfrac(), rat(), and rats() for continuous fractions.
* meshgrid() and magic().

## pracma 0.2-8

* quad() adaptive Simpson quadrature.
* Minimum finding with fibsearch() and golden_ratio().
* Root finding with newton(), secant(), and brentDekker().

## pracma 0.2-7

* Regular expression functions regexp(), regexpi(), regexprep() and
  refindall().

## pracma 0.2-6

* String functions blanks(), strtrim(), deblank(), strjust(), strrep().
* interp1() one-dimensional interpolation (incl. spline)

## pracma 0.2-5

* MATLAB functions mode(), clear() and beep().

## pracma 0.2-4

* primroot() finds the smallest primitive root modulo a given n;
  needed functions are modpower() and modorder().
* humps() and sinc(): MATLAB test functions.
* Root finding through bisection: bisect(), regulaFalsi().
* outlierMAD(), findpeaks(), and piecewise().
* polycnv() for polynomial multiplication.
* Functions extgcd(), gcd(), and lcm() have been renamed to extGCD(),
  GCD(), and LCM() respectively.

## pracma 0.2-3

* strfind(), strfindi(), and findstr().
* circlefit() fitting a circle to plane points.
* mldivide() and mrdivide(), emulating the MATLAB backslash operator.

## pracma 0.2-2

* vnorm() vector norm
* Warning about a nasty "non-ASCII input" in the savgol.RD file resolved.

## pracma 0.2-1 (2011-03-17)

* horner() implementing the horner scheme for evaluating a polynomial
  and its derivative.
* savgol() Savitzki-Golay smoothing and needed pseudoinverse pinv().

## pracma 0.2-0

* Package renamed to 'pracma' to avoid name clashes with packages
  such as 'matlab' that are sticking closer to the original.
* Added 'pracma-package' section to the manual.

## pracma 0.1-9 (2011-03-13)

* reshape(), repmat(), and blkdiag() matrix functions.
* combs() chooses all combinations of k elements out of n, and
  randcomb() generates a random selection.
* perms() generates all permutations, randperm() a random permutation.
* Pascal triangle as pascal(); nchoosek() returns binomial coefficients.
* Some string functions: strcmp(), strcmpi(), strcat().

## pracma 0.1-8 (2011-03-10)

* std() as refinement of the standard deviation function.
* ceil() and fix() as aliases for ceiling() and trunc().
  [floor() and round() already exist in R.]
* Modulo functions mod(), rem() and integer division idiv().
* Integer functions related to the Euclidean algorithm:
  extgcd(), gcd(), lcm(), coprime(), and modinv().
* distmat() and crossn(), the vector product in n-dimensional space.

## pracma 0.1-7 (2011-03-08)

* size(), numel(), ndims(), isempty(), and find().
* eye(), ones(), zeros().
* Functions returning random numbers: rand(), randn(), randi().
* linspace(), logspace(), and logseq() for linearly, logarithmically,
  and exponentially spaced sequences.

## pracma 0.1-6 (2011-03-06)

* Matrix functions mdiag() and mtrace() added.
  inv() is introduced as an alias for solve() in R.
* Generate special matrices hankel(), rosser(), and wilkinson().
  kron() is an alias for the R function kronecker().
* Renamed factors() to ifactor() to distinguish it more clearly
  from factors as used in R.

## pracma 0.1-5

* Added function for flipping or rotating numeric and complex
  matrices: flipdim(). flipud(), fliplr(), and rot90().

## pracma 0.1-4

* Added functions for generating sequences of (log-)linearly spaced
  numeric values: linspace() and logspace().
* Added basic complex functions real(), imag(), conj(), and angle()
  which are essentially only aliases of the R functions Re(), Im(),
  or Conj().
  angle() returns the angle of a complex number in radians. The R
  function Mod() is here only available as abs().

## pracma 0.1-3 (2011-02-20)

* Added compan() function for the `companion' matrix;
  the eig() function is an alias for the R eigen()values function.
* Added the polynomial functions poly(), polyder(), polyfit(),
  polyint(), and polyval().
* roots() returns real and complex roots of polynomials.
* Simplified the trapz() function.

## pracma 0.1-2

* Added functions from number theory: primes(), isprime() and factors().
  The corresponding function for factors() in MATLAB/Octave is called
  factor(), but that name should not be shadowed in R!
* Added the polyarea() and trapz() functions.

## pracma 0.1-1

* Added some simple functions such as nthroot(), pow2(), and nextpow2().
* dot() and cross() functions for scalar and vector product.
* Generate matrices through vander() and hilb().

## pracma 0.1-0 (Initial Version)

### INSTALLATION

* 'pracma' will be a pure R package without using any source code.
  Therefore, installation will be immediate on all platforms.

### INTENTION

* This package provides R implementations of more advanced math
  functions from MATLAB and Octave (and the Euler Math Toolbox)
  with a special view on optimization and time series routines.
