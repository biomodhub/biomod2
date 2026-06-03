The R package optimParallel 
===========================

The package provides a parallel versions of the L-BFGS-B optim method.
If the evaluation of the function fn takes more than 0.1 seconds,
optimParallel can significantly reduce the optimization time. For a p-parameter optimization,
the speed increase is about factor 1+2p when no analytic gradient is specified and
1+2p processor cores are available.

See the R Journal article https://doi.org/10.32614/RJ-2019-030 for more information.
It is also available as vignette.

R> vignette("optimParallel")