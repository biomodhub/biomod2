#include <RcppArmadillo.h>
using namespace arma;

#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

int matmodifcmCpp (arma::mat& tab, const arma::vec& poili);
int matmodifcnCpp (arma::mat& tab, const arma::vec& poili);
int matmodifcsCpp (arma::mat& tab, const arma::vec& poili);
int matmodifcpCpp (arma::mat& tab, const arma::vec& poili);
int matmodiffcCpp (arma::mat& tab, const arma::vec& poili);
int matcentrageCpp (arma::mat& A, const arma::vec& poili, const int typ);
double betweenvarCpp (const arma::mat& tab, const arma::vec& pl, Rcpp::IntegerVector fac);
double inerbetweenCpp (const arma::vec& pl, const arma::vec& pc, const int moda, Rcpp::IntegerVector indica, const arma::mat& tab);
