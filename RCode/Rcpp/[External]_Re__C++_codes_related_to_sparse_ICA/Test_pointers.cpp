#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


//[[Rcpp::export]]
arma::vec timesTwo(arma::vec x) {
  x = x * 2;
  return x;
}

//[[Rcpp::export]]
arma::vec timesTwo_pointer(arma::vec& x) {
  x = x * 2;
  return x;
}


//[[Rcpp::export]]
arma::vec timesTwo_pointer(const arma::vec& x)  {
   x = x * 2;
   return x;
}