// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::mat A(const arma::mat& X){
  arma::uvec indexX = arma::find (X > 1.5);
  return (X(indexX));
}

// [[Rcpp::export]]
arma::mat B(const arma::mat& X, const arma::vec& Y){
  arma::uvec indexY = arma::find(Y > 0);
  return(X.rows(indexY));
}