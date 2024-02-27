#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat matrix_mult(const arma::mat& X, const arma::mat& Y) {
  int m = X.n_rows;
  int n = Y.n_cols;
  arma::mat Z(m,n);
  Z = X * Y;
  return Z;
}

// Compute coefficients and their standard error
// during multiple linear regression given a
// design matrix X containing N observations with
// P regressors and a vector y containing of
// N responses
// [[Rcpp::export]]
Rcpp::List fastLm(const arma::mat& X,
                  const arma::colvec& y) { 
  // Dimension information
  int n = X.n_rows, p = X.n_cols;
  // Fit model y ~ X
  arma::colvec coef = arma::solve(X, y); 
  // Compute the residuals
  arma::colvec res  = y - X * coef;
  // Estimated variance of the random error
    double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0) / (n - p);
  // Standard error matrix of coefficients
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(X.t()*X)));
  // Create named list with the above quantities
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef, Rcpp::Named("stderr") = std_err, Rcpp::Named("df.residual") = n - p );
}

// [[Rcpp::export]]
double normArmaM(const arma::mat& X, int p = 2){
  return(arma::norm(X, p));
}

// [[Rcpp::export]]
double normArmaV(const arma::colvec& Y, int p = 2){
  return(arma::norm(Y, p));
}

// [[Rcpp::export]]
arma::mat procrustes(arma::mat& X, arma::mat &V){
  arma::mat U;
  
  arma::vec s;
  arma::mat Q, R;
  
  arma::svd_econ(Q, s, R, X * V);
  U = Q * R.t();
  return(U);
}


// [[Rcpp::export]]
double soft_I(double a, double lambda){
  if (a > lambda){
    return(a-lambda);
  }else if (a < -lambda){
    return(a + lambda);
  }else{
    return(0);
  }
}

// [[Rcpp::export]]
Rcpp::List sparsePCA(arma::mat& X, arma::mat& Vstart, double lambda, double eps = 0.0001){
  
  // initialize U
  int r = Vstart.n_cols;
  int p = Vstart.n_rows;
  arma::vec s;
  arma::mat Q, R;
  arma::mat U;
  arma::svd_econ(Q, s, R, X * Vstart);
  U = Q * R.t();
  arma::mat V(p, r);
  
  // Calculate current objective
  double fold, fnew;
  fold = accu(square(X - U * Vstart.t()))/2 + lambda * accu(abs(Vstart));
  
  // Additional intermediate
  arma::mat XU;
  
  // To store error and objective function difference
  double error = 1000;
  
  // Alternate updates of U with updates of V
  while (error > eps){
    XU = X.t() * U;
    // Update V
    for (int j = 0; j < p; j++){
      for (int k = 0; k < r; k++){
        V(j, k) = soft_I(XU(j, k), lambda);
      }
    }
    
    // Update U
    arma::svd_econ(Q, s, R, X * V);
    U = Q * R.t();
    
    // Calculate new objective
    fnew = accu(square(X - U * V.t()))/2 + lambda * accu(abs(V));
    
    // Calculate error
    error = std::fabs(fold - fnew);
    
    fold = fnew;
  }
  return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("V") = V, Rcpp::Named("error") = error);
}