#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// soft-thresholding function
// [[Rcpp::export]]
double soft_I_c(double a, double lambda){
  if (a > lambda){
    return(a-lambda);
  }else if (a < -lambda){
    return(a + lambda);
  }else{
    return(0);
  }
}

// Lasso objective function
// [[Rcpp::export]]
double lasso_I_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Compute n
  int n = Xtilde.n_rows;
  // Compute current residual value
  arma::colvec res  = Ytilde - Xtilde*beta;
  // Compute the loss part of the objective
  double loss = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(2*n);
  // Compute the penalty part of the objective
  double penalty = lambda * sum(abs(beta));
  return(loss + penalty);
}

// fit LASSO on standardized data
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_I_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){

  int p = Xtilde.n_cols;
  int n = Xtilde.n_rows;
  
  double error = 1000;
  // Calculate full residual
  arma::colvec r = Ytilde - Xtilde * beta_start; 
  arma::colvec beta = beta_start;
  arma::colvec beta_old = beta_start;
  while (error > eps){
    beta_old = beta;
    for (int j = 0; j < p; j++){
      // calculate inner product between Xtilde[,j] and r over n
      double inn = std::inner_product(r.begin(), r.end(), &Xtilde(0,j), 0.0)/n;
      // beta update
      beta[j] = soft_I_c(beta_old[j] + inn, lambda);
      // residual update
      r = r - Xtilde.col(j)*(beta[j]-beta_old[j]);
    }
    // Difference in objective function values
    error = lasso_I_c(Xtilde, Ytilde, beta_old, lambda) - lasso_I_c(Xtilde, Ytilde, beta, lambda);
  }
  return beta;
}

// Solve lasso on a sequence
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_I_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Initialize all the parameters
  int p = Xtilde.n_cols;
  int n_lambda = lambda_seq.n_rows;
  arma::mat beta_mat = arma::zeros(p, n_lambda);
  arma::colvec fmin_vec = arma::zeros(n_lambda);
  arma::colvec beta0 = arma::zeros(p);
  // Calculate solution for each lambda
  for (int i = 0; i < n_lambda; i++){
    beta_mat.col(i) = fitLASSOstandardized_I_c(Xtilde, Ytilde, lambda_seq[i], beta0, eps);
    // Update starting point to use "warm starts"
    beta0 = beta_mat.col(i);
  }
  return beta_mat;
}