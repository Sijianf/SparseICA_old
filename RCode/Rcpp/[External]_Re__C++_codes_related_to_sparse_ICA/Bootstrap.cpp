#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix bootstrap_cpp(NumericVector ds, int B = 1000) {
  // Preallocate storage for statistics
  NumericMatrix boot_stat(B, 2); 
  // Number of observations
  int n = ds.size();
  // Perform bootstrap
  for(int i = 0; i < B; i++) {
    // Sample initial data
    NumericVector gen_data = ds[floor(runif(n, 0, n))]; 
    // Calculate sample mean and std dev 
    boot_stat(i, 0) = mean(gen_data);
    boot_stat(i, 1) = sd(gen_data);
  }
  // Return bootstrap results
  return boot_stat; 
}

