#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
double soft_thresh(double x, double alpha){
  if (x > alpha){
    return(x-alpha);
  }else if (x < -alpha){
    return(x + alpha);
  }else{
    return(0);
  }
}


// [[Rcpp::export]]
arma::mat procrustes(arma::mat& X, arma::mat& V){
  arma::mat U;
  
  arma::vec s;
  arma::mat Q, R;
  
  arma::svd_econ(Q, s, R, X * V);
  U = Q * R.t();
  return(U);
}


// [[Rcpp::export]]
Rcpp::List relax_laplace(const arma::mat& xData,const arma::mat& newV, double nu, double lambda, double maxit, double eps = 0.000001){
  // initialize 
  int p = newV.n_rows;
  int r = newV.n_cols;
  int t = xData.n_cols;
  arma::mat newU(t,r);
  arma::vec converge(maxit);
  arma::mat XTV(t,r);

  // Intermediate V and XU
  arma::mat lagV = newV;
  arma::mat newV2 = newV;
  arma::mat XU(p,r);
  arma::mat IdenMat=arma::eye(r,r);

  for (int i = 0; i < maxit; i++)
  {
    // update U, use irina's function
    XTV = xData.t()*newV2;
    newU = procrustes(XTV,IdenMat);

    // update V
    XU = xData*newU;
    for (int j = 0; j < p; j++){
      for (int k = 0; k < r; k++){
        newV2(j, k) = soft_thresh(XU(j, k), nu/lambda);
      }
    }

    // check convergence
    converge[i] = arma::norm(lagV-newV2,"fro");
    //converge[i] = accu(square(lagV-newV2))/(p*r*r);
    if (converge[i] < eps){
      break;
    }
    lagV = newV2;
  }
  return Rcpp::List::create(Rcpp::Named("newU") = newU, Rcpp::Named("newV") = newV2, Rcpp::Named("converge") = converge);
}