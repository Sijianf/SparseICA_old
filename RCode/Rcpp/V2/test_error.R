rm(list = ls())
library(Rcpp)
library(RcppArmadillo)

source("0_jngcaFunctions.R")
#source("0_SparseICA_Rcpp.R")
#sourceCpp("[External]_Re__C++_codes_related_to_sparse_ICA/ArmadilloExamples.cpp")

##########################################################################################
# simulate data from noisyICA setting
set.seed(2023)
simData = SimFMRI123(var.inactive = 0,noisyICA = TRUE, snr=0.4) 
xmat = simData$X
smat = simData$S


xData=xmat
n.comp=3
nu = 1
U.list=NULL
whiten = 'eigenvec'
orth.method='svd'
restarts.pbyd = 10
lambda = sqrt(2)/2
irlba = FALSE
eps = 0.000001
maxit.laplace = 500


p = ncol(xData) # Extract the dimensions to initialize the V(0) for "relax_logistic" and "relax_soft" methods. p = number of column
d = n.comp # Extract the components need for V(0)

# Data whiten
xData = scale(xData, center=TRUE, scale=FALSE)
if (d > p) stop('d must be less than or equal to p')
if (whiten=='eigenvec') {
  # Use whitener=='eigenvec' so that restarts.dbyd initiates from the span of the first d eigenvectors.
  temp = whitener(X = xData,n.comp = p,irlba=irlba)
  xData = temp$Z
  whitener = temp$whitener
  rm(temp)
}  else if (whiten=='sqrtprec') {
  est.sigma = cov(xData)  ## Use eigenvalue decomposition rather than SVD.
  evd.sigma = svd(est.sigma)
  whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
  xData = xData%*%whitener
}else {
  whitener = diag(p)
}

# Randomly make an input for U(0), by default here used orth.method = "svd"
W.list = U.list
if (is.null(W.list)){
  W.list = gen.inits(p=p, d=d, runs = restarts.pbyd, orth.method=orth.method)
} else{
  W.list = list(U.list)
}
runs = length(W.list)

# Create a NULL list for storing outcomes
out.list = NULL

sourceCpp("Rcpp_func_sparseICA.cpp")

a=W.list[[1]]
b=runs_relax_laplace(xData,W.list,runs,n.comp,nu,lambda,maxit.laplace,eps)

c=xData%*%a

identical(c,b$newV)




k=1
newV = xData %*% W.list[[k]]  

#sourceCpp("Rcpp_func_sparseICA_test.cpp")
sourceCpp("Rcpp_func_sparseICA_test.cpp")
soft_thresh2 = function(x, nu = 1, lambda = sqrt(2)/2) {
  xmin = pmax(abs(x)-nu/lambda,0)*sign(x)
  return(xmin)
}

soft_thresh(3,2/sqrt(2))
soft_thresh2(3)


a=procrustes(t(xData)%*%newV,diag(3))
XU = xData%*%a
for (j in 1:1089){
  for (k in 1:3){
    XU[j, k] = soft_thresh(XU[j, k], nu/lambda);
  }
}




# add cpp
Rcpp_res=relax_laplace(xData,newV,nu,lambda,1,eps)

b=Rcpp_res$converge
b[80]





