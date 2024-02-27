##########################################################################################
########################### test CV on PCA+SparseICA ##############################
##########################################################################################

rm(list = ls())
library(irlba)
library(singR)

source('../References/0_jngcaFunctions.R')
source("../References/0_SparseICA.R")
source("../References/0_PCA_SparseICA.R")
source("../References/0_BIC_PCA_SparseICA.R")
source("../References/0_BIC_SparseICA.R")

source("../References/0_PCA_SparseICA_test.R")
##########################################################################################
# simulate data from noisyICA setting
set.seed(2023)
simData = SimFMRI123(var.inactive = 0,noisyICA = TRUE, snr=0.4) 
xmat = simData$X
smat = simData$S

# Check Data
par(mfrow=c(1,3))
image(matrix(smat[,1],33))
image(matrix(smat[,2],33))
image(matrix(smat[,3],33))
par(mfrow=c(1,1))

# examine the first three time points
par(mfrow=c(1,3))
image(matrix(xmat[,1],33))
image(matrix(xmat[,2],33))
image(matrix(xmat[,3],33))
par(mfrow=c(1,1))
# examine any time point
image(matrix(xmat[,35],33))

for (i in 1:50) {
  image(matrix(xmat[,i],33))
  Sys.sleep(0.5)
  cat(i,"\n")
}

################################################################################
# test BIC on PCA+sparseICA
BIC_results=BIC_PCA_sparseICA(xData = xmat,n.comp = 3,whiten = "eigenvec",method = "sparse_laplace",BIC_plot = T,show_message = T)

my_PCA_sparseICA = PCA_sparseICA(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", restarts.pbyd = 10, method = "sparse_laplace", lambda = sqrt(2)/2, nu = BIC_results$best_nu,eps = 1e-6,maxit.laplace = 500, maxit.logistic = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(my_PCA_sparseICA$estS[,1],33,33))
image(matrix(my_PCA_sparseICA$estS[,2],33,33))
image(matrix(my_PCA_sparseICA$estS[,3],33,33))
par(mfrow=c(1,1))

################################################################################
# test BIC on raw sparseICA
BIC_results2 = BIC_sparseICA(xData = xmat,n.comp = 3,whiten = "eigenvec",method = "sparse_laplace",BIC_plot = T,show_message = T)

best_nu=BIC_results2$best_nu
best_nu=1
my_sparseICA = sparseICA(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", restarts.pbyd = 10, method = "sparse_laplace", lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,maxit.laplace = 500, maxit.logistic = 500, converge_plot = T)
par(mfrow=c(1,3))
image(matrix(my_sparseICA$estS[,1],33,33))
image(matrix(my_sparseICA$estS[,2],33,33))
image(matrix(my_sparseICA$estS[,3],33,33))
par(mfrow=c(1,1))


################################################################################
rm(list = ls())

source("0_jngcaFunctions.R")

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

k=1
newV = xData %*% W.list[[k]]  

soft_thresh = function(x, nu = 1, lambda = sqrt(2)/2) {
  xmin = pmax(abs(x)-nu/lambda,0)*sign(x)
  return(xmin)
}


soft_thresh(3,1)

txv = t(xData)%*%newV
svd.txv = La.svd(txv)
newU = svd.txv$u%*%svd.txv$vt # define objective function for orthogonal U


