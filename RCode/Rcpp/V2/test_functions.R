##########################################################################################
########################### test CV on PCA+SparseICA ##############################
##########################################################################################

rm(list = ls())
#library(irlba)
#library(singR)

source("0_jngcaFunctions.R")
source("0_SparseICA_Rcpp.R")
source("0_PCA_SparseICA_Rcpp.R")

##########################################################################################
# simulate data from noisyICA setting
set.seed(2023)
simData = SimFMRI123(var.inactive = 0,noisyICA = TRUE, snr=0.4) 
xmat = simData$X
smat = simData$S

best_nu=1.24

my_sparseICA = sparseICA_Rcpp(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", restarts.pbyd = 10,lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,maxit.laplace = 500, converge_plot = T)
par(mfrow=c(1,3))
image(matrix(my_sparseICA$estS[,1],33,33))
image(matrix(my_sparseICA$estS[,2],33,33))
image(matrix(my_sparseICA$estS[,3],33,33))
par(mfrow=c(1,1))

best_nu=1.46

my_PCA_sparseICA = PCA_sparseICA_Rcpp(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", restarts.pbyd = 10, lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,maxit.laplace = 500, converge_plot = T)
par(mfrow=c(1,3))
image(matrix(my_PCA_sparseICA$estS[,1],33,33))
image(matrix(my_PCA_sparseICA$estS[,2],33,33))
image(matrix(my_PCA_sparseICA$estS[,3],33,33))
par(mfrow=c(1,1))


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

my_sparseICA = sparseICA(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", restarts.pbyd = 10, method = "sparse_laplace", lambda = sqrt(2)/2, nu = BIC_results2$best_nu,eps = 1e-6,maxit.laplace = 500, maxit.logistic = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(my_sparseICA$estS[,1],33,33))
image(matrix(my_sparseICA$estS[,2],33,33))
image(matrix(my_sparseICA$estS[,3],33,33))
par(mfrow=c(1,1))
