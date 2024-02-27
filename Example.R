##########################################################################################
########################### Toy Example of Using Sparse ICA ##############################
##########################################################################################

rm(list = ls())

source("0_Support_Functions.R")
source("1_SparseICA_R.R")
source("2_SparseICA_Rcpp.R")

##########################################################################################
# Read simulated data from Data folder
#xmat = read.csv("Data/Xmat.csv",header = F)
#smat = read.csv("Data/Smat.csv",header = F)

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

##########################################################################################
# use BIC on sparse LNGCA
# R version
BIC_results_R = BIC_sparseICA_R(xData = xmat,n.comp = 3,whiten = "eigenvec",BIC_plot = T,show_message = T)
best_nu=BIC_results_R$best_nu
my_sparseICA_R = sparseICA_R(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", restarts.pbyd = 10,lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,maxit.laplace = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(my_sparseICA_R$estS[,1],33,33))
image(matrix(my_sparseICA_R$estS[,2],33,33))
image(matrix(my_sparseICA_R$estS[,3],33,33))
par(mfrow=c(1,1))

# Rcpp version
BIC_results_Rcpp = BIC_sparseICA_Rcpp(xData = xmat,n.comp = 3,whiten = "eigenvec",BIC_plot = T,show_message = T)
best_nu=BIC_results_Rcpp$best_nu
my_sparseICA_Rcpp = sparseICA_Rcpp(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", restarts.pbyd = 10,lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,maxit.laplace = 500, converge_plot = F)

par(mfrow=c(1,3))
image(matrix(my_sparseICA_Rcpp$estS[,1],33,33))
image(matrix(my_sparseICA_Rcpp$estS[,2],33,33))
image(matrix(my_sparseICA_Rcpp$estS[,3],33,33))
par(mfrow=c(1,1))

###############################################################################################
# Use BIC on PCA+sparseICA
BIC_results2_R=BIC_PCA_sparseICA_R(xData = xmat,n.comp = 3,whiten = "eigenvec",BIC_plot = T,show_message = T)
best_nu=BIC_results2_R$best_nu
my_PCA_sparseICA_R = PCA_sparseICA_R(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", restarts.pbyd = 10, lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,maxit.laplace = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(my_PCA_sparseICA_R$estS[,1],33,33))
image(matrix(my_PCA_sparseICA_R$estS[,2],33,33))
image(matrix(my_PCA_sparseICA_R$estS[,3],33,33))
par(mfrow=c(1,1))


BIC_results2_Rcpp=BIC_PCA_sparseICA_Rcpp(xData = xmat,n.comp = 3,whiten = "eigenvec",BIC_plot = T,show_message = T)
best_nu=BIC_results2_Rcpp$best_nu
my_PCA_sparseICA_Rcpp = PCA_sparseICA_Rcpp(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", restarts.pbyd = 10, lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,maxit.laplace = 500, converge_plot = F)

par(mfrow=c(1,3))
image(matrix(my_PCA_sparseICA_Rcpp$estS[,1],33,33))
image(matrix(my_PCA_sparseICA_Rcpp$estS[,2],33,33))
image(matrix(my_PCA_sparseICA_Rcpp$estS[,3],33,33))
par(mfrow=c(1,1))

