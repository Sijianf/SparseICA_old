##########################################################################################
########################### Toy Example of Using Sparse ICA ##############################
##########################################################################################

rm(list = ls())
library(irlba)
library(singR)

source('../References/0_jngcaFunctions.R')
source("../References/0_SparseICA.R")
source("../References/0_Cross_Validation_SparseICA.R")
source("../References/0_Cross_Validation_single_PCA_SparseICA.R")

##########################################################################################
# Read simulated data from Data folder
# xmat = read.csv("../References/Xmat.csv",header = F)
# smat = read.csv("../References/Smat.csv",header = F)

# simulate data from noisyICA setting
set.seed(2022)
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

##########################################################################################
# Cross-Validation to find best nu
my_nulist=c(0.001,0.003,0.005,0.007,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7,0.8,0.9,1)
my_nulist=seq(0.001,1.8,0.01)

##########################################################################################
# crude
best_nu=CV_sparseICA(xmat,n.comp=3,fold=5,CV_method="projection",restarts.pbyd = 10, method = "sparse_laplace", lambda = sqrt(2)/2, eps = 1e-6,loop_laplace = 500, loop_logistic = 500,nu_list=c(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2),show_message=T,make_plot=T)

# Fit SparseICA using the best nu selected from CV
estX_sparse_laplace = sparseICA(xData = xmat, n.comp = 3, restarts.pbyd = 10, method = "sparse_laplace", lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,loop_laplace = 500, loop_logistic = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(estX_sparse_laplace$estV[,1],33,33))
image(matrix(estX_sparse_laplace$estV[,2],33,33))
image(matrix(estX_sparse_laplace$estV[,3],33,33))
par(mfrow=c(1,1))


#####################################################################################
# Compare with fastICA method
estX_fast_logistic = sparseICA(xData = xmat, n.comp = 3, restarts.pbyd = 10, method = "fast_logistic", lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,loop_laplace = 500, loop_logistic = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(estX_fast_logistic$estV[,1],33,33))
image(matrix(estX_fast_logistic$estV[,2],33,33))
image(matrix(estX_fast_logistic$estV[,3],33,33))
par(mfrow=c(1,1))


##########################################################################################
# PCA + sparseICA
best_nu=CV_PCA_sparseICA(random_seed=2021,xmat,n.comp=3,fold=5,CV_method="projection",restarts.pbyd = 10, method = "sparse_laplace", lambda = sqrt(2)/2, eps = 1e-6,loop_laplace = 500, loop_logistic = 500,nu_list=my_nulist,show_message=T,make_plot=T)

# PCA + sparse ICA
dat_all = scale(xmat)
dat_all = xmat
# perform PCA
subj_PCA=prcomp_irlba(dat_all,3,center = T,scale. = T)
myPC=subj_PCA$x
dimnames(myPC)=NULL
my_PCA_sparseICA = sparseICA(xData = myPC, n.comp = 3,whiten = "eigenvec", restarts.pbyd = 10, method = "sparse_laplace", lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,loop_laplace = 500, loop_logistic = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(my_PCA_sparseICA$estV[,1],33,33))
image(matrix(my_PCA_sparseICA$estV[,2],33,33))
image(matrix(my_PCA_sparseICA$estV[,3],33,33))
par(mfrow=c(1,1))

##########################################################################################
# sparse PCA
my_sparsePCA = ssvd(t(xmat),k=3,n=70)

par(mfrow=c(1,3))
image(matrix(my_sparsePCA$v[,1],33,33))
image(matrix(my_sparsePCA$v[,2],33,33))
image(matrix(my_sparsePCA$v[,3],33,33))
par(mfrow=c(1,1))


