##########################################################################################
################### Investigate what's going on when changing nTR ########################
##########################################################################################

rm(list = ls())

source("0_refinedICA_v2_03312022.R")
#source('jngcaFunctions_raw.r')
source("jngcaFunctions_raw_edited.r")
source("0_CV_refinedICA_04132022.R")

##########################################################################################
# Simulate data from SimFMRI123() function
simData = SimFMRI123(noisyICA=F, nTR=100,var.inactive = 0.0001, snr=0.4) 
#simData2 = SimFMRI123(nTR=50)
xmat = simData$X
smat = simData$S

# Check Data
par(mfrow=c(1,3))
image(matrix(smat[,1],33))
image(matrix(smat[,2],33))
image(matrix(smat[,3],33))
par(mfrow=c(1,1))

# Play the whole data
for (i in 1:ncol(xmat)) {
  image(matrix(xmat[,i],33))
  Sys.sleep(0.5)
  cat(i,"\n")
}

# examine any time point
image(matrix(xmat[,20],33))

##########################################################################################
# Cross-Validation to find best nu

best_nu=CV_refinedICA(xmat,Q=3,fold=5,CV_method="projection",restarts.pbyd = 10, method = "relax_laplace", lambda = sqrt(2)/2, eps = 1e-6,loop_laplace = 500, loop_logistic = 500,nu_list=c(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2),show_message=T,make_plot=T)

##########################################################################################
# Fit SparseICA using best nu
best_nu=0.1
Q=3
estX_relax_laplace = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 10, method = "relax_laplace", lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,loop_laplace = 500, loop_logistic = 500, converge_plot = T)
par(mfrow=c(1,3))
image(matrix(estX_relax_laplace$estV[,1],33,33))
image(matrix(estX_relax_laplace$estV[,2],33,33))
image(matrix(estX_relax_laplace$estV[,3],33,33))
par(mfrow=c(1,1))

estX_FastICA_tanh = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 10, method = "FastICA_tanh", lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,loop_laplace = 500, loop_logistic = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(estX_FastICA_tanh$estV[,1],33,33))
image(matrix(estX_FastICA_tanh$estV[,2],33,33))
image(matrix(estX_FastICA_tanh$estV[,3],33,33))
par(mfrow=c(1,1))








