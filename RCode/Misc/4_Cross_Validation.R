
#########################################################################################################
############################ Cross-Validation for tuning parameter nu ###################################
########################### Zihang Wang ####### Last updates: 2/17/2022 ##################################
#########################################################################################################

rm(list = ls())

###########################################################
# Source ICA functions
source("0_refinedICA.R")
source('0_jngcaFunctions.R')

##########################################################
# Read simulated data
# xmat = as.matrix(read.csv('../Data/X.csv'))
# smat = as.matrix(read.csv('../Data/S.csv'))

# Simulate data from SimFMRI123() function
# Seetings: noisyICA=TRUE, nTR=100
simData = SimFMRI123(var.inactive = 0, snr=0.2) 
xmat = simData$whitened.X
smat = simData$S
# Sample size (Number of columns of X)
n=100

###########################################################
# Cross Validation
# Version 1
# Compare the estimated V using Sparse ICA on training set and test set

## Step1, group columns into K equal-sized groups
# Number of equal-sized groups (folds)
K=5
# Number of components
Q=3
index=sample(1:n,size=n)
groups=matrix(index,nrow = K)

## Step2
CV_score1=c()
for (nu in c(0.01,0.03,0.05,0.07,0.1,0.2,0.3,0.5,0.7,1)) {
  sum1=c()
  for (k in 1:K) {
    V_negK_relax_laplace = refinedICA(xData = xmat[,c(groups[-k,])], n.comp = Q, restarts.pbyd = 7, method = "relax_laplace", lambda = sqrt(2)/2, nu = nu, converge_plot = F)
    V_K_relax_laplace = refinedICA(xData = xmat[,groups[k,]], n.comp = Q, restarts.pbyd = 5, method = "relax_laplace", lambda = sqrt(2)/2, nu = nu, converge_plot = F)
    
    sum1=c(sum1,frobICA(S1=V_negK_relax_laplace$estS,S2=V_K_relax_laplace$estS,standardize=TRUE))
  }
  CV_score1=c(CV_score1,mean(sum1))
}
CV_score1
plot(c(0.01,0.03,0.05,0.07,0.1,0.2,0.3,0.5,0.7,1),CV_score1,xlab = "nu")


###########################################################
# Cross Validation
# Version 2
# Compare X_k and and X^k projected by V^-k

## Step1, group columns into K equal-sized groups
K=5
Q=3
index=sample(1:n,size=n)
groups=matrix(index,nrow = K)

## Step2
CV_score2=c()
for (nu in c(0.001,0.005,0.01,0.05,0.1,0.4,0.7,1)) {
  score_k=c()
  for (k in 1:K) {
    X_negk = xmat[,c(groups[-k,])]
    V_negK_relax_laplace = refinedICA(xData = X_negk, n.comp = Q, restarts.pbyd = 7, method = "relax_laplace", lambda = sqrt(2)/2, nu = nu, converge_plot = F)
    V_negk = V_negK_relax_laplace$estS
    X_k = xmat[,c(groups[k,])]
    score1 = sum((X_k-V_negk%*%solve(crossprod(V_negk), t(V_negk))%*%X_k)^2)
    # score1 = (norm(X_k-V_negk%*%solve(t(V_negk)%*%V_negk)%*%t(V_negk)%*%X_k,type = "F"))^2
    score_k = c(score_k,score1)
  }
  CV_score2=c(CV_score2,mean(score_k))
}
CV_score2
plot(c(0.001,0.005,0.01,0.05,0.1,0.4,0.7,1),CV_score2,xlab = "nu")

image(matrix(V_negK_relax_laplace$estS[,1],33,33))
image(matrix(V_negK_relax_laplace$estS[,2],33,33))
image(matrix(V_negK_relax_laplace$estS[,3],33,33))


###########################################################
# Cross Validation
# Version 3
# Apply regular fast-ICA on X_k

## Step1, group columns into K equal-sized groups
K=5
Q=3
index=sample(1:n,size=n)
groups=matrix(index,nrow = K)

## Step2
CV_score3=c()
for (nu in c(0.01,0.03,0.05,0.07,0.1,0.2,0.3,0.5,0.7,1)) {
  sum1=c()
  for (k in 1:K) {
    V_negK_relax_laplace = refinedICA(xData = xmat[,c(groups[-k,])], n.comp = Q, restarts.pbyd = 7, method = "relax_laplace", lambda = sqrt(2)/2, nu = nu, converge_plot = F)
    V_k_fastICA = mlcaFP(xData = xmat[,groups[k,]], n.comp = Q, whiten = "eigenvec", restarts.pbyd = 30, distribution='logistic')
    sum1=c(sum1,frobICA(S1=V_negK_relax_laplace$estS,S2=V_k_fastICA$S,standardize=TRUE))
  }
  CV_score3=c(CV_score3,mean(sum1))
}
CV_score3
plot(c(0.01,0.03,0.05,0.07,0.1,0.2,0.3,0.5,0.7,1),CV_score3,xlab = "nu")

############################################################################
# Fast ICA tanh
CV_score4=c()
for (nu in c(0.01,0.03,0.05,0.07,0.1,0.2,0.3,0.5,0.7,1)) {
  sum1=c()
  for (k in 1:K) {
    V_negK_relax_laplace = refinedICA(xData = xmat[,c(groups[-k,])], n.comp = Q, restarts.pbyd = 7, method = "relax_laplace", lambda = sqrt(2)/2, nu = nu, converge_plot = F)
    V_k_fastICA = mlcaFP(xData = xmat[,groups[k,]], n.comp = Q, whiten = "eigenvec", restarts.pbyd = 1, distribution='tanh')
    sum1=c(sum1,frobICA(S1=V_negK_relax_laplace$estS,S2=V_k_fastICA$S,standardize=TRUE))
  }
  CV_score4=c(CV_score4,mean(sum1))
}
CV_score4
plot(c(0.01,0.03,0.05,0.07,0.1,0.2,0.3,0.5,0.7,1),CV_score4,xlab = "nu")





