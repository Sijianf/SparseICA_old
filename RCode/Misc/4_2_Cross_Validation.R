
#########################################################################################################
############################ Cross-Validation for tuning parameter nu ###################################
########################### Zihang Wang ####### Last updates: 3/31/2022 ##################################
#########################################################################################################

rm(list = ls())

#########################################################################################################
# Source ICA functions
source("0_refinedICA_v2_03312022.R")
#source("0_refinedICA_v3_03312022.R")
source('0_jngcaFunctions.R')

#########################################################################################################
# Read simulated data
# xmat = as.matrix(read.csv('../Data/X.csv'))
# smat = as.matrix(read.csv('../Data/S.csv'))

# Simulate data from SimFMRI123() function
# Seetings: noisyICA=TRUE, nTR=100
simData = SimFMRI123(noisyICA=FALSE, nTR=100,var.inactive = 0, snr=0.2) 
xmat = simData$X
smat = simData$S
# Check Data
image(matrix(xmat[,15],33,33))
image(matrix(smat[,3],33,33))
# Sample size (Number of columns of X)
n=100

#########################################################################################################
# Cross Validation
## Step1, group columns into K equal-sized groups
# Number of equal-sized groups (folds)
K=5
# Number of components
Q=3
index=sample(1:n,size=n)
groups=matrix(index,nrow = K)
#########################################################################################################

# CV -version2
nu_list = c(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2)
CV_score2=rep(NA,length(nu_list))
for (i in 1:length(nu_list)) {
  nu = nu_list[i]
  score_k=rep(NA,K)
  for (k in 1:K) {
    X_negk = xmat[,c(groups[-k,])]
    V_negK_relax_laplace = refinedICA(xData = X_negk, n.comp = Q, whiten = "eigenvec", restarts.pbyd = 10, method = "relax_laplace", lambda = sqrt(2)/2, nu = nu, eps = 1e-5,converge_plot = F)
    V_negk = V_negK_relax_laplace$estV
    X_k = xmat[,c(groups[k,])]
    score1 = sum((X_k-V_negk%*%solve(crossprod(V_negk), t(V_negk))%*%X_k)^2)
    score_k[k] = score1
  }
  CV_score2[i]=mean(score_k)
}
CV_score2
plot(nu_list,CV_score2,xlab = "nu",type = "l")
nu_list[which(CV_score2==min(CV_score2))]

Q=3
estX_relax_laplace = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 10, method = "relax_laplace", lambda = sqrt(2)/2, nu = 0.1,eps = 1e-5,converge_plot = F)
par(mfrow=c(1,3))
image(matrix(estX_relax_laplace$estV[,1],33,33))
image(matrix(estX_relax_laplace$estV[,2],33,33))
image(matrix(estX_relax_laplace$estV[,3],33,33))
par(mfrow=c(1,1))


# CV - version 3 - FastICA-Logistic
nu_list = c(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2)
CV_score3=rep(NA,length(nu_list))
for (i in 1:length(nu_list)) {
  nu = nu_list[i]
  sum1=rep(NA,K)
  for (k in 1:K) {
    V_negK_relax_laplace = refinedICA(xData = xmat[,c(groups[-k,])], n.comp = Q, whiten = "eigenvec", restarts.pbyd = 10, method = "relax_laplace", lambda = sqrt(2)/2, nu = nu,eps = 1e-5,loop_laplace = 500, loop_logistic = 500, converge_plot = F)
    V_k_fastICA = refinedICA(xData = xmat[,groups[k,]], n.comp = Q, whiten = "eigenvec", restarts.pbyd = 10, method = "FastICA_logistic")
    sum1[k]=frobICA(S1=V_negK_relax_laplace$estV,S2=V_k_fastICA$estV,standardize=TRUE)
  }
  CV_score3[i]=mean(sum1)
}
CV_score3
plot(nu_list,CV_score3,xlab = "nu",type = "l")

# CV - version 3 - FastICA-tanh
nu_list = c(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2)
CV_score4=rep(NA,length(nu_list))
for (i in 1:length(nu_list)) {
  nu=nu_list[i]
  sum1=rep(NA,K)
  for (k in 1:K) {
    V_negK_relax_laplace = refinedICA(xData = xmat[,c(groups[-k,])], n.comp = Q, whiten = "eigenvec", restarts.pbyd = 10, method = "relax_laplace", lambda = sqrt(2)/2, nu = nu,eps = 1e-5,loop_laplace = 500, loop_logistic = 500, converge_plot = F)
    V_k_fastICA = refinedICA(xData = xmat[,groups[k,]], n.comp = Q, whiten = "eigenvec", restarts.pbyd = 10, method = "FastICA_tanh")
    sum1[k]=frobICA(S1=V_negK_relax_laplace$estV,S2=V_k_fastICA$estV,standardize=TRUE)
  }
  CV_score4[i]=mean(sum1)
}
CV_score4
plot(nu_list,CV_score4,xlab = "nu",type = "l")


which(CV_score3==min(CV_score3))
which(CV_score4==min(CV_score4))
