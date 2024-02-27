# Other packages may be required: 
# install.packages("steadyICA")
# install.packages("neuRosim")
# install.packages("ProDenICA")

rm(list = ls())

#source("0_refinedICA_v2_03312022.R")
#source("0_refinedICA_v2_03242022.R")
source("0_refinedICA_v3_03312022.R")
source('0_jngcaFunctions.R')

#####################################################################################
# data showcase

## Data generation process: x and s matrices from function in jingcaFunctions.R
set.seed(2022)
simData = SimFMRI123(var.inactive = 0, snr=0.2)
xmat = simData$X #add random noise
smat = simData$S # true matrix
# each column contains pixels of 33x33 image, 50 time points in total

# Can directly use the provided datasets
# whitening data matrix + noise: X matrix
# true data matrix: S matrix
# xmat = as.matrix(read.csv('../Data/X.csv'))
# smat = as.matrix(read.csv('../Data/S.csv'))

# notice that the data matrix x is pre-whitened:
cov(xmat)[1:5,1:5]

# Examine the true components:
# There are 3 components, plot figures using image() function
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


############################################################################################################################
# check performance of four algorithms
#1. Check FastICA method with logistic density
estX_logis = mlcaFP(xData = xmat, n.comp = 3, whiten = "eigenvec", restarts.pbyd = 30, distribution='logistic')
par(mfrow=c(1,3))
image(matrix(estX_logis$S[,1],33,33))
image(matrix(estX_logis$S[,2],33,33))
image(matrix(estX_logis$S[,3],33,33))
par(mfrow=c(1,1))

#2. Check FastICA method with tanh density
estX_tanh = mlcaFP(xData = xmat, n.comp = 3, whiten = "eigenvec", restarts.pbyd = 1, distribution='tanh')
par(mfrow=c(1,3))
image(matrix(estX_tanh$S[,1],33,33))
image(matrix(estX_tanh$S[,2],33,33))
image(matrix(estX_tanh$S[,3],33,33))
par(mfrow=c(1,1))

#3. Check Relax-Logistics method
# Notice: nu=1 in Sijian's code
# Use nu=0.1 here
# Suggestion in Sijian's thesis: nu=0.1
Q=3
estX_relax_logistic = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 10, method = "relax_logistic", lambda = sqrt(2)/2, nu = 0.1,eps = 1e-5, loop_laplace = 500, loop_logistic = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(estX_relax_logistic$estV[,1],33,33))
image(matrix(estX_relax_logistic$estV[,2],33,33))
image(matrix(estX_relax_logistic$estV[,3],33,33))
par(mfrow=c(1,1))

Q=3
estX_relax_logistic = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 10, method = "relax_logistic", lambda = sqrt(2)/2, nu = 0.1,eps = 1e-5,converge_plot = F)
par(mfrow=c(1,3))
image(matrix(estX_relax_logistic$estV[,1],33,33))
image(matrix(estX_relax_logistic$estV[,2],33,33))
image(matrix(estX_relax_logistic$estV[,3],33,33))
par(mfrow=c(1,1))

#4. Check Relax-Softmax(Laplace) method
# Notice: nu=1 in Sijian's code
# Use nu=0.1 here, lambda=sqrt(2)/2
# Suggestion in Sijian's thesis: nu=1
Q=3
estX_relax_laplace = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 10, method = "relax_laplace", lambda = sqrt(2)/2, nu = 0.1,eps = 1e-5,loop_laplace = 500, loop_logistic = 500, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(estX_relax_laplace$estV[,1],33,33))
image(matrix(estX_relax_laplace$estV[,2],33,33))
image(matrix(estX_relax_laplace$estV[,3],33,33))
par(mfrow=c(1,1))

Q=3
estX_relax_laplace = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 10, method = "relax_laplace", lambda = sqrt(2)/2, nu = 0.1,eps = 1e-5,converge_plot = F)
par(mfrow=c(1,3))
image(matrix(estX_relax_laplace$estV[,1],33,33))
image(matrix(estX_relax_laplace$estV[,2],33,33))
image(matrix(estX_relax_laplace$estV[,3],33,33))
par(mfrow=c(1,1))


##################################################################################################################################################################
# Simulation study

# # Generate simulated data
# simData = SimFMRI123(var.inactive = 0, snr=0.2) 
# xmat = simData$whitened.X
# smat = simData$S

# number of components Q
Q = 3

# number of simulations itm
# e.g.: Do 10 Simulations 
itm = 30 

##1. FastICA logistic 
## (This is just to run the FastICA codes with logistic density)
FastICA_logistic_diff = c()
for (sim in 1:itm) {
  simData = SimFMRI123(var.inactive = 0, snr=0.2) 
  xmat = simData$whitened.X
  smat = simData$S
  Q = 3
  estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 1, method = "FastICA_logistic", lambda = 0.5, nu = 1)
  FastICA_logistic_diff[sim] = frobICA(S1 = estX$estS, S2 = smat, standardize=TRUE)
}
par(mfrow=c(1,1))
plot(1:itm, FastICA_logistic_diff, ylim = c(0,1.2),
     main = "FastICA-Logistic", xlab = paste0("Iteration (1~", itm,")"), ylab = "FastICA-Logistic Distance")


##2. FastICA tanh 
## (This is just to run the FastICA codes with tanh density)
FastICA_tanh_diff = c()
for (sim in 1:itm) {
  simData = SimFMRI123(var.inactive = 0, snr=0.2) 
  xmat = simData$whitened.X
  smat = simData$S
  Q = 3
  estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 1, method = "FastICA_tanh", lambda = 0.5, nu = 1)
  FastICA_tanh_diff[sim] = frobICA(S1 = estX$estS, S2 = smat, standardize=TRUE)
}
par(mfrow=c(1,1))
plot(1:itm, FastICA_tanh_diff, ylim = c(0,1.2),
     main = "FastICA-tanh", xlab = paste0("Iteration (1~", itm,")"), ylab = "FastICA-tanh Distance")


##3. Relax logistic 
## (Typically, this one is executed slowly, especially compared with others)
## Use nu=0.1
# Sijian's thesis suggestion: nu=0.1
relax_logistic_diff = c()
for (sim in 1:itm) {
  simData = SimFMRI123(var.inactive = 0, snr=0.2) 
  xmat = simData$whitened.X
  smat = simData$S
  Q = 3
  estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 2, method = "relax_logistic", lambda = 0.5, nu = 0.1, converge_plot = F)
  relax_logistic_diff[sim] = frobICA(S1 = estX$estS, S2 = smat, standardize=TRUE)
}
par(mfrow=c(1,1))
plot(1:itm, relax_logistic_diff, ylim = c(0,1.2),
     main = "Relax-Logistic", xlab = paste0("Iteration (1~", itm,")"), ylab = "Relax-Logistic Distance")


##4. Relax Softmax (Laplace) 
## (This runs fast, and it can allow sparsity results)
# Use nu=0.1
# Sijian's thesis suggestion: nu=1
relax_softmax_diff = c()
for (sim in 1:itm) {
  simData = SimFMRI123(var.inactive = 0, snr=0.2) 
  xmat = simData$whitened.X
  smat = simData$S
  Q = 3
  estX = refinedICA(xData = xmat, n.comp = Q, restarts.pbyd = 5, method = "relax_laplace", lambda = 0.5, nu = 0.1, converge_plot = F)
  relax_softmax_diff[sim] = frobICA(S1 = estX$estS, S2 = smat, standardize=TRUE)
}
par(mfrow=c(1,1))
plot(1:itm, relax_softmax_diff, ylim = c(0,1.2),
     main = "Relax-Softmax", xlab = paste0("Iteration (1~", itm,")"), ylab = "Relax-Softmax Distance")



## Combine the above plots into one. Need to change the margins for a good plot.
par(mfrow = c(2,2),mar=c(2,2,4,1))
plot(1:itm, relax_logistic_diff, ylim = c(0,1.2),
     main = "Relax-Logistic", xlab = paste0("Iteration (1~", itm,")"), ylab = "Relax-Logistic Distance")
plot(1:itm, relax_softmax_diff, ylim = c(0,1.2),
     main = "Relax-Softmax", xlab = paste0("Iteration (1~", itm,")"), ylab = "Relax-Softmax Distance")
plot(1:itm, FastICA_logistic_diff, ylim = c(0,1.2),
     main = "FastICA-Logistic", xlab = paste0("Iteration (1~", itm,")"), ylab = "FastICA-Logistic Distance")
plot(1:itm, FastICA_tanh_diff, ylim = c(0,1.2),
     main = "FastICA-tanh", xlab = paste0("Iteration (1~", itm,")"), ylab = "FastICA-tanh Distance")
par(mfrow=c(1,1))


## Make boxplot
library(ggplot2)
relax_logistic_diff_df = data.frame(Distance = relax_logistic_diff)
relax_logistic_diff_df$Method = "Relax-Logistic"
relax_softmax_diff_df = data.frame(Distance = relax_softmax_diff)
relax_softmax_diff_df$Method = "Relax-Softmax"
FastICA_logistic_diff_df = data.frame(Distance = FastICA_logistic_diff)
FastICA_logistic_diff_df$Method = "mlcaFP-Logistic"
FastICA_tanh_diff_df = data.frame(Distance = FastICA_tanh_diff)
FastICA_tanh_diff_df$Method = "mlcaFP-tanh"
diff = rbind(relax_logistic_diff_df,relax_softmax_diff_df,FastICA_logistic_diff_df,FastICA_tanh_diff_df)

ggplot(diff, aes(x = Method, y = Distance))+
  geom_boxplot(outlier.alpha = 0.5)




