
#########################################################################################################
################# Calculation of Blue and Red parts in optimization function ############################
########################### Zihang Wang ####### Last updates: 2/3/2022 ##################################
#########################################################################################################

###########################################################
# Source ICA functions
source("refinedICA_edited.R")
source('jngcaFunctions.R')

##########################################################
# Read simulated data
xmat = as.matrix(read.csv('../Data/X.csv'))
smat = as.matrix(read.csv('../Data/S.csv'))

######################################################################################################
#For Relax-Softmax(Laplace) method, check BLUE and RED parts in the optimization function
# Notice: nu=1 in Sijian's code
# Use nu=0.1 here, lambda=sqrt(2)/2
# Suggestion in Sijian's thesis: nu=1
######################################################################################################
Q=3
estX_relax_laplace = refinedICA(xData = xmat, n.comp = Q,whiten="eigenvec",orth.method="svd", restarts.pbyd = 5, method = "relax_laplace", lambda = sqrt(2)/2, nu = 0.1, converge_plot = F)
par(mfrow=c(1,3))
image(matrix(estX_relax_laplace$estS[,1],33,33))
image(matrix(estX_relax_laplace$estS[,2],33,33))
image(matrix(estX_relax_laplace$estS[,3],33,33))
par(mfrow=c(1,1))

dim(estX_relax_laplace$estS)
dim(estX_relax_laplace$U)
dim(xmat)
dim(estX_relax_laplace$X_whiten)

###################################################################################################
# lambda = sqrt(2)/2, nu = 0.1
###################################################################################################

blue01=c()
red01=c()
lam=sqrt(2)/2
nu=0.1
for (i in 1:10) {
  estX_relax_laplace = refinedICA(xData = xmat, n.comp = Q,whiten="eigenvec",orth.method="svd", restarts.pbyd = 5, method = "relax_laplace", lambda = sqrt(2)/2, nu = 0.1, converge_plot = F)
  X=estX_relax_laplace$X_whiten
  U=estX_relax_laplace$U
  V=estX_relax_laplace$estS
  
  #blue1=1/2*(norm(X-V%*%t(U),type = "F"))^2
  blue1=1/2*(norm(V-X%*%U,type = "F"))^2 #whitened
  red1=nu/lam*sum(abs(V))
  #red1=nu/2*lam*sum(exp(-abs(V)/lam))
  
  blue01=c(blue01,blue1)
  red01=c(red01,red1)
}
blue01
red01

###################################################################################################
# lambda = sqrt(2)/2, nu = 0.5
###################################################################################################
blue05=c()
red05=c()
lam=sqrt(2)/2
nu=0.5
for (i in 1:10) {
  estX_relax_laplace = refinedICA(xData = xmat, n.comp = Q,whiten="eigenvec",orth.method="svd", restarts.pbyd = 5, method = "relax_laplace", lambda = sqrt(2)/2, nu = 0.5, converge_plot = F)
  X=estX_relax_laplace$X_whiten
  U=estX_relax_laplace$U
  V=estX_relax_laplace$estS
  
  #blue1=1/2*(norm(X-V%*%t(U),type = "F"))^2
  blue1=1/2*(norm(V-X%*%U,type = "F"))^2
  red1=nu/lam*sum(abs(V))
  #red1=nu/2*lam*sum(exp(-abs(V)/lam))
  
  blue05=c(blue05,blue1)
  red05=c(red05,red1)
}
blue05
red05


###################################################################################################
# lambda = sqrt(2)/2, nu = 1
###################################################################################################
blue10=c()
red10=c()
lam=sqrt(2)/2
nu=1
for (i in 1:10) {
  estX_relax_laplace = refinedICA(xData = xmat, n.comp = Q,whiten="eigenvec",orth.method="svd", restarts.pbyd = 5, method = "relax_laplace", lambda = sqrt(2)/2, nu = 1, converge_plot = F)
  X=estX_relax_laplace$X_whiten
  U=estX_relax_laplace$U
  V=estX_relax_laplace$estS
  
  #blue1=1/2*(norm(X-V%*%t(U),type = "F"))^2
  blue1=1/2*(norm(V-X%*%U,type = "F"))^2
  red1=nu/lam*sum(abs(V))
  #red1=nu/2*lam*sum(exp(-abs(V)/lam))
  
  blue10=c(blue10,blue1)
  red10=c(red10,red1)
}
blue10
red10

