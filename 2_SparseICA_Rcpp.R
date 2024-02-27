##########################################################################################
################### Functions for sparse ICA method (Rcpp version) #######################
############################### Last Updates: 12/19/2022 #################################
##########################################################################################
require(irlba)
require(Rcpp)
require(RcppArmadillo)
sourceCpp("Rcpp_func_sparseICA.cpp")

##########################################################################################
##################### Function for generating random starting points #####################
##########################################################################################
gen.inits <- function(p,d,runs,orth.method=c('svd','givens')) {
  orth.method=match.arg(orth.method)
  W.list = list()
  for(i in 1:runs) {
    if(orth.method=='givens') {
      W.list[[i]] <- as.matrix(theta2W(runif(n=choose(p,2),min=0,max=2*pi)))[,1:d]
    } else {
      temp = matrix(rnorm(p*d),p,d)
      W.list[[i]] <- svd(temp)$u
    }
  }
  W.list
}

##########################################################################################
#################################### Whitening Function ##################################
##########################################################################################
whitener <- function(X,n.comp=ncol(X),center.row=FALSE,irlba=FALSE) {
  require(MASS)
  #X must be n x d
  if(ncol(X)>nrow(X)) warning('X is whitened with respect to columns')
  #Corresponds to model of the form X.center = S A, where S are orthogonal with covariance = identity.
  x.center=scale(X,center=TRUE,scale=FALSE)
  if(center.row==TRUE) x.center = x.center - rowMeans(x.center)
  n.rep=dim(x.center)[1]
  if(irlba==FALSE) svd.x=svd(x.center,nu=n.comp,nv=n.comp)
  if(irlba==TRUE) {
    require(irlba)
    svd.x=irlba(x.center,nu=n.comp,nv=n.comp)
  }
  #RETURNS PARAMETERIZATION AS IN fastICA (i.e., X is n x d)
  #NOTE: For whitened X, re-whitening leads to different X
  #The output for square A is equivalent to solve(K)
  return(list(whitener=t(ginv(svd.x$v%*%diag(svd.x$d[1:n.comp])/sqrt(n.rep-1))),Z=sqrt(n.rep-1)*svd.x$u,mean=apply(X,2,mean)))
}

##########################################################################################
############################# Sparse LNGCA Rcpp version ##################################
##########################################################################################

sparseICA_Rcpp = function(xData,n.comp,nu = 1,U.list=NULL,whiten = c('eigenvec','sqrtprec','none'), orth.method=c('svd','givens'), restarts.pbyd = 10, lambda = sqrt(2)/2, irlba = FALSE, eps = 0.000001, maxit.laplace = 500, show_message=T, converge_plot = F){
  
  start.time = Sys.time()

  ##########################################################################################
  ################ soft_thresh to update matrix in relax and split method ##################
  ##########################################################################################
  # Laplace distribution is used, the default nu is 1
  # Has been moved to Rcpp 

  
  ##################################################################################### 0
  xData = as.matrix(xData) # make sure the input data are in matrix format
  whiten=match.arg(whiten)
  orth.method= match.arg(orth.method)
  
  if(irlba) require(irlba)
  p = ncol(xData) # Extract the dimensions to initialize the V(0) for "relax_logistic" and "relax_soft" methods. p = number of column
  d = n.comp # Extract the components need for V(0)
  
  ####################################################################################
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
  
  ###########################################################################################
  # Randomly make an input for U(0), by default here used orth.method = "svd"
  W.list = U.list
  if (is.null(W.list)){
    W.list = gen.inits(p=p, d=d, runs = restarts.pbyd, orth.method=orth.method)
  } else{
    W.list = list(U.list)
  }
  runs = length(W.list)
  
  ###########################################################################################
  # Create a NULL list for storing outcomes
  out.list = NULL
  
  converge = c()
  #converge_plot = converge_plot
  
  # Store log likelihood
  loglik = c()
  
  final_Rcpp=runs_relax_laplace(xData,W.list,runs,n.comp,nu,lambda,maxit.laplace,eps)
  
  if(show_message){
    cat("MESSAGE: The algorithm for relax laplase converges (<",eps,") in",maxit.laplace,"iterations within",runs,"different start values!\n")
    #cat("Consider increasing the number of loop_laplase.\n")
  }
  
    out.list$loglik = final_Rcpp$loglik
    out.list$estS = final_Rcpp$newV
    out.list$estM = final_Rcpp$newU
    out.list$xData = xData
    out.list$converge = final_Rcpp$converge
    out.list$distribution = "Laplace"
    out.list$whitener = whitener
    
    if (converge_plot == TRUE) {
      converge = out.list$converge
      plot(converge, main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
    }
    

  
  # Print running time
  end.time = Sys.time()
  time.taken = end.time - start.time
  cat("Running time:",time.taken,"s.\n")
  
  return(out.list)
  
} # end of sparseICA function


##########################################################################################
########################## Sparse ICA (PCA+ICA) Rcpp version #############################
##########################################################################################

PCA_sparseICA_Rcpp = function(xData,n.comp,nu = 1,U.list=NULL,whiten = c('eigenvec','sqrtprec','none'), orth.method=c('svd','givens'), restarts.pbyd = 10, lambda = sqrt(2)/2, eps = 1e-06, maxit.laplace = 500, show_message=T, converge_plot = F){
  
  start.time = Sys.time()
  
  ##########################################################################################
  ################ soft_thresh to update matrix in relax and split method ##################
  ##########################################################################################
  # Laplace distribution is used, the default nu is 1
  # Has been moved to Rcpp module
  
  ##################################################################################### 0
  xData = as.matrix(xData) # make sure the input data are in matrix format
  whiten=match.arg(whiten)
  orth.method= match.arg(orth.method)
  
  # do PCA
  # perform PCA, center+scale included
  subj_PCA=prcomp_irlba(xData,n.comp,center = T,scale. = T)
  myPC=subj_PCA$x
  dimnames(myPC)=NULL
  
  p = ncol(myPC) # Extract the dimensions to initialize the V(0) for "relax_logistic" and "relax_soft" methods. p = number of column
  d = n.comp # Extract the components need for V(0)
  
  # Data whiten
  myPC = scale(myPC, center=TRUE, scale=FALSE)
  if (d > p) stop('d must be less than or equal to p')
  if (whiten=='eigenvec') {
    # Use whitener=='eigenvec' so that restarts.dbyd initiates from the span of the first d eigenvectors.
    temp = whitener(X = myPC,n.comp = p)
    myPC = temp$Z
    whitener = temp$whitener
    rm(temp)
  }  else if (whiten=='sqrtprec') {
    est.sigma = cov(myPC)  ## Use eigenvalue decomposition rather than SVD.
    evd.sigma = svd(est.sigma)
    whitener = evd.sigma$u%*%diag(evd.sigma$d^(-1/2))%*%t(evd.sigma$u)
    myPC = myPC%*%whitener
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
  
  converge = c()
  #converge_plot = converge_plot
  
  # Store log likelihood
  loglik = c()
  
  final_Rcpp=runs_relax_laplace(myPC,W.list,runs,n.comp,nu,lambda,maxit.laplace,eps)
  
  if(show_message){
    cat("MESSAGE: The algorithm for relax laplase converges (<",eps,") in",maxit.laplace,"iterations within",runs,"different start values!\n")
    #cat("Consider increasing the number of loop_laplase.\n")
  }
  
  out.list$loglik = final_Rcpp$loglik
  out.list$estS = final_Rcpp$newV
  out.list$estM = final_Rcpp$newU
  out.list$xData = xData
  out.list$converge = final_Rcpp$converge
  out.list$distribution = "Laplace"
  out.list$whitener = whitener
  
  if (converge_plot == TRUE) {
    converge = out.list$converge
    plot(converge, main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
  }
  
  # Print running time
  end.time = Sys.time()
  time.taken = end.time - start.time
  cat("Running time:",time.taken,"s.\n")
  
  return(out.list)
  
} # end of PCA + sparseICA function


##########################################################################################
############### Use BIC to select tunning parameter on Sparse LNGCA ######################
##########################################################################################

BIC_sparseICA_Rcpp = function(xData,n.comp,nu_list = seq(0.1,4,0.01),U.list=NULL,whiten = c('eigenvec','sqrtprec','none'), orth.method=c('svd','givens'), restarts.pbyd = 10, lambda = sqrt(2)/2, irlba = FALSE, eps = 1e-06, maxit.laplace = 500, show_message=T, BIC_plot = F){
  
  start.time = Sys.time()
  
  ##############################################################################
  xData = as.matrix(xData) # make sure the input data are in matrix format
  whiten=match.arg(whiten)
  orth.method= match.arg(orth.method)
  
  ##########################################################################################
  # reference Sparse ICA
  best_nu = 0.0000000001
  ref_sparseICA = sparseICA_Rcpp(xData = xmat, n.comp = n.comp,whiten = whiten, restarts.pbyd = restarts.pbyd, lambda = sqrt(2)/2, nu = best_nu,eps = eps,maxit.laplace = 500, show_message = show_message, converge_plot = F)
  
  # center + scale raw data matrix
  xmat_centered = scale(xmat,center = T,scale = T)
  
  v = dim(xData)[1]   # the number of voxels
  t = dim(xData)[2]   # the number of time points
  my_U_list = ref_sparseICA$estM   # initial U for next start
  
  ##############################################################################
  # start BIC process
  out_BIC = rep(NA,length(nu_list))
  out_ratio = rep(NA,length(nu_list))
  out_sparsity = rep(NA,length(nu_list))
  
  for (i in 1:length(nu_list)) {
    my_nu = nu_list[i]
    temp_sparseICA = sparseICA_Rcpp(xData = xmat,n.comp = n.comp,U.list = my_U_list,whiten = whiten, restarts.pbyd = restarts.pbyd, lambda = sqrt(2)/2, nu = my_nu,eps = eps,maxit.laplace = 500, show_message = show_message, converge_plot = F)
    S_nu = temp_sparseICA$estS
    temp_mat = solve(crossprod(S_nu), t(S_nu))%*%xmat_centered
    e_nu = sum((xmat_centered-S_nu%*%temp_mat)^2) # Calculate e_nu
    
    my_U_list = temp_sparseICA$estM   # Warm start for next iteration
    
    out_ratio[i] = e_nu/(v*t)
    out_sparsity[i] = sum(temp_sparseICA$estS!=0)*log(v*t)/(v*t)
    out_BIC[i] =e_nu/(v*t)+sum(temp_sparseICA$estS!=0)*log(v*t)/(v*t)
  }
  
  if (BIC_plot == TRUE) {
    plot(nu_list,out_BIC, main = "BIC Plot", xlab = "nu", ylab = "BIC", type = "l")
  }
  
  best_nu = nu_list[which(out_BIC==min(out_BIC))]
  cat("The best nu selected by BIC is ",best_nu,".\n")
  
  out.list=NULL
  out.list$ratio = out_ratio
  out.list$sparsity = out_sparsity
  out.list$BIC = out_BIC
  out.list$best_nu = best_nu
  
  # Print running time
  end.time = Sys.time()
  time.taken = end.time - start.time
  cat("Running time for BIC selection:",time.taken,"s.\n")
  
  return(out.list)
}

##########################################################################################
############### Use BIC to select tunning parameter on PCA+SparseICA #####################
##########################################################################################

BIC_PCA_sparseICA_Rcpp = function(xData,n.comp,nu_list = seq(0.1,4,0.01),U.list=NULL,whiten = c('eigenvec','sqrtprec','none'), orth.method=c('svd','givens'), restarts.pbyd = 10, lambda = sqrt(2)/2, irlba = FALSE, eps = 1e-06, maxit.laplace = 500, show_message=T, BIC_plot = F){
  start.time = Sys.time()
  
  library(irlba)
  #library(singR)
  
  ##############################################################################
  xData = as.matrix(xData) # make sure the input data are in matrix format
  whiten=match.arg(whiten)
  orth.method= match.arg(orth.method)
  
  ##########################################################################################
  # perform PCA, center+scale included
  subj_PCA=prcomp_irlba(xData,n.comp,center = T,scale. = T)
  myPC=subj_PCA$x
  dimnames(myPC)=NULL
  
  # reference Sparse ICA
  best_nu = 0.0000000001
  ref_sparseICA = sparseICA_Rcpp(xData = myPC, n.comp = n.comp,whiten = whiten, restarts.pbyd = restarts.pbyd,lambda = sqrt(2)/2, nu = best_nu,eps = eps,maxit.laplace = 500, show_message = show_message, converge_plot = F)
  
  # center + scale raw data matrix
  xmat_centered = scale(xmat,center = T,scale = T)
  
  v = dim(xData)[1]   # the number of voxels
  t = dim(xData)[2]   # the number of time points
  my_U_list = ref_sparseICA$estM   # initial U for next start
  
  ##############################################################################
  # start BIC process
  out_BIC = rep(NA,length(nu_list))
  out_ratio = rep(NA,length(nu_list))
  out_sparsity = rep(NA,length(nu_list))
  
  for (i in 1:length(nu_list)) {
    my_nu = nu_list[i]
    temp_sparseICA = sparseICA_Rcpp(xData = myPC,n.comp = n.comp,U.list = my_U_list,whiten = whiten, restarts.pbyd = restarts.pbyd, lambda = sqrt(2)/2, nu = my_nu,eps = eps,maxit.laplace = 500, show_message = show_message, converge_plot = F)
    S_nu = temp_sparseICA$estS
    temp_mat = solve(crossprod(S_nu), t(S_nu))%*%xmat_centered
    e_nu = sum((xmat_centered-S_nu%*%temp_mat)^2) # Calculate e_nu
    
    my_U_list = temp_sparseICA$estM   # Warm start for next iteration
    
    out_ratio[i] = e_nu/(v*t)
    out_sparsity[i] = sum(temp_sparseICA$estS!=0)*log(v*t)/(v*t)
    out_BIC[i] =e_nu/(v*t)+sum(temp_sparseICA$estS!=0)*log(v*t)/(v*t)
  }
  
  if (BIC_plot == TRUE) {
    plot(nu_list,out_BIC, main = "BIC Plot", xlab = "nu", ylab = "BIC", type = "l")
  }
  
  best_nu = nu_list[which(out_BIC==min(out_BIC))]
  cat("The best nu selected by BIC is ",best_nu,".\n")
  
  out.list=NULL
  out.list$ratio = out_ratio
  out.list$sparsity = out_sparsity
  out.list$BIC = out_BIC
  out.list$best_nu = best_nu
  
  # Print running time
  end.time = Sys.time()
  time.taken = end.time - start.time
  cat("Running time for BIC selection:",time.taken,"s.\n")
  
  return(out.list)
}

