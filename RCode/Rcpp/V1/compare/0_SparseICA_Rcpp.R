##########################################################################################
########################### Function for sparse ICA method ##############################
############################### Last Updates: 12/8/2022 ##################################
##########################################################################################

##########################################################################################
# Updates information:
# Updated 2022.12.8  for moving soft thresholding and ICA part to Rcpp modules; keep only Laplace sparse ICA algorithm
# Updated 2022.11.10 for add U.list as user-specified initial values for U, change the names of return values V->S, U->M
# Updated 2022.11.10 for change the function name, method names, information displayed, delete W.list from the arguments
# Updated 2022.03.31 for changing the norm of convergence of V matrix, defining soft_thresh and prox_newton inside
# Updated 2022.03.20 for optimization function, adding parameters to control number of iterations, adding warnings, adding running time
# Updated 2020.04.28 for (1) rename relax-softmax to relax-laplace; 
#                        (2) add the log likelihood value for all the methods; 
#                        (3) allow data whitening with choices of whiten = c('eigenvec','sqrtprec','none'); 
#                        (4) allow different initialization methods for the start up unmixing (W) matrix, orth.method=c('svd','givens');
#                        (5) allow restarts > 1 for relax-laplace and relax-logistic; (FastICAs are originally allowed by mlcaFP);
#                        (6) allow to plots the convergence trend for relax-laplace and relax-logistic, by default it is hidden; 
#                        (7) note: the output format is different from the previous one. 
# Updated 2020.04.06 for convergence rule
# Updated 2020.04.01 for eps
# Updated 2020.03.30 for whiten
##########################################################################################
require(Rcpp)
require(RcppArmadillo)
sourceCpp("Rcpp_func_sparseICA.cpp")

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
  
  converge = c()
  #converge_plot = converge_plot
  
  # Store log likelihood
  loglik = c()

  ################################################################################## 
  for (k in 1:runs) {  # make a list if restarts.pbyd > 1, we should update the index of W.list if restarts.pbyd > 1
    # Initial value of V(0)
    newV = xData %*% W.list[[k]]  
    
    # add cpp
    Rcpp_res=relax_laplace(xData,newV,nu,lambda,maxit.laplace,eps)
    
    # Store the results
    if (k != 1) {out.list[[k]] = out.list[[1]]}
    out.list[[k]]$loglik = sum(-abs(Rcpp_res$newV)/lambda-log(2*lambda)) - 1/(2*nu)*(norm(Rcpp_res$newV-xData%*%Rcpp_res$newU,type = "F"))^2
    loglik[k] = out.list[[k]]$loglik
    out.list[[k]]$estS = Rcpp_res$newV
    out.list[[k]]$estM = Rcpp_res$newU
    out.list[[k]]$xData = xData
    out.list[[k]]$converge = Rcpp_res$converge
    
  } #----------------------------------------------------------- end the loop to make a list
  
  if(show_message){
    cat("MESSAGE: The algorithm for relax laplase converges (<",eps,") in",maxit.laplace,"iterations within",runs,"different start values!\n")
    #cat("Consider increasing the number of loop_laplase.\n")
    
    out.list = out.list[[which.max(loglik)]]
    out.list$distribution = "Laplace"
    out.list$whitener = whitener
    
    if (converge_plot == TRUE) {
      converge = out.list$converge
      plot(converge, main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
    }
    
  }
  
  # Print running time
  end.time = Sys.time()
  time.taken = end.time - start.time
  cat("Running time:",time.taken,"s.\n")
  
  return(out.list)
  
} # end of sparseICA function


