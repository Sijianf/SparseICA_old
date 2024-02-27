##########################################################################################
############### Use BIC to select tunning parameter on PCA+SparseICA #####################
##########################################################################################

BIC_PCA_sparseICA = function(xData,n.comp,nu_list = seq(0.1,4,0.01),U.list=NULL,whiten = c('eigenvec','sqrtprec','none'), orth.method=c('svd','givens'), restarts.pbyd = 10, method = c( "sparse_laplace","sparse_logistic", "fast_logistic", "fast_tanh"), lambda = sqrt(2)/2, irlba = FALSE, eps = 1e-06, maxit.laplace = 500, maxit.logistic = 500, show_message=T, BIC_plot = F){
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
  ref_sparseICA = sparseICA(xData = myPC, n.comp = n.comp,whiten = whiten, restarts.pbyd = restarts.pbyd, method = method, lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,maxit.laplace = 500, maxit.logistic = 500,show_message = show_message, converge_plot = F)
  
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
    temp_sparseICA = sparseICA(xData = myPC,n.comp = n.comp,U.list = my_U_list,whiten = whiten, restarts.pbyd = restarts.pbyd, method = method, lambda = sqrt(2)/2, nu = my_nu,eps = 1e-6,maxit.laplace = 500, maxit.logistic = 500,show_message = show_message, converge_plot = F)
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

