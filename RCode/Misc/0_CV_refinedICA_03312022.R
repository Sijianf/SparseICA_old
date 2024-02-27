
CV_refinedICA = function(data,Q,fold=5,CV_method=c("projection","fastICA_logistic","fastICA_tanh"),restarts.pbyd = 10, method = "relax_laplace", lambda = sqrt(2)/2, eps = 1e-5,loop_laplace = 500, loop_logistic = 500,nu_list=c(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2),show_message=T,make_plot=T){
  start_time = Sys.time()
  
  data = as.matrix(data)
  CV_method = match.arg(CV_method)
  
  # Sample size (Number of columns of X)
  n = ncol(data)
  # Number of equal-sized groups (folds)
  K = fold
  # Number of components
  Q = Q
  
  p = ncol(data)
  d = Q
  
  # # Data whiten
  # xData <- scale(data, center=TRUE, scale=FALSE)
  # if (d > p) stop('d must be less than or equal to p')
  # temp = whitener(X = data,n.comp = p,irlba=FALSE)
  # data = temp$Z
  # rm(temp)
  
  ## Divide columns into K equal-sized groups
  index=sample(1:n,size=n)
  groups=matrix(index,nrow = K)
  
  if(CV_method=="projection"){
    # CV -version2
    CV_score2 = rep(NA,length(nu_list))
    for (i in 1:length(nu_list)) {
      score_k = rep(NA,K)
      nu = nu_list[i]
      for (j in 1:K) {
        X_negk = data[,c(groups[-j,])]
        V_negK_relax_laplace = refinedICA(xData = X_negk, n.comp = Q, whiten = "eigenvec", restarts.pbyd = restarts.pbyd, method = method, lambda = lambda, nu = nu, eps = eps,loop_laplace = loop_laplace, loop_logistic = loop_laplace, show_message=show_message,converge_plot = F)
        V_negk = V_negK_relax_laplace$estV
        X_k = data[,c(groups[j,])]
        score1 = sum((X_k-V_negk%*%solve(crossprod(V_negk), t(V_negk))%*%X_k)^2)
        score_k[j] = score1
      }
      CV_score2[i]=mean(score_k)
    }
    out = nu_list[which(CV_score2==min(CV_score2))]
    cat("The best nu is",out,".\n")
    
    if (make_plot) {
      plot(nu_list,CV_score2,xlab = "nu",type = "l")
    }
  }
  
  if(CV_method=="fastICA_logistic"){
    # CV - version 3 - FastICA-Logistic
    CV_score3 = rep(NA,length(nu_list))
    for (i in 1:length(nu_list)) {
      sum1 = rep(NA,K)
      nu = nu_list[i]
      for (j in 1:K) {
        V_negK_relax_laplace = refinedICA(xData = data[,c(groups[-j,])], n.comp = Q, whiten = "none", restarts.pbyd = restarts.pbyd, method = method, lambda = lambda, nu = nu,eps = eps,loop_laplace = loop_laplace, loop_logistic = loop_logistic, show_message=show_message,converge_plot = F)
        V_k_fastICA = refinedICA(xData = data[,groups[j,]], n.comp = Q, whiten = "none", restarts.pbyd = restarts.pbyd, method = "FastICA_logistic")
        sum1[j] = frobICA(S1=V_negK_relax_laplace$estV,S2=V_k_fastICA$estV,standardize=TRUE)
      }
      CV_score3[i]=mean(sum1)
    }
    out = nu_list[which(CV_score3==min(CV_score3))]
    cat("The best nu is",out,".\n")
    
    if (make_plot) {
      plot(nu_list,CV_score3,xlab = "nu",type = "l")
    }
  }
  
  if(CV_method=="fastICA_tanh"){
    # CV - version 3 - FastICA-tanh
    CV_score4 = rep(NA,length(nu_list))
    for (i in 1:length(nu_list)) {
      sum1 = rep(NA,K)
      nu = nu_list[i]
      for (j in 1:K) {
        V_negK_relax_laplace = refinedICA(xData = data[,c(groups[-j,])], n.comp = Q, whiten = "none", restarts.pbyd = restarts.pbyd, method = method, lambda = lambda, nu = nu,eps = eps,loop_laplace = loop_laplace, loop_logistic = loop_logistic, show_message=show_message,converge_plot = F)
        V_k_fastICA = refinedICA(xData = data[,groups[j,]], n.comp = Q, whiten = "none", restarts.pbyd = restarts.pbyd, method = "FastICA_tanh")
        sum1[j] = frobICA(S1=V_negK_relax_laplace$estV,S2=V_k_fastICA$estV,standardize=TRUE)
      }
      CV_score4[i]=mean(sum1)
    }
    out = nu_list[which(CV_score4==min(CV_score4))]
    cat("The best nu is",out,".\n")
    
    if (make_plot) {
      plot(nu_list,CV_score4,xlab = "nu",type = "l")
    }
  }
  
  end_time = Sys.time()
  my_time=end_time - start_time
  cat("Running time for CV:",my_time)
  
  return(out)
}

