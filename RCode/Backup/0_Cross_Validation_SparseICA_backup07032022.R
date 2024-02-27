##########################################################################################
#################### Cross-Validation for the Tuning Parameter nu ########################
###################### Zihang Wang ####### Last update: 5/10/2022 ########################
##########################################################################################

CV_sparseICA = function(xData,n.comp,fold=5,CV_method=c("projection","fastICA_logistic","fastICA_tanh"),restarts.pbyd = 10, method = c("sparse_laplace","sparse_logistic"), lambda = sqrt(2)/2, eps = 1e-5,loop_laplace = 500, loop_logistic = 500,nu_list=c(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2),show_message=T,make_plot=T){
  start_time = Sys.time()
  
  data = as.matrix(xData)
  CV_method = match.arg(CV_method)
  method=match.arg(method)
  
  # Sample size (Number of columns of X)
  n = ncol(data)
  # Number of equal-sized groups (folds)
  K = fold
  # Number of components
  Q = n.comp
  
  p = ncol(data)
  d = Q
  
  # # Data whiten
  # xData <- scale(data, center=TRUE, scale=FALSE)
  # if (d > p) stop('d must be less than or equal to p')
  # temp = whitener(X = data,n.comp = p,irlba=FALSE)
  # data = temp$Z
  # rm(temp)
  
  ## Divide columns into K equal-sized groups
  set.seed(2022)
  index=sample(1:n,size=n)
  groups=matrix(index,nrow = K)
  
  if(CV_method=="projection"){
    # CV -version2
    CV_score = matrix(NA,nrow=K,ncol=length(nu_list))
    for (j in 1:K) {
      X_negk = data[,c(groups[-j,])]
      X_k = data[,c(groups[j,])]
      
      # Whiten X_negk once
      X_negk <- scale(X_negk, center=TRUE, scale=FALSE)
      if (d > p) stop('d must be less than or equal to p')
      temp = whitener(X = X_negk,n.comp = ncol(X_negk),irlba=FALSE)
      X_negk = temp$Z
      rm(temp)
      
      for (i in 1:length(nu_list)) {
        nu = nu_list[i]
        V_negK_relax_laplace = sparseICA(xData = X_negk, n.comp = Q, whiten = "none", restarts.pbyd = restarts.pbyd, method = method, lambda = lambda, nu = nu, eps = eps,loop_laplace = loop_laplace, loop_logistic = loop_laplace, show_message=show_message,converge_plot = F)
        V_negk = V_negK_relax_laplace$estV
        CV_score[j,i] = sum((X_k-V_negk%*%solve(crossprod(V_negk), t(V_negk))%*%X_k)^2)
      }
    }
    CV_score2 = apply(CV_score, 2, mean)
    out = nu_list[which(CV_score2==min(CV_score2))]
    cat("The best nu is",out,".\n")
    
    if (make_plot) {
      plot(nu_list,CV_score2/sum(CV_score2),xlab = "nu",type = "l",ylab="CV Score",main="CV Score Plot")
    }
  }
  
  if(CV_method=="fastICA_logistic"){
    # CV - version 3 - FastICA-Logistic
    CV_score = matrix(NA,nrow=K,ncol=length(nu_list))
    for (j in 1:K) {
      X_negk = data[,c(groups[-j,])]
      X_k = data[,c(groups[j,])]
      
      # Whiten X_negk once
      X_negk <- scale(X_negk, center=TRUE, scale=FALSE)
      if (d > p) stop('d must be less than or equal to p')
      temp = whitener(X = X_negk,n.comp = ncol(X_negk),irlba=FALSE)
      X_negk = temp$Z
      rm(temp)
      
      for (i in 1:length(nu_list)) {
        nu = nu_list[i]
        V_negK_relax_laplace = sparseICA(xData = X_negk, n.comp = Q, whiten = "none", restarts.pbyd = restarts.pbyd, method = method, lambda = lambda, nu = nu, eps = eps,loop_laplace = loop_laplace, loop_logistic = loop_laplace, show_message=show_message,converge_plot = F)
        V_k_fastICA = sparseICA(xData = data[,groups[j,]], n.comp = Q, whiten = "none", restarts.pbyd = restarts.pbyd, method = "fast_logistic")
        CV_score[j,i] = frobICA(S1=V_negK_relax_laplace$estV,S2=V_k_fastICA$estV,standardize=TRUE)
      }
    }
    CV_score3 = apply(CV_score, 2, mean)
    out = nu_list[which(CV_score3==min(CV_score3))]
    cat("The best nu is",out,".\n")
    
    if (make_plot) {
      plot(nu_list,CV_score3/sum(CV_score3),xlab = "nu",type = "l",ylab="CV Score",main="CV Score Plot")
    }
  }
  
  if(CV_method=="fastICA_tanh"){
    # CV - version 3 - FastICA-tanh
    CV_score = matrix(NA,nrow=K,ncol=length(nu_list))
    for (j in 1:K) {
      X_negk = data[,c(groups[-j,])]
      X_k = data[,c(groups[j,])]
      
      # Whiten X_negk once
      X_negk <- scale(X_negk, center=TRUE, scale=FALSE)
      if (d > p) stop('d must be less than or equal to p')
      temp = whitener(X = X_negk,n.comp = ncol(X_negk),irlba=FALSE)
      X_negk = temp$Z
      rm(temp)
      
      for (i in 1:length(nu_list)) {
        nu = nu_list[i]
        V_negK_relax_laplace = sparseICA(xData = X_negk, n.comp = Q, whiten = "none", restarts.pbyd = restarts.pbyd, method = method, lambda = lambda, nu = nu, eps = eps,loop_laplace = loop_laplace, loop_logistic = loop_laplace, show_message=show_message,converge_plot = F)
        V_k_fastICA = sparseICA(xData = data[,groups[j,]], n.comp = Q, whiten = "none", restarts.pbyd = restarts.pbyd, method = "fast_tanh")
        CV_score[j,i] = frobICA(S1=V_negK_relax_laplace$estV,S2=V_k_fastICA$estV,standardize=TRUE)
      }
    }
    CV_score4 = apply(CV_score, 2, mean)
    out = nu_list[which(CV_score4==min(CV_score4))]
    cat("The best nu is",out,".\n")
    
    if (make_plot) {
      plot(nu_list,CV_score4/sum(CV_score4),xlab = "nu",type = "l",ylab="CV Score",main="CV Score Plot")
    }
  }
  
  end_time = Sys.time()
  my_time=end_time - start_time
  cat("Running time for CV:",my_time,"s.\n")
  
  return(out)
}

