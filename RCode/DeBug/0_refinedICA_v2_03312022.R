##########################################################################################
########################### Function for refined ICA method ##############################
##########################################################################################
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
refinedICA = function(xData = xmat, n.comp = 3, W.list = NULL, whiten = c('eigenvec','sqrtprec','none'), orth.method=c('svd','givens'), restarts.pbyd = 10, method = c( "relax_laplace","relax_logistic", "FastICA_logistic", "FastICA_tanh"), lambda = sqrt(2)/2, nu = 0.1, irlba = FALSE, eps = 1e-05, loop_laplace = 500, loop_logistic = 500,show_message=T, converge_plot = F){
  
  start.time <- Sys.time()

  ##########################################################################################
  ################ soft_thresh to update matrix in relax and split method ##################
  ##########################################################################################
  # Laplace distribution is used, the default nu is 1
  soft_thresh <- function(x, nu = 1, lambda = sqrt(2)/2) {
    xmin <- pmax(abs(x)-nu/lambda,0)*sign(x)
    return(xmin)
  }
  
  ##########################################################################################
  ####### Newton's method to find minimum/maximum (root for the first derivative) ##########
  ##########################################################################################
  # Updated 2020.04.27 Did not change the code, but the original objective function is equivalent to this one under pre-whitened data
  # Updated 2020.04.01 to change back the equation (the original objective function should be correct)
  # Updated 2020.03.31 to change the equation
  prox_newton <- function(v, nu = 1, itm = 50, tol= 1e-6) {
    
    # based on the function of logistic: G(y) = ln(1+exp(a*y))
    
    # initialization
    z = v
    a = -pi/sqrt(3.0)
    n = 0
    e = exp(a*v)
    g = z - v + 2*a*nu*e/(1.0 + e) - a*nu # it is the function that we need to find its root
    h = 1.0 + 4*a^2*nu*e/(1.0 + e)^2 # it is the first derivative of g
    # Additional note: the logistic density problem in Peng's ppt is to find the minimum of nu*log(1.0 + e) + 0.5*(z - v)^2
    
    while(abs(sum(g)) >= tol) {
      z = z - g/h # the Newton's Method equation
      
      # update gradient hessian
      n = n + 1;
      e = exp(a*z);
      g = z - v + 2*a*nu*e/(1.0 + e) - a*nu
      h = 1.0 + 4*a^2*nu*e/(1.0 + e)^2
      
      if(n >= itm) {
        message(paste0("Already achieved maximum iteration=",itm))
        break} # force to quit when more than the maximum iteration
      
    } # end of while loop
    
    return(z)
    
  } # end of prox_newton function
  
  ##################################################################################### 0
  xData = as.matrix(xData) # make sure the input data are in matrix format
  whiten=match.arg(whiten)
  orth.method= match.arg(orth.method)
  
  if(irlba) require(irlba)
  p = ncol(xData) # Extract the dimensions to initialize the V(0) for "relax_logistic" and "relax_soft" methods. p = number of column
  d = n.comp # Extract the components need for V(0)
  
  # Data whiten
  xData <- scale(xData, center=TRUE, scale=FALSE)
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
  }
  else {
    whitener = diag(p)
  }

  # Randomly make an input for U(0), by default here used orth.method = "svd"
  W.list = gen.inits(p=p, d=d, runs = restarts.pbyd, orth.method=orth.method) 
  runs = length(W.list)
  
  # Create a NULL list for storing outcomes
  out.list = NULL
  
  converge = c()
  converge_plot = converge_plot
  
  # Store log likelihood
  loglik = c()

  ################################################################################## 1
  if (method == "relax_laplace") {
    m=0
    for (k in 1:runs) {  # make a list if restarts.pbyd > 1, we should update the index of W.list if restarts.pbyd > 1
      # Initial value of V(0)
      newV = xData %*% W.list[[k]]  
      lagV = newV
      
      # loop for relax_laplace
      for (i in 1:loop_laplace) { 
        iteration = i
        # update U:
        txv = t(xData)%*%newV
        svd.txv = La.svd(txv)
        newU = svd.txv$u%*%svd.txv$vt # define objective function for orthogonal U
        # update V:
        newV = xData%*%newU
        newV = soft_thresh(newV, nu = nu, lambda = lambda)
        # if end after thresholding, ensures that the solution has exact zeros. 
        converge[iteration] = mean((lagV-newV)^2)/ncol(lagV) # Mean squared norm
        lagV = newV
        if (converge[iteration] < eps) {break}
      } 
      # End the loop for relax_laplace
      
      # Store Warning
      if(iteration==loop_laplace){
        #cat("Does not converage!\n")
        m=m+1
      }
      
      # Store the results
      if (k != 1) {out.list[[k]] = out.list[[1]]}
      out.list[[k]]$loglik = sum(-abs(lagV)/lambda-log(2*lambda)) - 1/(2*nu)*(norm(newV-xData%*%newU,type = "F"))^2
      loglik[k] = out.list[[k]]$loglik
      out.list[[k]]$estV = newV
      out.list[[k]]$estU = newU
      out.list[[k]]$xData = xData
      out.list[[k]]$converge = converge
      
    } #----------------------------------------------------------- end the loop to make a list
    
    if(show_message){
      cat("WARNING: The algorithm for relax laplase does not converge (<",eps,") in",loop_laplace,"iterations for",m,"/",runs,"different start values!\n")
      cat("Consider increasing the number of loop_laplase.\n")
    }
    
    out.list = out.list[[which.max(loglik)]]
    out.list$distribution = "Laplace"
    out.list$whitener = whitener
    
    if (converge_plot == TRUE) {
      converge = out.list$converge
      plot = recordPlot()
      plot(converge, main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
      out.list$converge = plot
    }
    
  }
  
  ################################################################################# 2
  if (method == "relax_logistic") {
    m=0
    for (k in 1:runs) { # Make a list if restarts.pbyd > 1
      newV = xData %*% W.list[[k]] # To be more accurate, we should update the index of W.list if restarts.pbyd > 1 
      lagV = newV
      # loop for relax_logistic
      for (i in 1:loop_logistic) { 
        iteration = i
        # update U:
        txv = t(xData)%*%newV
        svd.txv = La.svd(txv)
        newU = svd.txv$u%*%svd.txv$vt # define objective function for orthogonal U
        # update V:
        newV = xData%*%newU
        newV = prox_newton(v = newV, nu = nu)# here use the prox_newton function. 
        
        converge[iteration] = norm((lagV-newV)^2)/ncol(lagV) # Mean squared norm
        lagV = newV
        if (converge[iteration] < eps) {break}
      } 
      # End the loop for relax_logistic
      
      # Store Warning
      if(iteration==loop_logistic){
        #cat("Not converge!\n")
        m=m+1
      }
      
      if (k != 1) {out.list[[k]] = out.list[[1]]}
      out.list[[k]]$loglik = sum(-newV/(sqrt(3)/pi) - 2*log(1+exp(-newV/(sqrt(3)/pi)))) - 
        1/(2*nu)*(norm(newV-xData%*%newU,type = "F"))^2
      loglik[k] = out.list[[k]]$loglik
      out.list[[k]]$estV = newV
      out.list[[k]]$estU = newU
      out.list[[k]]$xData = xData
      out.list[[k]]$converge = converge
    } 
    # End the loop to make a list

    if(show_message){
      cat("WARNING: The algorithm for relax logistic does not converge (<",eps,")in",loop_logistic,"iterations for",m,"/",runs,"different start values!\n")
      cat("Consider increasing the number of loop_laplase.\n")
    }
    
    out.list = out.list[[which.max(loglik)]]
    out.list$distribution = "logistic"
    out.list$whitener = whitener
    
    if (converge_plot == TRUE) {
      converge = out.list$converge
      plot = recordPlot()
      plot(converge, main = "Convergence Plot", xlab = "Iteration", ylab = "Norm Difference", type = "o")
      out.list$converge = plot
    }
    
  }
  
  ########################### 3
  if (method == "FastICA_logistic") {
    estX_logis = mlcaFP(xData = xData, n.comp = n.comp, whiten = "none", restarts.pbyd = restarts.pbyd, distribution='logistic', eps = eps)
    out.list$loglik = estX_logis$loglik
    out.list$estV = estX_logis$S
    out.list$distribution = estX_logis$distribution
    out.list$whitener = whitener
    if (converge_plot == T) {
      out.list$converge = "No converge plot for FastICA_logistic and FastICA_tanh"
    }
  }
  
  
  ########################### 4
  if (method == "FastICA_tanh") {
    estX_tanh = mlcaFP(xData = xData, n.comp = n.comp, whiten = "none", restarts.pbyd = restarts.pbyd, distribution='tanh', eps = eps)
    out.list$loglik = estX_tanh$loglik
    out.list$estV = estX_tanh$S
    out.list$distribution = estX_tanh$distribution
    out.list$whitener = whitener
    if (converge_plot == T) {
      out.list$converge = "No converge plot for FastICA_logistic and FastICA_tanh"
    }
  }
  
  # Print running time
  end.time = Sys.time()
  time.taken = end.time - start.time
  cat("Running time:",time.taken,".\n")
  
  return(out.list)
  
} # end of refinedICA function


