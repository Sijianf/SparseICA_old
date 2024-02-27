# Copy wrapper for Procrustes
procrustesR <- function(X, V){
  svdXV <- svd(X %*% V)
  U <- tcrossprod(svdXV$u, svdXV$v)
  return(U)
}

# Wrapper for Sparse PCA
sparsePCAR <- function(X, Vstart, lambda, tol){
  # Evaluate U for given Vstart, and the value of objective
  U = procrustesR(X, Vstart)
  # Make V same as Vstart
  V = Vstart
  # objective function value at start
  obj_old = sum((X - tcrossprod(U, V))^2)/2 + lambda * sum(abs(V))
  # While not converged, repeat
  error = 100
  while (error > tol){
    ## Update V via soft-thresholding of X'U
    XtU = crossprod(X, U) # this is X'U
    V = sign(XtU) * pmax(abs(XtU) - lambda, 0) # soft-thresholding
    ## Update U via Procrustes
    U = procrustesR(X, V)
    # New objective function value after update
    obj_new = sum((X - tcrossprod(U, V))^2)/2 + lambda * sum(abs(V))
    # Difference in objectives as error to monitor convergence
    error = abs(obj_new - obj_old)
    # Make new old
    obj_old = obj_new
  }
  # Return a list of U, V and error on solution
  return(list(U = U, V = V, error = error))
}

set.seed(308723)
X <- matrix(rnorm(110), 11, 5)
V <- matrix(rnorm(30), 5, 3)
lambda = 1
eps = 1e-2
outR = sparsePCAR(X, V, lambda, eps)
outR$V