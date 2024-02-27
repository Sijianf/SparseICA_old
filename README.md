# Sparse Independent Component Analysis
**Sparse ICA** is a statistical framework to perform independent component analysis with sparse results. 

## 1. Requirements
The following packages are required in **R (version 4.1.3 recommended)**:   
1. **MASS**
2. **neuRosim**
3. **steadyICA**
4. **ProDenICA**
5. **clue**
6. **splines**
7. **foreach**
8. **gam**
9. **irlba**
10. **evd**
```
install.packages(c("MASS","neuRosim","steadyICA","ProDenICA","clue","splines","foreach","gam","irlba","evd"))
```

## 2. Tutorial
The `0_SparseICA.R` contains the main function `sparseICA()` of our Sparse ICA algorithm.    
The `0_Cross_Validation_SparseICA.R` contains the cross validation algorithm `CV_sparseICA()` for the tunning parameter nu.     
The `0_jngcaFunctions.R` contains functions supporting ICA.

### Explanation of Arguments  
#### 1. sparseICA() function
```
sparseICA(xData, n.comp , whiten = c('eigenvec','sqrtprec','none'), orth.method=c('svd','givens'), restarts.pbyd = 10, method = c( "sparse_laplace","sparse_logistic", "fast_logistic", "fast_tanh"), lambda = sqrt(2)/2, nu = 0.1, irlba = FALSE, eps = 1e-06, loop_laplace = 500, loop_logistic = 500, show_message=T, converge_plot = F)
```
- `xData`: Input data matrix with dimension p x n. p is the number of features. n is the number of samples.
- `n.comp`: The number of components. 
- `whiten`: The method for whitening input xData. Could take `eigenvec`,`sqrtprec`, and `none`. The default is `eigenvec`.
- `orth.method `: The method used for generating initial values of U matrix. The default is `svd`.
- `restarts.pbyd`: The number of initial points.
- `method`: ICA method. The default is `sparse_laplace`: the SparseICA with Laplace density. `sparse_logistic`: the SparseICA with logistic density. `fast_logistic`: the FastICA with logistic function. `fast_tanh`: the FastICA with tanh function.
- `lambda`: The scale parameter in Laplace density. The default is sqrt(2)/2 to make the default situation with unit variance.
- `nu`: the tunning parameter controlling the accuracy and sparsity of the results. Should be selected by the cross-validation algorithm `CV_sparseICA()` or expert knowledge. The default is 0.1.
- `irlba`: Whether to use the `irlba` method to perform fast truncated singular value decompositions in whitening step. The default is FALSE.
- `eps`: The convergence threshold. The default is 1e-6.
- `loop_laplace`: The maximum numer of iterations in the Sparese ICA method with laplace density. The default number is 500.
- `loop_logistic`: The maximum numer of iterations in the Sparese ICA method with logistic density. The default number is 500.
- `show_message`: Whether to print the information about convergence. The default is TRUE.
- `converge_plot`: Whether to make a convergence plot for `sparse_laplace` or `sparse_logistic` methods. The default is FALSE.

#### 2. CV_sparseICA() function
```
CV_sparseICA(xData,n.comp,fold=5,CV_method=c("projection","fastICA_logistic","fastICA_tanh"),restarts.pbyd = 10, method = c("sparse_laplace","sparse_logistic"), lambda = sqrt(2)/2, eps = 1e-5,loop_laplace = 500, loop_logistic = 500, nu_list=c(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2),show_message=T,make_plot=T){
```
- `xData`: Input data matrix with dimension p x n. p is the number of features. n is the number of samples.
- `n.comp`: The number of components. 
- `fold`: The number of folds in cross-validation. Should be a resonable integer divisible by the sample size. The default is 5.
- `CV_method`: The method to calculate the CV score. The default is `projection`, which uses prjection matrix to calculate CV score. Other methods include `fastICA_logistic` and `fastICA_tanh`, which compare the sparseICA results with FastICA results.
- `restarts.pbyd`: The number of initial points.
- `method`: The ICA method. The default is `sparse_laplace`: the SparseICA with Laplace density. `sparse_logistic`: the SparseICA with logistic density. 
- `lambda`: The scale parameter in Laplace density. The default is sqrt(2)/2 to make the default situation with unit variance.
- `eps`: The convergence threshold. The default is 1e-6.
- `loop_laplace`: The maximum numer of iterations in the Sparese ICA method with laplace density. The default number is 500.
- `loop_logistic`: The maximum numer of iterations in the Sparese ICA method with logistic density. The default number is 500.
- `nu_list`: The vector of possible values of the tunning parameter nu. The default is `(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2)`.
- `show_message`: Whether to print the information about convergence. The default is TRUE.
- `make_plot`: Whether to make a plot for CV score. The default is TRUE.

### Explanation of Output
#### 1. sparseICA() function
The output will be a list with 7 components as such:
- `loglik`: The log-likelihood.
- `estV`: The estimated V matrix.
- `estU`: The estimated U matrix.
- `xData`: The original input data matrix.
- `converge`: The mean squared norm quantifying the convergence of the V matrix.
- `distribution`: The density used in ICA method.
- `whitener`: The whitener matrix used to perform data whitening.

#### 2. CV_sparseICA() function
The outputs are the best nu selected from cross-validation and a plot for CV score.

## 3. Example
- First, import related functions.
```
source('RCode/0_jngcaFunctions.R')
source("RCode/0_SparseICA.R")
source("RCode/0_Cross_Validation_SparseICA.R")
```

- Directed read the simulated observed data matrix `Xmat.csv` and the true signal data matrix `Smat.csv`. They were simulated using `SimFMRI123(noisyICA=F, nTR=50,var.inactive = 0, snr=0.2)` from `0_jngcaFunctions.R`.
```
xmat = read.csv("Data/Xmat.csv",header = F)
smat = read.csv("Data/Smat.csv",header = F)
```

- Check the data with true signals.
```
par(mfrow=c(1,3))
image(matrix(smat[,1],33))
image(matrix(smat[,2],33))
image(matrix(smat[,3],33))
par(mfrow=c(1,1))
```

![sim_true](https://user-images.githubusercontent.com/43104137/190882927-21ee4e0f-a4eb-4f67-b282-1514f1c8cc45.png)

- Check the whole simulated data.
```
for (i in 1:50) {
  image(matrix(xmat[,i],33))
  Sys.sleep(0.5)
  cat(i,"\n")
}
```
Examine any time point
```
image(matrix(xmat[,35],33))
```

![anytime_35](https://user-images.githubusercontent.com/43104137/190882987-04da2886-a9aa-4d13-aad2-619713d9a70f.png)

- Perform cross-validation to find the best tunning parameter nu for Sparse ICA method with Laplace density 
```
best_nu=CV_sparseICA(xmat,n.comp=3,fold=5,CV_method="projection",restarts.pbyd = 10, method = "sparse_laplace", lambda = sqrt(2)/2, eps = 1e-6,loop_laplace = 500, loop_logistic = 500,nu_list=c(0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1,1.5,2),show_message=T,make_plot=T)
```
![cv_plot](https://user-images.githubusercontent.com/43104137/190882993-9983e55f-f3f1-4698-84c3-8c3e02b5dfae.png)

The best nu selected here is 0.07.


- Fit SparseICA using the best nu selected from CV
```
estX_sparse_laplace = sparseICA(xData = xmat, n.comp = 3, restarts.pbyd = 10, method = "sparse_laplace", lambda = sqrt(2)/2, nu = best_nu,eps = 1e-6,loop_laplace = 500, loop_logistic = 500, converge_plot = F)
```

- Examine recovered sparse components.
```
par(mfrow=c(1,3))
image(matrix(estX_sparse_laplace$estV[,1],33,33))
image(matrix(estX_sparse_laplace$estV[,2],33,33))
image(matrix(estX_sparse_laplace$estV[,3],33,33))
par(mfrow=c(1,1))
```
![sparse_results](https://user-images.githubusercontent.com/43104137/190883041-1f1c228b-39c9-4890-8583-7ab0609e7a51.png)



