install.packages(c("Rcpp", "RcppArmadillo"))

library(Rcpp)
library(microbenchmark)


cppFunction("
bool isOddCpp(int num){
bool result = (num%2 == 1);
return result;
}")
isOddCpp(10)
isOddCpp(11)

################################################################################
cumul_sumR <- function(x){
  p <- length(x)
  s <- x
  for (i in 2:p){
    s[i] <- s[i-1] + x[i]
  }
  return(s)
}
cumul_sumR(x = c(1, 3, 5))

# should not use exact copy
cppFunction("NumericVector cumul_sumCpp(NumericVector x){
int p = x.size();
NumericVector s = x;
for(int i = 1; i < p; i++){
s[i] = s[i-1] + x[i];
}
return(s);
}")
cumul_sumCpp(x = c(1, 3, 5))

x = c(1, 3, 5)
s = cumul_sumCpp(x)
print(s)
print(x)


cppFunction("NumericVector cumul_sumCpp(NumericVector x){
int p = x.size();
NumericVector s = clone(x); // cloning
for(int i = 1; i < p; i++){
s[i] = s[i-1] + x[i];
}
return(s);}")
x = c(1, 3, 5)
s = cumul_sumCpp(x)
print(s)
print(x)


################################################################################
# Compare
p = 1000
x = rnorm(p)
identical(cumul_sumR(x), cumul_sumCpp(x))

microbenchmark(
  cumul_sumR(x),
  cumul_sumCpp(x), times = 50
)

################################################################################
# Bootstrap example
# ds - vector of observations
# B - number of bootstrap samples
bootstrap_r <- function(ds, B = 1000){
  boot_stat <- matrix(NA, nrow = B, ncol = 2)
  n <- length(ds)
  # Perform bootstrap
  for(i in 1:B) {
    # Create a sample of size n with replacement
    gen_data <- ds[sample(n, n, replace=TRUE)]
    # Calculate sample data mean and SD
    boot_stat[i,] <- c(mean(gen_data),sd(gen_data))
  }
  return(boot_stat)
}


ds <- rnorm(1000, mean = 10, sd = 5)
out <- bootstrap_r(ds)
microbenchmark(bootstrap_r(ds), times = 10)

hist(out[, 1])
hist(out[, 2])

sourceCpp("[External]_Re__C++_codes_related_to_sparse_ICA/Bootstrap.cpp")
set.seed(2308)
ds <- rnorm(1000, mean = 10, sd = 5)
set.seed(34)
outR <- bootstrap_r(ds)
set.seed(34)
outCpp <- bootstrap_cpp(ds)

microbenchmark(bootstrap_cpp(ds), bootstrap_r(ds), times =50)


################################################################################
# RcppArmadillo

library(RcppArmadillo)
sourceCpp("[External]_Re__C++_codes_related_to_sparse_ICA/ArmadilloExamples.cpp")
X = matrix(rnorm(300), 30, 10)
Y = matrix(rnorm(200), 10, 20)
prodCpp = matrix_mult(X, Y)
prodR = X%*%Y
all.equal(prodCpp, prodR)

X = matrix(rnorm(30000), 300, 100)
Y = matrix(rnorm(20000), 100, 200)
microbenchmark(
  matrix_mult(X, Y),
  X%*%Y
)


# Linear model fit
set.seed(20386)
X = matrix(rnorm(100), 25, 4)
beta = rep(1, 4)
Y = X %*% beta + rnorm(25, sd = 0.5)
outC = fastLm(X, Y); names(outC)
cbind(outC$coefficients,
      solve(crossprod(X), crossprod(X, Y)))

X
a=soft_I(2,1)










