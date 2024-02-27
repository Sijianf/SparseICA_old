# Add Rcpp and Rcpp Armadillo library
library(Rcpp)
library(RcppArmadillo)



# Source two versions of timesTwo function
sourceCpp("Test_pointers.cpp")

# Copy by value
x = c(1, 3, 5)
y = timesTwo(x)
y
x

# Copy by pointer
x = c(1, 3, 5)
y = timesTwo_pointer(x)
y
x
