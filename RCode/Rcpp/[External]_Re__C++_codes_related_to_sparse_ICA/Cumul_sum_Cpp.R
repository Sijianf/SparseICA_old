# Source R functions
source("Rfunctions.R")

# Add Rcpp library
library(Rcpp)

# Odd example
##################
cppFunction("bool isOddCpp(int num){
  bool result = (num%2 == 1);
  return result;
}")
isOddCpp(10)


# Compile Version 1 of cumul sum
########################################3
cppFunction("NumericVector cumul_sumCpp(NumericVector x){
  int p = x.size(); 
  NumericVector s = x;
  for(int i = 1; i < p; i++){
    s[i] = s[i-1] + x[i];
  }
  return(s);
}")
cumul_sumR(x = c(1, 3, 5))
cumul_sumCpp(x = c(1, 3, 5))

# What happened to x?
x = c(1, 3, 5)
s = cumul_sumCpp(x)
print(s)
print(x)


# Compile Version 2
########################################3
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

# Compile Version 3
########################################3
cppFunction("NumericVector cumul_sumCpp(NumericVector x){
            int p = x.size(); 
            NumericVector s(p); // new vector of the same size
            s[0] = x[0]; // extra 1st element initialization
            for(int i = 1; i < p; i++){
            s[i] = s[i-1] + x[i];
            }
            return(s);}")
x = c(1, 3, 5)
s = cumul_sumCpp(x); print(s)
print(x)

# Timing comparison
####################################################
library(microbenchmark)
p = 1000
x = rnorm(p)
identical(cumul_sumR(x), cumul_sumCpp(x))
microbenchmark(
  cumul_sumR(x),
  cumul_sumCpp(x), times = 50
)