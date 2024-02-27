# Source R functions
source("Rfunctions.R")

# Add Rcpp library
library(Rcpp)

# R bootsrap
ds <- rnorm(1000, mean = 10, sd = 5)
out <- bootstrap_r(ds)
library(microbenchmark)
microbenchmark(
  bootstrap_r(ds), times = 10
)

hist(out[,1])
hist(out[,2])

# Cpp bootstrap
sourceCpp("Bootstrap.cpp")
set.seed(2308)
ds <- rnorm(1000, mean = 10, sd = 5)
set.seed(34)
outR <- bootstrap_r(ds)
set.seed(34)
outCpp <- bootstrap_cpp(ds)


# Compare the two in terms of correctness
summary(outR[ , 1])
summary(outCpp[ , 1])

summary(outR[ , 2])
summary(outCpp[ , 2])

# Compare the two in terms of speed
library(microbenchmark)
microbenchmark(
  bootstrap_cpp(ds),
  bootstrap_r(ds), times = 50
)