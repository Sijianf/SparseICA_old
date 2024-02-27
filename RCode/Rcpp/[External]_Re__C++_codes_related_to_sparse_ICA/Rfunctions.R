# Toy R functions

# Identifying whether the supplied number is odd
# num - scalar input
isOddR <- function(num){
  result <- num %% 2 == 1
  return(result)
}



# Cumulitive sums of elements of a vector
# x - vector
cumul_sumR <- function(x){
  p <- length(x)
  s <- x
  for (i in 2:p){
    s[i] <- s[i-1] + x[i]
  }
  return(s)
}

# Bootstrap for sample mean and sample standard deviation
# ds - vector of observations
# B - number of bootstrap samples
bootstrap_r <- function(ds, B = 1000){ 
  boot_stat <- matrix(NA, nrow = B, ncol = 2) 
  n <- length(ds)
  # Perform bootstrap
  for(i in 1:B) {
    # Create a sample of size n with replacement
    gen_data <- ds[sample(n, n, replace = TRUE)] 
    # Calculate sample data mean and SD 
    boot_stat[i,] <- c(mean(gen_data), sd(gen_data))
  }
  return(boot_stat)
}  