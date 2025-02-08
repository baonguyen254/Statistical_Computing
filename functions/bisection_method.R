# bisection_method.R
bisection_method <- function(func, a0, b0, error = 1e-6, max_iter = 40) {
  # func: the function to find the root
  # a0: the initial lower bound
  # b0: the initial upper bound
  # error: the error tolerance
  # max_iter: the maximum number of iterations
  start_time <- Sys.time()
  count <- 0
  for (i in 1:max_iter) {
    x0 <- (a0 + b0) / 2 # midpoint
    if (calculate_derivative(func, a0) * calculate_derivative(func, x0) < 0) {
      b0 <- x0
    } else {
      a0 <- x0
    }
    if (func(x0) < error) {
      break
    }
    count <- count + 1
  }
  # print the number of iterations and time taken
  cat("Number of iterations:", count, "\n")
  cat("Time taken:", format(Sys.time() - start_time), "\n")
  #   print(paste("Number of iterations:", count), "Time taken:", Sys.time() - start_time)
  return(x0)
}
