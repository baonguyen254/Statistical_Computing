# Newton-Raphson method
source("functions/daoham_bac2.R")
source("functions/daoham.R")
# newton_raphson <- function(func, x0, error = 1e-6, max_iter = 40) {
#     count <- 0
#     start_time <- Sys.time()
#     for (i in 1:max_iter) {
#         x1 <- x0 - calculate_derivative(func, x0)/calculate_derivative_level_2(func, x0)
#         if (func(x1) < error) {
#             break
#         }
#         x0 <- x1
#         count <- count + 1
#     }
#     # print the number of iterations and time taken
#     cat("Number of iterations:", count, "\n")
#     cat("Time taken:", format(Sys.time() - start_time), "\n")
#     return(x1)
# }


newton_raphson <- function(func, x0, data, error = 1e-6, max_iter = 40) {
  count <- 0
  start_time <- Sys.time()
  x <- x0
  for (i in 1:max_iter) {
    first_derivative <- calculate_derivative(func, x, data)
    second_derivative <- calculate_derivative_level_2(func, x, data)

    # Check for NA or Inf values in derivatives
    if (is.na(first_derivative) ||
      is.na(second_derivative) ||
      is.infinite(first_derivative) ||
      is.infinite(second_derivative)) {
      cat("Derivative calculation failed. Exiting.\n")
      return(NA)
    }

    x1 <- x - first_derivative / second_derivative
    if (abs(func(x1, data)) < error) {
      break
    }
    x <- x1
    count <- count + 1
  }
  cat("Number of iterations:", count, "\n")
  cat("Time taken:", format(Sys.time() - start_time), "\n")
  return(x)
}