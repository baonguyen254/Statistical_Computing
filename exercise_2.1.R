# Define the Cauchy log-likelihood function
# This function calculates the log-likelihood of the Cauchy distribution given data and a parameter theta
log_likelihood <- function(theta, data) {
  n <- length(data)
  return(sum(-log(1 + (data - theta)^2)))
}

logL_cauchy <- function(data, theta) {
  sum(dcauchy(x = data, location = theta, scale = 1, log = TRUE))
}

# Derivative of log-likelihood function
# Computes the gradient of the log-likelihood function
log_likelihood_derivative <- function(theta, data) {
  return(sum(2 * (data - theta) / (1 + (data - theta)^2)))
}

# Second derivative of log-likelihood function
# Computes the Hessian (second derivative) of the log-likelihood function
log_likelihood_second_derivative <- function(theta, data) {
  return(sum(-2 * ((1 - (data - theta)^2) / (1 + (data - theta)^2)^2)))
}

# Newton-Raphson method
# Iteratively updates theta using its derivative and second derivative to find the MLE
newton_raphson <- function(data, start_points, tol = 1e-6, max_iter = 100) {
  results <- list()
  for (start in start_points) {
    theta <- start
    iter <- 0
    repeat {
      iter <- iter + 1
      grad <- log_likelihood_derivative(theta, data)
      hess <- log_likelihood_second_derivative(theta, data)
      theta_new <- theta - grad / hess
      if (abs(theta_new - theta) < tol || iter >= max_iter) break
      theta <- theta_new
    }
    results[[as.character(start)]] <- list(theta = theta, iterations = iter)
  }
  return(results)
}

# Bisection method
# Uses a bracketing approach to find the root of the derivative of the log-likelihood function
bisection_method <- function(data, a, b, tol = 1e-6, max_iter = 100) {
  iter <- 0
  repeat {
    iter <- iter + 1
    c <- (a + b) / 2
    f_a <- log_likelihood_derivative(a, data)
    f_c <- log_likelihood_derivative(c, data)
    if (f_a * f_c < 0) {
      b <- c
    } else {
      a <- c
    }
    if (abs(b - a) < tol || iter >= max_iter) break
  }
  return(list(theta = c, iterations = iter))
}

# Fixed-point iteration
# Iteratively updates theta using a scaling factor alpha to find the root of the derivative
fixed_point_iteration <- function(data, start, alpha, tol = 1e-6, max_iter = 100) {
  theta <- start
  iter <- 0
  repeat {
    iter <- iter + 1
    grad <- log_likelihood_derivative(theta, data)
    theta_new <- theta + alpha * grad
    if (abs(theta_new - theta) < tol || iter >= max_iter) break
    theta <- theta_new
  }
  return(list(theta = theta, iterations = iter))
}

# Secant method
# Uses two initial guesses and approximates the derivative to find the root iteratively
secant_method <- function(data, theta_0, theta_1, tol = 1e-6, max_iter = 100) {
  iter <- 0
  repeat {
    iter <- iter + 1
    f_0 <- log_likelihood_derivative(theta_0, data)
    f_1 <- log_likelihood_derivative(theta_1, data)
    theta_new <- theta_1 - f_1 * (theta_1 - theta_0) / (f_1 - f_0)
    if (abs(theta_new - theta_1) < tol || iter >= max_iter) break
    theta_0 <- theta_1
    theta_1 <- theta_new
  }
  return(list(theta = theta_new, iterations = iter))
}

# Main script to run the methods
# Given data points, we perform the required analyses and methods
data <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29,
          3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)

# Question (a): Plot the log-likelihood function
theta_values <- seq(-15, 15, length.out = 1000)
log_likelihood_values <- sapply(theta_values, log_likelihood, data = data)
plot(theta_values, log_likelihood_values, type = "l",
     xlab = "Theta", ylab = "Log-Likelihood",
     main = "Log-Likelihood Function for Cauchy Distribution")

# Question (a): Newton-Raphson method
newton_results <- newton_raphson(data, c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38))
print(newton_results)
# Hãy nhận xét kết quả nghiệm thu được (hội tụ, ổn định). Giá trị trung bình của dữ liệu có phải là điểm khởi đầu tốt không?
# Nhận xét:
# - Newton-Raphson method hội tụ nhanh chóng và ổn định.
# - Giá trị trung bình của dữ liệu không phải là điểm khởi đầu tốt.


# Question (b): Bisection method
bisection_result <- bisection_method(data, -1, 1)
print(bisection_result)

# Question (c): Fixed-point iteration
fixed_point_results <- list(
  alpha_1 = fixed_point_iteration(data, -1, 1),
  alpha_064 = fixed_point_iteration(data, -1, 0.64),
  alpha_025 = fixed_point_iteration(data, -1, 0.25)
)
print(fixed_point_results)

# Question (d): Secant method
secant_result <- secant_method(data, -2, -1)
print(secant_result)

# Explanation and Comparison (e):
# The above methods allow us to investigate the speed and stability of optimization techniques.
# Newton-Raphson generally converges quickly but depends heavily on a good starting point.
# Bisection is robust but slower. Fixed-point iteration's performance depends on the scaling factor.
# Secant method is faster than bisection but may fail if initial guesses are poor.

