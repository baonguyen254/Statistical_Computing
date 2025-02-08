# derivative_function.R
# calculate_derivative <- function(func, x_values, h = 1e-6) {
#   # func: The function to differentiate (must be a function)
#   # x_values: A vector of x values at which to calculate the derivative
#   # h: Step size for numerical differentiation (smaller h is more accurate but slower)

#   if(!require(Deriv)){install.packages("Deriv")}

# #Using Deriv package for numerical derivative
# f_prime <- Deriv(func, "x")
# derivatives <- f_prime(x_values)
# return(derivatives)

# }
calculate_derivative <- function(func, x, data, h = 1e-6) {
  # Sử dụng phương pháp sai phân hữu hạn trung tâm
  (func(x + h, data) - func(x - h, data)) / (2 * h)
}