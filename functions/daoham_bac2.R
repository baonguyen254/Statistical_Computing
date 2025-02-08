# đạo hàm bậc 2
# calculate_derivative_level_2 <- function(func, x_values, h = 1e-6) {
#   f_prime <- Deriv(func, "x")
#   f_prime_prime <- Deriv(f_prime, "x")
#   derivatives <- f_prime_prime(x_values)
#   return(derivatives)
# }

calculate_derivative_level_2 <- function(func, x, data, h = 1e-6) {
  # Sử dụng phương pháp sai phân hữu hạn trung tâm
  (func(x + h, data) - 2 * func(x, data) + func(x - h, data)) / (h^2)
}