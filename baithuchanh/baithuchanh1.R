cauchy_data <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71,
                 -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75, 0.27, 43.21)

# hàm dcauchy() chỉ nhận 1 giá trị theta trong 1 lần tính, nên để tính nhiều giá trị của hàm logL_cauchy() tại nhiều giá trị theta cùng 1 lúc
logL_cauchy <- function(data, theta) {
  sum(dcauchy(x = data, location = theta, scale = 1, log = TRUE))
}


logL_cauchy(cauchy_data, 1)

logL_cauchy_vec <- Vectorize(FUN = logL_cauchy, vectorize.args = "theta")


# Với θ ∈ (−10, 10), ta tính các giá trị của hàm log-likelihood dựa trên dữ liệu.
theta_x <- seq(from = -10, to = 10, by = 0.1)
logL_cauchy_data <- logL_cauchy_vec(data = cauchy_data, theta = theta_x)

## Question a. Vẽ biểu đồ hàm log-likelihood của θ dựa trên dữ liệu.
# Vẽ đồ thị của hàm log-likelihood
plot(theta_x, logL_cauchy_data, type = "l", xlab = "theta", ylab = "logL_cauchy_data", main = "Log-likelihood function")

## Question b. Áp dụng thuật toán Newton-Raphson để tìm ước lượng hợp lý cực đại cho θ, với mỗi điểm bắt đầu sau: 
# -11, -1, 0, 1.5, 4, 4.7, 7, 8, và 38. Hãy nhận xét kết quả nghiệm thu được (hội tụ, ổn định). Giá trị trung
# bình của dữ liệu có phải là điểm khởi đầu tốt không?

# tìm đạo hàm của logL_cauchy
logL_cauchy_derivative <- function(data, theta) {
  sum(2 * (data - theta) / (1 + (data - theta)^2))
}

# tìm đạo hàm cấp 2 của logL_cauchy
logL_cauchy_second_derivative <- function(data, theta) {
  # sum(-2 * ((1 - (data - theta)^2) / (1 + (data - theta)^2)^2))
  sum(2 * ((data - theta)^2 - 1) / (1 + (data - theta)^2)^2)

}

newton_raphson <- function(data, theta0, tol = 1e-9, max_iter = 20) {
  theta <- theta0
  res_x_newt <- theta
  for (i in 1:max_iter) {
    grad <- logL_cauchy_derivative(data, theta)
    hess <- logL_cauchy_second_derivative(data, theta)
    theta <- theta - grad / hess
    res_x_newt[i + 1] <- theta
    if (abs(res_x_newt[i + 1] - res_x_newt[i]) < tol) {
      print(paste("Converged at iteration", i))
      return(theta)
    }
  }
  # print(head(cbind(res_x_newt), 20))
  return(theta)
}

# tìm ước lượng hợp lý cực đại cho θ, với mỗi điểm bắt đầu sau: -11, -1, 0, 1.5, 4, 4.7, 7, 8, và 38.
start_points <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)

for (start in start_points) {
  print(paste("Start point:", start))
  print(newton_raphson(cauchy_data, start))
}


## Question c. Áp dụng phương pháp bisection với điểm bắt đầu là [−1, 1]. Sử dụng các lần chạy bổ sung để minh họa
# cách thức mà phương pháp bisection có thể không tìm được giá trị cực đại toàn cục.

bisection <- function(data, a, b, tol = 1e-9, max_iter = 20) {
  x <- a + (b - a) / 2
  res_a <- a
  res_b <- b
  res_x <- x
  for (i in 1:max_iter) {
    f_a <- logL_cauchy_derivative(data, a)
    f_x <- logL_cauchy_derivative(data, x)
    if (f_a * f_x <= 0) {
      res_a[i + 1] <- a
      res_b[i + 1] <- b <- x
    } else {
      res_a[i + 1] <- a <- x
      res_b[i + 1] <- b
    }
    x <- a + (b - a) / 2
    res_x[i + 1] <- x
    if (abs(res_x[i + 1] - res_x[i]) < tol) {
      print(paste("Converged at iteration", i))
      break
    }
  }
  # print(head(cbind(res_a, res_b, res_x), 20))
  return(x)
}

bisection(cauchy_data, -1, 1, tol = 1e-6, max_iter = 100)
# "Converged at iteration 20" -0.1922865

## Question d. Từ các giá trị khởi đầu của (θ(0), θ(1)) = (−2, −1), áp dụng phương pháp secant 
# để tìm ước lượng hợp lý cực đại cho θ. 
# Điều gì xảy ra khi (θ(0), θ(1)) = (−3, 3), và đối với các lựa chọn khởi đầu khác?

secant <- function(data, theta0, theta1, tol = 1e-9, max_iter = 50) {
  theta <- theta1
  res_theta <- theta
  for (i in 1:max_iter) {
    f_theta0 <- logL_cauchy_derivative(data, theta0)
    f_theta1 <- logL_cauchy_derivative(data, theta1)
    theta <- theta1 - f_theta1 * (theta1 - theta0) / (f_theta1 - f_theta0)
    res_theta[i + 1] <- theta
    if (abs(res_theta[i + 1] - res_theta[i]) < tol) {
      print(paste("Converged at iteration", i))
      break
    }
  }
  # print(head(cbind(res_theta), 20))
  return(theta)
}

secant(cauchy_data, -2, -1)
#[1] "Converged at iteration 2"
#[1] -0.4412058
secant(cauchy_data, -3, 3)
#[1] "Converged at iteration 2"
#[1] 2.685789
secant(cauchy_data, 1, 5)
# [1] "Converged at iteration 2"
# [1] 0.4813908


# nhận xét:
# - Khi (θ(0), θ(1)) = (−2, −1), phương pháp secant tìm được giá trị cực đại toàn cục.
# - Khi (θ(0), θ(1)) = (−3, 3), phương pháp secant không tìm được giá trị cực đại toàn cục.
# - Khi (θ(0), θ(1)) = (1, 5), phương pháp secant tìm được giá trị cực đại toàn cục.

## Question e. Sử dụng ví dụ này để so sánh tốc độ và tính ổn định của phương pháp Newton–Raphson, bisection,
# fixed-point và phương pháp secant. Kết luận của bạn có thay đổi khi bạn áp dụng các phương pháp
# này cho một mẫu ngẫu nhiên có kích thước 20 từ phân phối chuẩn (θ, 1) không?

fixed_point <- function(data, theta0, alpha, tol = 1e-9, max_iter = 20) {
  theta <- theta0
  res_theta <- theta
  for (i in 1:max_iter) {
    res_theta[i + 1] <- theta <- theta + alpha * logL_cauchy_derivative(data, theta)
    if (abs(res_theta[i + 1] - res_theta[i]) < tol) {
      print(paste("Converged at iteration", i))
      return(theta)
    }
  }
  print(head(cbind(res_theta), 20))
  return(theta)
}

# alpha = 0.1
fixed_point(cauchy_data, -1, 0.1, tol = 1e-3, max_iter = 200) # -0.1939695
fixed_point(cauchy_data, -1, 1, tol = 1e-3, max_iter = 1000)
fixed_point(cauchy_data, -1, 0.64, tol = 1e-3, max_iter = 250) # -0.1927677
fixed_point(cauchy_data, -1, 0.25, tol = 1e-3, max_iter = 250) # -0.1923656

# tạo một mẫu ngẫu nhiên có kích thước 20 từ phân phối chuẩn (θ, 1)
normal_data <- rnorm(20, mean = 0, sd = 1)
plot(normal_data, type = "l", xlab = "theta", ylab = "logL_cauchy_data", main = "Log-likelihood function")

fixed_point(normal_data, 0, 0.1)

newton_raphson(normal_data, 0)

bisection(normal_data, -1, 1)

secant(normal_data, -2, -1)

# nhận xét:
# - Phương pháp Newton-Raphson và phương pháp secant tìm được giá trị cực đại toàn cục.
# - Phương pháp bisection không tìm được giá trị cực đại toàn cục.
# - Phương pháp fixed-point tìm được giá trị cực đại toàn cục.
