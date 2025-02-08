## Xét mật độ f(x) = 1 − cos(x − θ)/2π trên 0 ≤ x ≤ 2π, trong đó θ là tham số và θ ∈ [−π, π]. Dữ liệu i.i.d. 
# sau phát sinh từ mật độ này: 3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96, 2.53, 3.88, 2.22,
# 3.47, 4.82, 2.46, 2.99, 2.54, 0.52, 2.50. Chúng ta muốn ước tính θ.

f <- function(x, theta) {
  (1 - cos(x - theta)) / (2 * pi)
}


data <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96, 2.53, 3.88, 2.22,
          3.47, 4.82, 2.46, 2.99, 2.54, 0.52, 2.50)

theta_x <- seq(from = -pi, to = pi, by = 0.1)
f_vec <- Vectorize(f, vectorize.args = "theta")

f_data <- f_vec(data = data, theta = theta_x)

## Question a. Vẽ đồ thị hàm log-likelihood với θ ∈ [−π, π].
logL_f <- function(data, theta) {
  n <- length(data)
  ll <- -n * log(2 * pi) + sum(log(1 - cos(data - theta)))
  return(ll)
}

logL_f_vec <- Vectorize(logL_f, vectorize.args = "theta")

logL_f_data <- logL_f_vec(data = data, theta = theta_x)

plot(theta_x, logL_f_data, type = "l", xlab = "theta", ylab = "logL_f_data", main = "Log-likelihood function")

## Question b. Tìm ước lượng method-of-moment cho θ.

# method-of-moment thì phải đưa về  bậc 0 của moment, dùng tích phân để tính moment
# tích phân của hàm f(x) = 1 - cos(x - theta) / (2*pi) trên 0 ≤ x ≤ 2π = sin(theta) + pi 
# E_f <- function(theta){
#     sin(theta) + pi 
# }


X_mean <- mean(data)

theta_hat <- asin(X_mean - pi) # 0.05844061


## Question c. Tìm MLE cho θ bằng phương pháp Newton–Raphson, sử dụng kết quả từ (b) làm giá trị bắt đầu. Bạn
# tìm thấy những nghiệm nào khi bắt đầu ở -2.7 và 2.7?
logL_f_derivative <- function(theta, x) {
  grad <- sum(sin(x - theta) / (1 - cos(x - theta)))
  return(grad)
}

logL_f_second_derivative <- function(theta, x) {
  hess <- sum(-1 / (1 - cos(x - theta)))
  return(hess)
}

newton_raphson <- function(data, theta0, tol = 1e-9, max_iter = 40) {
  theta <- theta0
  res_x_newt <- theta
  for (i in 1:max_iter) {
    grad <- logL_f_derivative(theta, data)
    hess <- logL_f_second_derivative(theta, data)
    theta <- theta - grad / hess
    res_x_newt[i + 1] <- theta
    if (abs(res_x_newt[i + 1] - res_x_newt[i]) < tol) {
      print(paste("Converged at iteration", i))
      break
    }
  }
  return(theta)
}

newton_raphson(data = data, theta0 = 0.05844061, tol = 1e-4, max_iter = 100)

newton_raphson(data = data, theta0 = -2.7, tol = 1e-6, max_iter = 40)
newton_raphson(data = data, theta0 = 2.7, tol = 1e-6, max_iter = 40)

## Question d. Lặp lại phần (c) bằng cách sử dụng 200 giá trị bắt đầu cách đều nhau giữa −π và π. Nói cách khác,
# chia tập hợp các giá trị bắt đầu thành các nhóm riêng biệt, với mỗi nhóm tương ứng với một kết quả
# duy nhất riêng biệt của quá trình tối ưu hóa (một chế độ cục bộ). Thảo luận về kết quả của bạn.

theta_start <- seq(from = -pi, to = pi, length.out = 200)

for (theta in theta_start) {
  print(paste("theta =", theta))
  newton_raphson(data = data, theta0 = theta, max_iter = 40)
}

# Kết quả là một tập hợp các giá trị theta mà mỗi giá trị theta tương ứng với một kết quả duy nhất


## Question e. Tìm hai giá trị bắt đầu, gần bằng nhau nhất có thể, mà phương pháp Newton–Raphson hội tụ thành
# hai nghiệm khác nhau.

theta_start <- seq(from = -pi, to = pi, by = 0.1)

res_newton <- c()
for (i in 1:length(theta_start)) {
  print(paste("theta =", theta_start[i]))
  res_newton[i] <- newton_raphson(data = data, theta0 = theta_start[i], max_iter = 40)
}
print(head(cbind(res_newton), 40))

# nhận xét:
# "theta = 1.85840734641021" "Converged at iteration 7"
# "theta = 1.95840734641021" "Converged at iteration 5"
