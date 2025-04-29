f <- function(x) {
  (-1) * ((x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2)
}

x1 <- seq(-5, 5, by = 0.1)
x2 <- seq(-5, 5, by = 0.1)
ff <- matrix(0, nrow = length(x1), ncol = length(x2))
for (i in 1:length(x1)) {
  for (j in 1:length(x2)) {
    ff[i, j] <- f(x = c(x1[i], x2[j]))
  }
}

# hàm contour() để vẽ đường đồng mức của f(x1, x2) nhằm xác định vùng cực trị
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, 0),
        xlab = "x1", ylab = "x2")

# hàm image() để vẽ hình ảnh của f(x1, x2)
# image(x1, x2, ff, col = heat.colors(100))


# áp dụng pp newton để tìm cực trị
f_prime <- function(x) {
  f_prime_x1 <- (-1) * (4 * x[1]^3 + 4 * x[1] * x[2] - 42 * x[1] + 2 * x[2]^2 - 14)
  f_prime_x2 <- (-1) * (2 * x[1]^2 - 26 * x[2] - 22 +
    4 * x[1] * x[2] +
    4 * x[2]^3)
  out <- c(f_prime_x1, f_prime_x2)
  return(out)
}

f_2prime <- function(x) {
  f_2prime_x11 <- (-1) * (12 * x[1]^2 + 4 * x[2] - 42)
  f_2prime_x22 <- (-1) * (12 * x[2]^2 + 4 * x[1] - 26)
  f_2prime_x12 <- (-1) * (4 * (x[1] + x[2]))
  out <- matrix(c(f_2prime_x11, f_2prime_x12, f_2prime_x12, f_2prime_x22),
                nrow = 2, byrow = TRUE)
  return(out)
}

iter <- 40
x <- c(0, 0)
x_vals <- matrix(NA, nrow = iter + 1, ncol = 2)
x_vals[1,] <- x
for (i in 1:iter) {
  x <- x - as.numeric(solve(f_2prime(x)) %*% f_prime(x))
  x_vals[i + 1,] <- x
  # if (abs(x_vals[i + 1, ] - x_vals[i, ]) < 1e-6) {
  #     print(paste("Converged at iteration", i))
  #     break
  # }
}

head(x_vals, n = 4)

contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")

## function with x0
newton_method <- function(x0, tol = 1e-6, iter = 40) {
  x <- x0
  x_vals <- matrix(NA, nrow = iter + 1, ncol = 2)
  x_vals[1,] <- x
  for (i in 1:iter) {
    x <- x - as.numeric(solve(f_2prime(x)) %*% f_prime(x))
    x_vals[i + 1,] <- x
    # if (abs(x_vals[i + 1, ] - x_vals[i, ]) < tol) {
    #     print(paste("Converged at iteration", i))
    #     break
    # }
  }
  return(x_vals)
}

## Thực hành 1: Tiến hành thuật toán Newton với các điểm bắt đầu sau: 
# x(0) = (−1.5, 1.5), x(0) = (2.5, 3),
# x(0) = (−1.5, −2), x(0) = (2, −2) và x(0) = (−4, 0). Nhận xét về các kết quả đạt được.

x0 <- c(-1.5, 1.5)
x_vals <- newton_method(x0, tol = 1e-6, iter = 100)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")
# không hội tụ

x0 <- c(2.5, 3)
x_vals <- newton_method(x0, tol = 1e-6, iter = 100)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")
# hội tụ

x0 <- c(-1.5, -2)
x_vals <- newton_method(x0, tol = 1e-3)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")
# không hội tụ

x0 <- c(2, -2)
x_vals <- newton_method(x0, tol = 1e-9)
head(x_vals, n = 12)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")
# không hội tụ

x0 <- c(-4, 0)
x_vals <- newton_method(x0, tol = 1e-6)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")
# không hội tụ

# nhận xét:
# x0 = (-1.5, 1.5) không hội tụ
# x0 = (2.5, 3) hội tụ Converged at iteration 6
# x0 = (-1.5, -2) không hội tụ
# x0 = (2, -2) hội tụ Converged at iteration 12
# x0 = (-4, 0) không hội tụ

# Phương pháp Steepest Ascent
M <- diag(c(-1, -1), 2, 2)
alpha_default <- 1
alpha <- alpha_default
x0 <- c(-1.5, 1.5)
iter <- 40
x_vals_steep <- matrix(NA, nrow = iter + 1, ncol = 2)
x_vals_steep[1,] <- x0
for (i in 1:iter) {
  score <- f_prime(x0)
  hessian_inv <- solve(M)
  xt <- x0 - alpha * hessian_inv %*% score
  # REDUCE ALPHA UNTIL A CORRECT STEP IS REACHED
  f_t <- f(xt)
  f_t0 <- f(x0)
  while (f_t < f_t0) {
    alpha <- alpha / 2
    xt <- x0 - alpha * hessian_inv %*% score
    f_t <- f(xt)
  }
  ## RESET alpha and x0
  x_vals_steep[i + 1,] <- x0 <- xt
  alpha <- alpha_default
}
head(x_vals_steep, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, -2, -1, 0),
        xlab = "x1", ylab = "x2")
points(x_vals_steep[, 1], x_vals_steep[, 2], pch = 16, type = "b")
# với điểm bắt đầu x(0) = (0, 0), M(t), phương pháp Steepest Ascent đưa nghiệm lặp di chuyển 
# vào vùng cực trị (toàn cục).

steepest_ascent <- function(x0, M, alpha_default = 1, tol = 1e-6, iter = 40) {
  alpha <- alpha_default
  x_vals <- matrix(NA, nrow = iter + 1, ncol = 2)
  x_vals[1,] <- x0
  for (i in 1:iter) {
    score <- f_prime(x0)
    hessian_inv <- solve(M)
    xt = x0 - alpha * hessian_inv %*% score
    f_t <- f(xt)
    f_t0 <- f(x0)
    while (f_t < f_t0) {
      alpha <- alpha / 2
      xt <- x0 - alpha * hessian_inv %*% score
      f_t <- f(xt)
    }
    x_vals[i + 1,] <- x0 <- xt
    alpha <- alpha_default
    # condition to stop
    if (abs(x_vals[i + 1,] - x_vals[i,]) < tol) {
      print(paste("Converged at iteration", i))
      break
    }
  }
  return(x_vals)
}

## thực hành 2: Tiến hành thuật toán Steepest Ascent với các điểm bắt đầu sau: x(0) = (−1.5, 1.5),
#x(0) = (2.5, 3), x(0) = (−1.5, −2), x(0) = (2, −2) và x(0) = (−4, 0). M(t) = (−1 0 0 −1)
#, với mọi t. Nhận xét về các kết quả đạt được.

x0 <- c(-1.5, 1.5)
M <- matrix(c(-1, 0, 0, -1), nrow = 2, byrow = TRUE)
x_vals <- steepest_ascent(x0, M, tol = 1e-6, iter = 100)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, -2, -1, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")


x0 <- c(2.5, 3)
M <- matrix(c(-1, 0, 0, -1), nrow = 2, byrow = TRUE)
x_vals <- steepest_ascent(x0, M, tol = 1e-6, iter = 100)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, -2, -1, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")

x0 <- c(-1.5, -2)
M <- matrix(c(-1, 0, 0, -1), nrow = 2, byrow = TRUE)
x_vals <- steepest_ascent(x0, M, tol = 1e-6, iter = 100)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, -2, -1, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")

x0 <- c(2, -2)
M <- matrix(c(-1, 0, 0, -1), nrow = 2, byrow = TRUE)
x_vals <- steepest_ascent(x0, M, tol = 1e-6, iter = 100)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, -2, -1, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")

x0 <- c(-4, 0)
M <- matrix(c(-1, 0, 0, -1), nrow = 2, byrow = TRUE)
x_vals <- steepest_ascent(x0, M, tol = 1e-6, iter = 100)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, -2, -1, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")

# nhận xét:
# x0 = (-1.5, 1.5) hội tụ
# x0 = (2.5, 3) hội tụ 
# x0 = (-1.5, -2) hội tụ
# x0 = (2, -2) hội tụ 
# x0 = (-4, 0) hội tụ


## thực hành 3: Đối với các trường hợp không hội tụ về nghiệm toàn cục trong Thực hành 2, hãy cố gắng
# hiệu chỉnh ma trận M(t) để đạt được sự hội tụ.


# Phương pháp Quasi-Newton BFGS
M <- diag(c(-1, -1), 2, 2)
alpha_default <- 1
alpha <- alpha_default
x0 <- c(0, 0)
iter <- 40
x_vals_bfgs <- matrix(NA, nrow = iter + 1, ncol = 2)
x_vals_bfgs[1,] <- x0
epsilon <- 1e-10
for (i in 1:iter) {
  hessian_inv <- solve(M)
  score <- f_prime(x0)
  xt <- x0 - alpha * hessian_inv %*% score
  # REDUCE ALPHA UNTIL A CORRECT STEP IS REACHED
  f_t <- f(xt)
  f_t0 <- f(x0)
  while (f_t < f_t0) {
    alpha <- alpha / 2
    xt <- x0 - alpha * hessian_inv %*% score
    f_t <- f(xt)
  }
  x_vals_bfgs[i + 1,] <- xt
  ## UPDATE M(t) - BFGS METHOD
  ht <- xt - x0
  score_t <- f_prime(xt)
  ut <- score_t - score
  v <- ut - M %*% ht
  M_old <- M
  M <- M - ((M %*% ht %*% t(M %*% ht)) / (as.numeric(t(ht) %*% M %*% ht))) +
    ((ut %*% t(ut)) / (as.numeric(t(ht) %*% ut)))
  #CHECK CONDITION TO UPDATE M
  if (abs((t(v) %*% ht)[1]) < epsilon) {
    M <- M_old
  }
  ## RESET alpha and x0
  alpha <- alpha_default
  x0 <- xt
}
head(x_vals_bfgs, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, -2, -1, 0),
        xlab = "x1", ylab = "x2")
points(x_vals_bfgs[, 1], x_vals_bfgs[, 2], pch = 16, type = "b")


# Thực hành 4: Tiến hành thuật toán Quasi-Newton BFGS với các điểm bắt đầu sau: 
# x(0) = (−1.5, 1.5), x(0) = (2.5, 3), x(0) = (−1.5, −2), x(0) = (2, −2) và x(0) = (−4, 0). 
# Nhận xét về các kết quả đạt được. Việc hiệu chỉnh M(t) có giúp ích gì hay không?

Quasi_Newton <- function(x0, M, alpha_default = 1, tol = 1e-6, iter = 40) {
  alpha <- alpha_default
  x_vals <- matrix(NA, nrow = iter + 1, ncol = 2)
  x_vals[1,] <- x0
  for (i in 1:iter) {
    hessian_inv <- solve(M)
    score <- f_prime(x0)
    xt <- x0 - alpha * hessian_inv %*% score
    f_t <- f(xt)
    f_t0 <- f(x0)
    while (f_t < f_t0) {
      alpha <- alpha / 2
      xt <- x0 - alpha * hessian_inv %*% score
      f_t <- f(xt)
    }
    x_vals[i + 1,] <- xt
    ## UPDATE M(t) - BFGS METHOD
    ht <- xt - x0
    score_t <- f_prime(xt)
    ut <- score_t - score
    v <- ut - M %*% ht
    M_old <- M
    M <- M - ((M %*% ht %*% t(M %*% ht)) / (as.numeric(t(ht) %*% M %*% ht))) +
      ((ut %*% t(ut)) / (as.numeric(t(ht) %*% ut)))
    #CHECK CONDITION TO UPDATE M
    if (abs((t(v) %*% ht)[1]) < epsilon) {
      M <- M_old
    }
    # condition to stop
    if (abs(x_vals[i + 1,] - x_vals[i,]) < tol) {
      print(paste("Converged at iteration", i))
      break
    }
    ## RESET alpha and x0
    alpha <- alpha_default
    x0 <- xt
  }
  return(x_vals)
}

x0 <- c(-1.5, 1.5)
M <- matrix(c(-1, 0, 0, -1), nrow = 2, byrow = TRUE)
x_vals <- Quasi_Newton(x0, M, tol = 1e-6, iter = 100)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, -2, -1, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")

x0 <- c(2.5, 3)
M <- matrix(c(-2, 1, 3, -3), nrow = 2, byrow = TRUE)
x_vals <- Quasi_Newton(x0, M, tol = 1e-6, iter = 100)
head(x_vals, n = 4)
contour(x1, x2, ff, levels = c(seq(-600, -50, by = 50), -30, -20, -10, -5, -2, -1, 0),
        xlab = "x1", ylab = "x2")
points(x_vals[, 1], x_vals[, 2], pch = 16, type = "b")


# nhận xét:
# Việc hiệu chỉnh M(t) sẽ giúp đạt được sự hội tụ. 


## Thực hành 5: Áp dụng các tiêu chuẩn sai số để dừng vòng lặp, với các phương pháp đã làm ở trên.


# Bai tap 1: Áp dụng mô hình phân phối Weibull để mô hình hóa dữ liệu về thời gian hỏng của lò xo trong
# thí nghiệm với mức ứng suất 950 N/mm2: 225, 171, 198, 189, 189, 135, 162, 135, 117, 162

# hàm mật độ của phân phối Weibull:
data <- c(225, 171, 198, 189, 189, 135, 162, 135, 117, 162)

f <- function(y, theta, alpha) {
  # if (theta <= 0 || alpha <= 0 || y <=0) return(0) #Handle invalid inputs
  alpha / theta *
    (y / theta)^(alpha - 1) *
    exp(-(y / theta)^alpha)
}

logLikelihood_weibull <- function(data, theta, alpha) {
  n <- length(data)
  # log_likelihood <- n * log(alpha) - n * log(theta) + (alpha - 1) * sum(log(data/theta)) - sum((data/theta)^alpha)
  likelihood <- sapply(data, function(y) f(y, theta, alpha))
  log_likelihood <- sum(log(likelihood))
  return(log_likelihood)
}

# Do θ, α > 0, nên để đảm bảo nghiệm lặp thỏa tính chất, ta thực hiện phép đổi biến ψ = (log(θ), log(α)), sau
#đó với phép biến đổi ngược, exp(·), ta thu được kết quả θ^, α^ > 0.

# a. Vẽ đồ thị hàm log-likelihood với θ và α trên miền phù hợp.

# Use a range that avoids zero and very small values for theta and alpha
theta <- seq(150, 250, by = 0.1) #Avoid theta = 0
alpha <- seq(1, 10, by = 0.1) #Avoid alpha = 0
logL_matrix <- matrix(NA, nrow = length(theta), ncol = length(alpha))
for (i in 1:length(theta)) {
  for (j in 1:length(alpha)) {
    logL_matrix[i, j] <- logLikelihood_weibull(data, theta[i], alpha[j])
  }
}
# contour plot, nlevels là giá trị đường đồng mức với giá trị đó, levels số  lượng giá trị của hàm f được vẽ
contour(theta, alpha, logL_matrix, nlevels = 100,
        xlab = "theta", ylab = "alpha")

# b. Tìm MLE cho θ bằng phương pháp Newton, giá trị bắt đầu θ(0) = 160, α(0) = 1.

logL_prime <- function(par, data) {
  n <- length(data)
  res <- numeric(2)
  res[1] <- -n * par[2] / par[1] + par[2] * sum(data^par[2]) / (par[1]^(par[2] + 1))
  res[2] <- n / par[2] - n * log(par[1]) + sum(log(data)) - sum((data / par[1])^par[2] * log(data / par[1]))
  return(res)
}

logL_2prime <- function(par, data) {
  n <- length(data)
  res <- matrix(0, nrow = 2, ncol = 2)
  res[1, 1] <- n * par[2] / par[1]^2 - par[2] * (par[2] + 1) * sum(data^par[2]) / (par[1]^(par[2] + 2))
  res[1, 2] <- res[2, 1] <- -n / par[1] + (1 / par[1]^(par[2] + 1)) * (sum(data^par[2]) + par[2] * sum(data^par[2] * log(data / par[1])))
  res[2, 2] <- -n / par[2]^2 - sum((data / par[1])^par[2] * (log(data / par[1]))^2)
  return(res)
}


newton_method <- function(par, data, tol = 1e-6, iter = 100) {
  par <- par0
  par_vals <- matrix(NA, nrow = iter + 1, ncol = 2)
  par_vals[1,] <- par
  for (i in 1:iter) {
    par <- par - solve(logL_2prime(par, data)) %*% logL_prime(par, data)
    if (abs(par - par0) < tol) {
      print(paste("Converged at iteration", i))
      break
    }
    par_vals[i + 1,] <- par0 <- par
  }
  print(paste("Not converged after", iter, "iterations"))
  return(par_vals)
}

par <- log(c(163, 1))
logL_prime(par, data)
solve(logL_2prime(par, data))
par_ests <- matrix(NA, nrow = 100, ncol = 2)
par_ests[1,] <- par
par_ests <- newton_method(par, data, tol = 1e-6, iter = 100)
head(par_ests, n = 4)
contour(theta, alpha, logL_matrix, nlevels = 100,
        xlab = "theta", ylab = "alpha")
points(par_ests[, 1], par_ests[, 2], pch = 16, type = "b")


