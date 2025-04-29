theta_true <- 0.2
x <- rexp(50, rate = theta_true)
z <- rbinom(50, size = 1, prob = 0.6) # Size , prod là phần trăm giá trị bị khuyết
x[z == 0] <- NA

head(x, 6)

n <- length(x)
n2 <- sum(is.na(x))
n1 <- n - n2
x1 <- na.omit(x) # Loại bỏ giá trị NA
theta_0 <- 0.01
theta_t <- theta_0
iter <- 500
res <- numeric(iter)
res <- theta_t
for (i in 2:iter) {
  theta_t <- n * theta_t / (sum(x1) * theta_t + n2)
  res[i] <- theta_t
  if (abs(res[i] - res[i - 1]) < 1e-6)
    print(paste("Converged at iteration", i))
  break
}
head(res, 6)

res <- c(theta_0, res)
err_res <- abs(diff(res)) / res[1:20]
cbind(res, erorr = c(NA, err_res))

#thuc hành 1.1: Kết hợp vòng lặp thuật toán EM với tiêu chuẩn dừng dựa trên đánh giá sai lệch giữa
#nghiệm trong các bước lặp, cho ví dụ 1.1. Và xác định số vòng lặp cần thiết cho tới khi hội tụ.
EM <- function(x, theta_0, tol = 1e-6, max_iter = 100) {
  n <- length(x)
  n2 <- sum(is.na(x))
  n1 <- n - n2
  x1 <- na.omit(x)
  theta_t <- theta_0
  vals <- numeric(max_iter)
  iter <- 1
  while (iter < max_iter) {
    theta_t1 <- (n * theta_t) / (sum(x1) * theta_t + n2)
    # print(theta_t1)
    vals[iter] <- theta_t1
    if (abs(theta_t1 - theta_t) < tol) {
      print(paste("Converged at iteration", iter))
      break
    }
    theta_t <- theta_t1
    iter <- iter + 1

  }
  return(list(theta = theta_t, vals = vals))
}

theta_true <- 0.2
x <- rexp(50, rate = theta_true)
z <- rbinom(50, size = 1, prob = 0.6) # Size , prod là phần trăm giá trị bị khuyết
x[z == 0] <- NA

results <- EM(x, theta_0 = 0.1, tol = 1e-8, max_iter = 50)
print(results$theta)
print(results$vals)

#thuc hanh 1.2: Tiến hành một mô phỏng Monte Carlo với số lần lặp mô phỏng là 1000 lần, với thông tin
#như trong ví dụ 1.1. Trong mỗi lần mô phỏng dữ liệu, áp dụng thuật toán EM, tìm ước lượng θb. Vẽ biểu đồ
#histogram và lấy trung bình của 1000 ước lượng θb so sánh với θ0 = 0.2

monte_carlo_EM <- function(n_sim, n, theta_true, prob_missing) {
  theta_0 <- 0.1
  theta_true <- 0.2
  res <- numeric(n_sim)
  for (i in 1:n_sim) {
    x <- rexp(n, rate = theta_true)
    z <- rbinom(n, size = 1, prob = prob_missing)
    x[z == 0] <- NA
    results <- EM(x, theta_0 = theta_0, tol = 1e-8, max_iter = 50)
    res[i] <- results$theta
  }
  return(res)
}

monte_carlo_results <- monte_carlo_EM(1000, 50, 0.2, 0.6)
hist(monte_carlo_results, main = "Histogram of Monte Carlo Estimates", xlab = "Theta", breaks = 30)
abline(v = 0.2, col = "red", lwd = 2)

# Vi du 1.2:

data(geyser, package = "MASS")
head(geyser)

hist(geyser$waiting, xlab = "waiting time (mins)", ylab = "density",
     probability = TRUE, main = "")

y <- geyser$waiting
n <- length(y)
theta_0 <- c(0.3, 55, 80, 4, 7)
iter <- 20
theta_t <- theta_0
res <- matrix(0, nrow = iter, ncol = 5)
for (i in 1:iter) {
  f_1 <- dnorm(y, mean = theta_t[2], sd = theta_t[4])
  f_2 <- dnorm(y, mean = theta_t[3], sd = theta_t[5])
  p1 <- f_1 * theta_t[1] / (f_1 * theta_t[1] + f_2 * (1 - theta_t[1]))
  sum_p1 <- sum(p1)
  p1_t <- mean(p1)
  mu1_t <- sum(y * p1) / sum_p1
  mu2_t <- sum(y * (1 - p1)) / (n - sum_p1)
  sig1_t <- sqrt(sum((y - mu1_t)^2 * p1) / sum_p1)
  sig2_t <- sqrt(sum((y - mu2_t)^2 * (1 - p1)) / (n - sum_p1))
  theta_t <- c(p1_t, mu1_t, mu2_t, sig1_t, sig2_t)
  res[i,] <- theta_t
}
# res
theta_est <- res[iter,]
theta_est

dnorm_mix <- function(x, p, mu, sig) {
  return(p * dnorm(x, mu[1], sig[1]) + (1 - p) * dnorm(x, mu[2], sig[2]))
}

y_grid <- seq(40, max(y), length.out = 101)
dens_est <- dnorm_mix(y_grid, p = theta_est[1], mu = theta_est[2:3], sig = theta_est[4:5])

hist(geyser$waiting, xlab = "waiting time (mins)", ylab = "density",
     probability = TRUE, main = "")
lines(y_grid, dens_est, col = "blue")

#thuc hanh 1.3:  Kết hợp vòng lặp thuật toán EM với tiêu chuẩn dừng dựa trên đánh giá sai lệch giữa
#nghiệm trong các bước lặp, cho ví dụ 1.2. Và xác định số vòng lặp cần thiết cho tới khi hội tụ.
EM_new <- function(x, theta_0, tol = 1e-6, max_iter = 100) {
  n <- length(x)
  theta_t <- theta_0
  vals <- matrix(0, nrow = max_iter, ncol = 5)
  for (i in 1:max_iter) {
    f_1 <- dnorm(x, mean = theta_t[2], sd = theta_t[4])
    f_2 <- dnorm(x, mean = theta_t[3], sd = theta_t[5])
    p1 <- f_1 * theta_t[1] / (f_1 * theta_t[1] + f_2 * (1 - theta_t[1]))
    sum_p1 <- sum(p1)
    p1_t <- mean(p1)
    mu1_t <- sum(x * p1) / sum_p1
    mu2_t <- sum(x * (1 - p1)) / (n - sum_p1)
    sig1_t <- sqrt(sum((x - mu1_t)^2 * p1) / sum_p1)
    sig2_t <- sqrt(sum((x - mu2_t)^2 * (1 - p1)) / (n - sum_p1))
    theta_t1 <- c(p1_t, mu1_t, mu2_t, sig1_t, sig2_t)
    vals[i,] <- theta_t1
    # if (sqrt(sum((theta_t1-theta_t)^2)) < tol) {
    # L2 norm
    if (sqrt(sum((theta_t1 - theta_t)^2)) < tol) {
      print(paste("Converged at iteration", i))
      break
    }
    theta_t <- theta_t1
  }
  return(list(theta = theta_t, vals = vals))
}

theta_0 <- c(0.3, 55, 80, 4, 7)
results <- EM_new(geyser$waiting, theta_0, tol = 1e-6, max_iter = 100)
print(results$theta)
# print(results$vals)
#Thực hành 1.4: Thử với các điểm bắt đầu khác nhau cho bài toán trong ví dụ 1.2.
theta_0 <- c(0.5, 50, 70, 5, 6)
results <- EM_new(geyser$waiting, theta_0, tol = 1e-6, max_iter = 100)
print(results$theta)

# Ví dụ 2.1: dữ liệu theo phân phối mũ Exp(λ)
# read data datasets/censoreddata.dat

if (!require(readr)) { install.packages("readr") }
library(readr)
data_censored <- read_delim(file = "/home/baron/BaoBao/master_2024/Statistical_Computing/dataset/censoreddata.dat")
y <- data_censored$obsvi
n <- length(y)

# MCEM
lambda_0 <- 1 / mean(y)
n_cens <- sum(data_censored$deltai == 0)
iter <- 30
lambda_t <- lambda_0
m_t <- 10
y_t <- y
mean_y_t <- numeric(m_t)
res <- numeric(iter)
for (i in 1:iter) {
  for (j in 1:m_t) {
    z <- rexp(n_cens, lambda_t) + data_censored$ci[data_censored$deltai == 0]
    y_t[data_censored$deltai == 0] <- z
    mean_y_t[j] <- mean(y_t)
  }
  lambda_t <- 1 / mean(mean_y_t)
  res[i] <- lambda_t
}
res
# par(mar = c(5, 4, 4, 2) + 0.1)
plot(1:iter, res, type = "b", pch = 16)
print(lambda_t)

#Thực hành 2.1: Thi hành thuật toán EM thông thường và so sánh với kết quả của MCEM.

# EM thông thường
EM_method <- function(y, lambda_0, m_t = 10, tol = 1e-6, max_iter = 100) {
  n <- length(y)
  n_cens <- sum(data_censored$deltai == 0)
  lambda_t <- lambda_0
  iter <- 1
  res <- numeric(max_iter)
  while (iter < max_iter) {
    lambda_t1 <- n / (n_cens / lambda_t + sum(y * data_censored$deltai + data_censored$ci * (1 - data_censored$deltai)))
    res[iter] <- lambda_t1
    if (abs(lambda_t1 - lambda_t) < tol) {
      print(paste("Converged at iteration", iter))
      break
    }
    lambda_t <- lambda_t1
    iter <- iter + 1
  }
  return(list(lambda = lambda_t, vals = res))
}

results_exp <- EM_method(y, lambda_0, m_t = 10, tol = 1e-6, max_iter = 400)
print(results_exp$lambda)


# MCEM new
MCEM_new <- function(y, lambda_0, m_t = 10, tol = 1e-6, max_iter = 100) {
  n <- length(y)
  n_cens <- sum(data_censored$deltai == 0)
  lambda_t <- lambda_0
  y_t <- y
  mean_y_t <- numeric(m_t)
  res <- numeric(max_iter)
  for (i in 1:max_iter) {
    for (j in 1:m_t) {
      z <- rexp(n_cens, lambda_t) + data_censored$ci[data_censored$deltai == 0]
      y_t[data_censored$deltai == 0] <- z
      mean_y_t[j] <- mean(y_t)
    }
    lambda_t1 <- 1 / mean(mean_y_t)
    res[i] <- lambda_t1
    if (abs(lambda_t1 - lambda_t) < tol) {
      print(paste("Converged at iteration", i))
      break
    }
    lambda_t <- lambda_t1
  }
  return(list(lambda = lambda_t, vals = res))
}

results_MCEM <- MCEM_new(y, lambda_0, tol = 1e-6, max_iter = 400)
print(results_MCEM$lambda)
# compare MCEM and EM
# par(mfrow = c(2, 1))
# plot(1:length(results_exp$vals), results_exp$vals, type = "b", pch = 16, main = "EM")
# plot(1:length(results_MCEM$vals), results_MCEM$vals, type = "b", pch = 16, main = "MCEM")

# thực hành 2.2: tăng m_t, chẳng hạn, 20, 50, 100, tăng dần m_t theo bước lặp, chẳng hạn, m_t nhỏ (5, 10, 15) ở các lần lặp đầu nhưng lớn dần (50, 100,
# 200, ...) ở các lần lặp sau; hay m_t = 51+⌊t/10⌋ với = ⌊u⌋ là phần nguyên của u.
# Hãy áp dụng các giải pháp này cho đoạn chương trình MCEM ở trên và so sánh các kết quả.

results_MCEM_20 <- MCEM_new(y, lambda_0, m_t = 20, tol = 1e-6, max_iter = 400)
print(results_MCEM_20$lambda)
results_MCEM_50 <- MCEM_new(y, lambda_0, m_t = 50, tol = 1e-6, max_iter = 400)
print(results_MCEM_50$lambda)
results_MCEM_100 <- MCEM_new(y, lambda_0, m_t = 100, tol = 1e-6, max_iter = 400)
print(results_MCEM_100$lambda)

# Ví dụ 2.2: EM gradient algorithm
theta_true <- 0.2
x <- rexp(50, rate = theta_true)
z <- rbinom(50, size = 1, prob = 0.6)
x[z == 0] <- NA
head(x, 6)

n <- length(x)
n2 <- sum(is.na(x))
n1 <- n - n2
x1 <- na.omit(x)
sum_x1 <- sum(x1)
theta_0 <- 0.01
theta_t <- theta_0
iter <- 20
res_emg <- numeric(iter)

l_obs_exp <- function(x1, theta_t) {
  length(x1) * log(theta_t) - theta_t * sum(x1)
}

for (i in 1:iter) {
  # l_obs at thetaˆ(t)
  alpha_t <- 2
  l_obs_old <- l_obs_exp(x1, theta_t)
  theta_t_new <- theta_t + alpha_t * theta_tˆ2 / n * (n1 / theta_t - sum_x1)
  # l_obs at thetaˆ(t + 1)
  ## adjust alpha_t
  if (theta_t_new > 0) {
    l_obs_new <- l_obs_exp(x1, theta_t_new)
  }
  while (theta_t_new <= 0 || l_obs_old > l_obs_new) {
    alpha_t <- alpha_t / 2
    theta_t_new <- theta_t + alpha_t * theta_tˆ2 / n * (n1 / theta_t - sum_x1)
    if (theta_t_new > 0) {
      l_obs_new <- l_obs_exp(x1, theta_t_new)
    }
  }
  theta_t <- theta_t_new
  res_emg[i] <- theta_t_new
}

res_emg <- c(theta_0, res_emg)
err_res_emg <- abs(diff(res_emg)) / res_emg[1:20]
cbind(res_emg, erorr = c(NA, err_res_emg))