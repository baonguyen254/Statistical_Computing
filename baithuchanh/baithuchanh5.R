# mô phỏng, tạo mẫu ngẫu nhiên

I_est <- function(data) {
  return(mean(3 - data^2))
}

# rnom mẫu ngẫu nhiên
x1 <- rnorm(1000, mean = 0, sd = 1)
hist(y, breaks = 30, main = "Histogram of y", xlab = "y")
I <- I_est(x1)
I
hist(y, breaks = 30, main = paste("I estimate: ", I), xlab = "y")
qqnorm(x1)
qqline(x1)

# transform
my_rnorm <- function(n, mean = 0, sd = 1) {
  u1 <- runif(n = n)
  u2 <- runif(n = n)
  x1 <- mean + sd * sqrt(-2 * log(u1)) * cos(2 * pi * u2)
  x2 <- mean + sd * sqrt(-2 * log(u1)) * sin(2 * pi * u2)
  x <- c(x1, x2)
  return(sample(x, size = n))
}

x2 <- my_rnorm(1000, mean = 0, sd = 1)
hist(x2, breaks = 30, main = "Histogram of x2", xlab = "x2")

qqnorm(x2)
qqline(x2)

I_est(x2)

# nội suy tuyến tính
my_rnorm2 <- function(n, mean = 0, sd = 1) {
  # chia lưới
  x_grid <- seq(mean - 5 * sd, mean + 5 * sd, by = 0.01)
  u_grid <- pnorm(x_grid, mean, sd)
  U <- runif(n)
  X <- numeric(n)
  for (i in 1:n) {
    v1 <- which(U[i] > u_grid)
    u_i <- u_grid[length(v1)]
    u_j <- u_grid[length(v1) + 1]
    x_i <- x_grid[length(v1)]
    x_j <- x_grid[length(v1) + 1]
    X[i] <- (u_j - U[i]) * x_i / (u_j - u_i) + (U[i] - u_i) * x_j / (u_j - u_i)
  }
  return(X)
}

x3 <- my_rnorm2(1000, mean = 0, sd = 1)
hist(x3, breaks = 30, main = "Histogram of x3", xlab = "x3")

qqnorm(x3)
qqline(x3)

# mô phỏng monte carlo 10 lần
n_sim <- 100
total_1 <- numeric(n_sim)
total_2 <- numeric(n_sim)
total_3 <- numeric(n_sim)
for (i in 1:n_sim) {
  total_1[i] <- I_est(my_rnorm2(100, mean = 0, sd = 1))
  total_2[i] <- I_est(my_rnorm2(100, mean = 0, sd = 1))
  total_3[i] <- I_est(my_rnorm2(100, mean = 0, sd = 1))
}
print(mean(total_1))
print(mean(total_2))
print(mean(total_3))

out_MC <- sapply(1:1000, FUN = function(i, n) {
  x_data <- my_rnorm2(n, mean = 0, sd = 1)
  return(I_est(x_data))
}, n = 15)
mean(out_MC)
var(out_MC)
hist(out_MC, breaks = 30, main = "Histogram of I estimate", xlab = "I estimate")

# xấp xỉ số pi
pi_est <- function(x1, x2) {
  return(mean(x1^2 + x2^2 <= 1) * 4)
}

n_sim <- 5000
x1 <- runif(n_sim, min = -1, max = 1) # r uniform
x2 <- runif(n_sim, min = -1, max = 1)
pi_t <- pi_est(x1, x2)
pi_t

# mô phỏng monte carlo 1000 lần pi_est
n_sim <- 50000
total_pi <- numeric(n_sim)
for (i in 1:n_sim) {
  x1 <- runif(n_sim, min = -1, max = 1)
  x2 <- runif(n_sim, min = -1, max = 1)
  total_pi[i] <- pi_est(x1, x2)
}
mean(total_pi)

data <- sapply(1:5000, FUN = function(i, n) {
  x1 <- runif(n_sim, min = -1, max = 1)
  x2 <- runif(n_sim, min = -1, max = 1)
  return(pi_est(x1, x2))
}, n = 1000)

mean(data)
# plot circle est pi
plot(x1, x2, col = ifelse(x1^2 + x2^2 <= 1, "red", "blue"), pch = 19, asp = 1)
points(0, 0, col = "green", pch = 19)
points(mean(x1), mean(x2), col = "black", pch = 19)

data_circle <- data.frame(x = x1, y = x2, inside = (x1^2 + x2^2 <= 1))

library(ggplot2)
ggplot(data = data_circle, mapping = aes(x = x, y = y, colour = inside)) +
  geom_point() +
  stat_function(fun = function(x)sqrt(1 - x^2), inherit.aes = FALSE) +
  stat_function(fun = function(x)-sqrt(1 - x^2), inherit.aes = FALSE) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  scale_color_manual(breaks = c("TRUE", "FALSE"),
                     values = c("blue", "gray")) +
  theme_bw()


# location scale t-student
x_t <- rt(10, df = 3)

post_mu <- function(data, mu) {
  prod(3 + (data - mu)^2)^(-2)
}

post_mu <- Vectorize(FUN = post_mu, vectorize.args = "mu")
post_mu(x_t, 2)
post_mu(x_t, c(2, 2.1))

mu_grid <- seq(-3, 3, by = 0.01)
post_mu_est <- post_mu(x_t, mu_grid)
plot(mu_grid, post_mu_est, type = "l", xlab = "mu", ylab = "posterior", main = "Posterior of mu")

x_t
mu_sim <- x_t[3] + rt(10, df = 3)
num <- sum(mu_sim * post_mu(x_t[-3], mu_sim))
den <- sum(post_mu(x_t[-3], mu_sim))
num / den

plot(mu_grid, post_mu_est, type = "l", main = "Posterior of mu2")
abline(v = num / den, col = "red", lwd = 2)