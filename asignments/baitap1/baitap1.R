# Bài tập 3. Có 46 vụ tràn dầu thô với khối lượng ít nhất là 1000 thùng từ tàu chở dầu trên
# vùng biển Hoa Kỳ trong giai đoạn 1974-1999. Tập tin oilspills.dat chứa dữ liệu sau:
# year – năm thứ i;
# spills – số vụ tràn dầu Ni trong năm thứ i;
# importexport – lượng dầu ước tính được vận chuyển qua vùng biển Hoa Kỳ như một
#phần của hoạt động xuất nhập khẩu của Hoa Kỳ trong năm thứ i, được điều chỉnh cho
#sự cố tràn ở vùng biển quốc tế hoặc nước ngoài, bi1 ;
# domestic – lượng dầu được vận chuyển qua vùng biển Hoa Kỳ trong các chuyến hàng
#trong nước trong năm thứ i, bi2 .

# Khối lượng dầu được vận chuyển là thước đo mức độ tiếp xúc với rủi ro tràn dầu. Giả sử chúng ta sử dụng giả định quá trình Poisson được đưa ra bởi Ni |bi1 , bi2 ∼ P(λi ) trong đó
# λi = α1 bi1 + α2 bi2 . Các tham số của mô hình này là α1 và α2 , biểu thị tỷ lệ xảy ra tràn dầu
# trên mỗi Bbbl dầu được vận chuyển trong quá trình xuất nhập khẩu và trong nước.


dataset <- read.table("Statistical_Computing/dataset/oilspills.dat", header = TRUE)
data <- list(N = dataset$spills, b1 = dataset$importexport, b2 = dataset$domestic)


# (a) Xác hàm log-likelihood cho các hệ số α1 và α2 .
# lambda function
lambda <- function(alpha1, alpha2, data) {
  return(alpha1 * data$b1 + alpha2 * data$b2)
}

# log-likelihood function
log_likelihood <- function(alpha1, alpha2, data) {
  lambda_i <- lambda(alpha1, alpha2, data)
  return(sum(log(lambda_i) * data$N - lambda_i))
}

# plot log-likelihood
alpha1 <- seq(0, 5, by = 0.1)
alpha2 <- seq(0, 5, by = 0.1)

logL_matrix <- matrix(0, nrow = length(alpha1), ncol = length(alpha2))

for (i in 1:length(alpha1)) {
  for (j in 1:length(alpha2)) {
    logL_matrix[i, j] <- log_likelihood(alpha1[i], alpha2[j], data)
  }
}
contour(alpha1, alpha2, logL_matrix, nlevels = 10, xlab = "alpha1", ylab = "alpha2")


# (b) Viết các bước giải thuật để xác định ước lượng MLE của α1 và α2 bằng phương pháp
#Newton. Viết một chương trình thi hành giải thuật này.
log_likelihood_prime <- function(alpha1, alpha2, data) {
  lambda_i <- lambda(alpha1, alpha2, data)
  grad1 <- sum(data$N * data$b1 / lambda_i - data$b1)
  grad2 <- sum(data$N * data$b2 / lambda_i - data$b2)
  return(c(grad1, grad2))
}

log_likelihood_2prime <- function(alpha1, alpha2, data) {
  lambda_i <- lambda(alpha1, alpha2, data)
  h11 <- sum(-data$N * (data$b1^2) / lambda_i^2)
  h22 <- sum(-data$N * (data$b2^2) / lambda_i^2)
  h12 <- sum(-data$N * (data$b1 * data$b2) / lambda_i^2)
  return(matrix(c(h11, h12, h12, h22), nrow = 2, byrow = TRUE))
}

newton_raphson <- function(alpha_init, data, max_iter = 100, tol = 1e-3) {
  alpha <- alpha_init
  vals <- matrix(NA, nrow = max_iter + 1, ncol = 2)
  vals[1,] <- alpha
  for (i in 1:max_iter) {
    gradient <- log_likelihood_prime(alpha[1], alpha[2], data)
    hessian <- log_likelihood_2prime(alpha[1], alpha[2], data)
    update <- solve(hessian) %*% gradient
    alpha <- alpha - update
    vals[i + 1,] <- alpha
    if (sqrt(sum(update^2)) < tol)
      print(paste("Converged at iteration", i))
    break
  }
  return(vals)
}

alpha_init <- c(0.1, 0.1)
alpha_mle_newton <- newton_raphson(alpha_init, data)
contour(alpha1, alpha2, logL_matrix, nlevels = 10, xlab = "alpha1", ylab = "alpha2")
points(alpha_mle_newton[, 1], alpha_mle_newton[, 2], pch = 16, type = "b")
print(paste("alpha_mle_newton:", alpha_mle_newton[-1,]))


#  (c) Viết các bước giải thuật để xác định ước lượng MLE của α1 và α2 bằng phương pháp
#Fisher Scoring. Viết một chương trình thi hành giải thuật này. So sánh với kết quả câu (b).
fisher_information <- function(alpha1, alpha2, data) {
  # return(-log_likelihood_2prime(alpha1, alpha2, data))
  lambda_i <- lambda(alpha1, alpha2, data)
  I11 <- sum(data$b1^2 / lambda_i)
  I22 <- sum(data$b2^2 / lambda_i)
  I12 <- sum((data$b1 * data$b2) / lambda_i)
  return(matrix(c(I11, I12, I12, I22), nrow = 2, byrow = TRUE))
}

# Fisher Scoring method
fisher_scoring <- function(data, alpha_init, tol = 1e-3, max_iter = 100) {
  alpha <- alpha_init
  vals <- matrix(NA, nrow = max_iter + 1, ncol = 2)
  vals[1,] <- alpha
  for (iter in 1:max_iter) {
    grad <- log_likelihood_prime(alpha[1], alpha[2], data)
    fisher_info <- fisher_information(alpha[1], alpha[2], data)
    update <- solve(fisher_info) %*% grad
    alpha <- alpha + update
    vals[iter + 1,] <- alpha
    if (sqrt(sum(update^2)) < tol) break
  }
  return(vals)
}

alpha_init <- c(0.4, 0.6)
alpha_mle_fisher <- fisher_scoring(data, alpha_init)
contour(alpha1, alpha2, logL_matrix, nlevels = 10, xlab = "alpha1", ylab = "alpha2")
points(alpha_mle_fisher[1], alpha_mle_fisher[2], pch = 16, type = "b")
print(paste("alpha_mle_fisher:", alpha_mle_fisher))

difference <- abs(alpha_mle_newton - alpha_mle_fisher)
print(paste("Difference in alpha1:", difference[1]))
print(paste("Difference in alpha2:", difference[2]))


# (d) Ước lượng sai số chuẩn (standard errors) của ước lượng MLE.


# (e) Áp dụng phương pháp quasi-Newton với hai cách chọn M(0) : (1) ma trận đơn vị âm
#và (2) ma trận thông tin Fisher, −I α(0) . So sánh tính ổn định, tốc độ của hai cách
#chọn này.

quasi_newton <- function(data, alpha_init, M, tol = 1e-6, max_iter = 100) {
  iter_identity <- 0
  alpha <- alpha_init
  vals <- matrix(NA, nrow = max_iter + 1, ncol = 2)
  vals[1,] <- alpha
  for (iter in 1:max_iter) {
    grad <- log_likelihood_prime(alpha[1], alpha[2], data)
    update <- M %*% grad
    alpha <- alpha + as.vector(update)
    vals[iter + 1,] <- alpha
    if (sqrt(sum(update^2)) < tol)
      break
    iter_identity <- iter_identity + 1
  }
  print(paste("Converged at iter:", iter))
  return(list(alpha = alpha, vals = vals))
}

alpha_init <- c(0.3, 0.3)
M <- diag(c(-1, -1), 2, 2)
alpha_mle_quasi_newton <- quasi_newton(data, alpha_init, M, tol = 1e-6, max_iter = 1000)
contour(alpha1, alpha2, logL_matrix, nlevels = 10, xlab = "alpha1", ylab = "alpha2")
points(alpha_mle_quasi_newton$vals[, 1], alpha_mle_quasi_newton$vals[, 2], pch = 16, type = "b")
print(alpha_mle_quasi_newton$alpha)

# Quasi-Newton BFGS method
quasi_newton_bfgs <- function(data, alpha_init, M, tol = 1e-6, max_iter = 100) {
  alpha <- alpha_init
  learning_rate <- 1
  vals <- matrix(0, nrow = max_iter + 1, ncol = 2)
  vals[1,] <- alpha

  for (iter in 1:max_iter) {
    grad <- log_likelihood_prime(alpha[1], alpha[2], data)
    update <- solve(M) %*% grad
    alpha_new <- alpha - learning_rate - update
    # REDUCE ALPHA UNTIL A CORRECT STEP IS REACHED
    f_t <- log_likelihood(alpha_new[1], alpha_new[2], data)
    f_t0 <- log_likelihood(alpha[1], alpha[2], data)
    print(paste('check f_t and f_t0', f_t, f_t0))
    while (f_t < f_t0) {
      print(paste('check f_t and f_t0', f_t, f_t0))
      learning_rate <- learning_rate / 2
      alpha_new <- alpha - learning_rate * update
      f_t <- log_likelihood(alpha_new[1], alpha_new[2], data)
      # Check for NaN values
      if (is.nan(f_t)) {
        print("NaN encountered in log-likelihood calculation")
        break
      }
    }

    ## UPDATE M(t) - BFGS METHOD
    ht <- alpha_new - alpha
    ut <- log_likelihood_prime(alpha_new[1], alpha_new[2], data) - grad
    v <- ut - M %*% ht
    M_old <- M
    M <- M - ((M %*% ht %*% t(M %*% ht)) / (as.numeric(t(ht) %*% M %*% ht))) +
      ((ut %*% t(ut)) / (as.numeric(t(ht) %*% ut)))

    ## CHECK CONDITION TO UPDATE M
    if (abs((t(v) %*% ht)[1]) < tol) {
      M <- M_old
    }

    # Check for convergence
    if (sqrt(sum(update^2)) < tol) {
      print(paste("Converged at iteration", iter))
      break
    }
    ## RESET ALPHA
    alpha <- alpha_new
    vals[iter + 1,] <- alpha
  }
  if (iter == max_iter) {
    print(paste("Not converged after", max_iter, "iterations"))
  }
  return(list(alpha = alpha, vals = vals))
}


quasi_newton_bfgs <- function(data, alpha_init, M, tol = 1e-6, max_iter = 100) {
  alpha <- alpha_init
  learning_rate <- 1
  vals <- matrix(0, nrow = max_iter + 1, ncol = 2)
  vals[1,] <- alpha

  for (iter in 1:max_iter) {
    grad <- log_likelihood_prime(alpha[1], alpha[2], data)

    # Kiểm tra ma trận Fisher Information M có xác định dương không
    if (any(eigen(M)$values <= 0)) {
      M <- M + diag(rep(tol, 2))  # Thêm epsilon vào đường chéo nếu M không xác định dương
    }

    update <- solve(M) %*% grad
    alpha_new <- alpha - learning_rate * update

    # Kiểm tra alpha_new có hợp lệ không
    if (any(alpha_new <= 0)) {
      print("Warning: alpha_new có giá trị âm, điều chỉnh learning_rate")
      learning_rate <- learning_rate / 2
      next
    }

    f_t <- log_likelihood(alpha_new[1], alpha_new[2], data)
    f_t0 <- log_likelihood(alpha[1], alpha[2], data)

    # Kiểm tra NaN
    if (is.nan(f_t)) {
      print("NaN encountered in log-likelihood calculation, adjusting step size")
      learning_rate <- learning_rate / 2
      next
    }

    # Giảm learning rate nếu f_t < f_t0
    while (f_t < f_t0) {
      learning_rate <- learning_rate / 2
      alpha_new <- alpha - learning_rate * update
      f_t <- log_likelihood(alpha_new[1], alpha_new[2], data)

      if (is.nan(f_t) || learning_rate < 1e-6) {
        print("Stopping due to NaN or very small learning rate")
        return(list(alpha = alpha, vals = vals))
      }
    }

    # Cập nhật ma trận M theo phương pháp BFGS
    ht <- alpha_new - alpha
    ut <- log_likelihood_prime(alpha_new[1], alpha_new[2], data) - grad
    v <- ut - M %*% ht
    M_old <- M
    M <- M - ((M %*% ht %*% t(M %*% ht)) / as.numeric(t(ht) %*% M %*% ht)) +
      ((ut %*% t(ut)) / as.numeric(t(ht) %*% ut))

    # Kiểm tra điều kiện cập nhật M
    if (abs((t(v) %*% ht)[1]) < tol) {
      M <- M_old
    }

    # Kiểm tra hội tụ
    if (sqrt(sum(update^2)) < tol) {
      print(paste("Converged at iteration", iter))
      break
    }

    # Cập nhật alpha
    alpha <- alpha_new
    vals[iter + 1,] <- alpha
  }

  if (iter == max_iter) {
    print(paste("Not converged after", max_iter, "iterations"))
  }

  return(list(alpha = alpha, vals = vals))
}

# Initial parameters
alpha_init <- c(0.3, 0.3)
M <- -diag(length(alpha_init))
# M <- fisher_information(alpha_init[1], alpha_init[2], data)
# Run Quasi-Newton BFGS method
alpha_mle_quasi_newton_bfgs <- quasi_newton_bfgs(data, alpha_init, M, tol = 1e-6, max_iter = 1000)
print(alpha_mle_quasi_newton_bfgs$alpha)
# Plot results
contour(alpha1, alpha2, logL_matrix, nlevels = 10, xlab = "alpha1", ylab = "alpha2")
points(alpha_mle_quasi_newton_bfgs$vals[, 1], alpha_mle_quasi_newton_bfgs$vals[, 2], pch = 16, type = "b")
print(alpha_mle_quasi_newton_bfgs$alpha)

alpha_init <- c(0.1, 0.1)
M <- -fisher_information(alpha_init[1], alpha_init[2], data)
# Run Quasi-Newton BFGS method
alpha_mle_quasi_newton_bfgs <- quasi_newton_bfgs(data, alpha_init, M)
print(alpha_mle_quasi_newton_bfgs$alpha)
# Plot results
contour(alpha1, alpha2, logL_matrix, nlevels = 10, xlab = "alpha1", ylab = "alpha2")
points(alpha_mle_quasi_newton_bfgs$vals[, 1], alpha_mle_quasi_newton_bfgs$vals[, 2], pch = 16, type = "b")
print(alpha_mle_quasi_newton_bfgs$alpha)
# (f) Xây dựng một đồ thị biểu diễn vùng tối ưu của hàm log-likelihood và để so sánh các
#đường đi được thực hiện bởi các phương pháp được sử dụng trong (b), (c) và (e). Chọn
#vùng vẽ đồ thị và điểm bắt đầu để minh họa tốt nhất các tính năng về hiệu suất của
#thuật toán.