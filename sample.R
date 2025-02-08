f <- function(x) {
  log(x) / (1 + x)
}

f_prime <- function(x) {
  (1 + (1 / x) - log(x)) / ((1 + x)^2)
}

f_2prime <- function(x) {
  -(3 + 1 / x^2 + 4 / x - 2 * log(x)) / ((1 + x)^3)
}

## bisection method
bisection <- function(f, a, b, tol = 1e-6, max_iter = 40) {
  x <- a + (b - a) / 2
  res_a <- a
  res_b <- b
  res_x <- x
  for (i in 1:max_iter) {
    # x <- (a + b)/2
    if (f(x) < tol || (b - a) / 2 < tol) {
      print(paste("Converged at iteration", i))
      return(x)
    }
    if (f(a) * f(x) <= 0) {
      res_a[i + 1] <- a
      res_b[i + 1] <- b <- x
    } else {
      res_a[i + 1] <- a <- x
      res_b[i + 1] <- b
    }
    x <- a + (b - a) / 2
    res_x[i + 1] <- x
  }
  print(head(cbind(res_a, res_b, res_x), 20))
}

bisection(f_prime, 1, 5)


## newton-raphson method
newton_raphson <- function(f_prime, f_2prime, x0, tol = 1e-6, max_iter = 50) {
  x <- x0
  res_x_newt <- x
  for (i in 1:max_iter) {
    res_x_newt[i + 1] <- x <- x - f_prime(x) / f_2prime(x)
    if (abs(res_x_newt[i + 1] - res_x_newt[i]) < tol) {
      print(paste("Converged at iteration", i))
      # return(x)
      break
    }
  }
  print(head(cbind(res_x_newt), 20))
  return(x)
}

newton_raphson(f_prime, f_2prime, 3)

## secant method
secant <- function(f_prime, x0, x1, tol = 1e-9, max_iter = 40) {
  x <- x1
  res_x_sec <- x
  for (i in 1:max_iter) {
    res_x_sec[i + 1] <- x <- x - f_prime(x) * (x - x0) / (f_prime(x) - f_prime(x0))
    # update x0
    x0 <- res_x_sec[i]
    if (abs(res_x_sec[i + 1] - res_x_sec[i]) < tol) {
      print(paste("Converged at iteration", i))
      # return(x)
      break
    }
  }
  print(head(cbind(res_x_sec), 20))
  return(x)
}

secant(f_prime, 1, 5)

## fixed point iteration
fixed_point <- function(f_prime, x0, alpha, tol = 1e-9, max_iter = 100) {
  x <- x0
  res_x_fix <- x
  for (i in 1:max_iter) {
    # check f_prime(x) is NaN then break
    if (is.nan(f_prime(x))) {
      print(paste("f_prime(x) is NaN at iteration", i))
      break
    }
    res_x_fix[i + 1] <- x <- x + alpha * f_prime(x)
    # if (i == 11) {
    #     print(paste("alpha =", alpha))
    # }
    if (abs(res_x_fix[i + 1] - res_x_fix[i]) < tol) {
      print(paste("Converged at iteration", i))
      break
    }
  }
  print(head(cbind(res_x_fix), -10))
  return(x)
}

# alpha = 4 
fixed_point(f_prime, 3, 4, tol = 1e-6, max_iter = 200)

# alpha = 1 
fixed_point(f_prime, 1, 1, tol = 1e-4, max_iter = 200)


## Fisher Scoring method
log_likelihood <- function(x, theta) {
  sum(log(1 / (pi * (1 + (x - theta)^2))))
}

fisher_scoring <- function(log_likelihood, x, theta0, tol = 1e-6, max_iter = 40) {
  theta <- theta0
  res_theta_fisher <- theta
  for (i in 1:max_iter) {
    res_theta_fisher[i + 1] <- theta <- theta - f_prime(theta) / f_2prime(theta)
  }
  print(head(cbind(res_theta_fisher), 20))
}

fisher_scoring(log_likelihood, x, 1)