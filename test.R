---
  title:"test"
author:"Nguyen Van A"
date:"2025-01-11"
output:html_document
---


  # Cài đặt các gói cần thiết (nếu chưa có)
  if (!require(ggplot2)) { install.packages("ggplot2") }
if (!require(Deriv)) { install.packages("Deriv") }
if (!require(languageserver)) { install.packages("languageserver") }


# Load the function from derivative_function.R
source("functions/daoham.R")
# Define the function to differentiate
f <- function(x) log(x) / (1 + x)

# Generate x values
x_values <- seq(1, 8, 0.1)

f(x_values)
plot(x_values, f(x_values), type = "l")
image(x_values, f(x_values), type = "l")


#------ Bisection method ------
# fx_prime <- function(x) calculate_derivative(f, x)

# plot(x_values, fx_prime(x_values), type = "l")
# abline(h = 0, col = "red", lty = 2)

# step by step with bisection method
a0 <- 1
b0 <- 5
source("functions/bisection_method.R")
bisection_method(f, a0, b0)


#------ Newton-Raphson method ------
source("functions/Newton-Raphson.R")
newton_raphson(f, 1)

#------ Secant method ---------


#------ Fixed point method --------


# # Calculate the derivative
# derivatives <- calculate_derivative(my_function, x_values)

# #Now you can use the 'derivatives' vector.  For example, to plot:
# data <- data.frame(x = x_values, y = my_function(x_values), y_prime = derivatives)
# ggplot(data, aes(x = x)) +
#   geom_line(aes(y = y, color = "Hàm số gốc"), size = 1) +
#   geom_line(aes(y = y_prime, color = "Đạo hàm"), size = 1) +
#   labs(title = "Đồ thị hàm số và đạo hàm",
#        x = "x",
#        y = "y",
#        color = "Legend") +
#   scale_color_manual(values = c("Hàm số gốc" = "blue", "Đạo hàm" = "red")) +
#   theme_bw()

# hist(rnorm(100))
# plot(x <- sort(rnorm(47)), type = "s", main = "plot(x, type = \"s\")")
# points(x, cex = .5, col = "dark red")