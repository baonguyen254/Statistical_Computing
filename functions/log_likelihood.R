# Define the log-likelihood function for Cauchy(θ, 1)
log_likelihood <- function(theta, data) {
  likelihoods <- 1 / (pi * (1 + (data - theta)^2))
  #Check for invalid likelihoods
  likelihoods[likelihoods <= 0] <- 1e-10 #Handle cases where likelihood is zero or negative.  Replace with a small positive value.
  sum(log(likelihoods))
}