# Load required libraries
library(circular)
library(CircStats)
library(ggplot2)
library(fitdistrplus)

# Define the CDD PDF
dcdd <- function(theta, mu, kappa, lambda) {
  R <- sqrt(kappa^2 + lambda^2)
  denom <- 2 * pi * besselI(R, nu = 0)  # Normalizing constant
  exponent <- kappa * cos(theta - mu) + lambda * sin(theta - mu)
  out <- exp(exponent) / denom
  return(out)
}

# Negative log-likelihood for CDD
cdd_neg_loglik <- function(params, data) {
  mu <- params[1]
  kappa <- params[2]
  lambda <- params[3]
  -sum(log(dcdd(data, mu, kappa, lambda)))
}

# Fit CDD to data
fit_cdd <- function(data) {
  start_params <- c(mean.circular(data), 1, 0)  # Initial guesses
  fit <- optim(
    par = start_params,
    fn = cdd_neg_loglik,
    data = data,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  return(fit$par)
}

# Fit von Mises to data
fit_vm <- function(data) {
  mu <- mean.circular(data)
  rho <- rho.circular(data)
  kappa <- A1inv(rho)  # Approximation for kappa
  return(c(mu, kappa))
}

# Simulation and density evaluation for von Mises
dvm <- function(theta, mu, kappa) {
  1 / (2 * pi * besselI(kappa, nu = 0)) * exp(kappa * cos(theta - mu))
}

# Load and preprocess the wind dataset
wind_data <- as.numeric(wind)  # Convert to numeric
wind_data <- circular(wind_data, units = "degrees", template = "geographics")
wind_data_radians <- conversion.circular(wind_data, units = "radians")  # Convert to radians

# Fit both distributions to the dataset
cdd_params <- fit_cdd(wind_data_radians)
vm_params <- fit_vm(wind_data_radians)

# Create density grids for comparison
theta_grid <- seq(0, 2 * pi, length.out = 360)
cdd_density <- dcdd(theta_grid, cdd_params[1], cdd_params[2], cdd_params[3])
vm_density <- dvm(theta_grid, vm_params[1], vm_params[2])

# Normalize densities for plotting
cdd_density <- cdd_density / max(cdd_density)
vm_density <- vm_density / max(vm_density)

# Visualize the dataset and fitted densities
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))

# Rose diagram of the data
rose.diag(wind_data_radians, bins = 36, main = "Rose Diagram of Wind Data")

# Overlay fitted densities
plot(circular(theta_grid), cdd_density, type = "l", col = "red", lwd = 2, main = "Fitted Densities")
lines(circular(theta_grid), vm_density, col = "blue", lwd = 2)
legend("topright", legend = c("CDD", "von Mises"), col = c("red", "blue"), lwd = 2)

# Print fitted parameters
cat("Fitted parameters for CDD:\n")
print(cdd_params)
cat("\nFitted parameters for von Mises:\n")
print(vm_params)

# Evaluate and compare goodness of fit
loglik_cdd <- -cdd_neg_loglik(cdd_params, wind_data_radians)
loglik_vm <- sum(log(dvm(wind_data_radians, vm_params[1], vm_params[2])))
cat("\nLog-likelihood for CDD:", loglik_cdd, "\n")
cat("Log-likelihood for von Mises:", loglik_vm, "\n")

