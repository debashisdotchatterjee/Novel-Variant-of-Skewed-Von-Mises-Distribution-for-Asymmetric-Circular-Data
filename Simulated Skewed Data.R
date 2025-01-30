# Install required packages if not already installed
# install.packages(c("circular", "CircStats", "ggplot2"))

library(circular)   # For circular data handling and visualization
library(CircStats)  # For fitting von Mises distribution
library(ggplot2)    # For advanced plotting
set.seed(123)  # For reproducibility

# Generate data from two von Mises distributions
n1 <- 300
n2 <- 200
angles1 <- rvonmises(n1, mu = circular(pi/2), kappa = 4)         # Concentrated around 90° with high concentration
angles2 <- rvonmises(n2, mu = circular(pi/2 + 0.5), kappa = 2)   # Slightly shifted mean, lower concentration

# Combine and convert to numeric values
skewed_data <- c(angles1, angles2)
skewed_angles <- as.numeric(skewed_data)
rose.diag(circular(skewed_angles), bins = 36, main = "Rose Diagram of Skewed Data")
fit_vm_skewed <- mle.vonmises(circular(skewed_angles))
print(fit_vm_skewed)
