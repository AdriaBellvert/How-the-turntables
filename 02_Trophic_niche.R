# Load required libraries with error handling
required_packages <- c("factoextra", "tidyr", "SIBER")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Package ", package, " is not available.")
    }
  }
}

library(factoextra)
library(tidyr)
library(SIBER)

# Load data
data <- read.table("data_isotops.txt", header = TRUE, sep = "\t")

# Convert factors and character columns
data$sp_island <- as.character(data$sp_island)
data$island <- as.factor(data$island)

# Prepare data
datamt <- data.frame(
  iso1 = as.numeric(data$d13C),
  iso2 = as.numeric(data$d15N),
  group = data$sp_island,
  community = data$island
)

# Create SIBER object
siberdata <- createSiberObject(datamt)

# Group metrics
group_ML <- groupMetricsML(siberdata)
print(group_ML)

# Set parameters
parms <- list(
  n.iter = 2 * 10^4,
  n.burnin = 1 * 10^3,
  n.thin = 10,
  n.chains = 10
)

# Define priors
priors <- list(
  R = diag(2),
  k = 2,
  tau.mu = 1.0E-3
)

# Calculate ellipses
ellipses_posterior <- siberMVN(siberdata, parms, priors)

# Plot ellipses
SEA_B <- siberEllipses(ellipses_posterior)
siberDensityPlot(SEA_B, xticklabels = colnames(group_ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2)),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group",
                 ylims = c(0, 6))

# Calculate probability between groups
Pg1.1_lt_g1.2 <- sum(SEA_B[, 1] < SEA_B[, 5]) / nrow(SEA_B)
print(Pg1.1_lt_g1.2)
