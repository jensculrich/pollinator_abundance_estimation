library(rstan)
library(bayesplot)

#### mark_recapture to estimate population sizes in a meta population
#### with a site covariate (plant community/management)
# Building from text and examples provided in Kery & Schaub 2011: 
# Bayesian population analysis using WinBUGS: a hierarchical perspective

# JCU May, 2022

### Model 1: the simplest case with constant parameters
# both the ecological (population size) and observation (detection) processes
# described only by an intercept.

## Start by simulating data
# To generate count data y under this null model for 
# R = 200 spatial replicates (sites)
# and T = 3 temporal replicates:

# Determine sample sizes (spatial and temporal replication)
R <- 200
T <- 3
lambda <- 2 # mean abundance per site of 2
p <- 0.5 # mean detection per individual of 0.5

# Create structure to contain counts
y <- array(dim = c(R, T))
# Sample abundance from a Poisson(lambda = 2)
N <- rpois(n = R, lambda = lambda)
# Sample counts from a Binomial(N, p = 0.5)
for (j in 1:T){
  y[,j] <- rbinom(n = R, size = N, prob = p)
}

# Look at realization of biological and observation processes
cbind(N, y)

## Fit the model to the simulated data
# Bundle data
# need to define a maximum population size K
K <- 100
stan_data <- c("y", "R", "T", "K")

# Initial values
Nst <- apply(y, 1, max) + 1 # This line is important
inits <- function() list(N = Nst)

# Parameters monitored
params <- c("lambda", "p")

# MCMC settings
n_iterations <- 1200
n_thin <- 1
n_burnin <- 200
n_chains <- 4

# Call STAN model from R 
stan_model <- "./mark_recapture_binomial_mixture_model.stan"

## Call Stan from R
out <- stan(stan_model,
            data = stan_data, 
            init = inits, 
            pars = params,
            chains = n_chains, iter = n_iterations, 
            warmup = n_burnin, thin = n_thin,
            seed = 2,
            open_progress = FALSE)

## Summarize posteriors
print(out, digits = 3)
traceplot(out, pars = c("lambda", "p"), inc_warmup = TRUE, nrow = 2)

color_scheme_set("pink")
mcmc_dens_overlay(out, pars = c("lambda", "p"))

### Model 2: Introducing a covariate
# consider the effect of site covariates on lambda and on detection probability

# Define function for generating binom-mix model data
data.fn <- function(R = 200, T = 3, xmin = -1, xmax = 1, alpha0 = 1,
                    alpha1 = 3, beta0 = 0, beta1 = -5){
  # R: number of sites at which counts were made (= number of spatial reps)
  # T: number of times that counts were made at each site
  # (= number of temporal reps)
  # xmin, xmax: define range of the covariate X
  # alpha0 and alpha1: intercept and slope of log-linear regression
  # relating abundance to the site covariate A
  # beta0 and beta1: intercept and slope of logistic-linear regression
  # of detection probability on A
  
  y <- array(dim = c(R, T)) # Array for counts
  
  # Ecological process
  # Covariate values: sort for ease of presentation
  X <- sort(runif(n = R, min = xmin, max = xmax))
  # Relationship expected abundance – covariate
  lam <- exp(alpha0 + alpha1 * X)
  # Add Poisson noise: draw N from Poisson(lambda)
  N <- rpois(n = R, lambda = lam)
  table(N) # Distribution of abundances across sites
  sum(N > 0) / R # Empirical occupancy (proportion of sites w/ 1 or more individuals)
  totalN <- sum(N) ; totalN
  
  # Observation process
  # Relationship detection prob – covariate
  p <- plogis(beta0 + beta1 * X)
  # Make a 'census' (i.e., go out and count things)
  for (i in 1:T){
    y[,i] <- rbinom(n = R, size = N, prob = p)
  }
  
  # Naïve regression
  naive.pred <- exp(predict(glm(apply(y, 1, max) ~ X + I(X^2),
                                family = poisson)))
  
  # Plot features of the simulated system
  par(mfrow = c(2, 2))
  plot(X, lam, main = "Expected abundance", xlab = "Covariate",
       ylab = "lambda", las = 1, type = "l", col = "red", lwd = 3,
       frame.plot = FALSE)
  plot(X, N, main = "Realised abundance", xlab = "Covariate", ylab =
         "N", las = 1, frame.plot = FALSE, col = "red", cex = 1.2)
  plot(X, p, ylim = c(0, 1), main = "Detection probability", xlab =
         "Covariate", ylab = "p", type = "l", col = "red", lwd = 3, las = 1,
       frame.plot = FALSE)
  plot(X, naive.pred, main = "Actual counts \n and naïve regression",
       xlab = "Covariate", ylab = "Relative abundance", ylim = c(min(y),
                                                                 max(y)), type = "l", lty = 2, lwd = 4, col = "blue", las = 1,
       frame.plot = FALSE)
  points(rep(X, T), y, col = "black", cex = 1.2)
  
  # Return stuff
  return(list(R = R, T = T, X = X, alpha0 = alpha0, alpha1 = alpha1,
              beta0 = beta0, beta1 = beta1, lam = lam, N = N, totalN = totalN,
              p = p, y = y))
}

data <- data.fn()
str(data)

## Fit the model to the simulated data
# Bundle data
# need to define a maximum population size K
K <- 100

y <- data$y
X <- data$X
R <- data$R
T <- data$T
stan_data <- c("y", "X", "R", "T", "K")

# Parameters monitored
params <- c("totalN", "alpha0", "alpha1", "beta0", "beta1")
# Could try adding "N" to track site abundances?

# MCMC settings
n_iterations <- 1000
n_thin <- 1
n_burnin <- 500
n_chains <- 4

## Initial values
inits <- lapply(1:n_chains, function(i)
  list(alpha0 = runif(1, -1, 1),
       alpha1 = runif(1, -1, 1),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, -1, 1)))

# Call STAN model from R 
stan_model_covariate <- "./mark_recapture_binomial_mixture_model_covariate.stan"

## Call Stan from R
out_cov <- stan(stan_model_covariate,
            data = stan_data, 
            init = inits, 
            pars = params,
            chains = n_chains, iter = n_iterations, 
            warmup = n_burnin, thin = n_thin,
            seed = 2,
            open_progress = FALSE)

print(out_cov, digits = 3)

# sum of the maximum counts at each site, a conventional estimator of total 
# population size, is only 233 individuals 
observed_counts <- sum(apply(data$y, 1, max))
# mean totalN for this run of the model was 1848 (sd = 528).
# totalN in simulated data run was 1764. So does pretty good! even though the sd is large

# plot posterior distribution
p <- mcmc_hist(out_cov, pars = c("alpha0"))
p <- p + labs(x = "alpha0",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data$alpha0, linetype = "solid", size = 1)
p

mcmc_dens_overlay(out_cov, pars = c("alpha0", "alpha1", 
                                "beta0", "beta1"))
mcmc_pairs(out_cov, pars = c("alpha0", "alpha1", 
                             "beta0", "beta1"),
           off_diag_args = list(size = 1.5))

# Plot predicted covariate relationship with abundance
matrix <- as.matrix(out_cov)
str(matrix)
alpha0_mean = mean(matrix[,2])
alpha1_mean = mean(matrix[,3])

par(mfrow = c(1, 1))
plot(data$X, data$N, main = "", xlab = "Covariate", ylab = "Abundance",
     las = 1, ylim = c(0, max(data$N)), frame.plot = FALSE)
lines(data$X, data$lam, type = "l", col = "red", lwd = 3)
points(rep(X, T), y, col = "green", cex = 1.2)
GLM.pred <- exp(predict(glm(apply(data$y, 1, max) ~ X + I(X^2),
                            family = poisson, data = data)))
lines(data$X, GLM.pred, type = "l", lty = 2, col = "blue", lwd = 3)
Nmix.pred <- exp(alpha0_mean + alpha1_mean * data$X)
lines(data$X, Nmix.pred, type = "l", col = "blue", lwd = 3)

