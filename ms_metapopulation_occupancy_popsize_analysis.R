library(rstan)
library(bayesplot)

## estimating effect of habitat management on abundance 
## using a hierarchical multispecies, metapopulation occupancy model 
## applied to real data

### Model 1: a single covariate for binary habitat management (mowed v unmowed)
### as a predictor of actual abundance and detection ability;
### and partially-pooling variation in detection probability and abundance
### by species identity.

### Simulating data
# consider the effect of site covariates on lambda and on detection probability

# Define function for generating binom-mix model data
data.fn <- function(R = 10, T = 3, # R = n sites, T = number of repeat visits
                    alpha0 = 1.5, alpha1 = 2, # site-level effect on abundance
                    beta0 = 0, beta1 = -5, # site-level effect on detection
                    n_species = 7, # number of species to pool
                    u = 0, sigma_u = 1, # species-level effect on abundance
                    v = 0, sigma_v = 1, # species-level effect on detection
                    rho_uv = 0 # correlation between species abundance and detection
                    ){
  # R: number of sites at which counts were made (= number of spatial reps)
  # T: number of times that counts were made at each site
  # (= number of temporal reps)
  # prob = probability that a site is mowed or unmowed
  # alpha0 and alpha1: intercept and slope of log-linear regression
  # relating abundance to the site covariate A
  # beta0 and beta1: intercept and slope of logistic-linear regression
  # of detection probability on A
  
  y <- array(dim = c(R, T)) # Array for counts
  
  # Ecological process
  # Covariate values: sort for ease of presentation
  X <- sort(rep(c(0, 1), times = R/2)) # half sites in each category
  
  # set up a multivariate matrix for effect of species 
  # on variation in abundance and detection 
  mu_uv <- c(u, v)
  sigma_uv <- matrix(NA, nrow = 2, ncol = 2)
  sigma_uv[1,1] = (sigma_u)^2
  sigma_uv[2,2] = (sigma_v)^2
  sigma_uv[1,2] = sigma_u * sigma_v * rho_uv
  sigma_uv[2,1] = sigma_u * sigma_v * rho_uv
  
  # fill a vector of species effects on abundance and detection
  uv <- MASS::mvrnorm(n = n_species, mu = mu_uv, Sigma = sigma_uv, empirical = FALSE)
  uv
  (rho_uv_simmed <- cor(uv[,1], uv[,2])) # correlation of site occupancy and detection
  
  hist(uv[,1], xlim = c(-6,6)) # most species near the mean of abundance
  hist(uv[,2], xlim = c(-6,6)) # most species near the mean of detection
  plot(uv[,2] ~ uv[,1])
  lm2 <- lm(uv[,2] ~ uv[,1])
  abline(lm2)
  
  # Relationship expected abundance – covariate
  # determined in part by species-level effect on abundance
  Z <- matrix(NA, nrow = n_species, ncol = R)
  for(i in 1:n_species){
    for(j in 1:R){
      Z[i,j] <- exp(alpha0 + alpha1 * X[j] + uv[i, 1])
    }
  }
  # Add Poisson noise: draw N from Poisson(lambda)
  Z_with_sigma <- Z
  for(i in 1:n_species){
    for(j in 1:R){
      Z_with_sigma[i,j] <- rpois(n = 1, lambda = Z[i,j])
    }
  }
    
  Z_with_sigma # Distribution of abundances across sites
  sum(Z_with_sigma > 0) / (R*n_species) # Empirical occupancy 
  # (proportion of species X site combos w/ 1 or more individuals)
  Z_with_sigma <- sum(Z_with_sigma) ; Z_with_sigma
  
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

