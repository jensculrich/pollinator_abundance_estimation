library(rstan)
library(bayesplot)
library(tidyverse)
library(RColorBrewer)
library(lme4)

##########################
# Model M0 ###############
##########################

## highly simplified N-mix model to estimate the effect of habitat management on abundance
## used to test basic model output and goodness of fit (to see if problems in goodness in fit 
## are arising from multi-species model structure)

# in contrast to M2 includes NO random effects for species 

# As in M1 and M2, M0 considers abundance as arising from a simple Poisson distribution
# as opposed to M3, which uses a Zero inflated poisson distribution,
# or potential M4 which uses a negative-binomial dist to account for 
# overdispersion in abundance and/or detection.

##########################
### Simulate data ########
##########################

##########################
### Simulate data ########
##########################

# consider the effect of site covariates on lambda and on detection probability

# Define function for generating binom-mix model data
data.fn.2 <- function(n_sites = 100, #  R = sites, 
                      T = 3, # T = number of repeat visits (= number of temporal reps)
                      alpha0 = 0, # grand intercept abundance
                      alpha1 = 2, # site-level effect on abundance
                      beta0 = 0, # grand intercept detection probability
                      beta1 = -2 # site-level effect on detection
){
  # alpha0 and alpha1: intercept and slope of log-linear regression (abundance)
  # beta0 and beta1: intercept and slope of logistic-linear regression (detection probability)
  
  # number of unique data observation points (observations of unique species at unique site)
  R = n_sites
  
  # Covariate values: sort for ease of presentation
  X <- sort(rep(c(0, 1), times = (R/2))) # half sites in each category
  
  ## Ecological process
  
  # Relationship expected abundance – covariate + variation in species
  # determined in part by species-level effect on abundance
  Z <- vector(length=R)
  for(i in 1:length(Z)){
    Z[i] <- exp(alpha0 + alpha1 * X[i])
  }
  
  # Add Poisson noise: draw N from Poisson(lambda)
  Z_with_sigma <- Z
  for(i in 1:length(Z)){
    Z_with_sigma[i] <- rpois(n = 1, lambda = Z[i])
  }
  
  Z_with_sigma # Distribution of abundances across sites
  sum(Z_with_sigma > 0) / (R) # Empirical occupancy 
  totalN <- sum(Z_with_sigma)
  totalN
  
  ## Observation process
  # ilogit <- function(u) { return(1 / (1 + exp(-u))); }
  p <- vector(length=R)
  for(i in 1:length(p)){
    p[i] <- plogis(beta0 + beta1 * X[i])
  }
  p
  
  par(mfrow = c(2, 1))
  hist(p[1:(R/2)], xlim = c(0,1))
  hist(p[(R/2 + 1):(R)], xlim = c(0,1))
  
  # Make a 'census' (i.e., go out and count things)
  y <- matrix(NA, nrow = length(p), ncol = T) # Array for counts
  for(i in 1:nrow(y)){
    for(j in 1:T){
      y[i,j] <- rbinom(n = 1, size = Z_with_sigma[i], prob = p[i])
    }
  }
  
  # Plot features of the simulated system
  # need to work on this
  
  # Return stuff
  return(list(T = T, 
              X = X,
              n_sites = n_sites,
              alpha0 = alpha0, 
              alpha1 = alpha1,
              beta0 = beta0, 
              beta1 = beta1,
              Z = Z, 
              Z_with_sigma = Z_with_sigma, 
              totalN = totalN,
              p = p, 
              y = y
  ))
}

data_set <- data.fn.2()
str(data_set)

## Fit the model to the simulated data
# Bundle data
# need to define a maximum population size K
max(data_set$Z_with_sigma)
K <- as.integer(max(data_set$Z_with_sigma) * 1.5)

y <- data_set$y
X <- data_set$X
R <- data_set$n_sites
T <- data_set$T
totalN <- data_set$totalN

stan_data <- c("R", 
               "T", "y", "X", "K")

# Parameters monitored
params <- c("totalN",
            "alpha0",
            "alpha1",
            "beta0",
            "beta1",
            "fit",
            "fit_new"
)


# MCMC settings
n_iterations <- 1000
n_thin <- 1
n_burnin <- 500
n_chains <- 3
n_cores <- 3

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  list(alpha0 = runif(1, -1, 1),
       alpha1 = runif(1, -3, 3),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, -3, 3)
  )
)

# Call STAN model from R 
stan_model <- "./parks_pollinators_abundance_estimation/M0_Nmix.stan"

## Call Stan from R
M0_out_sim <- stan(stan_model,
                   data = stan_data, 
                   init = inits, 
                   pars = params,
                   chains = n_chains, iter = n_iterations, 
                   warmup = n_burnin, thin = n_thin,
                   seed = 1,
                   open_progress = FALSE,
                   cores = n_cores)

print(M0_out_sim, digits = 3)

saveRDS(M0_out_sim, "./parks_pollinators_abundance_estimation/M0_out_sim.rds")
M0_out_sim <- readRDS("./parks_pollinators_abundance_estimation/M0_out_sim.rds")

# plot posterior distribution
p <- mcmc_hist(M0_out_sim, pars = c("alpha0"))
p <- p + labs(x = "alpha0",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$alpha0, linetype = "solid", size = 1)
p

p <- mcmc_hist(M0_out_sim, pars = c("alpha1"))
p <- p + labs(x = "alpha1",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$alpha1, linetype = "solid", size = 1)
p

p <- mcmc_hist(M0_out_sim, pars = c("beta0"))
p <- p + labs(x = "beta0",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$beta0, linetype = "solid", size = 1)
p

p <- mcmc_hist(M0_out_sim, pars = c("beta1"))
p <- p + labs(x = "beta1",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$beta1, linetype = "solid", size = 1)
p

mcmc_dens_overlay(M0_out_sim, pars = c("alpha0", "alpha1", 
                                       "beta0", "beta1"))

traceplot(M0_out_sim, pars = c("alpha0", "alpha1", 
                               "beta0", "beta1"))

mcmc_pairs(M0_out_sim, pars = c("alpha0", "alpha1", 
                                "beta0", "beta1"),
           off_diag_args = list(size = 1.5))

mcmc_dens_overlay(M0_out_sim, pars = c("totalN"))


# plot posterior distribution totalN
q <- mcmc_hist(M0_out_sim, pars = c("totalN"))
q <- q + labs(x = "totalN",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN, linetype = "solid", size = 1)
q


# could set up a loop to draw plots for each species and then cowplot them together

# Evaluation of fit
list_of_draws_M0_out_sim <- as.data.frame(M0_out_sim)

par(mfrow = c(1, 1))

plot(list_of_draws_M0_out_sim$fit, list_of_draws_M0_out_sim$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 400),
     xlim = c(0, 400))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws_M0_out_sim$fit_new > list_of_draws_M0_out_sim$fit)
mean(list_of_draws_M0_out_sim$fit) / mean(list_of_draws_M0_out_sim$fit_new)

## let's see if we can do a posterior predictive check based on the parameter estimates from 
## the output. Maybe this will help to see if something is wrong with the posterior check 
## in the generated quantities block versus the model itself..


# inv logit transformation function
ilogit <- function(x){
  exp(x)/(1+exp(x))
}

# fit and fit new in one step of the mcmc search
fit = 0
fit_new = 0

eval <- matrix(nrow = R, ncol = T)

y <- data_set$y
E <- matrix(nrow = R, ncol = T)

y_new <- matrix(nrow = R, ncol = T)
E_new <-  matrix(nrow = R, ncol = T)

sim_num <- 52

for(i in 1:R){
  for(j in 1:T){
    # eval = N[i] * p[i]
    eval[i, 1:j] = (exp(list_of_draws_M0_out_sim[sim_num,]$alpha0 + 
                                 (list_of_draws_M0_out_sim[sim_num,]$alpha1 * X[i])
                                 ) * # N[i]
                   ilogit(list_of_draws_M0_out_sim[sim_num,]$beta0 +
                             list_of_draws_M0_out_sim[sim_num,]$beta1 * X[i]
                           ) # p[i]
    )
    
    E[i, j] <- (y[i, j] - eval[i, j])^2 / (eval[i, j] + 0.5)
    fit = sum(E)
    
    y_new[i,j] <- rbinom(1, 
                           rpois(1, exp(list_of_draws_M0_out_sim[sim_num,]$alpha0 + 
                                          (list_of_draws_M0_out_sim[sim_num,]$alpha1 * X[i]) 
                           )),
                           ilogit(list_of_draws_M0_out_sim[sim_num,]$beta0 +
                                    list_of_draws_M0_out_sim[sim_num,]$beta1 * X[i]
                           )
    )
    
    E_new[i, 1:j] = (y_new[i, j] - eval[i, j])^2 / (eval[i, j] + 0.5)
    fit_new = sum(E_new)
  }
}
    

par(mfrow = c(1, 1))
plot(fit, fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 500),
     xlim = c(0, 500))
abline(0, 1, lwd = 2, col = "black")


# Expected values for counts across all mcmc steps
# this seems to work as wanted to need to figure out how to write exactly this
# into the stan code
length <- 50
fit = vector(length = length)
fit_new = vector(length = length)

for(k in 1:length){
  eval <- matrix(nrow = R, ncol = T)
  
  y <- data_set$y
  E <- matrix(nrow = R, ncol = T)
  
  y_new <- matrix(nrow = R, ncol = T)
  E_new <-  matrix(nrow = R, ncol = T)
  
  for(i in 1:R){
    for(j in 1:T){
      # eval = N[i] * p[i]
      eval[i, j] = (rpois(1, exp(list_of_draws_M0_out_sim[1,]$alpha0 + 
                                   (list_of_draws_M0_out_sim[1,]$alpha1 * X[i])
      )) # N[i]
      * ilogit(list_of_draws_M0_out_sim[1,]$beta0 +
                 list_of_draws_M0_out_sim[1,]$beta1 * X[i]
      ) # p[i]
      )
      
      E[i, j] <- (y[i, j] - eval[i, j])^2 / (eval[i, j] + 0.5)
      fit[k] = sum(E)
      
      y_new[i,j] <- rbinom(1, 
                           rpois(1, exp(list_of_draws_M0_out_sim[1,]$alpha0 + 
                                          (list_of_draws_M0_out_sim[1,]$alpha1 * X[i]) 
                           )),
                           ilogit(list_of_draws_M0_out_sim[1,]$beta0 +
                                    list_of_draws_M0_out_sim[1,]$beta1 * X[i]
                           )
      )
      
      E_new[i, j] = (y_new[i, j] - eval[i, j])^2 / (eval[i, j] + 0.5)
      fit_new[k] = sum(E_new)
      
    }
  }
}

df <- cbind(fit[1:length], fit_new[1:length])
plot(df[,1], df[,2], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 500),
     xlim = c(0, 500))
abline(0, 1, lwd = 2, col = "black")

mean(df[,2] > df[,1])
mean(df[,1]) / mean(df[,2])

