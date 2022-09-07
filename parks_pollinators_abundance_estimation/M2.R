library(rstan)
library(bayesplot)
library(tidyverse)
library(RColorBrewer)
library(lme4)

##########################
# Model M2 ###############
##########################

## N-mix model to estimate the effect of habitat management on abundance
## using a hierarchical, multispecies, metapopulation design

# M2 includes random effects for species 
# (no random effect for site in M2, because too many parameters to estimate
# from the limited size of the real data set)

# As in M1, M2 considers abundance as arising from a simple Poisson distribution
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
data.fn.2 <- function(n_sites = 30, #  R = sites, 
                      T = 3, # T = number of repeat visits (= number of temporal reps)
                      alpha0 = 0, # grand intercept abundance
                      alpha1 = 2, # site-level effect on abundance
                      beta0 = 0, # grand intercept detection probability
                      beta1 = -2, # site-level effect on detection
                      n_species = 12, # number of species to pool
                      mu_eps = 0, sd_eps = .5, # species-level effect on abundance
                      mu_zeta = 0, sd_zeta = .5, # species-level effect on detection
                      mu_eta = 0, sd_eta = .5 # species-level effect on abundance slope
){
  # alpha0 and alpha1: intercept and slope of log-linear regression (abundance)
  # beta0 and beta1: intercept and slope of logistic-linear regression (detection probability)
  
  # number of unique data observation points (observations of unique species at unique site)
  R = n_sites*n_species
  
  # Covariate values: sort for ease of presentation
  X <- sort(rep(c(0, 1), times = (R/2))) # half sites in each category
  
  n_data_per_site <- n_species
  sites <- rep(1:n_sites, each = n_data_per_site)
  
  
  ## Ecological process
  
  ## now add a species-level random effect for abundance
  species_intercept <- rnorm(n_species, mu_eps, sd_eps)
  species_slope <- rnorm(n_species, mu_eta, sd_eta)
  
  n_data_per_species <- n_sites
  species <- rep(1:n_species, n_data_per_species)
  N <- length(species)
  
  species_intercept_data <- as.numeric(vector(length=N))
  species_intercept_data <- rep(species_intercept[1:n_species], 
                                n_data_per_species)
  
  species_slope_data <- as.numeric(vector(length=N))
  species_slope_data <- rep(species_slope[1:n_species], 
                            n_data_per_species)
  
  test <- cbind(species, species_intercept_data, species_slope_data)
  
  # Relationship expected abundance – covariate + variation in species
  # determined in part by species-level effect on abundance
  Z <- vector(length=R)
  for(i in 1:length(Z)){
    Z[i] <- exp(alpha0 + alpha1 * X[i] + 
                  species_intercept_data[i] + 
                  species_slope_data[i] * X[i])
  }
  
  # Add Poisson noise: draw N from Poisson(lambda)
  Z_with_sigma <- Z
  for(i in 1:length(Z)){
    Z_with_sigma[i] <- rpois(n = 1, lambda = Z[i])
  }
  
  Z_with_sigma # Distribution of abundances across sites
  sum(Z_with_sigma > 0) / (R) # Empirical occupancy 
  # (proportion of species X site combos w/ 1 or more individuals)
  totalN <- vector(length=n_species)
  totalN <- rowsum(Z_with_sigma, rep(1:n_species, times = n_sites))
  totalN
  
  ## Observation process
  # Relationship detection prob – covariate
  ## now add a species-level random effect for detection
  sd_zeta <- sd_zeta
  
  species_intercept_detection <- rnorm(n_species, mu_zeta, sd_zeta)
  
  species_intercept_data_detection <- as.numeric(vector(length=N))
  
  species_intercept_data_detection <- rep(species_intercept_detection[1:n_species], 
                                          n_data_per_species)
  test2 <- cbind(species, species_intercept_data_detection)
  
  # ilogit <- function(u) { return(1 / (1 + exp(-u))); }
  p <- vector(length=R)
  for(i in 1:length(p)){
    p[i] <- plogis(beta0 + beta1 * X[i] + 
                     species_intercept_data_detection[i])
  }
  p
  
  par(mfrow = c(2, 1))
  hist(p[1:(R/2)], xlim = c(0,1))
  hist(p[(R/2):(R)], xlim = c(0,1))
  
  # Make a 'census' (i.e., go out and count things)
  y <- matrix(NA, nrow = length(p), ncol = T) # Array for counts
  for(i in 1:nrow(y)){
    for(j in 1:T){
      y[i,j] <- rbinom(n = 1, size = Z_with_sigma[i], prob = p[i])
    }
  }
  
  # Naïve regression
  data <- vector(length = nrow(y))
  data <- apply(y, 1, max) # highest count for each row
  
  df <- as.data.frame(cbind(species, sites, X, data))
  df$species <- as.factor(df$species)
  
  naive.pred <- exp(predict(glmer(data = df, data ~ (X|species),
                                  family = poisson)))
  naive_fit <- glmer(data = df, data ~ (X|species),
                     family = poisson)
  coef(naive_fit)$species[1]
  ranef(naive_fit)$species # variation in species intercepts (around mean 0)
  mean((ranef(naive_fit)$species)[,1]) # should be 0
  mean((coef(naive_fit)$species)[,1]) # should be ~mu_u
  
  # Plot features of the simulated system
  # need to work on this
  
  # Return stuff
  return(list(T = T, 
              X = X,
              n_species = n_species, # number of species to pool
              species = species,
              n_sites = n_sites,
              sites = sites,
              alpha0 = alpha0, 
              alpha1 = alpha1,
              beta0 = beta0, 
              beta1 = beta1,
              Z = Z, 
              Z_with_sigma = Z_with_sigma, 
              totalN = totalN,
              p = p, 
              y = y,
              sd_eps = sd_eps,
              sd_zeta = sd_zeta,
              sd_eta = sd_eta,
              eps = species_intercept,
              zeta = species_intercept_detection,
              eta = species_slope
  ))
}

data_set <- data.fn.2()
str(data_set)

## Fit the model to the simulated data
# Bundle data
# need to define a maximum population size K
max(data_set$Z_with_sigma)
K <- as.integer(max(data_set$Z_with_sigma) * 1.25)

y <- data_set$y
X <- data_set$X
R <- data_set$n_sites*data_set$n_species
T <- data_set$T
sites <- data_set$sites
n_sites <- data_set$n_sites
n_species <- data_set$n_species
species <- data_set$species

eps <- data_set$eps
zeta <- data_set$zeta
eta <- data_set$eta

stan_data <- c("R", 
               "sites",
               "n_sites",
               "species",
               "n_species", 
               "T", "y", "X", "K")

# Parameters monitored
params <- c("totalN",
            "alpha0",
            "alpha1",
            "beta0",
            "beta1",
            "eps",
            "zeta",
            "eta",
            "sd_eps",
            "sd_eta",
            "sd_zeta",
            "fit",
            "fit_new"
)


# MCMC settings
n_iterations <- 1000
n_thin <- 1
n_burnin <- 500
n_chains <- 4
n_cores <- 2

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  list(alpha0 = runif(1, -1, 1),
       alpha1 = runif(1, -3, 3),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, -3, 3),
       eps = runif(data_set$n_species, -1, 1),
       sd_eps = runif(1, 0, 1),
       zeta = runif(data_set$n_species, -1, 1),
       sd_zeta = runif(1, 0, 1),
       eta = runif(data_set$n_species, -1, 1),
       sd_eta = runif(1, 0, 1)
  )
)

# Call STAN model from R 
stan_model_covariate <- "./parks_pollinators_abundance_estimation/M2_Nmix.stan"

## Call Stan from R
M2_out_sim <- stan(stan_model_covariate,
                   data = stan_data, 
                   init = inits, 
                   pars = params,
                   chains = n_chains, iter = n_iterations, 
                   warmup = n_burnin, thin = n_thin,
                   seed = 1,
                   open_progress = FALSE,
                   cores = n_cores)

print(M2_out_sim, digits = 3)

saveRDS(M2_out_sim, "./parks_pollinators_abundance_estimation/M2_out_sim.rds")
M2_out_sim <- readRDS("./parks_pollinators_abundance_estimation/M2_out_sim.rds")

# plot posterior distribution
p <- mcmc_hist(M2_out_sim, pars = c("alpha0"))
p <- p + labs(x = "alpha0",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$alpha0, linetype = "solid", size = 1)
p

p <- mcmc_hist(M2_out_sim, pars = c("alpha1"))
p <- p + labs(x = "alpha1",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$alpha1, linetype = "solid", size = 1)
p

p <- mcmc_hist(M2_out_sim, pars = c("beta0"))
p <- p + labs(x = "beta0",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$beta0, linetype = "solid", size = 1)
p

p <- mcmc_hist(M2_out_sim, pars = c("beta1"))
p <- p + labs(x = "beta1",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$beta1, linetype = "solid", size = 1)
p

mcmc_dens_overlay(M2_out_sim, pars = c("alpha0", "alpha1", 
                                       "beta0", "beta1"))

traceplot(M2_out_sim, pars = c("alpha0", "alpha1", 
                               "beta0", "beta1",
                               "eta[1]", "eps[1]", "zeta[1]"))

mcmc_pairs(M2_out_sim, pars = c("alpha0", "alpha1", 
                                "beta0", "beta1"),
           off_diag_args = list(size = 1.5))

mcmc_dens_overlay(M2_out_sim, pars = c("totalN[1]", "totalN[1]", 
                                       "totalN[3]", "totalN[4]"))

(p2 <- mcmc_hist(M2_out_sim, pars = c("totalN[1]", "totalN[2]", 
                                      "totalN[3]", "totalN[4]"))
)

# plot posterior distribution
q <- mcmc_hist(M2_out_sim, pars = c("totalN[1]"))
q <- q + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[1], linetype = "solid", size = 1)
q

# plot posterior distribution
s <- mcmc_hist(M2_out_sim, pars = c("totalN[2]"))
s <- s + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[2], linetype = "solid", size = 1)
s

# could set up a loop to draw plots for each species and then cowplot them together

# Evaluation of fit
list_of_draws_M2_out_sim <- as.data.frame(M2_out_sim)

par(mfrow = c(1, 1))

plot(list_of_draws_M2_out_sim$fit, list_of_draws_M2_out_sim$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 1200),
     xlim = c(0, 1200))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws_M2_out_sim$fit_new > list_of_draws_M2_out_sim$fit)
mean(list_of_draws_M2_out_sim$fit) / mean(list_of_draws_M2_out_sim$fit_new)

## let's see if we can do a posterior predictive check based on the parameter estimates from 
## the output. Maybe this will help to see if something is wrong with the posterior check 
## in the generated quantities block versus the model itself..


# inv logit transformation functio
ilogit <- function(x){
  exp(x)/(1+exp(x))
}

fit = vector(length = nrow(list_of_draws_M2_out_sim))
fit_new = vector(length = nrow(list_of_draws_M2_out_sim))

# Expected values for counts
for(k in 1:nrow(list_of_draws_M2_out_sim)){
  eval <- matrix(nrow = R, ncol = T)
  
  y <- data_set$y
  E <- matrix(nrow = R, ncol = T)
  
  y_new <- matrix(nrow = R, ncol = T)
  E_new <-  matrix(nrow = R, ncol = T)
  
  for(i in 1:R){
    for(j in 1:T){
      eval[i, j] = (rpois(1, exp(list_of_draws_M2_out_sim[1,]$alpha0 + 
                      (list_of_draws_M2_out_sim[k,]$alpha1 * X[i]) + 
                      list_of_draws_M2_out_sim[k,species[i]+16] +
                      (list_of_draws_M2_out_sim[k,species[i]+40] * X[i]))) # N[i]
      * ilogit(list_of_draws_M2_out_sim[k,]$beta0 +
                 list_of_draws_M2_out_sim[k,]$beta1 * X[i] +
                 list_of_draws_M2_out_sim[k,species[i]+28]) # p[i]
      )
      
      E[i, j] <- (y[i, j] - eval[i, j])^2 / (eval[i, j] + 0.5)
      fit[k] = sum(E)
      
      y_new[i,j] <- rbinom(1, 
                           rpois(1, exp(list_of_draws_M2_out_sim[k,]$alpha0 + 
                                          (list_of_draws_M2_out_sim[k,]$alpha1 * X[i]) + 
                                          list_of_draws_M2_out_sim[k,species[i]+16] +
                                          (list_of_draws_M2_out_sim[k,species[i]+40] * X[i]))),
                                 ilogit(list_of_draws_M2_out_sim[k,]$beta0 +
                                          list_of_draws_M2_out_sim[k,]$beta1 * X[i] +
                                          list_of_draws_M2_out_sim[k,species[i]+28]))
      
      E_new[i, j] = (y_new[i, j] - eval[i, j])^2 / (eval[i, j] + 0.5)
      fit_new[k] = sum(E_new)
      
    }
  }
}

df <- cbind(fit[1:100], fit_new[1:100])
plot(df[,1], df[,2], main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 5000),
     xlim = c(0, 5000))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws_M2_out_sim$fit_new > list_of_draws_M2_out_sim$fit)
mean(list_of_draws_M2_out_sim$fit) / mean(list_of_draws_M2_out_sim$fit_new)

####################################
# analysis of real data ############
####################################

df <- read.csv("./raw_data_2022.csv")
str(df)

## Clean and prep data for model fitting
# select needed columns
df <- select(df, sites, species, X, count)

n_species <- nrow(distinct(df, species)) # number of species
n_sites <- nrow(distinct(df, sites)) # number of sites

# long  to wide
visit <- rep(1:3, times = (n_species*n_sites) )
df <- cbind(df, visit)
df <- pivot_wider(df, names_from = visit, values_from = count)

species <- as.vector(df$species) # vector of species
sites <- as.vector(df$sites) # vector of sites

y <- as.matrix(df[,4:6]) # count columns
K <- as.integer(max(y)*2) # need to define a maximum population size K
X <- as.vector(df$X)
R <- nrow(y)
T <- ncol(y)

## Fit the model to the real data
# Bundle data
stan_data <- c("R", 
               "sites",
               "n_sites",
               "species",
               "n_species", 
               "T", "y", "X", "K")

# call model M2, no random effect for site, distes fully interchangeable
# Parameters monitored
params <- c("totalN",
            "alpha0",
            "alpha1",
            "beta0",
            "beta1",
            "eps",
            "zeta",
            "eta",
            "sd_eps",
            "sd_eta",
            "sd_zeta",
            "fit",
            "fit_new"
)

# MCMC settings
n_iterations <- 1600
n_thin <- 1
n_burnin <- 800
n_chains <- 4
n_cores <- 2

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  list(alpha0 = runif(1, -1, 1),
       alpha1 = runif(1, -3, 3),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, -3, 3),
       eps = runif(n_species, -1, 1),
       sd_eps = runif(1, 0, 1),
       sd_zeta = runif(1, 0, 1),
       sd_eta = runif(1, 0, 1))
)


# Call STAN model from R 
stan_model <- "./M2_Nmix.stan"

## Call Stan from R
out_real_M2 <- stan(stan_model,
                    data = stan_data, 
                    init = inits, 
                    pars = params,
                    chains = n_chains, iter = n_iterations, 
                    warmup = n_burnin, thin = n_thin,
                    seed = 1,
                    open_progress = FALSE,
                    cores = n_cores)

print(out_real_M2, digits = 3)

saveRDS(out_real_M2, "out_real_M2.rds")
out_real_M2 <- readRDS("out_real_M2.rds")

list_of_draws <- as.data.frame(out_real_M2)

# Evaluation of fit
plot(list_of_draws$fit, list_of_draws$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 2000),
     xlim = c(0, 2000))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws$fit_new > list_of_draws$fit)
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)

## model M2 doesn't seem like a great fit!!

# plot posterior distribution
p <- mcmc_hist(out_real_M2, pars = c("alpha0"))
p <- p + labs(x = "alpha0",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$alpha0), linetype = "solid", size = 1)
p

p <- mcmc_hist(out_real_M2, pars = c("alpha1"))
p <- p + labs(x = "alpha1",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$alpha1), linetype = "solid", size = 1)
p

p <- mcmc_hist(out_real_M2, pars = c("beta0"))
p <- p + labs(x = "beta0",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$beta0), linetype = "solid", size = 1)
p

p <- mcmc_hist(out_real_M2, pars = c("beta1"))
p <- p + labs(x = "beta1",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$beta1), linetype = "solid", size = 1)
p

mcmc_dens_overlay(out_real_M2, pars = c("alpha0", "alpha1", 
                                        "beta0", "beta1"))

traceplot(out_real_M2, pars = c("alpha0", "alpha1", 
                                "beta0", "beta1",
                                "eta"))

mcmc_pairs(out_real_M2, pars = c("alpha0", "alpha1", 
                                 "beta0", "beta1"),
           off_diag_args = list(size = 1.5))

mcmc_dens_overlay(out_real_M2, pars = c("totalN[1]", "totalN[2]", 
                                        "totalN[3]", "totalN[4]"))

(p2 <- mcmc_hist(out_real_M2, pars = c("totalN[1]", "totalN[2]", 
                                       "totalN[3]", "totalN[4]"))
)

# plot posterior distribution
q <- mcmc_hist(out_real, pars = c("totalN[1]"))
q <- q + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[1], linetype = "solid", size = 1)
q

# plot posterior distribution
s <- mcmc_hist(out_real, pars = c("totalN[2]"))
s <- s + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[2], linetype = "solid", size = 1)
s

# could set up a loop to draw plots for each species and then cowplot them together


