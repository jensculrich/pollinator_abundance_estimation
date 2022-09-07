library(rstan)
library(bayesplot)
library(tidyverse)
library(RColorBrewer)
library(lme4)

##########################
# Model M1 ###############
##########################

## N-mix model to estimate the effect of habitat management on abundance
## using a hierarchical, multispecies, metapopulation design

# M1 includes random effects for both species AND site 
# (no random effect for site in M2, becuase too many parameters to estimate
# from the limited size of the real data set)

# M1 considers abundance as arising from a simple Poisson distribution
# as opposed to M3, which uses a Zero inflated poisson distribution,
# or potential M4 which uses a negative-binomial dist to account for 
# overdispersion in abundance and/or detection.

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
                      mu_nu = 0, sd_nu = .25, # site-level effect on abundance
                      mu_xi = 0, sd_xi = .25, # site-level effect on detection
                      mu_eta = 0, sd_eta = 0.5 # species-level effect on abundance slope,
){
  # alpha0 and alpha1: intercept and slope of log-linear regression (abundance)
  # beta0 and beta1: intercept and slope of logistic-linear regression (detection probability)
  
  # number of unique data observation points (observations of unique species at unique site)
  R = n_sites*n_species
  
  # Covariate values: sort for ease of presentation
  X <- sort(rep(c(0, 1), times = (R/2))) # half sites in each category
  
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
  
  ## now add a site-level random effect for abundance
  site_intercept <- rnorm(n_sites, mu_nu, sd_nu)
  
  n_data_per_site <- n_species
  sites <- rep(1:n_sites, each = n_data_per_site)
  N <- length(sites)
  
  site_intercept_data <- as.numeric(vector(length=N))
  
  site_intercept_data <- rep(site_intercept[1:n_sites], 
                             each = n_data_per_site)
  
  test <- cbind(species, sites, species_intercept_data, site_intercept_data)
  
  # Relationship expected abundance – covariate + variation in species
  # determined in part by species-level effect on abundance
  Z <- vector(length=R)
  for(i in 1:length(Z)){
    Z[i] <- exp(alpha0 + alpha1 * X[i] + 
                  species_intercept_data[i] + site_intercept_data[i] +
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
  
  ## now add a site-level random effect for abundance
  site_intercept_detection <- rnorm(n_sites, mu_xi, sd_xi)
  
  site_intercept_data_detection <- as.numeric(vector(length=N))
  
  site_intercept_data_detection <- rep(site_intercept_detection[1:n_sites], 
                                       each = n_data_per_site)
  
  test2 <- cbind(species, sites, species_intercept_data_detection, site_intercept_data_detection)
  
  # ilogit <- function(u) { return(1 / (1 + exp(-u))); }
  p <- vector(length=R)
  for(i in 1:length(p)){
    p[i] <- plogis(beta0 + beta1 * X[i] + 
                     species_intercept_data_detection[i] + site_intercept_data_detection[i])
  }
  p
  
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
  
  df <- as.data.frame(cbind(species, X, data))
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
  plotcolors <- brewer.pal(n = n_species, name = "Dark2")
  
  # take means of species in X = 0 or X = 1 - 
  plot_data <- as.data.frame(
    cbind(X, species, sites, Z, Z_with_sigma, naive.pred))
  
  plot_data_expected_abundance <- plot_data %>%
    group_by(species, X) %>%
    mutate(mean_lambda = mean(Z)) %>%
    distinct(species, X, mean_lambda, .keep_all = FALSE) %>%
    pivot_wider(names_from = species, values_from = mean_lambda)
  
  plot_data_expected_abundance <- as.data.frame(plot_data_expected_abundance)
  
  colnames(plot_data_expected_abundance) <- c(
    "X", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8"
  )
  
  # and now start plotting...
  par(mfrow = c(2, 2))
  
  # expected abundance
  plot(plot_data_expected_abundance$X, plot_data_expected_abundance$s1, 
       main = "Expected abundance", xlab = NA,
       ylab = "lambda", las = 1, type = "l", col = plotcolors[1], lwd = 1.5,
       ylim = c(0, max(plot_data_expected_abundance)*1.5),
       xaxt="n",
       frame.plot = FALSE)
  for(i in 1:n_species){
    lines(plot_data_expected_abundance$X, as.vector(plot_data_expected_abundance[,i]), 
          type = "l", col = plotcolors[i], lwd = 1.5)
  }
  
  axis(1, at = c(0, 1), las = 1, labels = c("mowed", "unmowed"))
  
  # realized abundance (expected plus sigma)
  # needs more work
  plot_data_realized_abundance <- plot_data
  
  data_subset <- rbind(plot_data_realized_abundance[1,],
                       plot_data_realized_abundance[R/2+1,]
  )
  
  plot(data_subset$X, data_subset$Z_with_sigma, 
       main = "Realized abundance", xlab = NA, ylab =
         "N", las = 1, 
       xlim = c(0, 1),
       ylim = c(0, max(plot_data_expected_abundance)*1.5),
       xaxt="n",
       frame.plot = FALSE, col = plotcolors[1], cex = 1.2)
  
  
  
  # Return stuff
  return(list(T = T, X = X,
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
              p = p, y = y,
              sd_eps = sd_eps,
              sd_zeta = sd_zeta,
              sd_nu = sd_nu, # species-level effect on abundance
              sd_xi = sd_xi, # species-level effect on detection))
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
            "sd_nu",
            "sd_xi",
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
       sd_eta = runif(1, 0, 1),
       nu = runif(data_set$n_sites, -1, 1),
       sd_nu = runif(1, 0, 1),
       xi = runif(data_set$n_sites, -1, 1),
       sd_xi = runif(1, 0, 1))
)


# Call STAN model from R 
stan_model_covariate <- "./parks_pollinators_abundance_estimation/M1_Nmix_w_site_random_effects.stan"

## Call Stan from R
M1_out_sim <- stan(stan_model_covariate,
                data = stan_data, 
                init = inits, 
                pars = params,
                chains = n_chains, iter = n_iterations, 
                warmup = n_burnin, thin = n_thin,
                seed = 1,
                open_progress = FALSE,
                cores = n_cores)

print(M1_out_sim, digits = 3)

saveRDS(M1_out_sim, "./parks_pollinators_abundance_estimation/M1_out_sim.rds")
M1_out_sim <- readRDS("./parks_pollinators_abundance_estimation/M1_out_sim.rds")

# plot posterior distribution
p <- mcmc_hist(M1_out_sim, pars = c("alpha0"))
p <- p + labs(x = "alpha0",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$alpha0, linetype = "solid", size = 1)
p

p <- mcmc_hist(M1_out_sim, pars = c("alpha1"))
p <- p + labs(x = "alpha1",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$alpha1, linetype = "solid", size = 1)
p

p <- mcmc_hist(M1_out_sim, pars = c("beta0"))
p <- p + labs(x = "beta0",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$beta0, linetype = "solid", size = 1)
p

p <- mcmc_hist(M1_out_sim, pars = c("beta1"))
p <- p + labs(x = "beta1",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$beta1, linetype = "solid", size = 1)
p

mcmc_dens_overlay(M1_out_sim, pars = c("alpha0", "alpha1", 
                                    "beta0", "beta1"))

traceplot(M1_out_sim, pars = c("alpha0", "alpha1", 
                            "beta0", "beta1",
                            "eta[1]", "eps[1]", "zeta[1]"))

mcmc_pairs(M1_out_sim, pars = c("alpha0", "alpha1", 
                             "beta0", "beta1"),
           off_diag_args = list(size = 1.5))

mcmc_dens_overlay(M1_out_sim, pars = c("totalN[1]", "totalN[1]", 
                                    "totalN[3]", "totalN[4]"))

(p2 <- mcmc_hist(M1_out_sim, pars = c("totalN[1]", "totalN[2]", 
                                   "totalN[3]", "totalN[4]"))
)

# plot posterior distribution
q <- mcmc_hist(M1_out_sim, pars = c("totalN[1]"))
q <- q + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[1], linetype = "solid", size = 1)
q

# plot posterior distribution
s <- mcmc_hist(M1_out_sim, pars = c("totalN[2]"))
s <- s + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[2], linetype = "solid", size = 1)
s

# could set up a loop to draw plots for each species and then cowplot them together

# Plot predicted covariate relationship with abundance
matrix <- as.matrix(M1_out_sim)
str(matrix)
alpha0_mean = mean(matrix[,13])
alpha1_mean = mean(matrix[,14])


list_of_draws_M1_out_sim <- as.data.frame(M1_out_sim)

# Evaluation of fit
plot(list_of_draws_M1_out_sim$fit, list_of_draws_M1_out_sim$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 1500),
     xlim = c(0, 1500))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws_M1_out_sim$fit_new > list_of_draws_M1_out_sim$fit)
mean(list_of_draws_M1_out_sim$fit) / mean(list_of_draws_M1_out_sim$fit_new)

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
            "sd_nu",
            "sd_xi",
            "fit",
            "fit_new"
)


# MCMC settings
n_iterations <- 400
n_thin <- 1
n_burnin <- 200
n_chains <- 3
n_cores <- 3

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
       sd_eta = runif(1, 0, 1),
       nu = runif(data_set$n_sites, -1, 1),
       sd_nu = runif(1, 0, 1),
       xi = runif(data_set$n_sites, -1, 1),
       sd_xi = runif(1, 0, 1))
)


# Call STAN model from R 
stan_model_covariate <- "./individual_based_occupancy_model_M1.stan"

## Call Stan from R
M1_out_real <- stan(stan_model_covariate,
                 data = stan_data, 
                 init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 seed = 1,
                 open_progress = FALSE,
                 cores = n_cores)

print(M1_out_real, digits = 3)

saveRDS(M1_out_real, "./parks_pollinators_abundance_estimation/M1_out_real.rds")
M1_out_real <- readRDS("./parks_pollinators_abundance_estimation/M1_out_real.rds")

list_of_draws <- as.data.frame(M1_out_real)

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
p <- mcmc_hist(M1_out_real, pars = c("alpha0"))
p <- p + labs(x = "alpha0",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$alpha0), linetype = "solid", size = 1)
p

p <- mcmc_hist(M1_out_real, pars = c("alpha1"))
p <- p + labs(x = "alpha1",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$alpha1), linetype = "solid", size = 1)
p

p <- mcmc_hist(M1_out_real, pars = c("beta0"))
p <- p + labs(x = "beta0",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$beta0), linetype = "solid", size = 1)
p

p <- mcmc_hist(M1_out_real, pars = c("beta1"))
p <- p + labs(x = "beta1",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$beta1), linetype = "solid", size = 1)
p

mcmc_dens_overlay(M1_out_real, pars = c("alpha0", "alpha1", 
                                        "beta0", "beta1"))

traceplot(M1_out_real, pars = c("alpha0", "alpha1", 
                                "beta0", "beta1",
                                "eta"))

mcmc_pairs(M1_out_real, pars = c("alpha0", "alpha1", 
                                 "beta0", "beta1"),
           off_diag_args = list(size = 1.5))

mcmc_dens_overlay(M1_out_real, pars = c("totalN[1]", "totalN[2]", 
                                        "totalN[3]", "totalN[4]"))

(p2 <- mcmc_hist(M1_out_real, pars = c("totalN[1]", "totalN[2]", 
                                       "totalN[3]", "totalN[4]"))
)

# plot posterior distribution
q <- mcmc_hist(M1_out_real, pars = c("totalN[1]"))
q <- q + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[1], linetype = "solid", size = 1)
q

# plot posterior distribution
s <- mcmc_hist(M1_out_real, pars = c("totalN[2]"))
s <- s + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[2], linetype = "solid", size = 1)
s

# could set up a loop to draw plots for each species and then cowplot them together
