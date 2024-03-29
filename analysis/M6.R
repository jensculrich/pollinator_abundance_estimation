library(rstan)
library(tidyverse)
library(bayesplot)


##########################
# Model M4 ###############
##########################

## N-mix model to estimate the effect of habitat management on abundance
## using a hierarchical, multispecies, metapopulation design

# M4 includes random effects for species 
# (no random effect for site in M4, because too many parameters to estimate
# from the limited size of the real data set)

# M4 will consider abundance as arising from a NEGATIVE BINOMIAL distribution
# rather than a poisson or zero inflated poisson distribution

##########################
### Simulate data ########
##########################

##########################
### Simulate data ########
##########################

# consider the effect of site covariates on lambda and on detection probability

# Define function for generating binom-mix model data
data.fn.2 <- function(n_sites = 30, #  R = sites, 
                      n_species = 20, # number of species to pool
                      T = 3, # T = number of repeat visits (= number of temporal reps)
                      
                      alpha0 = 0, # grand intercept abundance
                      sigma_alpha0_species = 1.5, # community variation in abundance intercept
                      mu_alpha1 = 1, # community mean in abundance response to management
                      sigma_alpha1_species = 0.5, # community variation in abundance intercept
                      
                      beta0 = 0, # grand intercept abundance
                      sigma_beta0_species = 0.25, # community variation in abundance intercept
                      mu_beta1 = -1, # community mean in abundance response to management
                      sigma_beta1_species = 0.25, # community variation in abundance intercept

                      phi = 1 # negbin overdispersion control
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
  
  ## add a species-level random effect for abundance
  # both on the intercept (some species more abundant than others) centered on 0
  species_intercept <- rnorm(n_species, 0, sigma_alpha0_species)
  (mean(species_intercept)) # should be near 0
  # and on the slope (change in abundance of different species responds differently to management)
  # centered on community mean
  species_slope <- rnorm(n_species, mu_alpha1, sigma_alpha1_species)
  (mean(species_slope)) # should be near mu_alpha1
  
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
  # species effect on management parameter (alpha1) should be akin to something like this: 
  # psi.meanmaxt[species]   * meanmaxt[site,era] 
  Z <- vector(length=R)
  for(i in 1:length(Z)){ # rather than alpha1*X + species slope * X I've now combined them
    Z[i] <- exp(alpha0 +
                  species_intercept_data[i] +
                  species_slope_data[i] * X[i] 
                  )
  }
  
  # Add noise: draw N from neg_bin(mu, phi)
  Z_with_sigma <- Z
  for(i in 1:length(Z)){
    Z_with_sigma[i] <- rnbinom(n = 1, mu = Z[i], size = phi)
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
  ## add a species-level random effect for abundance
  # both on the intercept (some species more abundant than others)
  species_intercept_detection <- rnorm(n_species, 0, sigma_beta0_species)
  (mean(species_intercept_detection)) # should be near 0
  # and on the slope (change in abundance of different species responds differently to management)
  species_slope_detection <- rnorm(n_species, mu_beta1, sigma_beta1_species)
  (mean(species_slope_detection)) # should be near mu_alpha1
  
  species_intercept_data_detection <- as.numeric(vector(length=N))
  species_intercept_data_detection <- rep(species_intercept_detection[1:n_species], 
                                          n_data_per_species)
  
  species_slope_data_detection <- as.numeric(vector(length=N))
  species_slope_data_detection <- rep(species_slope_detection[1:n_species], 
                            n_data_per_species)
  
  test2 <- cbind(species, species_intercept_data_detection, species_slope_data_detection)
  
  # ilogit <- function(u) { return(1 / (1 + exp(-u))); }
  p <- vector(length=R)
  for(i in 1:length(p)){
    p[i] <- plogis(beta0 + 
                     species_intercept_data_detection[i] +
                     species_slope_data_detection[i] * X[i])
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
  # Plot features of the simulated system
  # need to work on this
  
  # Return stuff
  return(list(
              # study design properties
              T = T, 
              X = X,
              n_species = n_species, # number of species to pool
              species = species,
              n_sites = n_sites,
              sites = sites,
              
              # parameters
              # ecological
              alpha0 = alpha0, 
              sigma_alpha0_species = sigma_alpha0_species,
              mu_alpha1 = mu_alpha1,
              sigma_alpha1_species = sigma_alpha1_species,
              phi = phi,
              # observation
              beta0 = beta0, 
              sigma_beta0_species = sigma_beta0_species,
              mu_beta1 = mu_beta1,
              sigma_beta1_species = sigma_beta1_species,
              
              # simulated data outcomes 
              Z = Z, 
              Z_with_sigma = Z_with_sigma, 
              totalN = totalN,
              p = p, 
              y = y # the counts themselves that will be fed to the model
  ))
}

data_set <- data.fn.2()
str(data_set)

## Fit the model to the simulated data
# Bundle data
# need to define a maximum population size K
View(as.data.frame(data_set$Z_with_sigma))
max(data_set$Z_with_sigma)
K <- as.integer(max(data_set$Z_with_sigma) * 1.25)

totalN <- data_set$totalN

y <- data_set$y
X <- data_set$X
R <- data_set$n_sites*data_set$n_species
T <- data_set$T
sites <- data_set$sites
n_sites <- data_set$n_sites
n_species <- data_set$n_species
species <- data_set$species

stan_data <- c("R", 
               "sites",
               "n_sites",
               "species",
               "n_species", 
               "T", "y", "X", "K")

# Parameters monitored
params <- c("totalN",
            "alpha0",
            "alpha0_species",
            "sigma_alpha0_species",
            "mu_alpha1",
            "sigma_alpha1_species",
            "alpha1",
            "beta0",
            "beta0_species",
            "sigma_beta0_species",
            "mu_beta1",
            "sigma_beta1_species",
            "beta1",
            "phi",
            "fit",
            "fit_new"
)


# MCMC settings
n_iterations <- 800
n_thin <- 1
n_burnin <- 400
n_chains <- 3
n_cores <- 3

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  list(alpha0 = runif(1, -1, 1),
       sigma_alpha0_species = runif(1, 0, 1),
       mu_alpha1 = runif(1, -3, 3),
       sigma_alpha1_species = runif(1, 0, 1),
       beta0 = runif(1, -1, 1),
       sigma_beta0_species = runif(1, 0, 1),
       mu_beta1 = runif(1, -3, 3),
       sigma_beta1_species = runif(1, 0, 1),
       phi = runif(1, 0, 4)
  )
)

# Call STAN model from R 
stan_model_covariate <- "./parks_pollinators_abundance_estimation/M5_Nmix_OD_2.stan"

## Call Stan from R
M5_out_sim <- stan(stan_model_covariate,
                   data = stan_data, 
                   init = inits, 
                   pars = params,
                   chains = n_chains, iter = n_iterations, 
                   warmup = n_burnin, thin = n_thin,
                   seed = 1,
                   open_progress = FALSE,
                   cores = n_cores)

print(M5_out_sim, digits = 3)

saveRDS(M5_out_sim, "./parks_pollinators_abundance_estimation/M5_out_sim.rds")
M5_out_sim <- readRDS("./parks_pollinators_abundance_estimation/M5_out_sim.rds")

# plot posterior distribution
p <- mcmc_hist(M5_out_sim, pars = c("alpha0"))
p <- p + labs(x = "alpha0",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$alpha0, linetype = "solid", size = 1)
p

p <- mcmc_hist(M5_out_sim, pars = c("mu_alpha1"))
p <- p + labs(x = "alpha1",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$mu_alpha1, linetype = "solid", size = 1)
p

p <- mcmc_hist(M5_out_sim, pars = c("beta0"))
p <- p + labs(x = "beta0",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$beta0, linetype = "solid", size = 1)
p

p <- mcmc_hist(M5_out_sim, pars = c("mu_beta1"))
p <- p + labs(x = "beta1",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$mu_beta1, linetype = "solid", size = 1)
p

mcmc_dens_overlay(M5_out_sim, pars = c("alpha0", "mu_alpha1", 
                                       "beta0", "mu_beta1"))

traceplot(M5_out_sim, pars = c("alpha0", "mu_alpha1", 
                               "beta0", "mu_beta1"))

mcmc_pairs(M5_out_sim, pars = c("alpha0", "mu_alpha1", 
                                "beta0", "mu_beta1"),
           off_diag_args = list(size = 1.5))

(p2 <- mcmc_hist(M4_out_sim, pars = c("totalN[1]", "totalN[2]", 
                                      "totalN[3]", "totalN[4]"))
)

# plot posterior distribution
q <- mcmc_hist(M4_out_sim, pars = c("totalN[1]"))
q <- q + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[1], linetype = "solid", size = 1)
q

# plot posterior distribution
s <- mcmc_hist(M4_out_sim, pars = c("totalN[2]"))
s <- s + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = data_set$totalN[2], linetype = "solid", size = 1)
s

# could set up a loop to draw plots for each species and then cowplot them together

# Evaluation of fit
list_of_draws_M5_out_sim <- as.data.frame(M5_out_sim)

par(mfrow = c(1, 1))

plot(list_of_draws_M5_out_sim$fit, list_of_draws_M5_out_sim$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(3000, 9000),
     xlim = c(3000, 9000))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws_M5_out_sim$fit_new > list_of_draws_M5_out_sim$fit)
mean(list_of_draws_M5_out_sim$fit) / mean(list_of_draws_M5_out_sim$fit_new)


####################################
# analysis of real data ############
####################################

df <- read.csv("./parks_pollinators_abundance_estimation/raw_data_2022.csv")
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

y <- pull(df[,4]) # ONLY THE FIRST COUNT COLUMN as opposed to [,4:6]
K <- as.integer(max(y)*2.5) # need to define a maximum population size K
X <- as.vector(df$X)
R <- length(y)
T <- ncol(y)

names <- rbind("Bombus mixtus", "Bombus flavifrons", "Agapostemon texanus",
          "Anthidium oblongatum", "Halictus rubicundus",
          "Melissodes microsticus", "Megachile montivaga")

species_names <- cbind(cbind(1:7), names)

## Fit the model to the real data
# Bundle data
stan_data <- c("R", 
               "sites",
               "n_sites",
               "species",
               "n_species",
               "y", "X", "K")

# call model M2, no random effect for site, distes fully interchangeable
# Parameters monitored
# Parameters monitored
params <- c("totalN",
            "alpha0",
            "alpha0_species",
            "sigma_alpha0_species",
            "mu_alpha1",
            "sigma_alpha1_species",
            "alpha1",
            "phi",
            "fit",
            "fit_new"
)


# MCMC settings
n_iterations <- 6000
n_thin <- 1
n_burnin <- 2000
n_chains <- 4
n_cores <- 4
adapt_delta <- 0.99

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  list(alpha0 = runif(1, -1, 1),
       sigma_alpha0_species = runif(1, 0, 1),
       mu_alpha1 = runif(1, -3, 3),
       sigma_alpha1_species = runif(1, 0, 1),
       phi = runif(1, 0, 4)
  )
)

# Call STAN model from R 
stan_model <- "./parks_pollinators_abundance_estimation/M6_perfect_detection.stan"

## Call Stan from R
out_real_M6 <- stan(stan_model,
                    data = stan_data, 
                    init = inits, 
                    pars = params,
                    chains = n_chains, iter = n_iterations, 
                    warmup = n_burnin, thin = n_thin,
                    seed = 1,
                    open_progress = FALSE,
                    cores = n_cores,
                    control=list(adapt_delta = adapt_delta)
)

print(out_real_M6, digits = 3)

saveRDS(out_real_M6, "out_real_M6.rds")
out_real_M6 <- readRDS("out_real_M6.rds")
# perhaps revisit: adapt delta, priors, inits, iterations seems ok.

list_of_draws <- as.data.frame(out_real_M6)

# Evaluation of fit
plot(list_of_draws$fit, list_of_draws$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 10000),
     xlim = c(0, 10000))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws$fit_new > list_of_draws$fit)
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)
# decreased model fit! versus model M5

# plot posterior distribution
p <- mcmc_hist(out_real_M6, pars = c("alpha0"))
p <- p + labs(x = "alpha0",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$alpha0), linetype = "solid", size = 1)
p

p <- mcmc_hist(out_real_M6, pars = c("mu_alpha1"))
p <- p + labs(x = "Community mean effect of no-mow habitat on abundance",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1), linetype = "solid", size = 1) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1) + 
               sd(list_of_draws$mu_alpha1), linetype = "dashed", size = .75) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1) - 
               sd(list_of_draws$mu_alpha1), linetype = "dashed", size = .75) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1) + 
               (2 * sd(list_of_draws$mu_alpha1)), color = "red", linetype = "dashed", size = .75) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1) - 
               (2 * sd(list_of_draws$mu_alpha1)), color = "red", linetype = "dashed", size = .75) +
  theme(axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16))    
p

# would be good to add a g linear model plot with mean and species specific response lines

mcmc_dens_overlay(out_real_M6, pars = c("alpha0", "mu_alpha1"))

traceplot(out_real_M6, pars = c("alpha0", "mu_alpha1"))

mcmc_pairs(out_real_M6, pars = c("alpha0", "mu_alpha1"),
           off_diag_args = list(size = 1.5))



#################
# Now let's try with max_y per site (most individuals detected on any visit to site)
#############
y <- apply(as.matrix(df[,4:6]) , 1, max) # max value from all three visits
K <- as.integer(max(y)*2.5) # need to define a maximum population size K
X <- as.vector(df$X)
R <- length(y)
T <- ncol(y)

names <- rbind("Bombus mixtus", "Bombus flavifrons", "Agapostemon texanus",
               "Anthidium oblongatum", "Halictus rubicundus",
               "Melissodes microsticus", "Megachile montivaga")

species_names <- cbind(cbind(1:7), names)

## Fit the model to the real data
# Bundle data
stan_data <- c("R", 
               "sites",
               "n_sites",
               "species",
               "n_species",
               "y", "X", "K")

# call model M2, no random effect for site, distes fully interchangeable
# Parameters monitored
# Parameters monitored
params <- c("totalN",
            "alpha0",
            "alpha0_species",
            "sigma_alpha0_species",
            "mu_alpha1",
            "sigma_alpha1_species",
            "alpha1",
            "phi",
            "fit",
            "fit_new"
)


# MCMC settings
n_iterations <- 6000
n_thin <- 1
n_burnin <- 2000
n_chains <- 4
n_cores <- 4
adapt_delta <- 0.99

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  list(alpha0 = runif(1, -1, 1),
       sigma_alpha0_species = runif(1, 0, 1),
       mu_alpha1 = runif(1, -3, 3),
       sigma_alpha1_species = runif(1, 0, 1),
       phi = runif(1, 0, 4)
  )
)

# Call STAN model from R 
stan_model <- "./parks_pollinators_abundance_estimation/M6_perfect_detection.stan"

## Call Stan from R
out_real_M6_ymax <- stan(stan_model,
                    data = stan_data, 
                    init = inits, 
                    pars = params,
                    chains = n_chains, iter = n_iterations, 
                    warmup = n_burnin, thin = n_thin,
                    seed = 1,
                    open_progress = FALSE,
                    cores = n_cores,
                    control=list(adapt_delta = adapt_delta)
)

print(out_real_M6_ymax, digits = 3)

saveRDS(out_real_M6_ymax, "out_real_M6_ymax.rds")
out_real_M6_ymax <- readRDS("out_real_M6_ymax.rds")
# perhaps revisit: adapt delta, priors, inits, iterations seems ok.

list_of_draws <- as.data.frame(out_real_M6_ymax)

# Evaluation of fit
plot(list_of_draws$fit, list_of_draws$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 10000),
     xlim = c(0, 10000))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws$fit_new > list_of_draws$fit)
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)
# decreased model fit! versus model M5

# plot posterior distribution
p <- mcmc_hist(out_real_M6_ymax, pars = c("alpha0"))
p <- p + labs(x = "alpha0",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$alpha0), linetype = "solid", size = 1)
p

p <- mcmc_hist(out_real_M6_ymax, pars = c("mu_alpha1"))
p <- p + labs(x = "Community mean effect of no-mow habitat on abundance",
              y = "Frequency in 3,200 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1), linetype = "solid", size = 1) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1) + 
               sd(list_of_draws$mu_alpha1), linetype = "dashed", size = .75) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1) - 
               sd(list_of_draws$mu_alpha1), linetype = "dashed", size = .75) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1) + 
               (2 * sd(list_of_draws$mu_alpha1)), color = "red", linetype = "dashed", size = .75) +
  geom_vline(xintercept = median(list_of_draws$mu_alpha1) - 
               (2 * sd(list_of_draws$mu_alpha1)), color = "red", linetype = "dashed", size = .75) +
  theme(axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16))    
p

# would be good to add a g linear model plot with mean and species specific response lines

mcmc_dens_overlay(out_real_M6_ymax, pars = c("alpha0", "mu_alpha1"))

traceplot(out_real_M6_ymax, pars = c("alpha0", "mu_alpha1"))

mcmc_pairs(out_real_M6_ymax, pars = c("alpha0", "mu_alpha1"),
           off_diag_args = list(size = 1.5))
