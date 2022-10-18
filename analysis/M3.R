library(rstan)
library(bayesplot)
library(tidyverse)
library(lme4)

##########################
# Model M3 ###############
##########################

## N-mix model to estimate the effect of habitat management on abundance
## using a hierarchical, multispecies, metapopulation design

# M1 includes random effects for both species AND site 
# (no random effect for site in M2, becuase too many parameters to estimate
# from the limited size of the real data set)

# in contrast to M1/M2, M3 consdiers abundance as arising from a Zero inflated
# poisson distribution, which divides speciesXsite combinations into those that
# are suitable versus unsuitable. Only at suitable sites is abundance drawn from a poisson distr.
# potential M4 will instead use a negative-binomial dist to account for 
# overdispersion in abundance and/or detection.

##########################
### Simulate data ########
##########################

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
            "fit_new",
            "omega"
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
       eps = runif(n_species, -1, 1),
       sd_eps = runif(1, 0, 1),
       sd_zeta = runif(1, 0, 1),
       sd_eta = runif(1, 0, 1))
)


# Call STAN model from R 
stan_model <- "./M3_Nmix_ZIP.stan"

## Call Stan from R
out_real_M3 <- stan(stan_model,
                    data = stan_data, 
                    init = inits, 
                    pars = params,
                    chains = n_chains, iter = n_iterations, 
                    warmup = n_burnin, thin = n_thin,
                    seed = 1,
                    open_progress = FALSE,
                    cores = n_cores)

print(out_real_M3, digits = 3)

saveRDS(out_real_M3, "out_real_M3.rds")
out_real_M3 <- readRDS("out_real_M3.rds")

list_of_draws <- as.data.frame(out_real_M3)

# Evaluation of fit
plot(list_of_draws$fit, list_of_draws$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(0, 2000),
     xlim = c(0, 2000))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws$fit_new > list_of_draws$fit)
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)
