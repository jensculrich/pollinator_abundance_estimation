library(rstan)
library(bayesplot)
library(tidyverse)

## Mark-recapture populations size estimates
# JCU May, 2022

# This R script first calls a Stan program defining a Lincoln-Peterson mark-recapture model, 
# using data from captured and recaptured individuals to estimate population size,
# as outlined in the STAN manual "Chapter 7.3 Mark-Recapture Models". Here, 
# population size is a probabilistic maximum likelihood estimate, and population size
# is treated as continuous.

# The second model called below (an M0 type model, Otis et al., 1978) 
# will be used on simulated data and the real data to estimate
# population size of bees in urban parks using mark-recapture detection history data.

### Model 1: Lincoln-Peterson model
# The L-P mathematical model is defined as: Nhat = (M*C)/R
# This population estimate would arise from a probabilistic model in 
# which the number of animals recaptured is binomially distributed:
# R ~ binom(C, M/N)

# define some reasonable 
M <- 20 # Initially captured, marked individuals 
C <- 25 # Individuals captured on a revisit
R <- 18 # Number of individuals captured on the revisit that were marked

my_dat <- c("M", "C", "R")

stan_model <- "./mark_recapture_simple_LincolnPeterson.stan"

# keep iter and chains small just for purposes of making sure the model runs
n_iter <- 4000
n_chains <- 4
n_cores <- 4
sim_fit <- stan(stan_model, data = my_dat, 
                    iter = n_iter, warmup = n_iter*0.5, 
                    chains = n_chains, cores = n_cores, 
                    control=list(adapt_delta=0.9))

options(width="120")
print(sim_fit, c("N"))



### Model 2: M0 model for population size estimation (Otis et al. 1978)
# estimating population size using detection histories following examples from 
# Bayesian population analysis using WinBUGS: a hierarchical perspective
# Kery & Schaub 2011

## Define function to simulate data under M0 data.
# M0, model where detection is constant across individuals and through time
# N = population size
# p = detection probability
# and T = number of surveys

# simulate some detection histories, given N, p and T:
fn <- function(N = 100, p = 0.5, T = 3){
  yfull <- yobs <- array(NA, dim = c(N, T))
  
  for (j in 1:T){ 
    yfull[,j] <- rbinom(n = N, size = 1, prob = p) 
  }
  
  ever.detected <- apply(yfull, 1, max) 
  C <- sum(ever.detected) 
  
  yobs <- yfull[ever.detected == 1,] 
  cat(C, "out of", N, "individuals present were detected.\n") 
  return(list(N = N, p = p, C = C, T = T, yfull = yfull, yobs = yobs)) 
} 

# identify number of individuals ever detected (out of maximum N):
data = fn() # note, not all individuals are detected and thus detection is imperfect!

str(data)
# here, yfull is the full capture history matrix of all N individuals, including those
# that were never detected and which have an all-zero capture history. 
data$yfull
# the observed data are called yobs, and they have C rows only, where C is less than or 
# equal to N (i.e., N includes all observed individuals PLUS any that we never detected)
data$yobs

# We will apply a model to use these data to make an inference, based on the rules of 
# probability, of how great N might likely be.
# When estimating population size using Bayesian MCMC approaches, a technical
# challenge is that the dimension of the parameter N may change at every iteration of
# the MCMC algorithm. Kery & Schaub use a data augmentation approach to circumvent this issue.

# Data augmentation works by adding a large number of potential unobserved individuals.
# The augmented data set has dimensions of M x T, where M >> N. 
# To this augmented data set we fit a zero inflated model.
# we add to the model a binary indicator variable, z, that represents whether a row
# in the augmented data set is a 'real' individual or not. These indicators are given
# a bernoulli prior distribution, and the parameter of that distribution, Omega or,
# the inclusion probability, is estimable from the data.
# Augmentation translates the problem of estimating N into the problem of estimating
# Omega, since the expectation of N is Omega*M

# essentially we are now fitting an occupancy model where we will sum up the latent
# indicators, z, which represent whether potential individuals are 
# 'occupied' by real individuals despite them potentially not having ever
# been detected.

# Augment the observed dataset by 150 potential individuals:
nz <- 150  # choose an arbitrarily large number of new individuals
# can later check that it was large enough by confirming that the posterior distr
# of N is not right truncated. Too large slows down computation time.
yaug <- rbind(data$yobs, array(0, dim = c(nz, data$T)))
M <- nrow(yaug) # number of all potential individuals
T <- 3 # number of survey events (max detections = max number of surveys)
N = 100 # "True" population size to estimate
omega_point <- N / M # augmentation mediating parameter to estimate.

## Bundle data to feed the model
stan_data <- c("yaug", "M", "T")

## Initial values
inits <- function() {
  list(p = runif(1, 0, 1), omega = 0.5)}

## Parameters monitored
params <- c("N", "p", "omega")

## MCMC settings
n_iterations <- 2000
n_thin <- 1
n_burnin <- 1000
n_chains <- 4

stan_model <- "./mark_recapture_data_augmentation_M0.stan"

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
traceplot(out, pars = c("N", "p", "omega"), inc_warmup = TRUE, nrow = 2)

nobs_in_sim <- nrow(data$yobs) # number of individuals observed in the surveys
# i.e., the minimum population size or size if detection is perfect
color_scheme_set("pink")

# gather the mean pop size for plot from posterior distribution of N
matrix <- as.matrix(out)
N_mean = mean(matrix[,1])

# plot posterior distribution
p <- mcmc_hist(out, pars = c("N"))
p <- p + labs(x = "Population Size",
              y = "Frequency in 4000 Draws") +
  xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = nobs_in_sim, linetype = "dashed", size = 1) +
  geom_vline(xintercept = N_mean, linetype = "solid", size = 1)
p
    
mcmc_dens_overlay(out, pars = c("N", "p", "omega"))
mcmc_pairs(out, pars = c("N", "p", "omega"),
           off_diag_args = list(size = 1.5))
