// Binomial mixture model with covariates

// adapted from:
// https://discourse.mc-stan.org/t/calculating-abundance-from-marginalized-n-mixture-model/12751
// and from Kery and Schaub: population analysis using WinBUGS

// same as model M1 but without random site effect

data {
  int<lower=0> R;       // Number of observation units
  int<lower=0> T;       // Number of temporal replications ( = number of visits = 3)
  
  int<lower=0> y[R, T]; // Counts at observation 1:R on visits 1:T
  
  vector[R] X;          // Covariate ( = X: mowed or unmowed)
  
  int<lower=0> n_sites;       // Number of species
  int<lower=1> sites[R];      // vector of species
  
  int<lower=0> n_species;       // Number of species
  int<lower=1> species[R];      // vector of species
  
  int<lower=0> K;       // Upper bound of population size
}

transformed data {
  int<lower=0> max_y[R];
  int<lower=0> N_ll;
  int tmp[R];
  
  for (i in 1:R) {
    max_y[i] = max(y[i]);
    tmp[i] = K - max_y[i] + 1;
  }
  
  N_ll = sum(tmp);
}

parameters {
  real alpha0;
  real alpha1;

  real beta0;
  real beta1;
  
  vector[n_species] eps;         // Random (species) effects on abundance
  real<lower=0> sd_eps;          // Variation (scale) in (species) effects on abundance
  
  vector[n_species] zeta;        // Random (species) effects on detection
  real<lower=0> sd_zeta;         // Variation (scale) in (species) effects on detection
  
  vector[n_species] eta;     // Random (species) effects on abundance slope
  real<lower=0> sd_eta;   // Variation (scale) in (species) effects on abundance slope
  
}

transformed parameters {
  vector[R] log_lambda; // Log population size
  vector[R] logit_p; // Logit detection probability
  
  vector[n_species] eps_star; // scaled species random effects
  vector[n_species] zeta_star; // scaled species random effects
  
  vector[n_species] eta_star; // scaled species random effects

  eps_star = eps * sd_eps;
  zeta_star = zeta * sd_zeta;
  
  eta_star = eta * sd_eta;
  
  for (i in 1:R) { // for each site
    log_lambda[i] = alpha0 + (alpha1 * X[i]) + 
        eps_star[species[i]] +
        (eta_star[species[i]] * X[i]); // non-centered formulation of random effect (see Monnahan et al. 2017)
    logit_p[i] = beta0 + beta1 * X[i] + 
        zeta_star[species[i]]; // non-centered formulation 
  }
}

model {
  // Priors
  alpha0 ~ normal(0, 5);
  alpha1 ~ normal(0, 5);
  
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  
  eps ~ normal(0, 1);
  sd_eps ~ cauchy(0, 2.5);
  zeta ~ normal(0, 1);
  sd_zeta ~ cauchy(0, 2.5);
  
  eta ~ normal(0, 1);
  sd_eta ~ cauchy(0, 2.5);
  
  // Likelihood
  for (i in 1:R) { // for each siteXspecies
    vector[K - max_y[i] + 1] lp; // lp vector of length of possible abundances 
    // (from max observed to K) 
    
    for (j in 1:(K - max_y[i] + 1)) // for each possible abundance:
    // lp of abundance given ecological model and observational model
    lp[j] = poisson_log_lpmf(max_y[i] + j - 1 | log_lambda[i])
      + binomial_logit_lpmf(y[i] | max_y[i] + j - 1, logit_p[i]); // vectorized over T
    target += log_sum_exp(lp);
  }
  
}

//generated quantities {
//  int<lower=0> N[R]; // predicted abundance at each site for each species
//  vector[R] p; // detection probability of speciesXsite combos
  
//  matrix[R, T] eval; // Expected values
  
//  int y_new[R, T]; // new data for counts generated from eval
    
//  matrix[R, T] E; // squared scaled distance of real data from expected value
//  matrix[R, T] E_new; // squared scaled distance of new data from expected value
  
//  real fit = 0; // sum squared distances of real data across all observation intervals
//  real fit_new = 0; // sum squared distances of new data across all observation intervals

//  vector[n_species] totalN; // total pop size PER SPECIES
  
  // predict abundance given log_lambda
//  for (i in 1:R) {
//    N[i] = poisson_log_rng(log_lambda[i]);
//    p[i] = inv_logit(logit_p[i]);
//  }
  
  // Bayesian p-value fit. 
    
  // Initialize E and E_new
//  for (i in 1:1) {
//    for(j in 1:T) {
//      E[i, j] = 0;
//      E_new[i, j] = 0;
//    }
//  }
  
//  for (i in 2:R) {
//    E[i] = E[i - 1];
//    E_new[i] = E_new[i - 1];
//  }
  
//  for (i in 1:R) {
//    for (j in 1:T) {
      // Assess model fit using Chi-squared discrepancy
      // Compute fit statistic E for observed data
//      eval[i, j] = p[i] * N[i]; // expected value at observation i for visit j 
        // (probabilty across visits is fixed) is = expected detection prob * expected abundance
      // Compute fit statistic E_new for real data (y)
//      E[i, j] = square(y[i, j] - eval[i, j]) / (eval[i, j] + 0.5);
      // Generate new replicate data and
//      y_new[i, j] = binomial_rng(N[i], p[i]);
      // Compute fit statistic E_new for replicate data
//      E_new[i, j] = square(y_new[i, j] - eval[i, j]) / (eval[i, j] + 0.5);
//    }
    
//    fit = fit + sum(E[i]); // descrepancies for each siteXspecies combo (across 1:T visits)
//    fit_new = fit_new + sum(E_new[i]); // descrepancies for generated data for each 
                                      // siteXspecies combos (across 1:T visits)
//  }
  
  // to find total abundance PER SPECIES:
  // sum abundance by rows for each species i
//  for (i in 1:n_species){
//    vector[n_sites] totalSpecies;
//    for (j in 1:n_sites){
//      totalSpecies[j] = N[(i+((n_species*j)-n_species))];
//    }
//    totalN[i] = sum(totalSpecies);
//  }
//}
