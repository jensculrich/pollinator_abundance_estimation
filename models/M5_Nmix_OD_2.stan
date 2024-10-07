// Binomial mixture model with covariates
// including random effect of abundance to account for overdispersion in data

// adapted from:
// https://discourse.mc-stan.org/t/calculating-abundance-from-marginalized-n-mixture-model/12751
// and from Kery and Schaub: population analysis using WinBUGS

// same as model M4 but with random effects redefined 
// as in model M4, M5 includes an additional latent effect in the linear predictor of
// abundance (phi) for scaling of abundance error. Expected to correct for extra-Poisson dispersion in the data, 
// but at the cost of decreased certainty in the parameter estimates.

data {
  int<lower=0> R;       // Number of observation units
  int<lower=0> T;       // Number of temporal replications ( = number of visits = 3)
  
  int<lower=0> y[R, T]; // Counts at observation 1:R on visits 1:T
  
  vector[R] X;          // Covariate ( = X: mowed or unmowed)
  
  int<lower=0> n_sites;       // Number of sites
  int<lower=1> sites[R];      // vector of sites
  
  int<lower=0> n_species;       // Number of species
  int<lower=1> species[R];      // vector of species
  
  int<lower=0> K;       // Upper bound of population size
}

transformed data {
  int<lower=0> max_y[R];
  //int<lower=0> N_ll;
  //int tmp[R];
  
  for (i in 1:R) {
    max_y[i] = max(y[i]);
    //tmp[i] = K - max_y[i] + 1;
  }
  
  //N_ll = sum(tmp);
}

parameters {
  real alpha0; // global interecept for abundance
  
  // species specific intercept allows some species to be more abundant than others, 
  // but with overall estimates for abundances partially informed by the data pooled across all species.
  vector[n_species] alpha0_species; // species specific intercept for abundance
  real<lower=0> sigma_alpha0_species; // variance in species intercepts

  // random slope for the effect of management on abundance for each species
  vector[n_species] alpha1; // vector of species specific slope estimates
  real mu_alpha1; // community mean of species specific slopes
  real<lower=0> sigma_alpha1_species; // variance in species slopes

  real beta0; // global interecept for detection probability
  
  // species specific intercept allows some species to be easier to detect than others, 
  // but with overall estimates for detection partially informed by the data pooled across all species.
  vector[n_species] beta0_species; // species specific intercept for detection
  real<lower=0> sigma_beta0_species; // variance in species intercepts

  // random slope for the effect of management for each species
  vector[n_species] beta1; // vector of species specific slope estimates
  real mu_beta1; // community mean of species specific slopes
  real<lower=0> sigma_beta1_species; // variance in species slopes

  real<lower=0> phi; // abundance overdispersion parameter
}

transformed parameters {
  vector[R] log_mu; // Log population size (mu as the centering parameter
    // for the negative binomial count distribution).
  vector[R] logit_p; // Logit detection probability

  
  for (i in 1:R) { // for each site*species combination
    log_mu[i] = // log abundance is equal to..
      alpha0 + // a global intercept plus
      alpha0_species[species[i]] + // a species specific intercept plus 
      (alpha1[species[i]] * X[i]); // a species specific slope
    logit_p[i] = // logit (detection probability) is equal to..
      beta0 + // a global intercept plus
      beta0_species[species[i]] + // a species specific intercept plus 
      (beta1[species[i]] * X[i]); // a species specific slope
  }
}

model {
  // Priors
  
  // Abundance (Ecological Process)
  alpha0 ~ normal(0, 5); // global intercept for abundance
  
  alpha0_species ~ normal(0, sigma_alpha0_species); 
  // abundance intercept for each species drawn from the community
  // distribution (variance defined by 1/sigma^2), centered at 0. 
  sigma_alpha0_species ~ cauchy(0, 2.5);
  
  alpha1 ~ normal(mu_alpha1, sigma_alpha1_species);
  // abundance slope (effect of management) for each species drawn from the community
  // distribution (variance defined by 1/sigma^2), centered at mu_alpha1. 
  // centering on mu (rather than 0) allows us to estimate the average effect of
  // the management on abundance across all species.
  mu_alpha1 ~ normal(0, 5); // community mean
  sigma_alpha1_species ~ cauchy(0, 2.5); // community variance
  
  // Detection (Observation Process)
  beta0 ~ normal(0, 5); // global intercept for abundance
  
  beta0_species ~ normal(0, sigma_beta0_species); 
  // abundance intercept for each species drawn from the community
  // distribution (variance defined by 1/sigma^2), centered at 0. 
  sigma_beta0_species ~ cauchy(0, 2.5);
  
  beta1 ~ normal(mu_beta1, sigma_beta1_species);
  // abundance slope (effect of management) for each species drawn from the community
  // distribution (variance defined by 1/sigma^2), centered at mu_alpha1. 
  // centering on mu (rather than 0) allows us to estimate the average effect of
  // the management on abundance across all species.
  mu_beta1 ~ normal(0, 5); // community mean
  sigma_beta1_species ~ cauchy(0, 2.5); // community variance
  
  // abundance overdispersion scale parameter
  phi ~ cauchy(0, 2.5);
  
  // Likelihood
  for (i in 1:R) { // for each siteXspecies
    vector[K - max_y[i] + 1] lp; // lp vector of length of possible abundances 
    // (from max observed to K) 
    
    for (j in 1:(K - max_y[i] + 1)) // for each possible abundance:
    // lp of abundance given ecological model and observational model
    lp[j] = neg_binomial_2_lpmf(max_y[i] + j - 1 | exp(log_mu[i]), phi)
      + binomial_logit_lpmf(y[i] | max_y[i] + j - 1, logit_p[i]); // vectorized over T
    target += log_sum_exp(lp);
  }
  
}

generated quantities {
  
  // Posterior Predictive Check
  
  int<lower=0> N[R]; // predicted abundance at each site for each species
  vector[R] p; // detection probability of speciesXsite combos
  
  matrix[R, T] eval; // Expected values
  
  int y_new[R, T]; // new data for counts generated from eval
    
  matrix[R, T] E; // squared scaled distance of real data from expected value
  matrix[R, T] E_new; // squared scaled distance of new data from expected value
  
  real fit = 0; // sum squared distances of real data across all observation intervals
  real fit_new = 0; // sum squared distances of new data across all observation intervals

  vector[n_species] totalN; // total pop size PER SPECIES
 
  // predict abundance given log_mu
  for (i in 1:R) {
    N[i] = neg_binomial_2_rng(exp(log_mu[i]), phi);
    p[i] = inv_logit(logit_p[i]);
  }
    
  // Initialize E and E_new
  for (i in 1:1) {
    for(j in 1:T) {
      E[i, j] = 0;
      E_new[i, j] = 0;
    }
  }
  
  for (i in 2:R) {
    E[i] = E[i - 1];
    E_new[i] = E_new[i - 1];
  }
  
  for (i in 1:R) {
    for (j in 1:T) {
      // Assess model fit using Chi-squared discrepancy
      // Compute fit statistic E for observed data
      eval[i, j] = inv_logit(logit_p[i]) * neg_binomial_2_rng(exp(log_mu[i]), phi); // expected value at observation i for visit j 
        // (probabilty across visits is fixed) is = expected detection prob * expected abundance
      // Compute fit statistic E_new for real data (y)
      E[i, j] = square(y[i, j] - eval[i, j]) / (eval[i, j] + 0.5);
      // Generate new replicate data and
      y_new[i, j] = binomial_rng(N[i], p[i]);
      // Compute fit statistic E_new for replicate data
      E_new[i, j] = square(y_new[i, j] - eval[i, j]) / (eval[i, j] + 0.5);
    }
    
    fit = fit + sum(E[i]); // descrepancies for each siteXspecies combo (across 1:T visits)
    fit_new = fit_new + sum(E_new[i]); // descrepancies for generated data for each 
                                      // siteXspecies combos (across 1:T visits)
  }
  
  // to find total abundance PER SPECIES:
  // sum abundance by rows for each species i
  for (i in 1:n_species){
    vector[n_sites] totalSpecies;
    for (j in 1:n_sites){
      totalSpecies[j] = N[(i+((n_species*j)-n_species))];
    }
    totalN[i] = sum(totalSpecies);
  }
}
