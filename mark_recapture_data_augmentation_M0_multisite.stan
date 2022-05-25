data {
  int<lower=0> M; // Size of augumented data set
  int<lower=0> T; // Number of sampling occasions
  int<lower=0, upper=T> yaug[M, T]; // Capture-history matrix
}

transformed data {
  int s[M]; // Totals in each row
  int<lower=0> C; // Size of observed data set
  
  C = 0;
  for (i in 1:M) { // for each of the individuals in the augmented dataset
    s[i] = sum(yaug[i]); // s[i] is then the number of times each individ detected (0-3)
    if (s[i] > 0) { // if individual i is ever observed,
      C = C + 1; // we have another observed individual
    }
  }
}

parameters {
  real<lower=0, upper=1> omega; // Inclusion probability
  real<lower=0, upper=1> p; // Detection probability
}

model {
  // Priors are implicitly defined in the parameters block;
  //  omega ~ uniform(0, 1);
  //  p ~ uniform(0, 1);
  
  // Likelihood
  for (i in 1:M) { 
    // for each of the individuals in M (whether observed, not observed, or not 'real') 
    if (s[i] > 0) { // if detected at least once 
      // z[i] == 1, that is the individual must exist
      target += bernoulli_lpmf(1 | omega) + binomial_lpmf(s[i] | T, p);
      // likelihood of existing given estimate for omega
      // likelihood of s[i] detections (1-3) out of 3, given estimate for p
    } else // s[i] == 0, that is individual not detected
    { // could be a missed detection, or non-real individual
      target += log_sum_exp(bernoulli_lpmf(1 | omega)
                            // z[i] == 1 (the individual exists)
                            // likelihood of existing given estimate for omega
                            + binomial_lpmf(0 | T, p), 
                            // + likelihood was never detected, 
                            // given T visits and estimate for detection prob p
                            bernoulli_lpmf(0 | omega) // z[i] == 0
                            // was never detected because i is non-real, i.e., 
                            // position i in a theoretical population of size N is 
                            // not occupied by any individuals. 
                            );
    }
  }
}

generated quantities {
  // the expectation of N is M*omega
  // the generated value can't be less that C, since we did observe C individuals
  // so is instead represented by C plus the expected number of 
  // undetected individuals in the population (M-C individuals with prob. omega_nd)
  
  // prob present given never detected
  real omega_nd = (omega * (1 - p) ^ T) / (omega * (1 - p) ^ T + (1 - omega));
  int<lower=C, upper=M> N = C + binomial_rng(M - C, omega_nd);
}
