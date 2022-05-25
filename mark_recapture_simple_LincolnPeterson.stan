//
// This Stan program defines a Lincoln-Peterson mark-recapture model, 
// using data from captured and recaptured individuals to estimate
// population size.
// 
// M: number of animals marked in first capture
// C: number animals in second capture, and
// R: number of marked animals in second capture
// and the estimand of interest
// N: number of animals in the population.

// The L-P mathematical model is defined as: Nhat = (M*C)/R
// This population estimate would arise from a probabilistic model in 
// which the number of animals recaptured is binomially distributed:
// R ~ binom(C, M/N)

// The input data are integer values of M, C, and R 
// (will need to be transformed into vectors when multiple sites included)
data {
  int<lower = 0 > M;
  int<lower = 0 > C;
  int<lower = 0, upper= min(M, C) > R;
}

// The parameters accepted by the model. Our model
// accepts one parameter 'N'.
parameters {
  real<lower = (C - R + M) > N;
}

// The model to be estimated. We model the output
// Recapture should be composed of a binomial distribution 
// for each capture, with probability M / N
model {
  // priors
  N ~ cauchy(0, 1000);
  
  // likelihood
  R ~ binomial(C, M / N);
}
