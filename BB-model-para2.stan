data{
  int<lower=0> J;     // number of experiments
  int<lower=0> N[J];  // sample size in each experiment
  int<lower=0> y[J];  // observation in each experiment
}


parameters{
  real<lower=0> alpha;
  real<lower=0> beta;
  vector<lower=0, upper=1>[J] prob;
}



model{
  y ~ binomial(N, prob);
  prob ~ beta(alpha, beta);
  alpha ~ gamma(1, 1);
  beta ~ gamma(1, 1);
}
