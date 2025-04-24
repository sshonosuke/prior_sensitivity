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


transformed parameters{
  real<lower=0> trans1;
  real<lower=0> trans2;
  trans1 = - log(alpha / (alpha + beta));
  trans2 = 1 / sqrt(alpha + beta);
}


model{
  y ~ binomial(N, prob);
  prob ~ beta(alpha, beta);
  trans1 ~ gamma(1, 1);
  trans2 ~ gamma(1, 1);
}
