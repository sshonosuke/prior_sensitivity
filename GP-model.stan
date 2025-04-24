data {
  int<lower=1> n;           // Number of data points
  real x[n];                // Input data as an array
  vector[n] y;              // Output data (target)
}

parameters {
  real<lower=0> tau_mu;     // Kernel amplitude parameter for mean GP
  real<lower=0> psi_mu;     // Kernel length-scale parameter for mean GP
  vector[n] mu;             // Latent mean vector for Gaussian Process
  real<lower=0> sig;        
}

model {
  matrix[n, n] cov_mu;      // Covariance matrix for mean GP
  vector[n] GP_mean;        // Zero mean vector
  
  // Covariance matrices for mean GP and variance GP using the squared exponential kernel
  cov_mu = cov_exp_quad(x, tau_mu, psi_mu);
  
  // Cholesky decomposition of covariance matrices with jitter for numerical stability
  matrix[n, n] L_mu = cholesky_decompose(cov_mu + diag_matrix(rep_vector(1e-6, n)));
  
  // Zero mean vector
  GP_mean = rep_vector(0, n);
  
  // Gaussian Process models for mean and log standard deviation
  mu ~ multi_normal_cholesky(GP_mean, L_mu);
  
  // Likelihood of observed data
  y ~ normal(mu, sig);
  
  // Priors for hyperparameters
  //tau_mu ~ cauchy(0, 2);
  //psi_mu ~ cauchy(0, 2);
  tau_mu ~ gamma(1, 1);
  psi_mu ~ gamma(1, 1);
}
