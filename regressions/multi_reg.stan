data {
  int<lower=1> K;
  int<lower=1> J;
  int<lower=1> n_obs;
  matrix[n_obs,J] x;
  vector[K] y[n_obs];
  vector[J] b0;
  matrix[J,J] inv_XX;
}

parameters {
  matrix[K, J] delta;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;
  real<lower=0> sigma;
}

transformed parameters {
  vector[K] mu[n_obs];
  matrix[K, K] L_Sigma;
  real<lower=0> sigma_sq;
  cov_matrix[J] B0;
  
  for (i in 1:n_obs)
    mu[i] = delta * (x[i,])';

  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  
  sigma_sq = sigma^2;
  
  B0 = n_obs * sigma_sq * inv_XX;
}

model {
  for (i in 1:K){
    delta[i,] ~ multi_normal(b0, B0);
  }
  L_Omega ~ lkj_corr_cholesky(0.5);
  L_sigma ~ normal(0, 10);
  sigma ~ normal(0, 10);
  
  y ~ multi_normal_cholesky(mu, L_Sigma);
}

generated quantities {
  real dev;                    // deviance
  vector[n_obs] log_lik;       // log-likelihood for loo package
  matrix[K,K] Omega;
  matrix[K,K] Sigma;
  
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  Sigma = quad_form_diag(Omega, L_sigma);
  
  for (i in 1:n_obs) {
    log_lik[i] = multi_normal_cholesky_lpdf(y[i]| mu[i], L_Sigma);
 }
    dev = (-2) * sum(log_lik);
}