// =========================================================================
// 2nd Order DLM Model
// y_t      = theta_t1 + nu_t, nu_t ~ N (0, V)
// theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1,  omega_t1 ~ N(0, W1)
// theta_t2 =               + theta_{t-1,2} + omega_t2,  omega_t2 ~ N(0, W2)
// theta_01 ~ N(mu_01, sigma2_01)
// theta_02 ~ N(mu_02, sigma2_02)
// sigma_V  ~ Exp(lambda_V)
// sigma_W1 ~ Exp(lambda_W1)
// sigma_W2 ~ Exp(lambda_W2)
// =========================================================================

data {
  int<lower=1> T;          // number of observations
  array[T] real y; // observed data y_1, ..., y_T

  // theta_01 hyperparameters
  real mu_01;               // theta_01 ~ N(mu_01, sigma2_01)
  real<lower=0> sigma2_01;

  // theta_02 hyperparameters
  real mu_02;               // theta_02 ~ N(mu_02, sigma2_02)
  real<lower=0> sigma2_02;

  // PC Prior hyperparameters - P(sigma > U) = alpha
  real<lower=0,upper=1> alpha;
  real<lower=0> U_V;
  real<lower=0> U_W1;
  real<lower=0> U_W2;
}

transformed data {
  // PC Prior hyperparameters
  real lambda_V  = -log(alpha) / U_V;
  real lambda_W1 = -log(alpha) / U_W1;
  real lambda_W2 = -log(alpha) / U_W2;
}

parameters {
  array[T] real theta1;     // {theta_11, ..., theta_T1}
  array[T] real theta2;     // {theta_12, ..., theta_T2}

  real theta_01;
  real theta_02;

  real<lower=0> sigma_V;
  real<lower=0> sigma_W1;
  real<lower=0> sigma_W2;
}

transformed parameters {
  real<lower=0> V  = sigma_V^2;
  real<lower=0> W1 = sigma_W1^2;
  real<lower=0> W2 = sigma_W2^2;
}

model {
  // ---- Priors ----

  // theta_01 ~ N(mu_01, sigma2_01)
  // theta_02 ~ N(mu_02, sigma2_02)
  theta_01 ~ normal(mu_01, sqrt(sigma2_01));
  theta_02 ~ normal(mu_02, sqrt(sigma2_02));

  // PC Prior
  sigma_V  ~ exponential(lambda_V);
  sigma_W1 ~ exponential(lambda_W1);
  sigma_W2 ~ exponential(lambda_W2);

  // ---- System equation ----
  // theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1, omega_t1 ~ N(0,W1)
  // theta_t2 =                 theta_{t-1,2} + omega_t2, omega_t1 ~ N(0,W2)
  theta1[1] ~ normal(theta_01 + theta_02, sigma_W1);
  for (t in 2:T) {
    theta1[t] ~ normal(theta1[t-1] + theta2[t-1], sigma_W1);
  }

  theta2[1] ~ normal(theta_02, sigma_W2);
  for (t in 2:T) {
    theta2[t] ~ normal(theta2[t-1], sigma_W2);
  }

  // ---- Likelihood ----
  // y_t | theta_t1 ~ N(theta_t1, V)
  for (t in 1:T) {
    y[t] ~ normal(theta1[t], sigma_V);
  }
}
