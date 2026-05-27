// =========================================================================
// 2nd Order DLM Model
// y_t      = theta_t1 + nu_t, nu_t ~ N (0, V)
// theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1,  omega_t1 ~ N(0, W1)
// theta_t2 =               + theta_{t-1,2} + omega_t2,  omega_t2 ~ N(0, W2)
// theta_01 ~ N(mu_01, sigma2_01)
// theta_02 ~ N(mu_02, sigma2_02)
// phi_V = 1/V ~ Gamma(nu_V, eta_V)
// phi1 = 1/W1 ~ Gamma(nu_01, eta_01)
// phi2 = 1/W2 ~ Gamma(nu_02, eta_02)
// =========================================================================

data {
  int<lower=1> T;          // number of observations
  array[T] real y; // observed data y_1, ..., y_T

  // Prior hyperparameters
  real mu_01;               // theta_01 ~ N(mu_01, sigma2_01)
  real<lower=0> sigma2_01;

  real mu_02;               // theta_02 ~ N(mu_02, sigma2_02)
  real<lower=0> sigma2_02;

  real<lower=0> nu_V;       // V  ~ Gamma(nu_V, eta_V)
  real<lower=0> eta_V;

  real<lower=0> nu_01;      // W1 ~ Gamma(nu_01, eta_01)
  real<lower=0> eta_01;

  real<lower=0> nu_02;      // W2 ~ Gamma(nu_02, eta_02)
  real<lower=0> eta_02;
}

parameters {
  array[T] real theta1;     // {theta_11, ..., theta_T1}
  array[T] real theta2;     // {theta_12, ..., theta_T2}

  real theta_01;
  real theta_02;

  real<lower=0> phi_V;      // phi_V = 1/V
  real<lower=0> phi1;       // phi1 = 1/W1
  real<lower=0> phi2;       // phi2 = 1/W2
}

transformed parameters {
  real<lower=0> V;
  V = 1.0 / phi_V;

  real<lower=0> W1;
  W1 = 1.0 / phi1;

  real<lower=0> W2;
  W2 = 1.0 / phi2;
}

model {
  // ---- Priors ----

  // theta_01 ~ N(mu_01, sigma2_01)
  // theta_02 ~ N(mu_02, sigma2_02)
  theta_01 ~ normal(mu_01, sqrt(sigma2_01));
  theta_02 ~ normal(mu_02, sqrt(sigma2_02));

  // phi_V = 1/V ~ Gamma(nu_V, eta_V)
  // phi1 = 1/W1 ~ Gamma(nu_01, eta_01)
  // phi2 = 1/W2 ~ Gamma(nu_02, eta_02)
  phi_V ~ gamma(nu_V, eta_V);
  phi1 ~ gamma(nu_01, eta_01);
  phi2 ~ gamma(nu_02, eta_02);

  // ---- System equation ----
  // theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1, omega_t1 ~ N(0,W1)
  // theta_t2 =                 theta_{t-1,2} + omega_t2, omega_t1 ~ N(0,W2)
  theta1[1] ~ normal(theta_01 + theta_02, sqrt(W1));
  for (t in 2:T) {
    theta1[t] ~ normal(theta1[t-1] + theta2[t-1], sqrt(W1));
  }

  theta2[1] ~ normal(theta_02, sqrt(W2));
  for (t in 2:T) {
    theta2[t] ~ normal(theta2[t-1], sqrt(W2));
  }

  // ---- Likelihood ----
  // y_t | theta_t1 ~ N(theta_t1, V)
  for (t in 1:T) {
    y[t] ~ normal(theta1[t], sqrt(V));
  }
}
