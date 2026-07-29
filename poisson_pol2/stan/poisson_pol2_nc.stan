// Poisson Local Trend Model
//
// Model:
//   y_t ~ Poisson(exp(theta_t1))
//   theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1,  omega_t1 ~ N(0, W1)
//   theta_t2 =               + theta_{t-1,2} + omega_t2,  omega_t2 ~ N(0, W2)
//
// Priors:
//   theta_01 ~ N(mu_01, sigma2_01)
//   theta_02 ~ N(mu_02, sigma2_02)
//   phi1 = 1/W1 ~ Gamma(nu_01, eta_01)
//   phi2 = 1/W2 ~ Gamma(nu_02, eta_02)
//
// Parameterization: NON-CENTERED for theta1/theta2. W1, W2 are themselves
// estimated parameters that scale the innovations of theta1/theta2 (a
// classic Neal's funnel setup - HMC step size struggles to adapt to both
// the "funnel neck" (small W) and the "funnel mouth" (large W) under a
// single global step size). theta_01/theta_02 are NOT reparameterized: their
// prior scales (sigma2_01, sigma2_02) are fixed data, not estimated
// parameters, so there is no funnel there. z1/z2 (raw N(0,1) innovations)
// are the actual sampled parameters; theta1/theta2 are rebuilt
// deterministically in transformed parameters.
//
// Author: Cleiton Moya de Ameida

data {
  int<lower=1> Tt; // number of observations
  array[Tt] int<lower=0> y; // observed data y_1, ..., y_T

  // Prior hyperparameters
  real mu_01;
  real<lower=0> sigma2_01;

  real mu_02;
  real<lower=0> sigma2_02;

  real<lower=0> nu_01;
  real<lower=0> eta_01;

  real<lower=0> nu_02;
  real<lower=0> eta_02;
}

parameters {
  array[Tt] real z1; // raw (non-centered) innovations for theta1, z1_t ~ N(0,1)
  array[Tt] real z2; // raw (non-centered) innovations for theta2, z2_t ~ N(0,1)
  real theta_01;
  real theta_02;
  real<lower=0> phi1;
  real<lower=0> phi2;
}

transformed parameters {
  array[Tt] real theta1;
  array[Tt] real theta2;

  real<lower=0> W1;
  W1 = 1.0 / phi1;

  real<lower=0> W2;
  W2 = 1.0 / phi2;

  {
    real sd_W1 = sqrt(W1);
    real sd_W2 = sqrt(W2);

    theta2[1] = theta_02 + sd_W2 * z2[1];
    theta1[1] = theta_01 + theta_02 + sd_W1 * z1[1];
    for (t in 2:Tt) {
      theta2[t] = theta2[t-1] + sd_W2 * z2[t];
      theta1[t] = theta1[t-1] + theta2[t-1] + sd_W1 * z1[t];
    }
  }
}

model {
  // Priors
  theta_01 ~ normal(mu_01, sqrt(sigma2_01));
  theta_02 ~ normal(mu_02, sqrt(sigma2_02));

  phi1 ~ gamma(nu_01, eta_01);
  phi2 ~ gamma(nu_02, eta_02);

  // Non-centered innovations (replace the previous direct system-equation
  // sampling statements on theta1/theta2)
  z1 ~ std_normal();
  z2 ~ std_normal();

  // Observation equation
  // y_t | theta_t1 ~ Poisson(exp(theta_t1))
  for (t in 1:Tt) {
    y[t] ~ poisson_log(theta1[t]);
  }
}
