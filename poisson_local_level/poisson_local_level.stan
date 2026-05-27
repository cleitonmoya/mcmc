// ============================================================
// Poisson Local Level Modelo
// y_t ~ Poisson(exp(theta_t))
// theta_t = theta_{t-1} + omega_t,  omega_t ~ N(0, W)
// theta_0 ~ N(mu_0, sigma2_0)
// phi = 1/W ~ Gamma(nu_0, eta_0)  [parametrization: shape, rate]
// ============================================================

data {
  int<lower=1> T;          // number of observations
  array[T] int<lower=0> y; // observed data y_1, ..., y_T

  // Prior hyperparameters
  real mu_0;               // média prior de theta_0
  real<lower=0> sigma2_0;  // variância prior de theta_0
  real<lower=0> nu_0;      // shape da Gamma para phi = 1/W
  real<lower=0> eta_0;     // rate  da Gamma para phi = 1/W
}

parameters {
  array[T] real theta;     // estados latentes theta_1, ..., theta_T
  real theta_init;         // theta_0 (estado inicial)
  real<lower=0> phi;       // phi = 1/W (precisão do sistema)
}

transformed parameters {
  real<lower=0> W;
  W = 1.0 / phi;           // variância do sistema
}

model {
  // ---- Priori ----

  // theta_0 ~ N(mu_0, sigma2_0)
  theta_init ~ normal(mu_0, sqrt(sigma2_0));

  // phi = 1/W ~ Gamma(nu_0, eta_0)
  phi ~ gamma(nu_0, eta_0);

  // ---- Equação do sistema ----
  // theta_t | theta_{t-1} ~ N(theta_{t-1}, W)

  theta[1] ~ normal(theta_init, sqrt(W));
  for (t in 2:T) {
    theta[t] ~ normal(theta[t-1], sqrt(W));
  }

  // ---- Verossimilhança ----
  // y_t | theta_t ~ Poisson(exp(theta_t))

  for (t in 1:T) {
    y[t] ~ poisson_log(theta[t]);
    // poisson_log(theta) é equivalente a Poisson(exp(theta))
    // mas numericamente mais estável
  }
}

generated quantities {
  // Médias ajustadas (para diagnóstico e visualização)
  array[T] real lambda_hat;
  for (t in 1:T) {
    lambda_hat[t] = exp(theta[t]);
  }
}
