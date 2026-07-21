# Poisson - 2nd Order Polynomial Dynamic Model
# Particle Gibbs (PG) with Backward Sampling + Component-wise Gibbs for theta_t2
# Strategy: SMC only for theta_t1 (scalar state); theta_t2 sampled via
#           component-wise Gibbs (Normal conjugate full conditionals)
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269-342.
# pg4: theta_t1 marginalized
# Author: Cleiton Moya de Almeida

library(invgamma)

rm(list = ls())
set.seed(42)
options(error = function() traceback(2))
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source_file <- "poisson_sin_2000"
datA <- readRDS(paste("../data/", source_file, ".rds", sep=""))
t_observed <- c(200, 500, 1000, 1500)

#y <- as.data.frame(Seatbelts)$DriversKilled
# t_observed  <- c(50, 75, 100, Tt)
y <- data$y
Tt <- length(y)

theta1_present <- TRUE
theta2_present <- FALSE

if (theta1_present) {
    theta1_true <- data$theta1
    lambda_true <- exp(theta1_true)
}

if (theta2_present) theta2_true <- data$theta2

# Auxiliary functions
printf <- function(...) cat(paste(sprintf(...), "\n"))

logsumexp <- function(x) {
  cc <- max(x)
  return(cc + log(sum(exp(x - cc))))
}

log_p_yt <- function(yt, theta_t1) {
  res <- yt * theta_t1 - exp(theta_t1)
  res[!is.finite(res)] <- -Inf
  return(res)
}

# Prior Hyperparameters

# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- log(y[1] + 0.5)
sigma2_01 <- 10

# theta_02 ~ N(mu_02, sigma2_02)
mu_02     <- 0
sigma2_02 <- 1

# W1 ~ InvGamma(alpha_W1, beta_W1)
alpha_W1 <- 2
beta_W1  <- 0.01

# W2 ~ InvGamma(alpha_W2, beta_W2)
alpha_W2 <- 2
beta_W2  <- 0.001

N      <- 1000
K      <- 200
burnin <- 200

W1_hist          <- numeric(N)
W2_hist          <- numeric(N)
theta_01_hist    <- numeric(N)
theta_02_hist    <- numeric(N)
theta1_star_hist <- matrix(0, N, Tt)
theta2_star_hist <- matrix(0, N, Tt)

# Initial Values
theta1_star <- stats::filter(log(y + 0.5), rep(1/5, 5), sides=2)
theta1_star[is.na(theta1_star)] <- log(y[is.na(theta1_star)] + 0.5)
theta1_star <- as.numeric(theta1_star)
theta2_star <- c(0, diff(theta1_star))

theta_01 <- log(y[1] + 0.5)
theta_02 <- 0
W1 <- 0.01
W2 <- 0.001

#####
start_time <- proc.time()

ess_smc <- numeric(N)

for (n in 1:N) {

  # ------------------------------------------------------------------
  # STEP 1: Sample theta_02 | theta_01, theta1_star[1], theta2_star[1], W1, W2
  # Termos: omega_11 = theta1_star[1] - theta_01 - theta_02 ~ N(0,W1)
  #         omega_12 = theta2_star[1] - theta_02            ~ N(0,W2)
  #         priori: theta_02 ~ N(mu_02, sigma2_02)
  # ------------------------------------------------------------------
  sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
  mu_02_bar <- sigma2_02_bar * (
    (theta1_star[1] - theta_01) / W1 +
    theta2_star[1]              / W2 +
    mu_02 / sigma2_02
  )
  theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))

  # ------------------------------------------------------------------
  # STEP 2: Sample theta_01 | theta_02, theta1_star[1], W1
  # Termos: omega_11 = theta1_star[1] - theta_01 - theta_02 ~ N(0,W1)
  #         priori: theta_01 ~ N(mu_01, sigma2_01)
  # ------------------------------------------------------------------
  sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
  mu_01_bar <- sigma2_01_bar * (
    mu_01 / sigma2_01 +
    (theta1_star[1] - theta_02) / W1
  )
  theta_01 <- rnorm(1, mean=mu_01_bar, sd=sqrt(sigma2_01_bar))

  # ------------------------------------------------------------------
  # STEP 3: Sample W1 | theta_01, theta_02, theta1_star, theta2_star
  # omega_t1 = theta1_star[t] - theta1_star[t-1] - theta2_star[t-1]
  # ------------------------------------------------------------------
  dif1   <- theta1_star - c(theta_01, theta1_star[-Tt])
  diffs1 <- dif1 - c(theta_02, theta2_star[-Tt])
  alpha_W1_bar <- alpha_W1 + Tt/2
  beta_W1_bar  <- beta_W1 + 0.5 * sum(diffs1^2)
  W1 <- rinvgamma(1, shape=alpha_W1_bar, rate=beta_W1_bar)

  # ------------------------------------------------------------------
  # STEP 4: Sample W2 | theta_02, theta2_star
  # omega_t2 = theta2_star[t] - theta2_star[t-1]
  # ------------------------------------------------------------------
  diffs2 <- theta2_star - c(theta_02, theta2_star[-Tt])
  alpha_W2_bar <- alpha_W2 + Tt/2
  beta_W2_bar  <- beta_W2 + 0.5 * sum(diffs2^2)
  W2 <- rinvgamma(1, shape=alpha_W2_bar, rate=beta_W2_bar)

  sd_W1 <- sqrt(W1)
  sd_W2 <- sqrt(W2)
  sd = sqrt(W1 + W2)
  # ------------------------------------------------------------------
  # STEP 5: Conditional SMC para theta_t1 apenas
  # Proposta bootstrap: theta_t1^k ~ N(theta_{t-1,1}^k + theta2_star[t-1], W1)
  # theta2_star é tratado como fixo (da iteração anterior)
  # Trajetória condicionada na coluna K
  # ------------------------------------------------------------------
  theta_1_k   <- matrix(0, Tt, K)
  log_w_tilde <- matrix(0, Tt, K)

  # t = 1
  theta_1_k[1, ] <- rnorm(K, mean=theta_01 + theta_02, sd=sd_W1)
  theta_1_k[1, K] <- theta1_star[1]

  log_w_1 <- log_p_yt(y[1], theta_1_k[1, ])
  log_w_tilde[1, ] <- log_w_1 - logsumexp(log_w_1) # normalizing

  # t = 2, ..., T
  for (t in 2:Tt) {

    A <- sample(1:K, K, replace=TRUE,
                prob=exp(log_w_tilde[t-1, ] - max(log_w_tilde[t-1, ])))
    A[K] <- K   # <-- preserva a trajetoria de referencia (coluna K)

    # Proposta: theta_t1^k ~ N(theta_{t-1,1}^k + theta2_star[t-1], W1)
    # theta2_star[t-1] é escalar fixo — da iteração MCMC anterior
    drift <- if(t == 2) theta_02 else theta2_star[t-2]
    theta_t1_k <- rnorm(K, mean=theta_1_k[t-1, A] + drift, sd=sd)
    theta_t1_k[K] <- theta1_star[t]
    theta_1_k[t, ] <- theta_t1_k

    log_w_t <- log_p_yt(y[t], theta_1_k[t, ])
    log_w_tilde[t, ] <- log_w_t - logsumexp(log_w_t)
  }


  ess_smc[n] <- 1 / sum(exp(2*(log_w_tilde[Tt, ] - logsumexp(log_w_tilde[Tt, ]))))

  # ------------------------------------------------------------------
  # STEP 6: Backward sampling para theta1_star
  # Pesos backward: w_t^k * N(theta1_star[t+1] | theta_t1^k + theta2_star[t], W1)
  # theta2_star[t] é escalar fixo — da iteração MCMC anterior
  # ------------------------------------------------------------------
  k_final <- sample(1:K, 1,
                    prob=exp(log_w_tilde[Tt, ] - max(log_w_tilde[Tt, ])))
  theta1_star[Tt] <- theta_1_k[Tt, k_final]

  for (t in (Tt-1):1) {

    drift_bw <- if(t == 1) theta_02 else theta2_star[t-1]
    log_bw <- log_w_tilde[t, ] +
              dnorm(theta1_star[t+1],
                    mean=theta_1_k[t, ] + drift_bw, sd=sd, log=TRUE)
    log_bw <- log_bw - max(log_bw)
    bw <- exp(log_bw)
    bw <- bw / sum(bw)
    b  <- sample(1:K, 1, prob=bw)
    theta1_star[t] <- theta_1_k[t, b]
  }

  # ------------------------------------------------------------------
  # STEP 7: Component-wise Gibbs para theta_t2 | theta1_star, W1, W2
  #
  # Para t = 1, ..., T-1:
  #   Termos que contêm theta_t2:
  #     omega_{t+1,1} = theta1_star[t+1] - theta1_star[t] - theta_t2  ~ N(0,W1)
  #     omega_{t+1,2} = theta2_star[t+1] - theta_t2                   ~ N(0,W2)
  #     omega_{t2}    = theta_t2 - theta2_star[t-1]                   ~ N(0,W2)
  #
  #   sigma2_t2_bar = (1/W1 + 2/W2)^(-1)
  #   mu_t2_bar = sigma2_t2_bar * (
  #     (theta1_star[t+1] - theta1_star[t]) / W1 +
  #     theta2_star[t+1]                    / W2 +
  #     theta2_star[t-1]                    / W2
  #   )
  #
  # Para t = T:
  #   Só omega_{T2} = theta_T2 - theta2_star[T-1] ~ N(0,W2)
  #   (theta_T2 não influencia nenhum tempo futuro)
  #
  #   sigma2_T2_bar = W2
  #   mu_T2_bar     = theta2_star[T-1]
  # ------------------------------------------------------------------

  # Variâncias condicionais (constantes para todos os t interiores)
  sigma2_t2_bar_interior <- (1/W1 + 2/W2)^(-1)
  sigma2_t2_bar_last     <- W2

  # t = 1: theta2_star[-1] = theta_02 (estado inicial)
  # Tratar separadamente pois theta2_star[t-1] = theta_02
  sigma2_bar <- (1/W1 + 2/W2)^(-1)
  mu_bar <- sigma2_bar * (
    (theta1_star[2] - theta1_star[1]) / W1 +
    theta2_star[2]                    / W2 +
    theta_02                          / W2
  )
  theta2_star[1] <- rnorm(1, mean=mu_bar, sd=sqrt(sigma2_bar))

  # t = 2, ..., T-1
  for (t in 2:(Tt-1)) {
    mu_bar <- sigma2_t2_bar_interior * (
      (theta1_star[t+1] - theta1_star[t]) / W1 +
      theta2_star[t+1]                    / W2 +
      theta2_star[t-1]                    / W2
    )
    theta2_star[t] <- rnorm(1, mean=mu_bar, sd=sqrt(sigma2_t2_bar_interior))
  }

  # t = T: só evolução própria, sem termos futuros
  theta2_star[Tt] <- rnorm(1, mean=theta2_star[Tt-1], sd=sd_W2)

  # Store the results
  theta_01_hist[n]      <- theta_01
  theta_02_hist[n]      <- theta_02
  W1_hist[n]            <- W1
  W2_hist[n]            <- W2
  theta1_star_hist[n, ] <- theta1_star
  theta2_star_hist[n, ] <- theta2_star

  if (n %% 100 == 0) {
    elapsed <- (proc.time() - start_time)[[3]]
    printf("n=%d/%d | t=%.0fs | W1=%.6f | W2=%.6f",
           n, N, elapsed, W1, W2)
  }
}

elapsed_time <- (proc.time() - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)

# Results
theta1_mean <- colMeans(theta1_star_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_star_hist[-(1:burnin), ])
W1_mean <- mean(W1_hist[-(1:burnin)])
W2_mean <- mean(W2_hist[-(1:burnin)])
printf("W1 posterior mean: %.6f", W1_mean)
printf("W2 posterior mean: %.6f", W2_mean)

# Plot
x <- 1:Tt

# lambda ####
lambda_true <- exp(theta1_true)
lambda_mean <- exp(theta1_mean)
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(x, y, type="l", col="gray", xlab="t", ylab="",
     main="Poisson 2nd Order Polynomial Model")
points(x, y, pch=20)
lines(x, lambda_mean, col="red",  lwd=2)
lines(x, lambda_true, col="blue", lwd=2)
legend("topright",
       legend=expression(y[t], lambda[t], hat(lambda)[t]),
       col=c("black","blue","red"),
       lty=c(NA,1,1), lwd=c(NA,2,2), pch=c(20,NA,NA), bty="n")

# theta_t1 ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
ylim_range <- range(theta1_mean, theta1_true)
plot(x, theta1_mean, type="l", col="red", lwd=2, ylim=ylim_range,
     xlab="t", ylab="", main="theta_t1")
lines(x, theta1_true, col="blue", lwd=2)
legend("topright", legend=expression(hat(theta)[t1], theta[t1]),
       col=c("red","blue"), lwd=2, bty="n")

# theta_t2 ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(x, theta2_mean, type="l", col="red", lwd=2,
     xlab="t", ylab="", main="theta_t2")
legend("topright", legend=expression(hat(theta)[t2], theta[t2]),
       col=c("red","blue"), lwd=2, bty="n")

# Traceplot W1 ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(W1_hist[-(1:burnin)], type="l", xlab="n", ylab="W1",
     main="Traceplot of W1")

# Traceplot W2 ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(W2_hist[-(1:burnin)], type="l", xlab="n", ylab="W2",
     main="Traceplot of W2")

# Traceplots theta1_star ####
par(mfrow=c(2,2))
for (t in t_observed) {
  plot(theta1_star_hist[, t], type="l",
       main=bquote(theta[list(.(t),1)]), xlab="n", ylab="")
}

# Traceplots theta2_star
par(mfrow=c(2,2))
for (t in t_observed) {
  plot(theta2_star_hist[, t], type="l",
       main=bquote(theta[list(.(t),2)]), xlab="n", ylab="")
}

# Diagnóstico
library(coda)
idx <- (burnin + 1):N
t1  <- theta1_star_hist[idx, ]
t2  <- theta2_star_hist[idx, ]

rho  <- sapply(1:Tt, function(t) cor(t1[, t], t2[, t]))
ess1 <- apply(t1, 2, effectiveSize)
ess2 <- apply(t2, 2, effectiveSize)

cat(sprintf("mediana |rho_t| : %.3f\n", median(abs(rho))))
cat(sprintf("ESS theta1 (min/mediana): %.0f / %.0f\n", min(ess1), median(ess1)))
cat(sprintf("ESS theta2 (min/mediana): %.0f / %.0f\n", min(ess2), median(ess2)))

# ESS
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_smc, type="l")

upd1 <- colMeans(theta1_star_hist[(burnin+2):N, ] !=
                     theta1_star_hist[(burnin+1):(N-1), ])
cat(sprintf("upd1 (min/mediana): %.3f / %.3f | t do min: %d\n",
            min(upd1), median(upd1), which.min(upd1)))


cat(sprintf("ESS W1: %.0f | ESS W2: %.0f\n",
            +             effectiveSize(W1_hist[(burnin+1):N]),
            +             effectiveSize(W2_hist[(burnin+1):N])))

