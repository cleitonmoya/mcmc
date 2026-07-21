# Poisson - 2nd Order Polynomial Dynamic Model
# Particle Gibbs (PG) with Backward Sampling + Component-wise Gibbs for theta_t2
# Strategy: SMC only for theta_t1 (scalar state); theta_t2 sampled via
#           component-wise Gibbs (Normal conjugate full conditionals)
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269-342.
# Author: Cleiton Moya de Almeida

library(invgamma)
library(coda)

rm(list = ls())
set.seed(42)
options(error = function() traceback(2))
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source_file <- "poisson_pol2_200"
data <- readRDS(paste("../data/", source_file, ".rds", sep=""))
y <- data$y

#theta2_true <- data$theta2
Tt          <- length(y)
t_observed <- c(50, 100, 150, 175)
#t_observed  <- c(250, 500, 750, 1000)

theta1_present <- TRUE
theta2_present <- FALSE

if (theta1_present) {
    theta1_true <- data$theta
    lambda_true <- exp(theta1_true)
}

if (theta2_present) theta2_true <- data$theta2

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
sigma2_02 <- 10

# W1 ~ InvGamma(alpha_W1, beta_W1)
alpha_W1 <- 2
beta_W1  <- 0.1

# W2 ~ InvGamma(alpha_W2, beta_W2)
alpha_W2 <- 2
beta_W2  <- 0.1

N      <- 10000
K      <- 200
burnin <- 1000

W1_hist          <- numeric(N)
W2_hist          <- numeric(N)
theta_01_hist    <- numeric(N)
theta_02_hist    <- numeric(N)
theta1_hist <- matrix(0, N, Tt)
theta2_hist <- matrix(0, N, Tt)

# Initial Values
theta1 <- log(y + 0.5)
theta2 <- c(diff(theta1), 0)
theta_01 <- log(y[1] + 0.5)
theta_02 <- y[2]-y[1]
W1 <- 0.01
W2 <- 0.01

#####
start_time <- proc.time()

ess_smc <- numeric(N)

for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[3]]
        printf("Iteration %d / %d | Elapsed time: %.0f s", n, N, elapsed_time)
    }

  # ------------------------------------------------------------------
  # STEP 1: Sample theta_02 | theta_01, theta1[1], theta2[1], W1, W2
  # Termos: omega_11 = theta1[1] - theta_01 - theta_02 ~ N(0,W1)
  #         omega_12 = theta2[1] - theta_02            ~ N(0,W2)
  #         priori: theta_02 ~ N(mu_02, sigma2_02)
  # ------------------------------------------------------------------
  sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
  mu_02_bar <- sigma2_02_bar * (
    (theta1[1] - theta_01) / W1 +
    theta2[1]              / W2 +
    mu_02 / sigma2_02
  )
  theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))

  # ------------------------------------------------------------------
  # STEP 2: Sample theta_01 | theta_02, theta1[1], W1
  # Termos: omega_11 = theta1[1] - theta_01 - theta_02 ~ N(0,W1)
  #         priori: theta_01 ~ N(mu_01, sigma2_01)
  # ------------------------------------------------------------------
  sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
  mu_01_bar <- sigma2_01_bar * (
    mu_01 / sigma2_01 +
    (theta1[1] - theta_02) / W1
  )
  theta_01 <- rnorm(1, mean=mu_01_bar, sd=sqrt(sigma2_01_bar))

  # ------------------------------------------------------------------
  # STEP 3: Sample W1 | theta_01, theta_02, theta1, theta2
  # omega_t1 = theta1[t] - theta1[t-1] - theta2[t-1]
  # ------------------------------------------------------------------
  dif1   <- theta1 - c(theta_01, theta1[-Tt])
  diffs1 <- dif1 - c(theta_02, theta2[-Tt])
  alpha_W1_bar <- alpha_W1 + Tt/2
  beta_W1_bar  <- beta_W1 + 0.5 * sum(diffs1^2)
  W1 <- rinvgamma(1, shape=alpha_W1_bar, rate=beta_W1_bar)

  # ------------------------------------------------------------------
  # STEP 4: Sample W2 | theta_02, theta2
  # omega_t2 = theta2[t] - theta2[t-1]
  # ------------------------------------------------------------------
  diffs2 <- theta2 - c(theta_02, theta2[-Tt])
  alpha_W2_bar <- alpha_W2 + Tt/2
  beta_W2_bar  <- beta_W2 + 0.5 * sum(diffs2^2)
  W2 <- rinvgamma(1, shape=alpha_W2_bar, rate=beta_W2_bar)

  sd_W1 <- sqrt(W1)
  sd_W2 <- sqrt(W2)

  # ------------------------------------------------------------------
  # STEP 5: Conditional SMC para theta_t1 apenas
  # Proposta bootstrap: theta_t1^k ~ N(theta_{t-1,1}^k + theta2[t-1], W1)
  # theta2 é tratado como fixo (da iteração anterior)
  # Trajetória condicionada na coluna K
  # ------------------------------------------------------------------
  theta_1_k   <- matrix(0, Tt, K)
  log_w_tilde <- matrix(0, Tt, K)

  # t = 1
  theta_1_k[1, ] <- rnorm(K, mean=theta_01 + theta_02, sd=sd_W1)
  theta_1_k[1, K] <- theta1[1]

  log_w_1 <- log_p_yt(y[1], theta_1_k[1, ])
  log_w_tilde[1, ] <- log_w_1 - logsumexp(log_w_1) # normalizing

  # t = 2, ..., T
  for (t in 2:Tt) {

    A <- sample(1:K, K, replace=TRUE,
                prob=exp(log_w_tilde[t-1, ] - max(log_w_tilde[t-1, ])))
    A[K] <- K   # <-- preserva a trajetoria de referencia (coluna K)

    # Proposta: theta_t1^k ~ N(theta_{t-1,1}^k + theta2[t-1], W1)
    # theta2[t-1] é escalar fixo — da iteração MCMC anterior
    theta_t1_k <- rnorm(K, mean=theta_1_k[t-1, A] + theta2[t-1],
                        sd=sd_W1)
    theta_t1_k[K] <- theta1[t]
    theta_1_k[t, ] <- theta_t1_k

    log_w_t <- log_p_yt(y[t], theta_1_k[t, ])
    log_w_tilde[t, ] <- log_w_t - logsumexp(log_w_t)
  }


  ess_smc[n] <- 1 / sum(exp(2*(log_w_tilde[Tt, ] - logsumexp(log_w_tilde[Tt, ]))))

  # ------------------------------------------------------------------
  # STEP 6: Backward sampling para theta1
  # Pesos backward: w_t^k * N(theta1[t+1] | theta_t1^k + theta2[t], W1)
  # theta2[t] é escalar fixo — da iteração MCMC anterior
  # ------------------------------------------------------------------
  k_final <- sample(1:K, 1,
                    prob=exp(log_w_tilde[Tt, ] - max(log_w_tilde[Tt, ])))
  theta1[Tt] <- theta_1_k[Tt, k_final]

  for (t in (Tt-1):1) {
    log_bw <- log_w_tilde[t, ] +
              dnorm(theta1[t+1],
                    mean=theta_1_k[t, ] + theta2[t],
                    sd=sd_W1, log=TRUE)
    log_bw <- log_bw - max(log_bw)
    bw <- exp(log_bw)
    bw <- bw / sum(bw)
    b  <- sample(1:K, 1, prob=bw)
    theta1[t] <- theta_1_k[t, b]
  }

  # ------------------------------------------------------------------
  # STEP 7: Component-wise Gibbs para theta_t2 | theta1, W1, W2
  #
  # Para t = 1, ..., T-1:
  #   Termos que contêm theta_t2:
  #     omega_{t+1,1} = theta1[t+1] - theta1[t] - theta_t2  ~ N(0,W1)
  #     omega_{t+1,2} = theta2[t+1] - theta_t2                   ~ N(0,W2)
  #     omega_{t2}    = theta_t2 - theta2[t-1]                   ~ N(0,W2)
  #
  #   sigma2_t2_bar = (1/W1 + 2/W2)^(-1)
  #   mu_t2_bar = sigma2_t2_bar * (
  #     (theta1[t+1] - theta1[t]) / W1 +
  #     theta2[t+1]                    / W2 +
  #     theta2[t-1]                    / W2
  #   )
  #
  # Para t = T:
  #   Só omega_{T2} = theta_T2 - theta2[T-1] ~ N(0,W2)
  #   (theta_T2 não influencia nenhum tempo futuro)
  #
  #   sigma2_T2_bar = W2
  #   mu_T2_bar     = theta2[T-1]
  # ------------------------------------------------------------------

  # Variâncias condicionais (constantes para todos os t interiores)
  sigma2_t2_bar_interior <- (1/W1 + 2/W2)^(-1)
  sigma2_t2_bar_last     <- W2

  # t = 1: theta2[-1] = theta_02 (estado inicial)
  # Tratar separadamente pois theta2[t-1] = theta_02
  sigma2_bar <- (1/W1 + 2/W2)^(-1)
  mu_bar <- sigma2_bar * (
    (theta1[2] - theta1[1]) / W1 +
    theta2[2]                    / W2 +
    theta_02                          / W2
  )
  theta2[1] <- rnorm(1, mean=mu_bar, sd=sqrt(sigma2_bar))

  # t = 2, ..., T-1
  for (t in 2:(Tt-1)) {
    mu_bar <- sigma2_t2_bar_interior * (
      (theta1[t+1] - theta1[t]) / W1 +
      theta2[t+1]                    / W2 +
      theta2[t-1]                    / W2
    )
    theta2[t] <- rnorm(1, mean=mu_bar, sd=sqrt(sigma2_t2_bar_interior))
  }

  # t = T: só evolução própria, sem termos futuros
  theta2[Tt] <- rnorm(1, mean=theta2[Tt-1], sd=sd_W2)

  # Store the results
  theta_01_hist[n]      <- theta_01
  theta_02_hist[n]      <- theta_02
  W1_hist[n]            <- W1
  W2_hist[n]            <- W2
  theta1_hist[n, ] <- theta1
  theta2_hist[n, ] <- theta2

}

# Simulation summary ####
elapsed_time <- (proc.time() - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)

# Results
theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
W1_mean <- mean(W1_hist[-(1:burnin)])
W2_mean <- mean(W2_hist[-(1:burnin)])
printf("W1 posterior mean: %.6f", W1_mean)
printf("W2 posterior mean: %.6f", W2_mean)


# Posterior mean
theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
theta2_median <- colMedians(theta2_hist[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)
printf("W1 mean: %.3f", mean(W1_hist[-(1:burnin)]))
printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))

# Log-likelihood
loglik <- sum(dpois(y, lambda_mean, log=TRUE))
printf("Log-likelihood: %.2f", loglik)

# Effective sample size ####
printf("Effective Sample Size:")
ess_w1 <- effectiveSize(mcmc(W1_hist[-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W2_hist[-(1:burnin)]))
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)

ess_theta1 <- effectiveSize(mcmc(theta1_hist[-(1:burnin),]))
ess_theta2 <- effectiveSize(mcmc(theta2_hist[-(1:burnin),]))
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ess_theta1, type="l", main=expression("Effective sample of " * theta[t1]), xlab="t")
plot(ess_theta2, type="l", main=expression("Effective sample of " * theta[t2]), xlab="t")

par(mfrow = c(1, 2), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
hist(ess_theta1)
hist(ess_theta2)
printf("\ttheta1 (mean): %.0f", mean(ess_theta1))
printf("\ttheta2 (mean): %.0f", mean(ess_theta2))

# Effective sample size per second ####
printf("Effective Sample Size / second:")
printf("\tW1: %.2f", ess_w1/elapsed_time)
printf("\tW2: %.2f", ess_w2/elapsed_time)

ess_sec_theta1 <- ess_theta1/elapsed_time
ess_sec_theta2 <- ess_theta2/elapsed_time
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ess_sec_theta1, type="l", main=expression("Effective sample size per second of " * theta[t1]), xlab="t")
plot(ess_sec_theta2, type="l", main=expression("Effective sample size per second of " * theta[t2]), xlab="t")

par(mfrow = c(1, 2), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
hist(ess_sec_theta1)
hist(ess_sec_theta2)
printf("\ttheta1 (mean): %.0f", mean(ess_sec_theta1))
printf("\ttheta2 (mean): %.0f", mean(ess_sec_theta2))


# Geweke diagnostic: Z test for two mean difference
#   H0: segments same means -> chain has converged
printf("Geweke convergence diagnostic")
z_w1 <- unname(geweke.diag(W1_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
z_w2 <- unname(geweke.diag(W2_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
printf("\tz_w1: %.2f", z_w1)
printf("\tz_w2: %.2f", z_w2)

z_theta1 <- unname(geweke.diag(theta1_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z_theta2 <- unname(geweke.diag(theta2_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z_theta1_mean <- mean(z_theta1)
z_theta2_mean <- mean(z_theta2)
printf("\tz_theta1 (mean): %.2f", z_theta1_mean)
printf("\tz_theta2 (mean): %.2f", z_theta2_mean)

# Percent of instants in the H_0 rejection region:
z1_out <- sum((z_theta1 < -1.96) | (z_theta1 > 1.96))/Tt
z2_out <- sum((z_theta2 < -1.96) | (z_theta2 > 1.96))/Tt
printf("\tPercent of theta1 out: %.3f", z1_out)
printf("\tPercent of theta2 out: %.3f", z2_out)

par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(z_theta1, type="l", main=expression("Geweke diagnostic for " * theta[t1]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")
abline(h=z_theta1_mean, col="blue")

plot(z_theta2, type="l", main=expression("Geweke diagnostic for " * theta[t2]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")
abline(h=z_theta1_mean, col="blue")



####
# Plots

x <- 1:Tt

# lambda ####
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

# Traceplots theta1 ####
par(mfrow=c(2,2))
for (t in t_observed) {
  plot(theta1_hist[, t], type="l",
       main=bquote(theta[list(.(t),1)]), xlab="n", ylab="")
}

# Traceplots theta2
par(mfrow=c(2,2))
for (t in t_observed) {
  plot(theta2_hist[, t], type="l",
       main=bquote(theta[list(.(t),2)]), xlab="n", ylab="")
}
