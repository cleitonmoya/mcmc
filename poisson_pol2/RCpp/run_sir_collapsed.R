# Driver R para o amostrador SIR-Collapsed implementado em C++/Rcpp
# (sir_collapsed.cpp). sir_collapsed.cpp inclui utils.h, onde ficam as
# funções compartilhadas com os demais algoritmos (logsumexp, log_p_yt,
# log_dnorm, gibbs_sample_theta01/phi1/phi2, sample_index_from_logw, método
# de Chan tridiagonal) — os dois arquivos precisam estar no mesmo diretório.

library(Rcpp)
library(coda)

rm(list = ls())
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

sourceCpp("sir_collapsed.cpp")

printf <- function(...) cat(paste(sprintf(...), "\n"))

# Load the data
functions_grid <- c("constant", "linear", "quadratic", "sinusoidal")
f <- 3
Tt <- 1600
replica <- 1
source_name <- sprintf("%s_%s_%s", functions_grid[f], Tt, replica)
data <- readRDS(paste("../../../cobalebeb2027/data/simulated/", source_name, ".rds", sep = ""))
y <- data$y
Tt_grid <- c(200, 400, 800, 1600)
method <- 2
tau <- match(Tt, Tt_grid)
seed <- method * 1e5 + f * 1e4 + tau * 1e3 + replica
set.seed(seed)
printf("Running for %s, seed=%d", source_name, seed)
if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 400) t_obs <- c(75, 100, 200, 300)
if (Tt == 800) t_obs <- c(100, 300, 500, 700)
if (Tt == 1600) t_obs <- c(400, 800, 1200, 1600)
theta1_present <- TRUE
theta2_present <- FALSE
if (theta1_present) {
    theta1_true <- data$theta
    lambda_true <- exp(theta1_true)
}
if (theta2_present) theta2_true <- data$theta2

# ---- Hiperparâmetros ----
mu_01 <- 0
sigma2_01 <- 100
mu_02 <- 0
sigma2_02 <- 100
nu_01 <- 2
eta_01 <- 0.01
nu_02 <- 2
eta_02 <- 0.0001

# ---- Parâmetros da simulação ----
R_prerun <- 3000
N <- 10000
burnin <- 1000
M_is_lik <- 3
M_sir_theta1 <- 3
M_irls_max <- 20
tol <- 1e-4

# ---- Valores iniciais ----
theta1 <- rep(0, Tt)
theta2 <- rep(0, Tt)
theta1_tilde <- rep(0, Tt)
theta_01 <- 0
theta_02 <- 0
W1 <- 0.01
W2 <- 0.01

# ---- Execução (Gibbs em C++) ----
start_time <- proc.time()
out <- sir_collapsed_cpp(
    y = y,
    R_prerun = R_prerun,
    N = N,
    mu_01 = mu_01,
    sigma2_01 = sigma2_01,
    mu_02 = mu_02,
    sigma2_02 = sigma2_02,
    nu_01 = nu_01,
    eta_01 = eta_01,
    nu_02 = nu_02,
    eta_02 = eta_02,
    theta1 = theta1,
    theta2 = theta2,
    theta1_tilde = theta1_tilde,
    theta_01 = theta_01,
    theta_02 = theta_02,
    W1 = W1,
    W2 = W2,
    M_is_lik = M_is_lik,
    M_sir_theta1 = M_sir_theta1,
    M_irls_max = M_irls_max,
    tol = tol,
    verbose = TRUE,
    print_every = 1000
)
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[1]]
printf("Total CPU time: %.2f s", elapsed_time)

theta1_hist <- out$theta1_hist
theta2_hist <- out$theta2_hist
theta_01_hist <- out$theta_01_hist
theta_02_hist <- out$theta_02_hist
W1_hist <- out$W1_hist
W2_hist <- out$W2_hist
accepted_hist <- out$accepted_hist
accepted2_hist <- out$accepted2_hist
itr_irls <- out$itr_irls
ess_is_hist <- out$ess_is_hist
ess_sir_hist <- out$ess_sir_hist

printf("CE Gamma proposal for phi1: shape = %.4f, rate = %.4f (mean phi1 = %.6f, mean W1 = %.6f)",
       out$ce1_shape, out$ce1_rate, out$ce1_shape / out$ce1_rate, out$ce1_rate / out$ce1_shape)
printf("CE Gamma proposal for phi2: shape = %.4f, rate = %.4f (mean phi2 = %.6f, mean W2 = %.6f)",
       out$ce2_shape, out$ce2_rate, out$ce2_shape / out$ce2_rate, out$ce2_rate / out$ce2_shape)

#### Simulation summary

printf("W1 MH acceptance rate: %.3f", mean(accepted_hist))
printf("W2 MH acceptance rate: %.3f", mean(accepted2_hist))

theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)

printf("W1 mean: %.5f", mean(W1_hist[-(1:burnin)]))
printf("W1 median: %.5f", median(W1_hist[-(1:burnin)]))
printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))

loglik <- sum(dpois(y, lambda_mean, log = TRUE))
printf("Log-likelihood: %.2f", loglik)

# Effective sample size
ess_theta01 <- effectiveSize(mcmc(theta_01_hist[-(1:burnin)]))
ess_theta02 <- effectiveSize(mcmc(theta_02_hist[-(1:burnin)]))
ess_w1 <- effectiveSize(mcmc(W1_hist[-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W2_hist[-(1:burnin)]))
ess_theta1 <- effectiveSize(mcmc(theta1_hist[-(1:burnin), ]))
ess_theta2 <- effectiveSize(mcmc(theta2_hist[-(1:burnin), ]))
printf("Effective Sample Size:")
printf("\ttheta_01: %.2f", ess_theta01)
printf("\ttheta_02: %.2f", ess_theta02)
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)
printf("\ttheta1 (mean): %.2f", mean(ess_theta1))
printf("\ttheta_11 %.2f", ess_theta1[1])
printf("\ttheta2 (mean): %.2f", mean(ess_theta2))

printf("Effective Sample Size / second:")
printf("\tW1: %.2f", ess_w1 / elapsed_time)
printf("\tW2: %.2f", ess_w2 / elapsed_time)
printf("\ttheta_01: %.2f", ess_theta01 / elapsed_time)
printf("\ttheta_02: %.2f", ess_theta02 / elapsed_time)

ess_sec_theta1 <- ess_theta1 / elapsed_time
ess_sec_theta2 <- ess_theta2 / elapsed_time
printf("\ttheta1 (mean): %.2f", mean(ess_sec_theta1))
printf("\ttheta2 (mean): %.2f", mean(ess_sec_theta2))

printf("Geweke convergence diagnostic")
z_w1 <- unname(geweke.diag(W1_hist[-(1:burnin)], frac1 = 0.1, frac2 = 0.5)[[1]])
z_w2 <- unname(geweke.diag(W2_hist[-(1:burnin)], frac1 = 0.1, frac2 = 0.5)[[1]])
printf("\tz_w1: %.2f", z_w1)
printf("\tz_w2: %.2f", z_w2)

z_theta1 <- unname(geweke.diag(theta1_hist[-(1:burnin), ], frac1 = 0.1, frac2 = 0.5)[[1]])
z_theta2 <- unname(geweke.diag(theta2_hist[-(1:burnin), ], frac1 = 0.1, frac2 = 0.5)[[1]])
z1_out <- sum((z_theta1 < -1.96) | (z_theta1 > 1.96)) / Tt
z2_out <- sum((z_theta2 < -1.96) | (z_theta2 > 1.96)) / Tt
printf("\tPercent of theta1 out: %.3f", z1_out)
printf("\tPercent of theta2 out: %.3f", z2_out)


#####
# Plots
# y, lambda_true, lambda_estimated
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(x, y, type = "l", xlab = "t", ylab = "", col = "gray",
     main = "Poisson Local Trend Polynomial Model")
points(x, y, pch = 20)
lines(x, lambda_mean, col = "red", lwd = 2)
if (theta1_present) {
    lines(x, lambda_true, col = "blue", lwd = 2)
    legend("topright",
           legend = expression(y[t], lambda[t], hat(lambda)[t]),
           col = c("black", "blue", "red"),
           lty = c(NA, 1, 1),
           lwd = c(NA, 2, 2),
           pch = c(20, NA, NA),
           bty = "n")
} else {
    legend("topright",
           legend = expression(y[t], hat(lambda)[t]),
           col = c("black", "red"),
           lty = c(NA, 1),
           lwd = c(NA, 2),
           pch = c(20, NA),
           bty = "n")
}


# theta1_true, theta1_mean
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
if (theta1_present) {
    ylim_range <- range(theta1_mean, theta1_true)
} else {
    ylim_range <- range(theta1_mean)
}
plot(x, theta1_mean, type = "l", col = "red", lwd = 2, ylim = ylim_range,
     xlab = "t", ylab = "", main = "theta_t1")
if (theta1_present) {
    lines(x, theta1_true, col = "blue", lwd = 2)
    legend("topright", legend = expression(hat(theta)[t1], theta[t1]),
           col = c("red", "blue"), lwd = 2, bty = "n")
} else {
    legend("topright", legend = expression(hat(theta)[t1]),
           col = "red", lwd = 2, bty = "n")
}


# theta2_true, theta2_mean
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
if (theta2_present) {
    ylim_range <- range(theta2_mean, theta2_true)
} else {
    ylim_range <- range(theta2_mean)
}
plot(x, theta2_mean, type = "l", col = "red", lwd = 2, ylim = ylim_range,
     xlab = "t", ylab = "", main = "theta_t2")
if (theta2_present) {
    lines(x, theta2_true, col = "blue", lwd = 2)
    legend("topright", legend = expression(hat(theta)[t2], theta[t2]),
           col = c("red", "blue"), lwd = 2, bty = "n")
} else {
    legend("topright", legend = expression(hat(theta)[t2]),
           col = "red", lwd = 2, bty = "n")
}


# Posterior distribution of theta_t1
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta1_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Posterior of " * theta[.(t) * "," * 1]))
    lines(density(theta1_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of theta_t2
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta2_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 2]),
         main = bquote("Posterior of " * theta[.(t) * "," * 2]))
    lines(density(theta2_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of W1
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_hist[-(1:burnin)], breaks = 50, freq = FALSE, main = "Posterior of W1")
lines(density(W1_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_hist[-(1:burnin)], breaks = 50, freq = FALSE, main = "Posterior of W2")
lines(density(W2_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot for W1 and W2
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W1_hist[-(1:burnin)], type = "l", xlab = "n", ylab = "W", main = "Traceplot of W1")
abline(v = burnin, col = "red")
plot(W2_hist[-(1:burnin)], type = "l", xlab = "n", ylab = "W", main = "Traceplot of W2")
abline(v = burnin, col = "red")


# Traceplot for theta_01 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(theta_01_hist, type = "l", main = "Traceplot of theta01", xlab = "", ylab = "")

# Traceplot for theta_02 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(theta_02_hist, type = "l", main = "Traceplot of theta01", xlab = "", ylab = "")


# Traceplots for theta_t1
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta1_hist[, t], type = "l", main = bquote(theta[.(t) * "," * 1]), xlab = "", ylab = "")
    abline(v = burnin, col = "red")
}


# Traceplots for theta_t2
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta2_hist[, t], type = "l", main = bquote(theta[.(t) * "," * 2]), xlab = "", ylab = "")
    abline(v = burnin, col = "red")
}


# Effective sample size ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(ess_theta1, type = "l", main = expression("Effective sample of " * theta[t1]), xlab = "t")
plot(ess_theta2, type = "l", main = expression("Effective sample of " * theta[t2]), xlab = "t")


# Geweke diagnostic ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(z_theta1, type = "l", main = expression("Geweke diagnostic for " * theta[t1]),
     xlab = "t", ylab = "Z score")
abline(h = c(-1.96, 1.96), col = "red")

plot(z_theta2, type = "l", main = expression("Geweke diagnostic for " * theta[t2]),
     xlab = "t", ylab = "Z score")
abline(h = c(-1.96, 1.96), col = "red")


# ACF for theta1 e theta2 ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
for (t in t_obs) {
    acf(theta1_hist[-(1:burnin), t], main = bquote(theta[.(t) * "," * 1]))
    acf(theta2_hist[-(1:burnin), t], main = bquote(theta[.(t) * "," * 2]))
}


# Prior vs posterior for phi1 and phi2 (comparing against the CE proposals)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
curve(dgamma(x, shape = nu_01, rate = eta_01), from = 1e-6, to = max(1 / W1_hist[-(1:burnin)]),
      main = "phi1 prior vs. posterior", col = "red", lwd = 2)
lines(density(1 / W1_hist[-(1:burnin)]), col = "blue", lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = c("red", "blue"), lwd = 2)

par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
curve(dgamma(x, shape = nu_02, rate = eta_02), from = 1e-6, to = max(1 / W2_hist[-(1:burnin)]),
      main = "phi2 prior vs. posterior", col = "red", lwd = 2)
lines(density(1 / W2_hist[-(1:burnin)]), col = "blue", lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = c("red", "blue"), lwd = 2)


# Effective Sample Size - IS of W1's Integrated Likelihood
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(ess_is_hist, type = "l", main = "ESS - IS Integrated Likelihood (W1)", xlab = "n")


# Effective Sample Size - SIR of theta1
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(ess_sir_hist, type = "l", main = "ESS - SIR theta1", xlab = "n")


# W1 and W2 acceptance rate over time (rolling window)
roll_mean <- function(x, k) {
    n <- length(x)
    out <- rep(NA, n)
    for (i in k:n) out[i] <- mean(x[(i - k + 1):i])
    out
}
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(roll_mean(as.numeric(accepted_hist), 200), type = "l",
     main = "W1 rolling acceptance rate (window=200)", xlab = "n", ylab = "rate")
plot(roll_mean(as.numeric(accepted2_hist), 200), type = "l",
     main = "W2 rolling acceptance rate (window=200)", xlab = "n", ylab = "rate")
