# Poisson - 2nd Order Polynomial Dynamic Model
# Particle Gibbs with Ancestral Sampling (PG-AS) + Component-wise Gibbs for theta_t2
# Strategy: - Bootstrap Particle Filter with Ancestral Sampling for theta_t1
#             (scalar state), reconstructed by simple Backward Tracking;
#           - theta2 sampled via Preicison Matrix (Chan Method)
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269-342.
#
# PERFORMANCE NOTE: the resampling of the K-1 non-reference particles now
# uses SYSTEMATIC resampling (Kitagawa, 1996; Doucet, de Freitas & Gordon,
# 2001) instead of multinomial resampling: a single runif() draw generates
# all K-1 indices via one monotone sweep through the cumulative weights
# (findInterval), instead of K-1 independent draws inside sample(). This is
# a standard, unbiased SMC resampling scheme with LOWER variance than
# multinomial resampling (not merely a speed shortcut). The single-index
# Ancestral Sampling draw (a_ref) is unchanged -- it was already a single
# draw, so there is nothing to gain there. Because systematic resampling
# consumes the RNG stream differently, results are NOT bit-identical to the
# multinomial version for the same seed -- validated by posterior mean/ESS
# equivalence instead.
# Author: Cleiton Moya de Almeida

library(Matrix)
library(coda)

#graphics.off()     # close the plots
#cat("\014")        # clear the console
rm(list = ls())     # clear the environment

options(error = function() traceback(2))  # more informative traceback

# Change the directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Auxiliary functions
printf <- function(...) cat(paste(sprintf(...), "\n"))

# Load the data
functions_grid <- c("constant","linear", "quadratic", "sinusoidal")

f <- 3
Tt <- 1600
replica <- 1

source <- sprintf("%s_%s_%s", functions_grid[f], Tt, replica)
data <- readRDS(paste("../../cobalebeb2027/data/simulated/", source, ".rds", sep=""))
y <- data$y

#
# Set the seed according to the task_grid
#
Tt_grid <- c(200, 400, 800, 1600)
method <- 2    # grid: (montoril, pg_apf, sir_laplace, sir_collapsed, stan)

tau <- match(Tt, Tt_grid)
seed = method*1e5 + f*1e4 + tau*1e3 + replica

set.seed(seed)

printf("Executing %s, seed=%d" , source, seed)

if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 400) t_obs <- c(75, 100, 200, 300)
if (Tt == 800) t_obs <- c(200, 300, 500, 700)
if (Tt == 1600) t_obs <- c(400, 800, 1200, 1500)

theta1_present <- TRUE
theta2_present <- FALSE


if (theta1_present) {
    theta1_true <- data$theta
    lambda_true <- exp(theta1_true)
}

if (theta2_present)
    theta2_true <- data$theta2

logsumexp <- function(x) {
  cc <- max(x)
  return(cc + log(sum(exp(x - cc))))
}

# Systematic resampling (Kitagawa, 1996): draws `ndraws` indices in 1:K from
# (unnormalized) log-weights using a SINGLE uniform draw, instead of the
# `ndraws` independent draws that sample(..., prob=) uses internally. A
# single u0 ~ U(0, total/ndraws) determines `ndraws` equally-spaced points
# u_j = u0 + total*j/ndraws, located via one monotone sweep through the
# cumulative weights (findInterval). Unbiased, lower resampling variance
# than multinomial (Doucet, de Freitas & Gordon, 2001).
systematic_resample <- function(logw, ndraws) {
  K <- length(logw)
  m <- max(logw)
  w <- exp(logw - m)
  cum_w <- cumsum(w)
  total <- cum_w[K]
  u0 <- runif(1, 0, total / ndraws)
  u <- u0 + total * (0:(ndraws - 1)) / ndraws
  findInterval(u, cum_w, left.open = TRUE) + 1
}

# Single-index draw from (unnormalized) log-weights via inverse-CDF, avoiding
# the fixed per-call overhead of sample() (S3 dispatch, argument validation,
# alias-table setup) -- R's sample() is disproportionately costly for a
# single draw, called here T times per Gibbs iteration (backward
# simulation / final-index / ancestral-sampling draws).
sample_one_from_logw <- function(logw) {
  m <- max(logw)
  w <- exp(logw - m)
  cum_w <- cumsum(w)
  findInterval(runif(1, 0, cum_w[length(cum_w)]), cum_w, left.open = TRUE) + 1
}

# Log-likelihood
log_p_yt <- function(yt, theta_t1) {
  res <- yt * theta_t1 - exp(theta_t1)
  res[!is.finite(res)] <- -Inf
  return(res)
}

# Prior Hyperparameters

# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- 0
sigma2_01 <- 100

# theta_02 ~ N(mu_02, sigma2_02)
mu_02     <- 0
sigma2_02 <- 100

# phi1 = 1/W1 ~ Gamma(nu_01, eta_01)
nu_01 <- 2
eta_01  <- 0.1

# phi2 = 1/W2 ~ Gamma(nu_02, eta_02)
nu_02 <- 2
eta_02  <- 0.01


# Simulation parameters
N <- 3000      # number of Gibbs iterations
K <- 30         # number of particles
burnin <- 1000

# Initial Values
theta1 <- log(y + 0.5)
theta2 <- c(diff(theta1), 0)
theta_01 <- theta1[1]
theta_02 <- theta2[1]
W1 <- var(diff(theta1))
W2 <- var(diff(theta2))


# Auxiliary variables
W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <- numeric(N)
theta_02_hist <- numeric(N)
theta1_hist <- matrix(0, N, Tt)
theta2_hist <- matrix(0, N, Tt)
ess_smc <- numeric(N)


#####
# FIXED SPARSE STRUCTURES FOR CHAN METHOD ####
start_time = proc.time() # execution time

# theta2 (EXTENDED, (T+1)-dimensional: theta_02 is node "0")
Ttp1 <- Tt + 1
sub_diag_ext <- rep(-1, Ttp1 - 1)
main_diag_rw <- c(1, rep(2, Ttp1 - 2), 1)   # pure RW skeleton, grade 1 at extremes

K0_ext_symbolic <- bandSparse(n = Ttp1, k = c(0, -1),
                              diagonals = list(main_diag_rw + 1, sub_diag_ext),
                              symmetric = TRUE)

diag_pattern_e <- bandSparse(n = Ttp1, k = c(0, -1),
                             diagonals = list(rep(TRUE, Ttp1), rep(FALSE, Ttp1-1)),
                             symmetric = TRUE)
idx_diag_e <- which(diag_pattern_e@x)

sub_pattern_e <- bandSparse(n = Ttp1, k = c(0, -1),
                            diagonals = list(rep(FALSE, Ttp1), rep(TRUE, Ttp1-1)),
                            symmetric = TRUE)
idx_sub_e <- which(sub_pattern_e@x)

Ch02_factor <- Cholesky(K0_ext_symbolic, perm = FALSE, LDL = TRUE)
P2_matrix <- K0_ext_symbolic

time1 <- proc.time()
building_time <- (time1 - start_time)[[1]]
printf("Sparse structures building: %.4f s", building_time)

#####
# CHAN METHOD
chan_smoothing_theta2 <- function(theta1, phi1, phi2, mu_02, sigma2_02, theta_01) {
    z <- diff(theta1)   # z_t = theta1[t+1]-theta1[t], t=1,...,T-1

    extra_diag <- c(1/sigma2_02 + phi1, rep(phi1, Tt-1), 0)
    P2_matrix@x[idx_diag_e] <- (main_diag_rw * phi2) + extra_diag
    P2_matrix@x[idx_sub_e]  <- -phi2
    ch <- update(Ch02_factor, P2_matrix)

    b <- numeric(Ttp1)
    b[1] <- mu_02 / sigma2_02 + phi1 * (theta1[1] - theta_01)
    b[2:Tt] <- z * phi1
    b[Ttp1] <- 0

    list(theta_hat = as.numeric(Matrix::solve(ch, b, system = "A")), ch = ch, z = z)
}


chan_sample_from_build <- function(build, Tt) {

    ch <- build$ch
    theta_hat <- build[[1]] # theta1_hat or theta2_hat, always the first element

    d <- Matrix::diag(ch)
    u <- rnorm(Tt)
    w <- u / sqrt(d)
    x <- as.vector(Matrix::solve(ch, w, system="Lt"))

    return(theta_hat + x)
}

# Gibbs sampling
start_time = proc.time() # execution time
for (n in 1:N) {

    if (n %% 100 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[1]]
        printf("Iteration %d / %d | Elapsed CPU time: %.0f s", n, N, elapsed_time)
    }

    # Sample W1
    dif1   <- theta1 - c(theta_01, theta1[-Tt])
    diffs1 <- dif1 - c(theta_02, theta2[-Tt])
    nu_01_bar <- nu_01 + Tt / 2
    eta_01_bar  <- eta_01 + 0.5 * sum(diffs1^2)
    phi1 <- rgamma(1, shape = nu_01_bar, rate = eta_01_bar)
    W1 <- 1 / phi1
    sd_W1 <- sqrt(W1)

    # Sample W2
    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    nu_02_bar <- nu_02 + Tt / 2
    eta_02_bar  <- eta_02 + 0.5 * sum(diffs2^2)
    phi2 <- rgamma(1, shape = nu_02_bar, rate = eta_02_bar)
    W2 <- 1 / phi2
    sd_W2 <- sqrt(W2)


    # Sample theta_01 (conjugated Normal)
    # theta_01 is NOT part of the particle system: it has no likelihood term
    # (no y_0), so its exact full conditional is available in closed form,
    # exactly as in the bootstrap sampler and mathematically equivalent to
    # the informed precision (1/sigma2_01 + phi1) used for node 0 in the
    # SIR-Laplace/Collapsed extended block. Sampling it as a K=50-particle
    # cloud from the raw N(mu_01, sigma2_01) prior (previous version) is
    # uninformed by phi1 and collapses the first-stage APF ESS to ~4/50.
    sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
    mu_01_bar <- sigma2_01_bar * (mu_01/sigma2_01 + (theta1[1] - theta_02)/W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))

    #
    # Conditional SMC for theta_t1 (Bootstrap PF + Ancestral Sampling)
    #
    theta_1_k   <- matrix(0, Tt, K)
    log_w_tilde <- matrix(0, Tt, K)
    A_hist      <- matrix(0, Tt, K)   # ancestor indices, for Backward Tracking

    # t = 1
    # theta_01 is now a fixed scalar (no particle diversity at t = 0), so
    # there is no first-stage (auxiliary) resampling to do here -- exactly
    # as in the bootstrap sampler's t = 1 step.
    theta_1_k[1, ] <- rnorm(K, mean = theta_01 + theta_02, sd = sd_W1)
    theta_1_k[1, K] <- theta1[1]

    # Updated weights
    log_w_1 <- log_p_yt(y[1], theta_1_k[1, ])
    log_w_tilde[1, ] <- log_w_1 - logsumexp(log_w_1) # normalizing

    # t = 2, ..., T
    for (t in 2:Tt) {
        # Resampling + propagation for the K - 1 non-reference particles
        A <- systematic_resample(log_w_tilde[t - 1, ], K - 1)
        theta_t1_k <- numeric(K)
        theta_t1_k[1:(K - 1)] <- rnorm(K - 1,
                                        mean = theta_1_k[t - 1, A] + theta2[t - 1],
                                        sd = sd_W1)

        # Reference particle: fixed value, ancestor drawn via Ancestral Sampling
        theta_t1_k[K] <- theta1[t]
        log_as <- log_w_tilde[t - 1, ] +
            dnorm(theta1[t],
                mean = theta_1_k[t - 1, ] + theta2[t - 1],
                sd = sd_W1,
                log = TRUE
            )
        log_as <- log_as - max(log_as)
        a_ref  <- sample_one_from_logw(log_as)

        A_hist[t, ] <- c(A, a_ref)
        theta_1_k[t, ] <- theta_t1_k

        # Updated weights (bootstrap: likelihood only)
        log_w_t <- log_p_yt(y[t], theta_1_k[t, ])
        log_w_tilde[t, ] <- log_w_t - logsumexp(log_w_t)
    }

    ess_smc[n] <- 1 / sum(exp(2 * (log_w_tilde[Tt, ])))


    # backward tracking for theta1
    # PG-AS: the genealogy already encodes the backward path, so we only
    # need to look up the stored ancestor indices -- no Backward Simulation.
    k_final <- sample_one_from_logw(log_w_tilde[Tt, ])
    theta1[Tt] <- theta_1_k[Tt, k_final]

    b <- k_final
    for (t in (Tt - 1):1) {
        b <- A_hist[t + 1, b]
        theta1[t] <- theta_1_k[t, b]
    }

    # theta_01 was already sampled exactly (conjugated Normal, above) --
    # no backward sampling needed for it.

    # (theta_02, theta2) jointly via extended block
    build2 <- chan_smoothing_theta2(theta1, phi1, phi2, mu_02, sigma2_02, theta_01)
    draw2 <- chan_sample_from_build(build2, Ttp1)
    theta_02 <- draw2[1]
    theta2 <- draw2[-1]

    # Store the results
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_hist[n, ] <- theta1
    theta2_hist[n, ] <- theta2
}

#### Simulation summary

# Execution time
end_time <- proc.time()
sampling_time <- (end_time - time1)[[1]]
elapsed_time <- (end_time - start_time)[[1]]
printf("Sampling: %.2f s", sampling_time)
printf("Total CPU time: %.0f s", elapsed_time)

# Posterior mean
theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)

printf("W1 mean: %.5f", mean(W1_hist[-(1:burnin)]))
printf("W1 median: %.5f", median(W1_hist[-(1:burnin)]))
printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))


# Log-likelihood
loglik <- sum(dpois(y, lambda_mean, log=TRUE))
printf("Log-likelihood: %.2f", loglik)


# Effective sample size
ess_theta01 <- effectiveSize(mcmc(theta_01_hist[-(1:burnin)]))
ess_theta02 <- effectiveSize(mcmc(theta_02_hist[-(1:burnin)]))
ess_w1 <- effectiveSize(mcmc(W1_hist[-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W2_hist[-(1:burnin)]))
ess_theta1 <- effectiveSize(mcmc(theta1_hist[-(1:burnin),]))
ess_theta2 <- effectiveSize(mcmc(theta2_hist[-(1:burnin),]))
printf("Effective Sample Size:")
printf("\ttheta_01: %.2f", ess_theta01)
printf("\ttheta_02: %.2f", ess_theta02)
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)
printf("\ttheta1 (mean): %.2f", mean(ess_theta1))
printf("\ttheta_11 %.2f", ess_theta1[1])
printf("\ttheta2 (mean): %.2f", mean(ess_theta2))

# Effective sample size per second
printf("Effective Sample Size / second:")
printf("\tW1: %.2f", ess_w1/elapsed_time)
printf("\tW2: %.2f", ess_w2/elapsed_time)

ess_sec_theta1 <- ess_theta1/elapsed_time
ess_sec_theta2 <- ess_theta2/elapsed_time
printf("\ttheta1 (mean): %.2f", mean(ess_sec_theta1))
printf("\ttheta2 (mean): %.2f", mean(ess_sec_theta2))


# Geweke diagnostic: Z test for two mean difference
#   H0: segments same means -> chain has converged
printf("Geweke convergence diagnostic")
z_w1 <- unname(geweke.diag(W1_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
z_w2 <- unname(geweke.diag(W2_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
printf("\tz_w1: %.2f", z_w1)
printf("\tz_w2: %.2f", z_w2)

# Percent of instants in the H_0 rejection region:
z_theta1 <- unname(geweke.diag(theta1_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z_theta2 <- unname(geweke.diag(theta2_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z1_out <- sum((z_theta1 < -1.96) | (z_theta1 > 1.96))/Tt
z2_out <- sum((z_theta2 < -1.96) | (z_theta2 > 1.96))/Tt
printf("\tPercent of theta1 out: %.3f", z1_out)
printf("\tPercent of theta2 out: %.3f", z2_out)


#####
# Plots
# y, lambda_true, lambda_estimated
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson Local Trend Polynomial Model")
points(x, y, pch = 20)
lines(x, lambda_mean, col="red", lwd=2)
if (theta1_present) {
    lines(x, lambda_true, col="blue", lwd=2)
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


# y, theta1_true, theta1_mean ####
x <- 1:Tt
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
if (theta1_present) {
    ylim_range <- range(theta1_mean, theta1_true)
} else {
    ylim_range <- range(theta1_mean)
}
plot(x, theta1_mean, type="l", col="red", lwd=2, ylim=ylim_range,
     xlab="t", ylab="", main="theta_t1")
if (theta1_present) {
    lines(x, theta1_true, col="blue", lwd=2)
    legend("topright", legend=expression(hat(theta)[t1], theta[t1]),
           col=c("red","blue"), lwd=2, bty="n")
} else {
    legend("topright", legend=expression(hat(theta)[t1]),
           col="red", lwd=2, bty="n")
}


# y, theta2_true, theta2_mean ####
#y_range <- range(theta2_true, theta2_mean)
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
if (theta2_present) {
    ylim_range <- range(theta2_mean, theta2_true)
} else {
    ylim_range <- range(theta2_mean)
}
plot(x, theta2_mean, type="l", col="red", lwd=2, ylim=ylim_range,
     xlab="t", ylab="", main="theta_t2")
if (theta2_present) {
    lines(x, theta2_true, col="blue", lwd=2)
    legend("topright", legend=expression(hat(theta)[t2], theta[t2]),
           col=c("red","blue"), lwd=2, bty="n")
} else {
    legend("topright", legend=expression(hat(theta)[t2]),
           col="red", lwd=2, bty="n")
}


# Posterior distribution of theta_t1 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta1_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Posterior of " * theta[.(t) * "," * 1]))
    lines(density(theta1_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of theta_t2 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta2_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 2]),
         main = bquote("Posterior of " * theta[.(t) * "," * 2]))
    lines(density(theta2_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of W1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W1")
lines(density(W1_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W2")
lines(density(W2_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot for W1 and W2 ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W1_hist, type="l", xlab="n", ylab="W", main="Traceplot of W1")
abline(v=burnin, col="red")
plot(W2_hist, type="l", xlab="n", ylab="W", main="Traceplot of W2")
abline(v=burnin, col="red")


# Traceplots for theta_01 and theta_02 ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(theta_01_hist, type="l", xlab="n", ylab="W", main="Traceplot of theta_01")
abline(v=burnin, col="red")
plot(theta_02_hist, type="l", xlab="n", ylab="W", main="Traceplot of theta_02")
abline(v=burnin, col="red")

# Traceplots for theta_t1 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
    abline(v=burnin, col="red")
}


# Traceplots for theta_t2 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta2_hist[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
    abline(v=burnin, col="red")
}


# Effective sample size ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ess_theta1, type="l", main=expression("Effective sample of " * theta[t1]), xlab="t")
plot(ess_theta2, type="l", main=expression("Effective sample of " * theta[t2]), xlab="t")

par(mfrow = c(1, 2), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
hist(ess_theta1)
hist(ess_theta2)


# Geweke diagnostic ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(z_theta1, type="l", main=expression("Geweke diagnostic for " * theta[t1]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")

plot(z_theta2, type="l", main=expression("Geweke diagnostic for " * theta[t2]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")


# ACF for theta1 e theta2
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
for (t in t_obs) {
    acf(theta1_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*1]))
    acf(theta2_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*2]))
}


# Prior vs posterior for phi2
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
curve(dgamma(x, shape=nu_02, rate=eta_02), from=0, to=max(1/W2_hist[-(1:burnin)]),
      main="phi2 prior vs. posterior", col="red", lwd=2)
lines(density(1/W2_hist[-(1:burnin)]), col="blue", lwd=2)
legend("topright", legend=c("Prior","Posterior"), col=c("red","blue"), lwd=2)


# Effective sample size ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_smc, type="l", main="Effective Sample Size - SMC")
abline(v=burnin, col="red")
