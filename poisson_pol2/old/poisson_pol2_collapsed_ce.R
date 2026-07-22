# Poisson - 2nd Order Polynomial Dynamic Model
# Gibbs with Collapsed Sampler with Cross-Entropy for W1
# Strategy:
#   - Pre-run: standard Gibbs (R_prerun iterations) to collect W1 samples
#   - CE calibration: fit Gamma(c_hat, d_hat) to pre-run W1 samples
#   - Main loop:
#       * theta_02, theta_01 : conjugated Normal (as before)
#       * W1                 : collapsed MH with CE proposal (marginalizes over theta1)
#       * theta1*            : IRLS + Chan sampler given accepted W1
#       * 1/W2               : conjugated Gamma
#       * theta2             : precision sampler (Chan)
#
# Author: Cleiton Moya de Almeida

library(invgamma)
library(Matrix)
library(coda)

graphics.off()
rm(list = ls())
options(error = function() traceback(2))
tp <- base::t
set.seed(42)

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source_file <- "poisson_pol2_200"
data <- readRDS(paste("../data/", source_file, ".rds", sep = ""))
y <- data$y

Tt <- length(y)
if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 2000) t_obs <- c(500, 1000, 1500, 1750)

theta1_present <- TRUE
theta2_present <- FALSE

if (theta1_present) {
    theta1_true <- data$theta
    lambda_true <- exp(theta1_true)
}
if (theta2_present)
    theta2_true <- data$theta2

#
# Auxiliary functions
#

printf <- function(...) cat(paste(sprintf(...), "\n"))

logsumexp <- function(x) {
    cc <- max(x)
    return(cc + log(sum(exp(x - cc))))
}

# Chan precision sampler
diag_base <- c(rep(2, Tt - 1), 1)
off_base  <- rep(-1, Tt - 1)
K0 <- bandSparse(Tt, k = c(0, 1), diagonals = list(diag_base, off_base), symmetric = TRUE)

chan_build <- function(z, f, theta2, theta_01, theta_02, W1) {
    P  <- K0 / W1 + Diagonal(x = 1 / f)
    ch <- Cholesky(P, LDL = FALSE, perm = FALSE)
    drift <- (c(theta_02, theta2) - c(theta2, 0)) / W1
    b    <- z / f + drift[1:Tt]
    b[1] <- b[1] + (theta_01 + theta_02) / W1
    eta_hat <- as.vector(solve(ch, b))
    list(eta_hat = eta_hat, ch = ch, b = b, f = f, z = z)
}

chan_sample <- function(res) {
    u <- rnorm(Tt)
    x <- as.vector(solve(res$ch, u, system = "Lt"))
    res$eta_hat + x
}

# IRLS: build Laplace approximation around theta1_tilde
# Returns the chan_build result at convergence, plus updated theta1_tilde
run_irls <- function(theta1_tilde, theta2, theta_01, theta_02, W1,
                     y, tol, M_irls_max) {
    for (j in 1:M_irls_max) {
        f_t <- exp(-theta1_tilde)
        z_t <- theta1_tilde + f_t * y - 1
        res <- chan_build(z_t, f_t, theta2, theta_01, theta_02, W1)
        theta1_tilde_new <- res$eta_hat
        if (max(abs(theta1_tilde_new - theta1_tilde)) < tol) {
            theta1_tilde <- theta1_tilde_new
            break
        }
        theta1_tilde <- theta1_tilde_new
    }
    list(res = res, theta1_tilde = theta1_tilde, itr = j)
}

# ---------------------------------------------------------------------------
# IS estimator of log p(y | W1, theta2)  [integrated over theta1]
#
# log w^(i) = log p(y | theta1^(i))
#           + log p(theta1^(i) | W1, theta2)   <- term missing before
#           - log q(theta1^(i) | W1, theta2, y)
#
# p(theta1 | W1, theta2): Gaussian via evolution equation
#   theta_t1 = theta_{t-1,1} + theta_{t-1,2} + w_t,  w_t ~ N(0, W1)
#   innovations: eps_t = theta_t1 - theta_{t-1,1} - theta_{t-1,2}
#   log p = -T/2 * log(2*pi*W1) - 1/(2*W1) * sum(eps_t^2)
#
# q(theta1 | ...): N(eta_hat, H^{-1}) built by chan_build
#   log q = -T/2*log(2*pi) + 1/2*log|H| - 1/2*(th - eta_hat)' H (th - eta_hat)
#   Since chan_sample gives th = eta_hat + C_H^{-T} u, u~N(0,I), we have
#   (th - eta_hat) = C_H^{-T} u, so H(th-eta_hat) = C_H^T C_H C_H^{-T} u = C_H u
#   => (th-eta_hat)' H (th-eta_hat) = u'u  and  log|H| = 2*sum(log(diag(C_H)))
# ---------------------------------------------------------------------------
is_log_lik <- function(irls_res, y, W1, theta2, theta_01, theta_02, M_is) {
    res     <- irls_res$res
    f_t     <- res$f
    z_t     <- res$z
    eta_hat <- res$eta_hat
    ch      <- res$ch

    # log|H| = 2 * log|C_H|, where C_H is the Cholesky factor (C_H' C_H = H)
    # determinant() on a CHMfactor object returns log|C_H| directly
    log_det_H <- 2 * as.numeric(determinant(ch, logarithm = TRUE)$modulus)

    # Precompute prior lag vectors (fixed across IS draws)
    th_lag2_fixed <- c(theta_02, theta2[-Tt])   # theta_{t-1,2} for all t

    log_norm_H <- -Tt / 2 * log(2 * pi) + 0.5 * log_det_H

    log_w <- numeric(M_is)
    u_mat <- matrix(rnorm(M_is * Tt), nrow = M_is)       # pre-draw standard normals

    for (i in 1:M_is) {
        u  <- u_mat[i, ]
        x  <- as.vector(solve(ch, u, system = "Lt"))
        th <- eta_hat + x   # draw from q

        # log p(y | theta1)
        log_py <- sum(y * th - exp(th))

        # log p(theta1 | W1, theta2):
        # eps_t = theta_t1 - theta_{t-1,1} - theta_{t-1,2}
        # theta_{0,1} = theta_01,  theta_{0,2} = theta_02
        th_lag1 <- c(theta_01, th[-Tt])
        eps     <- th - th_lag1 - th_lag2_fixed
        log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)

        # log q(theta1): N(eta_hat, H^{-1})
        # = -T/2*log(2*pi) + 1/2*log|H| - 1/2 * u'u
        log_q <- log_norm_H - 0.5 * sum(u^2)

        log_w[i] <- log_py + log_prior_th - log_q
    }
    # Diagnóstico: verificar ESS dos pesos IS internos
    w_norm <- exp(log_w - logsumexp(log_w))
    ess_is <- 1 / sum(w_norm^2)
    return(logsumexp(log_w) - log(M_is))
}

# ---------------------------------------------------------------------------
# Cross-entropy calibration: fit Gamma(c, d) to samples via moment matching
# Solves:  log(c) - digamma(c) = log(mean(s)) - mean(log(s))
# ---------------------------------------------------------------------------
calibrate_ce_gamma <- function(samples) {
    m1    <- mean(samples)
    m_log <- mean(log(samples))
    rhs   <- log(m1) - m_log   # always > 0 by Jensen

    # Initial guess (Choi & Wette 1969 approximation)
    c_hat <- if (rhs <= 0.5772) {
        (3 - rhs + sqrt((rhs - 3)^2 + 24 * rhs)) / (12 * rhs)
    } else {
        1 / rhs
    }

    # Newton-Raphson
    for (k in 1:100) {
        f_val <- log(c_hat) - digamma(c_hat) - rhs
        fp    <- 1 / c_hat - trigamma(c_hat)
        c_new <- c_hat - f_val / fp
        if (!is.finite(c_new) || c_new <= 0) break
        if (abs(c_new - c_hat) < 1e-10) { c_hat <- c_new; break }
        c_hat <- c_new
    }
    d_hat <- c_hat / m1
    list(shape = c_hat, rate = d_hat)
}

# ---------------------------------------------------------------------------
# Collapsed MH step for W1 (independence chain with CE proposal)
# ---------------------------------------------------------------------------
mh_w1_collapsed <- function(W1_cur, log_lik_cur, irls_cur,
                             ce_params, alpha_W1, beta_W1,
                             theta2, theta_01, theta_02,
                             theta1_tilde, y, tol, M_irls_max, M_is) {

    # Propose W1* from CE Gamma proposal
    W1_prop <- rgamma(1, shape = ce_params$shape, rate = ce_params$rate)

    # Build IRLS approximation for W1_prop (warm-start from current theta1_tilde)
    irls_prop <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                          W1_prop, y, tol, M_irls_max)

    # IS estimate of log p(y | W1_prop, theta2)
    log_lik_prop <- is_log_lik(irls_prop, y, W1_prop, theta2, theta_01, theta_02, M_is)

    # Log InvGamma prior:  log p(W) = alpha*log(beta) - lgamma(alpha)
    #                                 - (alpha+1)*log(W) - beta/W
    log_prior <- function(w)
        dinvgamma(w, shape = alpha_W1, rate = beta_W1, log = TRUE)

    # Log CE Gamma proposal
    log_prop <- function(w)
        dgamma(w, shape = ce_params$shape, rate = ce_params$rate, log = TRUE)

    log_alpha <- (log_lik_prop + log_prior(W1_prop) + log_prop(W1_cur)) -
                 (log_lik_cur  + log_prior(W1_cur)  + log_prop(W1_prop))

    accepted <- log(runif(1)) < log_alpha

    if (accepted) {
        list(W1 = W1_prop, log_lik = log_lik_prop,
             irls_res = irls_prop, accepted = TRUE)
    } else {
        list(W1 = W1_cur, log_lik = log_lik_cur,
             irls_res = irls_cur, accepted = FALSE)
    }
}

# ---------------------------------------------------------------------------
# Prior hyperparameters
# ---------------------------------------------------------------------------

mu_01     <- log(y[1] + 0.5)
sigma2_01 <- 10

mu_02     <- 0
sigma2_02 <- 10

alpha_W1 <- 2
beta_W1  <- 0.1

alpha_W2 <- 2
beta_W2  <- 0.1

# ---------------------------------------------------------------------------
# Algorithm settings
# ---------------------------------------------------------------------------

R_prerun   <- 1000   # pre-run iterations to calibrate CE
N          <- 10000  # main MCMC iterations
burnin     <- 1000
M_is       <- 3     # IS draws for integrated likelihood
M_irls_max <- 20
tol        <- 1e-4

# ---------------------------------------------------------------------------
# Storage
# ---------------------------------------------------------------------------

W1_hist       <- numeric(N)
W2_hist       <- numeric(N)
theta_01_hist <- numeric(N)
theta_02_hist <- numeric(N)
theta1_hist <- matrix(0, N, Tt)
theta2_hist      <- matrix(0, N, Tt)
accepted_hist    <- logical(N)
itr_irls         <- numeric(N)

# ---------------------------------------------------------------------------
# Initial values (shared for both pre-run and main run)
# ---------------------------------------------------------------------------

theta1_star <- stats::filter(log(y + 0.5), rep(1 / 5, 5), sides = 2)
theta1_star[is.na(theta1_star)] <- log(y[is.na(theta1_star)] + 0.5)
theta1_star <- as.numeric(theta1_star)

theta2 <- c(0, diff(theta1_star))
theta_01 <- log(y[1] + 0.5)
theta_02 <- 0
W1       <- 0.01
W2       <- 0.01

theta1_tilde <- log(y + 0.5)

# ===========================================================================
# PRE-RUN: standard Gibbs to collect W1 samples for CE calibration
# ===========================================================================

printf("Starting pre-run (%d iterations)...", R_prerun)
W1_prerun <- numeric(R_prerun)
start_time <- proc.time()

for (n in 1:R_prerun) {

    # theta_02
    sigma2_02_bar <- (1 / sigma2_02 + 1 / W1 + 1 / W2)^(-1)
    mu_02_bar <- sigma2_02_bar * ((theta1_star[1] - theta_01) / W1 +
                                      theta2[1] / W2 + mu_02 / sigma2_02)
    theta_02 <- rnorm(1, mean = mu_02_bar, sd = sqrt(sigma2_02_bar))

    # theta_01
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar * (mu_01 / sigma2_01 +
                                      (theta1_star[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))

    # W1 (conjugated InvGamma)
    dif1   <- theta1_star - c(theta_01, theta1_star[-Tt])
    diffs1 <- dif1 - c(theta_02, theta2[-Tt])
    alpha_W1_bar <- alpha_W1 + Tt / 2
    beta_W1_bar  <- beta_W1 + 0.5 * sum(diffs1^2)
    W1 <- rinvgamma(1, shape = alpha_W1_bar, rate = beta_W1_bar)

    # W2 (conjugated InvGamma)
    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    alpha_W2_bar <- alpha_W2 + Tt / 2
    beta_W2_bar  <- beta_W2 + 0.5 * sum(diffs2^2)
    W2 <- rinvgamma(1, shape = alpha_W2_bar, rate = beta_W2_bar)

    # IRLS + Chan for theta1*
    irls_out <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                         W1, y, tol, M_irls_max)
    theta1_tilde <- irls_out$theta1_tilde
    theta1_star  <- chan_sample(irls_out$res)

    # CW-Gibbs for theta2
    sigma2_t2_bar_interior <- (1 / W1 + 2 / W2)^(-1)
    sd_t2 <- sqrt(sigma2_t2_bar_interior)
    sd_W2 <- sqrt(W2)

    sigma2_bar <- (1 / W1 + 2 / W2)^(-1)
    mu_bar <- sigma2_bar * ((theta1_star[2] - theta1_star[1]) / W1 +
                                theta2[2] / W2 + theta_02 / W2)
    theta2[1] <- rnorm(1, mean = mu_bar, sd = sqrt(sigma2_bar))

    d1 <- diff(theta1_star) / W1
    for (t in 2:(Tt - 1)) {
        mu_bar   <- sigma2_t2_bar_interior * (d1[t] + (theta2[t + 1] + theta2[t - 1]) / W2)
        theta2[t] <- rnorm(1, mu_bar, sd_t2)
    }
    theta2[Tt] <- rnorm(1, mean = theta2[Tt - 1], sd = sd_W2)

    W1_prerun[n] <- W1
}

elapsed_prerun <- (proc.time() - start_time)[[3]]
printf("Pre-run done in %.0f s", elapsed_prerun)

# ===========================================================================
# CE CALIBRATION
# ===========================================================================

ce_params <- calibrate_ce_gamma(W1_prerun)
printf("CE Gamma proposal: shape = %.4f, rate = %.4f (mean = %.6f)",
       ce_params$shape, ce_params$rate, ce_params$shape / ce_params$rate)

# ===========================================================================
# MAIN MCMC LOOP
# ===========================================================================

# Compute initial IRLS and log integrated likelihood
irls_cur     <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                         W1, y, tol, M_irls_max)
theta1_tilde <- irls_cur$theta1_tilde
log_lik_cur  <- is_log_lik(irls_cur, y, W1, theta2, theta_01, theta_02, M_is)

printf("Starting main MCMC (%d iterations)...", N)
start_time <- proc.time()

for (n in 1:N) {

    if (n %% 1000 == 0) {
        elapsed <- (proc.time() - start_time)[[3]]
        acc_rate <- mean(accepted_hist[1:(n - 1)])
        printf("Iter %d / %d | Elapsed: %.0f s | W1 accept rate: %.2f",
               n, N, elapsed, acc_rate)
    }

    # ------------------------------------------------------------------
    # 1. Sample theta_02 (conjugated Normal)
    # ------------------------------------------------------------------
    sigma2_02_bar <- (1 / sigma2_02 + 1 / W1 + 1 / W2)^(-1)
    mu_02_bar <- sigma2_02_bar * ((theta1_star[1] - theta_01) / W1 +
                                      theta2[1] / W2 + mu_02 / sigma2_02)
    theta_02 <- rnorm(1, mean = mu_02_bar, sd = sqrt(sigma2_02_bar))

    # ------------------------------------------------------------------
    # 2. Sample theta_01 (conjugated Normal)
    # ------------------------------------------------------------------
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar * (mu_01 / sigma2_01 +
                                      (theta1_star[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))

    # ------------------------------------------------------------------
    # 3. Collapsed MH for W1 (marginalised over theta1)
    # ------------------------------------------------------------------
    mh_res <- mh_w1_collapsed(
        W1_cur       = W1,
        log_lik_cur  = log_lik_cur,
        irls_cur     = irls_cur,
        ce_params    = ce_params,
        alpha_W1     = alpha_W1,
        beta_W1      = beta_W1,
        theta2       = theta2,
        theta_01     = theta_01,
        theta_02     = theta_02,
        theta1_tilde = theta1_tilde,
        y            = y,
        tol          = tol,
        M_irls_max   = M_irls_max,
        M_is         = M_is
    )

    W1          <- mh_res$W1
    log_lik_cur <- mh_res$log_lik
    irls_cur    <- mh_res$irls_res
    accepted_hist[n] <- mh_res$accepted

    # ------------------------------------------------------------------
    # 4. Sample W2 (conjugated InvGamma)
    # ------------------------------------------------------------------
    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    alpha_W2_bar <- alpha_W2 + Tt / 2
    beta_W2_bar  <- beta_W2 + 0.5 * sum(diffs2^2)
    W2 <- rinvgamma(1, shape = alpha_W2_bar, rate = beta_W2_bar)

    sd_W2 <- sqrt(W2)

    # ------------------------------------------------------------------
    # 5. Sample theta1* via Chan given accepted W1
    #    irls_cur already holds the IRLS result for the accepted W1
    # ------------------------------------------------------------------
    theta1_tilde <- irls_cur$theta1_tilde
    itr_irls[n]  <- irls_cur$itr
    theta1_star  <- chan_sample(irls_cur$res)

    # ------------------------------------------------------------------
    # 6. Component-wise Gibbs for theta2 | theta1*, W1, W2
    # ------------------------------------------------------------------
    sigma2_t2_bar_interior <- (1 / W1 + 2 / W2)^(-1)
    sd_t2 <- sqrt(sigma2_t2_bar_interior)

    # t = 1
    sigma2_bar <- (1 / W1 + 2 / W2)^(-1)
    mu_bar <- sigma2_bar * ((theta1_star[2] - theta1_star[1]) / W1 +
                                theta2[2] / W2 + theta_02 / W2)
    theta2[1] <- rnorm(1, mean = mu_bar, sd = sqrt(sigma2_bar))

    # t = 2, ..., T-1
    d1 <- diff(theta1_star) / W1
    for (t in 2:(Tt - 1)) {
        mu_bar   <- sigma2_t2_bar_interior * (d1[t] + (theta2[t + 1] + theta2[t - 1]) / W2)
        theta2[t] <- rnorm(1, mu_bar, sd_t2)
    }

    # t = T
    theta2[Tt] <- rnorm(1, mean = theta2[Tt - 1], sd = sd_W2)

    # ------------------------------------------------------------------
    # Store
    # ------------------------------------------------------------------
    theta_01_hist[n]      <- theta_01
    theta_02_hist[n]      <- theta_02
    W1_hist[n]            <- W1
    W2_hist[n]            <- W2
    theta1_hist[n, ] <- theta1_star
    theta2_hist[n, ]      <- theta2

    # ------------------------------------------------------------------
    # Update irls_cur and log_lik_cur for next iteration
    # (theta2, theta_01, theta_02 changed in this iteration)
    # ------------------------------------------------------------------
    irls_cur     <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                             W1, y, tol, M_irls_max)
    theta1_tilde <- irls_cur$theta1_tilde
    log_lik_cur  <- is_log_lik(irls_cur, y, W1, theta2, theta_01, theta_02, M_is)

}


# ===========================================================================
# Results
# ===========================================================================

# Simulation summary ####
elapsed_time <- (proc.time() - start_time)[[3]]
printf("Main MCMC done in %.0f s", elapsed_time)
printf("W1 MH acceptance rate: %.3f", mean(accepted_hist))

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

plot(z_theta2, type="l", main=expression("Geweke diagnostic for " * theta[t2]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")


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
