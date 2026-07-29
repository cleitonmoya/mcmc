# Poisson - 2nd Order Polynomial Dynamic Model
# Gibbs with Collapsed Samplers for W1 AND W2
#
# theta_01/theta1: EXTENDED-state joint block ((T+1)-dimensional), sampled
#                  together via chan_smoothing_theta1/chan_sample_theta1 -
#                  still needs SIR (the Poisson likelihood on theta1[1..T]
#                  makes it only a Laplace/IRLS approximation; theta_01
#                  itself needs no approximation, no likelihood on it).
# theta_02/theta2: EXTENDED-state joint block ((T+1)-dimensional), sampled
#                  together via chan_smoothing_theta2/chan_sample_theta2 -
#                  exact, no SIR needed (no Poisson likelihood involved).
#
# Author: Cleiton Moya de Almeida

library(Matrix)
library(coda)

rm(list = ls())
options(error = function() traceback(2))
tp <- base::t
set.seed(40)

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "quadratic_2000_1"
data <- readRDS(paste("../../cobalebeb2027/data/simulated/", source, ".rds", sep=""))
y <- data$y

Tt <- length(y)
Ttp1 <- Tt + 1
if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 400) t_obs <- c(75, 100, 200, 300)
if (Tt == 800) t_obs <- c(200, 300, 500, 700)
if (Tt == 2000) t_obs <- c(500, 1000, 1500, 1750)

theta1_present <- TRUE
theta2_present <- FALSE

if (theta1_present) {
    theta1_true <- data$theta1
    lambda_true <- exp(theta1_true)
}
if (theta2_present)
    theta2_true <- data$theta2


#####
# Auxiliary functions

printf <- function(...) cat(paste(sprintf(...), "\n"))

logsumexp <- function(x) {
    cc <- max(x)
    return(cc + log(sum(exp(x - cc))))
}


#####
# Static sparse matrices for the Chan Method
#
# theta1 AND theta2 both use (T+1)-dimensional EXTENDED structures: theta_01
# (resp. theta_02) is folded in as node "0" of the chain (its own Normal
# prior), so (theta_01, theta1) and (theta_02, theta2[1:T]) are each sampled
# JOINTLY in one block. This breaks the slow-mixing centered-parameterization
# correlation between node "0" and node "1" that a separate conjugate-Normal
# + block-sampler alternation would have. For theta2 this joint block is
# EXACT (no Poisson likelihood involved). For theta1 it is still only a
# Laplace/IRLS approximation (theta1[1..T] carries the Poisson likelihood),
# so it still needs SIR/importance sampling on top - see run_irls/is_log_lik/
# sir_theta1 - theta_01 itself needs no linearization (no likelihood), only
# its exact Gaussian prior and its exact process link into theta1[1].
#
# The pure random-walk skeleton on T+1 nodes (main_diag_rw below) is
# singular on its own (a RW with no anchor has a free translation
# direction), so the initial symbolic Cholesky() call needs a
# positive-definite placeholder (main_diag_rw + 1) just to extract the
# sparsity pattern; update() is always called afterwards with the real
# (always PD, thanks to the prior + data terms on the diagonal) values.

start_time = proc.time()

# --- theta1 (T-dimensional) ---
sub_diag_base <- rep(-1, Tt-1)
main_diag_base <- c(rep(2, Tt-1), 1)
K0 <- bandSparse(n=Tt, k=c(0, -1),
                 diagonals=list(main_diag_base, sub_diag_base),
                 symmetric = TRUE)

diag_pattern <- bandSparse(n=Tt, k=c(0, -1),
                           diagonals=list(rep(TRUE, Tt), rep(FALSE, Tt-1)),
                           symmetric=TRUE)
idx_diag <- which(diag_pattern@x)

sub_pattern <- bandSparse(n=Tt, k=c(0, -1),
                          diagonals=list(rep(FALSE, Tt), rep(TRUE, Tt-1)),
                          symmetric=TRUE)
idx_sub <- which(sub_pattern@x)

Ch01_factor <- Cholesky(K0, perm = FALSE, LDL = TRUE)

# log|K0| (constant, precomputed once - used in the exact W2 marginal likelihood)
log_det_K0 <- 2 * as.numeric(determinant(Ch01_factor, logarithm = TRUE)$modulus)


# --- theta2 (EXTENDED, (T+1)-dimensional: theta_02 is node "0") ---
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

# --- theta1 (EXTENDED, (T+1)-dimensional: theta_01 is node "0") ---
# Same sparsity pattern as theta2's extended structure (both are first-order
# Gauss-Markov chains of length T+1) - reuses K0_ext_symbolic, but needs its
# OWN Cholesky object (Cholesky() called again, never assigned/shared) and
# its own working precision matrix, since P1_matrix/Ch01_factor get updated
# every Gibbs iteration with different numeric values (phi1, IRLS precisions
# - not phi2). This reassigns the T-dimensional Ch01_factor/P1_matrix above;
# log_det_K0 was already computed from the original values, so it is
# unaffected.
Ch01_factor <- Cholesky(K0_ext_symbolic, perm = FALSE, LDL = TRUE)
P1_matrix <- K0_ext_symbolic

# --- dedicated T-dimensional structure for the W2 integrated likelihood ---
# log_marginal_lik_w2 needs theta_02 held FIXED (not jointly marginalized) -
# a mathematically distinct computation from sampling (theta_02, theta2),
# so it gets its own small, clearly-scoped Cholesky factor (never touched
# by chan_smoothing_theta2/chan_sample_theta2 above).
Ch_w2lik_factor <- Cholesky(K0, perm = FALSE, LDL = TRUE)
P_w2lik_matrix <- K0

building_time <- (proc.time() - start_time)[[1]]
printf("Sparse structures building: %.4f s", building_time)


#####
# Chan Method functions

# theta1 | theta2 (fixed), phi1 - EXTENDED, (T+1)-dimensional: theta_01 is
# node "0", jointly sampled with theta1[1..T]. theta_01 needs no
# linearization itself (no likelihood) - only its exact Gaussian prior and
# the exact process link to theta1[1], both already captured by
# main_diag_rw*phi1 + extra_diag (no extra "+phi1" term at node 0, unlike
# theta_02, since theta_01 has no external channel into another block).
chan_smoothing_theta1 <- function(z_t, phi_V, phi1, theta_02, theta2) {
    extra_diag <- c(1/sigma2_01, phi_V)
    P1_matrix@x[idx_diag_e] <- (main_diag_rw * phi1) + extra_diag
    P1_matrix@x[idx_sub_e]  <- -phi1
    ch <- update(Ch01_factor, P1_matrix)

    RHS_ext <- c(theta_02, theta2[-Tt])
    Hb_ext <- numeric(Ttp1)
    Hb_ext[1] <- -RHS_ext[1]
    Hb_ext[2:Tt] <- RHS_ext[1:(Tt-1)] - RHS_ext[2:Tt]
    Hb_ext[Ttp1] <- RHS_ext[Tt]
    b <- c(mu_01/sigma2_01, phi_V*z_t) + phi1*Hb_ext

    theta1_hat <- as.numeric(Matrix::solve(ch, b, system="A"))
    list(theta1_hat=theta1_hat, ch=ch)
}


chan_sample_theta1 <- function(build_res) {
    d <- Matrix::diag(build_res$ch)
    u <- rnorm(Ttp1)
    w <- u / sqrt(d)
    x <- as.vector(Matrix::solve(build_res$ch, w, system="Lt"))
    build_res$theta1_hat + x
}


# IRLS: build Laplace approximation around theta1_tilde (T-dimensional)
run_irls <- function(theta1_tilde, theta2, theta_02, phi1,
                     y, tol, M_irls_max) {
    for (j in 1:M_irls_max) {
        f_t   <- exp(-theta1_tilde)
        phi_V <- 1 / f_t
        z_t   <- theta1_tilde + f_t * y - 1   # Poisson pseudo-observation

        res <- chan_smoothing_theta1(z_t, phi_V, phi1, theta_02, theta2)
        theta1_tilde_new <- res$theta1_hat[-1]   # drop node 0 (theta_01)
        if (any(!is.finite(theta1_tilde_new))) break

        if (max(abs(theta1_tilde_new - theta1_tilde)) < tol) {
            theta1_tilde <- theta1_tilde_new
            break
        }
        theta1_tilde <- theta1_tilde_new
    }
    return(list(res = res, theta1_tilde = theta1_tilde, itr = j,
         f_t = f_t, phi_V = phi_V))
}


# ---------------------------------------------------------------------------
# IS estimator of log p(y | phi1, theta2)  [integrated over (theta_01, theta1)
# JOINTLY]. Approximate: the Poisson likelihood makes p(theta1|y,phi1,theta2)
# non-Gaussian. theta_01 itself needs no approximation (no likelihood), but
# once the proposal q integrates it jointly (via the extended block), the
# target here must integrate it jointly too - holding it fixed would break
# the target/proposal cancellation, not just lose efficiency.
# ---------------------------------------------------------------------------
is_log_lik <- function(irls_res, y, phi1, theta2, theta_02, M_is_lik) {
    res     <- irls_res$res
    eta_hat <- res$theta1_hat
    ch      <- res$ch
    W1 <- 1/phi1

    log_det_H <- 2 * as.numeric(determinant(ch, logarithm = TRUE)$modulus)
    th_lag2_fixed <- c(theta_02, theta2[-Tt])
    log_norm_H <- -Ttp1 / 2 * log(2 * pi) + 0.5 * log_det_H

    log_w <- numeric(M_is_lik)
    u_mat <- matrix(rnorm(M_is_lik * Ttp1), nrow = M_is_lik)
    d_ch  <- Matrix::diag(ch)

    for (i in 1:M_is_lik) {
        u <- u_mat[i, ]
        w <- u / sqrt(d_ch)
        x <- as.vector(solve(ch, w, system = "Lt"))
        draw_i <- eta_hat + x
        theta_01_i <- draw_i[1]
        th         <- draw_i[-1]

        log_py <- sum(y * th - exp(th))
        log_prior_01 <- -0.5 * log(2 * pi * sigma2_01) - (theta_01_i - mu_01)^2 / (2 * sigma2_01)
        th_lag1 <- c(theta_01_i, th[-Tt])
        eps     <- th - th_lag1 - th_lag2_fixed
        log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)
        log_q <- log_norm_H - 0.5 * sum(u^2)

        log_w[i] <- log_py + log_prior_01 + log_prior_th - log_q
    }

    log_lik <- logsumexp(log_w) - log(M_is_lik)
    w_norm  <- exp(log_w - logsumexp(log_w))
    ess_is  <- 1 / sum(w_norm^2)

    return(list(log_lik = log_lik, ess_is = ess_is))
}


# SIR for (theta_01, theta1) JOINTLY (T+1-dimensional; corrects the
# Laplace/IRLS Gaussian approximation via importance resampling)
sir_theta1 <- function(irls_res, y, phi1, theta2, theta_02, M_sir_theta1) {
    eta_hat <- irls_res$res$theta1_hat
    ch      <- irls_res$res$ch
    W1 <- 1/phi1

    log_det_H     <- 2 * as.numeric(determinant(ch, logarithm = TRUE)$modulus)
    th_lag2_fixed <- c(theta_02, theta2[-Tt])
    log_norm_H    <- -Ttp1 / 2 * log(2 * pi) + 0.5 * log_det_H

    log_w  <- numeric(M_sir_theta1)
    draws  <- matrix(0, M_sir_theta1, Ttp1)
    d_ch   <- Matrix::diag(ch)

    for (i in 1:M_sir_theta1) {
        u  <- rnorm(Ttp1)
        w <- u / sqrt(d_ch)
        x  <- as.vector(Matrix::solve(ch, w, system = "Lt"))
        draw_i <- eta_hat + x
        draws[i, ] <- draw_i
        theta_01_i <- draw_i[1]
        th         <- draw_i[-1]

        log_py       <- sum(y * th - exp(th))
        log_prior_01 <- -0.5 * log(2 * pi * sigma2_01) - (theta_01_i - mu_01)^2 / (2 * sigma2_01)
        th_lag1      <- c(theta_01_i, th[-Tt])
        eps          <- th - th_lag1 - th_lag2_fixed
        log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)
        log_q        <- log_norm_H - 0.5 * sum(u^2)

        log_w[i] <- log_py + log_prior_01 + log_prior_th - log_q
    }

    w_norm  <- exp(log_w - logsumexp(log_w))
    ess_sir <- 1 / sum(w_norm^2)
    idx     <- sample.int(M_sir_theta1, size = 1, prob = w_norm)

    list(theta_01 = draws[idx, 1], theta1 = draws[idx, -1], ess = ess_sir)
}


# ---------------------------------------------------------------------------
# EXTENDED-STATE block for (theta_02, theta2[1:T]), given theta1 (fixed, all
# of it) and theta_01 (fixed). No Poisson likelihood involved here, so this
# block draw is EXACT (no SIR/IS needed).
#
# theta_02 appears in BOTH theta1's evolution (theta1[1] = theta_01+theta_02+w1)
# and theta2's evolution (theta2[1] = theta_02+w2), so its own diagonal entry
# combines both channels: 1/sigma2_02 [prior] + phi1 [from theta1[1]] + the
# usual phi2 RW contribution (handled via main_diag_rw*phi2, same as theta2[t]
# connecting to theta2[t-1]).
# ---------------------------------------------------------------------------
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


chan_sample_theta2 <- function(build_res) {
    d <- Matrix::diag(build_res$ch)
    u <- rnorm(Ttp1)
    w <- u / sqrt(d)
    x <- as.vector(Matrix::solve(build_res$ch, w, system = "Lt"))
    build_res$theta_hat + x
}


# ---------------------------------------------------------------------------
# EXACT integrated likelihood log p(z | phi2, phi1, theta_02) [integrated over
# theta2, theta_02 held FIXED]. Mathematically distinct from
# chan_smoothing_theta2 above (which marginalizes/samples theta_02 JOINTLY) -
# the W2 collapsed MH needs the likelihood conditional on the CURRENT theta_02,
# so this uses its own dedicated T-dimensional sparse structure
# (Ch_w2lik_factor / P_w2lik_matrix), never touched by the sampling functions.
# ---------------------------------------------------------------------------
log_marginal_lik_w2 <- function(theta1, phi1, phi2, theta_02) {
    z <- diff(theta1)
    diag_obs <- c(rep(phi1, Tt-1), 0)
    P_w2lik_matrix@x[idx_diag] <- (main_diag_base*phi2) + diag_obs
    P_w2lik_matrix@x[idx_sub]  <- -phi2
    ch <- update(Ch_w2lik_factor, P_w2lik_matrix)

    b <- numeric(Tt)
    b[1:(Tt-1)] <- z * phi1
    b[1] <- b[1] + theta_02 * phi2
    theta2_hat <- as.numeric(Matrix::solve(ch, b, system = "A"))

    log_pz <- -0.5*(Tt-1)*log(2*pi/phi1) - 0.5*phi1*sum((z - theta2_hat[1:(Tt-1)])^2)

    diffs2 <- theta2_hat - c(theta_02, theta2_hat[-Tt])
    log_p_theta2 <- -0.5*Tt*log(2*pi/phi2) + 0.5*log_det_K0 - 0.5*phi2*sum(diffs2^2)

    log_det_H <- 2 * as.numeric(determinant(ch, logarithm = TRUE)$modulus)
    log_q <- -0.5*Tt*log(2*pi) + 0.5*log_det_H

    list(log_lik = log_pz + log_p_theta2 - log_q)
}


calibrate_ce_gamma <- function(samples) {
    m1    <- mean(samples)
    m_log <- mean(log(samples))
    rhs   <- log(m1) - m_log

    c_hat <- if (rhs <= 0.5772) {
        (3 - rhs + sqrt((rhs - 3)^2 + 24 * rhs)) / (12 * rhs)
    } else {
        1 / rhs
    }

    for (k in 1:100) {
        f_val <- log(c_hat) - digamma(c_hat) - rhs
        fp    <- 1 / c_hat - trigamma(c_hat)
        c_new <- c_hat - f_val / fp
        if (!is.finite(c_new) || c_new <= 0) break
        if (abs(c_new - c_hat) < 1e-10) { c_hat <- c_new; break }
        c_hat <- c_new
    }
    d_hat <- c_hat / m1
    return(list(shape = c_hat, rate = d_hat))
}


mh_w1_collapsed <- function(phi1_cur, log_lik_cur, irls_cur,
                            ce_params, nu_01, eta_01,
                            theta2, theta_02,
                            theta1_tilde, y, tol, M_irls_max, M_is_lik) {

    phi1_prop <- rgamma(1, shape = ce_params$shape, rate = ce_params$rate)

    irls_prop    <- run_irls(theta1_tilde, theta2, theta_02,
                             phi1_prop, y, tol, M_irls_max)
    res_lik_prop <- is_log_lik(irls_prop, y, phi1_prop, theta2, theta_02, M_is_lik)
    log_lik_prop <- res_lik_prop$log_lik

    log_prior <- function(phi) dgamma(phi, shape = nu_01, rate = eta_01, log = TRUE)
    log_prop  <- function(phi) dgamma(phi, shape = ce_params$shape, rate = ce_params$rate, log = TRUE)

    log_alpha <- (log_lik_prop + log_prior(phi1_prop) + log_prop(phi1_cur)) -
        (log_lik_cur  + log_prior(phi1_cur)  + log_prop(phi1_prop))

    if (log(runif(1)) < log_alpha) {
        return(list(phi1 = phi1_prop, log_lik = log_lik_prop,
                    irls_res = irls_prop, accepted = TRUE))
    } else {
        return(list(phi1 = phi1_cur, log_lik = log_lik_cur,
                    irls_res = irls_cur, accepted = FALSE))
    }
}


mh_w2_collapsed <- function(phi2_cur, log_lik2_cur,
                            ce2_params, nu_02, eta_02,
                            theta1_star, phi1, theta_02) {

    phi2_prop <- rgamma(1, shape = ce2_params$shape, rate = ce2_params$rate)

    res_prop <- log_marginal_lik_w2(theta1_star, phi1, phi2_prop, theta_02)
    log_lik2_prop <- res_prop$log_lik

    log_prior <- function(phi) dgamma(phi, shape = nu_02, rate = eta_02, log = TRUE)
    log_prop  <- function(phi) dgamma(phi, shape = ce2_params$shape, rate = ce2_params$rate, log = TRUE)

    log_alpha <- (log_lik2_prop + log_prior(phi2_prop) + log_prop(phi2_cur)) -
        (log_lik2_cur  + log_prior(phi2_cur)  + log_prop(phi2_prop))

    if (!is.finite(log_alpha)) {
        return(list(phi2 = phi2_cur, log_lik2 = log_lik2_cur, accepted = FALSE))
    }

    if (log(runif(1)) < log_alpha) {
        return(list(phi2 = phi2_prop, log_lik2 = log_lik2_prop, accepted = TRUE))
    } else {
        return(list(phi2 = phi2_cur, log_lik2 = log_lik2_cur, accepted = FALSE))
    }
}


# Prior hyperparameters
mu_01     <- 0
sigma2_01 <- 100
mu_02     <- 0
sigma2_02 <- 100
nu_01 <- 2
eta_01  <- 0.01
nu_02 <- 2
eta_02  <- 0.0001


# Simulation parameters
R_prerun     <- 3000
N            <- 10000
burnin       <- 1000
M_is_lik     <- 3
M_sir_theta1 <- 3
M_irls_max   <- 20
tol          <- 1e-4


# Initial values
theta1_star  <- rep(0, Tt)
theta2       <- rep(0, Tt)
theta1_tilde <- rep(0, Tt)
theta_01 <- 0
theta_02 <- 0
W1 <- 0.01
phi1 <- 1/W1
W2 <- 0.01
phi2 <- 1/W2

W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <- numeric(N)
theta_02_hist <- numeric(N)
theta1_hist <- matrix(0, N, Tt)
theta2_hist <- matrix(0, N, Tt)
accepted_hist  <- logical(N)
accepted2_hist <- logical(N)
itr_irls <- numeric(N)
ess_is_hist <- numeric(N)
ess_sir_hist <- numeric(N)


#####
# PRE-RUN

printf("Starting pre-run (%d iterations)...", R_prerun)
phi1_prerun <- numeric(R_prerun)
phi2_prerun <- numeric(R_prerun)
time_prerun <- proc.time()

for (r in 1:R_prerun) {

    diffs1 <- theta1_star - c(theta_01, theta1_star[-Tt]) - c(theta_02, theta2[-Tt])
    nu_01_bar  <- nu_01 + Tt / 2
    eta_01_bar <- eta_01 + 0.5 * sum(diffs1^2)
    phi1 <- rgamma(1, shape = nu_01_bar, rate = eta_01_bar)
    W1 <- 1/phi1
    phi1_prerun[r] <- phi1

    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    nu_02_bar  <- nu_02 + Tt / 2
    eta_02_bar <- eta_02 + 0.5 * sum(diffs2^2)
    phi2 <- rgamma(1, shape = nu_02_bar, rate = eta_02_bar)
    W2 <- 1/phi2
    phi2_prerun[r] <- phi2

    # (theta_01, theta1) jointly (SIR, extended)
    irls_out <- run_irls(theta1_tilde, theta2, theta_02, phi1, y, tol, M_irls_max)
    theta1_tilde <- irls_out$theta1_tilde
    res_sir <- sir_theta1(irls_out, y, phi1, theta2, theta_02, M_sir_theta1)
    theta_01    <- res_sir$theta_01
    theta1_star <- res_sir$theta1

    # (theta_02, theta2) jointly via extended block
    build2 <- chan_smoothing_theta2(theta1_star, phi1, phi2, mu_02, sigma2_02, theta_01)
    draw2  <- chan_sample_theta2(build2)
    theta_02 <- draw2[1]
    theta2   <- draw2[-1]
}

time1 <- proc.time()
elapsed_prerun <- (time1 - time_prerun)[[1]]
printf("Pre-run done in %.0f s", elapsed_prerun)


#####
# CE Calibration
ce_params  <- calibrate_ce_gamma(phi1_prerun)
ce2_params <- calibrate_ce_gamma(phi2_prerun)

printf("CE Gamma proposal for phi1: shape = %.4f, rate = %.4f (mean phi1 = %.6f, mean W1 = %.6f)",
       ce_params$shape, ce_params$rate,
       ce_params$shape / ce_params$rate,
       ce_params$rate / ce_params$shape)

printf("CE Gamma proposal for phi2: shape = %.4f, rate = %.4f (mean phi2 = %.6f, mean W2 = %.6f)",
       ce2_params$shape, ce2_params$rate,
       ce2_params$shape / ce2_params$rate,
       ce2_params$rate / ce2_params$shape)


#####
# Gibbs Loop

irls_cur <- run_irls(theta1_tilde, theta2, theta_02, phi1, y, tol, M_irls_max)
theta1_tilde <- irls_cur$theta1_tilde

res_lik <- is_log_lik(irls_cur, y, phi1, theta2, theta_02, M_is_lik)
log_lik_cur <- res_lik$log_lik

res_lik2 <- log_marginal_lik_w2(theta1_star, phi1, phi2, theta_02)
log_lik2_cur <- res_lik2$log_lik

printf("Starting main MCMC (%d iterations)...", N)
start_time <- proc.time()

for (n in 1:N) {

    if (n %% 1000 == 0) {
        elapsed <- (proc.time() - start_time)[[1]]
        acc_rate  <- mean(accepted_hist[1:(n - 1)])
        acc_rate2 <- mean(accepted2_hist[1:(n - 1)])
        printf("Iter %d / %d | Elapsed: %.0f s | W1 accept: %.2f | W2 accept: %.2f",
               n, N, elapsed, acc_rate, acc_rate2)
    }

    # Collapsed MH for W1
    mh_res <- mh_w1_collapsed(
        phi1_cur     = phi1,
        log_lik_cur  = log_lik_cur,
        irls_cur     = irls_cur,
        ce_params    = ce_params,
        nu_01        = nu_01,
        eta_01       = eta_01,
        theta2       = theta2,
        theta_02     = theta_02,
        theta1_tilde = theta1_tilde,
        y            = y,
        tol          = tol,
        M_irls_max   = M_irls_max,
        M_is_lik     = M_is_lik
    )

    phi1 <- mh_res$phi1
    W1 <- 1 / phi1
    log_lik_cur <- mh_res$log_lik
    irls_cur <- mh_res$irls_res
    itr_irls[n] <- irls_cur$itr
    accepted_hist[n] <- mh_res$accepted

    # Collapsed MH for W2
    # theta1_star here is still the value from the previous iteration -
    # (theta_01, theta1)'s SIR step hasn't run yet this sweep, only phi1 is
    # fresh at this point (matches the reference ordering below).
    res_lik2_cur <- log_marginal_lik_w2(theta1_star, phi1, phi2, theta_02)
    log_lik2_cur <- res_lik2_cur$log_lik

    mh2_res <- mh_w2_collapsed(
        phi2_cur      = phi2,
        log_lik2_cur  = log_lik2_cur,
        ce2_params    = ce2_params,
        nu_02         = nu_02,
        eta_02        = eta_02,
        theta1_star   = theta1_star,
        phi1          = phi1,
        theta_02      = theta_02
    )

    phi2 <- mh2_res$phi2
    W2 <- 1 / phi2
    accepted2_hist[n] <- mh2_res$accepted

    # (theta_02, theta2) jointly via extended block given accepted phi2.
    # Sampled BEFORE (theta_01, theta1): theta1's SIR step needs the FRESH
    # theta2 to keep the phi2-theta1 coupling tight - doing it in the
    # opposite order (theta1 first) was found to degrade ESS(W2)
    # substantially (see poisson_pol2_gibbs_sir_collapsed_W1W2-separated.R,
    # ESS(W2) ~4500 vs ~350).
    build2 <- chan_smoothing_theta2(theta1_star, phi1, phi2, mu_02, sigma2_02, theta_01)
    draw2  <- chan_sample_theta2(build2)
    theta_02 <- draw2[1]
    theta2   <- draw2[-1]

    # (theta_01, theta1) jointly (SIR, extended, given accepted phi1 and the
    # fresh theta2)
    res_sir <- sir_theta1(irls_cur, y, phi1, theta2, theta_02, M_sir_theta1)
    theta_01     <- res_sir$theta_01
    theta1_star  <- res_sir$theta1
    theta1_tilde <- irls_cur$theta1_tilde
    ess_sir_hist[n] <- res_sir$ess

    # Update irls_cur and log_lik_cur for next iteration
    irls_cur <- run_irls(theta1_tilde, theta2, theta_02, phi1, y, tol, M_irls_max)
    theta1_tilde <- irls_cur$theta1_tilde

    res_lik <- is_log_lik(irls_cur, y, phi1, theta2, theta_02, M_is_lik)
    log_lik_cur <- res_lik$log_lik
    ess_is_hist[n] <- res_lik$ess_is

    # Store the results
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_hist[n, ] <- theta1_star
    theta2_hist[n, ] <- theta2
}


#####
# Simulation summary

end_time <- proc.time()
sampling_time <- (end_time - time1)[[1]]
elapsed_time <- (end_time - start_time)[[1]]
printf("Sampling: %.2f s", sampling_time)
printf("Total CPU time: %.0f s", elapsed_time)

printf("W1 MH acceptance rate: %.3f", mean(accepted_hist))
printf("W2 MH acceptance rate: %.3f", mean(accepted2_hist))

theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)

printf("W1 mean: %.5f", mean(W1_hist[-(1:burnin)]))
printf("W1 median: %.5f", median(W1_hist[-(1:burnin)]))
printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))

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


printf("Effective Sample Size / second:")
printf("\tW1: %.2f", ess_w1/elapsed_time)
printf("\tW2: %.2f", ess_w2/elapsed_time)
printf("\ttheta_01: %.2f", ess_theta01/elapsed_time)
printf("\ttheta_02: %.2f", ess_theta02/elapsed_time)
printf("\ttheta1 (mean): %.2f", mean(ess_theta1)/elapsed_time)
printf("\ttheta2 (mean): %.2f", mean(ess_theta2)/elapsed_time)

printf("Geweke convergence diagnostic")
z_w1 <- unname(geweke.diag(W1_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
z_w2 <- unname(geweke.diag(W2_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
printf("\tz_w1: %.2f", z_w1)
printf("\tz_w2: %.2f", z_w2)

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
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8)
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


# theta1_true, theta1_mean
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


# theta2_true, theta2_mean
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
hist(W1_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W1")
lines(density(W1_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W2")
lines(density(W2_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot for W1 and W2
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W1_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W1")
abline(v=burnin, col="red")
plot(W2_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W2")
abline(v=burnin, col="red")


# Traceplot for theta_01 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(theta_01_hist, type="l", main="Traceplot of theta01", xlab="", ylab="")

# Traceplot for theta_02 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(theta_02_hist, type="l", main="Traceplot of theta01", xlab="", ylab="")


# Traceplots for theta_t1
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
    abline(v=burnin, col="red")
}


# Traceplots for theta_t2
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta2_hist[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
    abline(v=burnin, col="red")
}


# Effective sample size ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(ess_theta1, type="l", main=expression("Effective sample of " * theta[t1]), xlab="t")
plot(ess_theta2, type="l", main=expression("Effective sample of " * theta[t2]), xlab="t")


# Geweke diagnostic ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(z_theta1, type="l", main=expression("Geweke diagnostic for " * theta[t1]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")

plot(z_theta2, type="l", main=expression("Geweke diagnostic for " * theta[t2]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")


# ACF for theta1 e theta2 ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8)
for (t in t_obs) {
    acf(theta1_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*1]))
    acf(theta2_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*2]))
}


# Prior vs posterior for phi1 and phi2 (comparing against the CE proposals)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8)
curve(dgamma(x, shape=nu_01, rate=eta_01), from=1e-6, to=max(1/W1_hist[-(1:burnin)]),
      main="phi1 prior vs. posterior", col="red", lwd=2)
lines(density(1/W1_hist[-(1:burnin)]), col="blue", lwd=2)
legend("topright", legend=c("Prior","Posterior"), col=c("red","blue"), lwd=2)

par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8)
curve(dgamma(x, shape=nu_02, rate=eta_02), from=1e-6, to=max(1/W2_hist[-(1:burnin)]),
      main="phi2 prior vs. posterior", col="red", lwd=2)
lines(density(1/W2_hist[-(1:burnin)]), col="blue", lwd=2)
legend("topright", legend=c("Prior","Posterior"), col=c("red","blue"), lwd=2)


# Effective Sample Size - IS of W1's Integrated Likelihood
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_is_hist, type="l", main="ESS - IS Integrated Likelihood (W1)", xlab="n")


# Effective Sample Size - SIR of theta1
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_sir_hist, type="l", main="ESS - SIR theta1", xlab="n")


# W1 and W2 acceptance rate over time (rolling window)
roll_mean <- function(x, k) {
    n <- length(x)
    out <- rep(NA, n)
    for (i in k:n) out[i] <- mean(x[(i-k+1):i])
    out
}
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(roll_mean(as.numeric(accepted_hist), 200), type="l",
     main="W1 rolling acceptance rate (window=200)", xlab="n", ylab="rate")
plot(roll_mean(as.numeric(accepted2_hist), 200), type="l",
     main="W2 rolling acceptance rate (window=200)", xlab="n", ylab="rate")
