# Poisson - 2nd Order Polynomial Dynamic Model (Local Trend)
# Liu & West (2001) Combined Parameter and State Estimation Filter
#
# Model (same as poisson_pol2_gibbs_mh_cw.R):
#  y_t ~ Poisson(exp{theta_t1})
#  theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1, omega_t1 ~ N(0, W1)
#  theta_t2 =                 theta_{t-1,2} + omega_t2, omega_t2 ~ N(0, W2)
#
# Priors (t=0):
#  theta_01 | D_0 ~ N(mu_01, sigma2_01)
#  theta_02 | D_0 ~ N(mu_02, sigma2_02)
#  1/W1 ~ Gamma(shape=nu_01, rate=eta_01)
#  1/W2 ~ Gamma(shape=nu_02, rate=eta_02)
#
# Filtering: Liu & West (2001) auxiliary particle filter with kernel
#            shrinkage smoothing of fixed parameters theta = (W1, W2),
#            worked on the log-scale phi = (log W1, log W2) so that the
#            normal kernel is appropriate (Liu & West, Sec. 10.4, end).
#
# Notation map (Liu & West <-> our usual W&H notation):
#   x_t     <-> (theta_t1, theta_t2)'
#   theta   <-> (W1, W2)              [fixed, unknown]
#   D_t     <-> D_t
#   Section/eq. references (10.2.x, 10.3.x, 10.4) refer to Liu & West (2001),
#   chapter 10 of "Sequential Monte Carlo Methods in Practice".
#
# Reference: Liu, J., & West, M. (2001). Combined parameter and state
#    estimation in simulation-based filtering. In Sequential Monte Carlo
#    methods in practice (pp. 197-223). Springer.
#
# Author: Cleiton Moya de Almeida


#graphics.off()     # close the plots
rm(list = ls())     # clear the environment
#cat("\014")         # clear the console
options(error = function() traceback(2)) # more informative traceback
set.seed(42)

library(coda)        # only used for effectiveSize on stored series if needed
library(mvtnorm)      # rmvnorm() for the parameter kernel

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data (same dataset as poisson_pol2_gibbs_mh_cw.R)
source <- "poisson_pol2_200"
data <- readRDS(paste("../data/", source, ".rds", sep=""))
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
if (theta2_present) theta2_true <- data$theta2

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))

#####
# AUXILIARY FUNCTIONS

# Systematic resampling (Kitagawa 1996) - lower variance than multinomial
systematic_resample <- function(w) {
    N <- length(w)
    positions <- (runif(1) + 0:(N-1)) / N
    cumw <- cumsum(w)
    idx <- numeric(N)
    i <- 1; j <- 1
    while (i <= N) {
        if (positions[i] < cumw[j]) {
            idx[i] <- j
            i <- i + 1
        } else {
            j <- j + 1
        }
    }
    return(idx)
}

# Weighted mean/covariance of a 2-column matrix
weighted_mean_cov <- function(X, w) {
    m <- c(sum(w * X[, 1]), sum(w * X[, 2]))
    Xc <- sweep(X, 2, m)
    V <- crossprod(Xc * sqrt(w))  # sum_j w_j (x_j - m)(x_j - m)'
    return(list(m = m, V = V))
}

#####
# SIMULATION MAIN PARAMETERS

# Prior hyperparameters (identical to poisson_pol2_gibbs_mh_cw.R)
mu_01     <- 0
sigma2_01 <- 100
mu_02     <- 0
sigma2_02 <- 100
nu_01  <- 2
eta_01 <- 0.01
nu_02  <- 2
eta_02 <- 0.0001

N <- 10000          # number of particles
delta <- 0.8        # discount factor (Liu & West, Sec 10.3.3, typ. 0.95-0.99)
a <- (3*delta - 1) / (2*delta)
h2 <- 1 - a^2
ess_threshold <- N/2 # adaptive resampling threshold

#####
# INITIALIZATION (t = 0): sample from the priors

theta1 <- rnorm(N, mean = mu_01, sd = sqrt(sigma2_01))   # theta_01^(j)
theta2 <- rnorm(N, mean = mu_02, sd = sqrt(sigma2_02))   # theta_02^(j)

phi1_prec <- rgamma(N, shape = nu_01, rate = eta_01)      # 1/W1
phi2_prec <- rgamma(N, shape = nu_02, rate = eta_02)      # 1/W2
logW1 <- -log(phi1_prec)
logW2 <- -log(phi2_prec)

w <- rep(1/N, N)   # particle weights (uniform at t = 0)

# Storage for history (filtered posterior at each t)
theta1_hist  <- matrix(nrow = Tt, ncol = N)
theta2_hist  <- matrix(nrow = Tt, ncol = N)
W1_hist      <- matrix(nrow = Tt, ncol = N)
W2_hist      <- matrix(nrow = Tt, ncol = N)
w_hist       <- matrix(nrow = Tt, ncol = N)
ess_hist     <- numeric(Tt)
resampled_at <- logical(Tt)

theta1_mean <- numeric(Tt)
theta2_mean <- numeric(Tt)
W1_mean     <- numeric(Tt)
W2_mean     <- numeric(Tt)

#####
# LIU & WEST FILTER MAIN LOOP

start_time <- proc.time()
for (t in 1:Tt) {

    if (t %% 200 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[1]]
        printf("t = %d / %d | Elapsed CPU time: %.0f s", t, Tt, elapsed_time)
    }

    # --- Step 0: weighted moments of phi = (logW1, logW2)  (eq. 10.3.8) ---
    Phi <- cbind(logW1, logW2)
    mv <- weighted_mean_cov(Phi, w)
    phi_bar <- mv$m
    V <- mv$V

    # --- kernel locations (shrinkage, eq. 10.3.9) ---
    m1 <- a * logW1 + (1 - a) * phi_bar[1]
    m2 <- a * logW2 + (1 - a) * phi_bar[2]

    # --- Step 1: prior point estimates of the state (mu_{t+1}) ---
    mu1 <- theta1 + theta2   # E(theta_t1 | theta_{t-1})
    mu2 <- theta2            # E(theta_t2 | theta_{t-1})

    # --- Step 2: auxiliary weights and index sampling ---
    log_g <- log(w) + dpois(y[t], lambda = exp(mu1), log = TRUE)
    g <- exp(log_g - max(log_g))
    g <- g / sum(g)
    k <- systematic_resample(g)

    theta1_k <- theta1[k]; theta2_k <- theta2[k]
    m1_k <- m1[k]; m2_k <- m2[k]
    mu1_k <- mu1[k]

    # --- Step 3: sample new parameter vector phi_{t}^(k) ~ N(m^(k), h^2 V) ---
    M <- cbind(m1_k, m2_k)
    Phi_new <- M + rmvnorm(N, mean = c(0, 0), sigma = h2 * V)
    logW1_new <- Phi_new[, 1]
    logW2_new <- Phi_new[, 2]
    W1_new <- exp(logW1_new)
    W2_new <- exp(logW2_new)

    # --- Step 4: sample new state from the system equation ---
    theta1_new <- rnorm(N, mean = theta1_k + theta2_k, sd = sqrt(W1_new))
    theta2_new <- rnorm(N, mean = theta2_k,            sd = sqrt(W2_new))

    # --- Step 5: importance weights ---
    log_w_new <- dpois(y[t], lambda = exp(theta1_new), log = TRUE) -
                 dpois(y[t], lambda = exp(mu1_k),       log = TRUE)
    w_new <- exp(log_w_new - max(log_w_new))
    w_new <- w_new / sum(w_new)

    # --- adaptive resampling ---
    ess <- 1 / sum(w_new^2)
    if (ess < ess_threshold) {
        idx <- systematic_resample(w_new)
        theta1_new <- theta1_new[idx]; theta2_new <- theta2_new[idx]
        W1_new <- W1_new[idx]; W2_new <- W2_new[idx]
        logW1_new <- logW1_new[idx]; logW2_new <- logW2_new[idx]
        w_new <- rep(1/N, N)
        resampled_at[t] <- TRUE
    } else {
        resampled_at[t] <- FALSE
    }

    # update current particle set
    theta1 <- theta1_new; theta2 <- theta2_new
    logW1 <- logW1_new; logW2 <- logW2_new
    w <- w_new

    # store
    theta1_hist[t, ] <- theta1
    theta2_hist[t, ] <- theta2
    W1_hist[t, ] <- W1_new
    W2_hist[t, ] <- W2_new
    w_hist[t, ] <- w
    ess_hist[t] <- ess

    theta1_mean[t] <- sum(w * theta1)
    theta2_mean[t] <- sum(w * theta2)
    W1_mean[t] <- sum(w * W1_new)
    W2_mean[t] <- sum(w * W2_new)
}

end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[1]]
printf("Total elapsed CPU time: %.0f s", elapsed_time)

#####
# SUMMARY

lambda_mean <- exp(theta1_mean)
loglik <- sum(dpois(y, lambda_mean, log = TRUE))
printf("Log-likelihood (filtered plug-in): %.2f", loglik)

printf("W1 filtered mean (t=T): %.5f", W1_mean[Tt])
printf("W2 filtered mean (t=T): %.5f", W2_mean[Tt])

printf("Mean particle ESS: %.0f / %d", mean(ess_hist), N)
printf("Percent of steps resampled: %.1f%%", 100 * mean(resampled_at))


#####
# Plots

# y, lambda_true, lambda_filtered ####
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(x, y, type = "l", xlab = "t", ylab = "", col = "gray",
     main = "Poisson Local Trend Model - Liu & West Filter")
points(x, y, pch = 20)
lines(x, lambda_mean, col = "red", lwd = 2)
if (theta1_present) {
    lines(x, lambda_true, col = "blue", lwd = 2)
    legend("topright",
           legend = expression(y[t], lambda[t], hat(lambda)[t]),
           col = c("black", "blue", "red"),
           lty = c(NA, 1, 1), lwd = c(NA, 2, 2), pch = c(20, NA, NA),
           bty = "n")
} else {
    legend("topright",
           legend = expression(y[t], hat(lambda)[t]),
           col = c("black", "red"),
           lty = c(NA, 1), lwd = c(NA, 2), pch = c(20, NA),
           bty = "n")
}

# theta_t1: true vs filtered mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
ylim_range <- if (theta1_present) range(theta1_mean, theta1_true) else range(theta1_mean)
plot(x, theta1_mean, type = "l", col = "red", lwd = 2, ylim = ylim_range,
     xlab = "t", ylab = "", main = "theta_t1 (filtered)")
if (theta1_present) {
    lines(x, theta1_true, col = "blue", lwd = 2)
    legend("topright", legend = expression(hat(theta)[t1], theta[t1]),
           col = c("red", "blue"), lwd = 2, bty = "n")
} else {
    legend("topright", legend = expression(hat(theta)[t1]), col = "red", lwd = 2, bty = "n")
}

# theta_t2: true vs filtered mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
ylim_range <- if (theta2_present) range(theta2_mean, theta2_true) else range(theta2_mean)
plot(20:200, theta2_mean[20:200], type = "l", col = "red", lwd = 2,
     xlab = "t", ylab = "", main = "theta_t2 (filtered)")
if (theta2_present) {
    lines(x, theta2_true, col = "blue", lwd = 2)
    legend("topright", legend = expression(hat(theta)[t2], theta[t2]),
           col = c("red", "blue"), lwd = 2, bty = "n")
} else {
    legend("topright", legend = expression(hat(theta)[t2]), col = "red", lwd = 2, bty = "n")
}

# Filtered posterior distribution of theta_t1 at selected times ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta1_hist[t, ], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Filtered posterior of " * theta[.(t) * "," * 1]))
    lines(density(theta1_hist[t, ]), col = "blue", lwd = 2)
}

# Filtered posterior distribution of theta_t2 at selected times ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta2_hist[t, ], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 2]),
         main = bquote("Filtered posterior of " * theta[.(t) * "," * 2]))
    lines(density(theta2_hist[t, ]), col = "blue", lwd = 2)
}

# Filtered posterior of W1, W2 (final time T) ####
par(mfrow = c(1, 2), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_hist[Tt, ], breaks = 50, freq = FALSE, main = "Filtered posterior of W1 (t=T)")
lines(density(W1_hist[Tt, ]), col = "blue", lwd = 2)
hist(W2_hist[Tt, ], breaks = 50, freq = FALSE, main = "Filtered posterior of W2 (t=T)")
lines(density(W2_hist[Tt, ]), col = "blue", lwd = 2)

# Time trajectory of filtered W1, W2 quantiles ####
qtl <- function(H) t(apply(H, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
W1_q <- qtl(W1_hist); W2_q <- qtl(W2_hist)

par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
matplot(x, W1_q, type = "l", lty = c(2, 2, 1, 2, 2), col = "red",
        xlab = "t", ylab = "W1", main = "Filtered quantiles of W1 over time")
matplot(x, W2_q, type = "l", lty = c(2, 2, 1, 2, 2), col = "red",
        xlab = "t", ylab = "W2", main = "Filtered quantiles of W2 over time")

# Particle ESS trajectory ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(ess_hist, type = "l", xlab = "t", ylab = "ESS",
     main = "Particle Effective Sample Size (pre-resampling)")
abline(h = ess_threshold, col = "red", lty = 2)
points(which(resampled_at), ess_hist[resampled_at], pch = 20, col = "blue")
legend("bottomright", legend = c("ESS", "threshold N/2", "resampled"),
       col = c("black", "red", "blue"), lty = c(1, 2, NA), pch = c(NA, NA, 20), bty = "n")

# Weight distribution at final time ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(w_hist[Tt, ], breaks = 50, main = "Distribution of weights at t = T",
     xlab = "weight")
