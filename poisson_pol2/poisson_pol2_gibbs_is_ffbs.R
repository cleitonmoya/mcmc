# Poisson - 2nd Order Polynomial Dynamic Model
#
# Strategy: IS for theta_t1, theta_t2 sampled via FFBS
# Author: Cleiton Moya de Almeida

library(Rfast)
library(coda)


graphics.off()      # close the plots
rm(list = ls())     # clear the environment
#cat("\014")         # clear the console
options(error = function() traceback(2))
tp <- base::t       # alias to transpose function
set.seed(42)

# Change de directory to the same of the current file
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

forward_filter <- function(y, m0, C0, V, W, drift=NULL) {
    Tt <- length(y)

    if (is.null(drift)) drift <- rep(0, Tt)
    if (length(V) == 1) V <- rep(V, Tt)

    a <- numeric(Tt)
    R <- numeric(Tt)
    m <- numeric(Tt)
    C <- numeric(Tt)

    mt <- m0
    Ct <- C0

    for (t in 1:Tt) {
        at <- mt + drift[t]
        Rt <- Ct + W

        Qt <- Rt + V[t]
        At <- Rt / Qt

        mt <- at + At*(y[t] - at)
        Ct <- Rt*(1 - At)

        a[t] <- at
        R[t] <- Rt
        m[t] <- mt
        C[t] <- Ct
    }
    return(list(a=a, R=R, m=m, C=C))
}

kalman_smoother <- function(kf) {
    Tt <- length(kf$m)
    m_s <- kf$m
    B   <- kf$C[-Tt] / kf$R[-1]
    for (t in seq(Tt-1, 1)) {
        m_s[t] <- kf$m[t] + B[t] * (m_s[t+1] - kf$a[t+1])
    }
    return(list(m_s=m_s, B=B))
}

ffbs <- function(kf, ks) {
    Tt <- length(kf$m)
    theta <- numeric(Tt)
    theta[Tt] <- rnorm(1, kf$m[Tt], sqrt(kf$C[Tt]))
    B <- ks$B

    H   <- kf$C[-Tt] - B^2 * kf$R[-1]
    sH  <- sqrt(H)
    mC  <- kf$m[-Tt]
    aC1 <- kf$a[-1]

    for (t in seq(Tt-1, 1)) {
        theta[t] <- rnorm(1, mC[t] + B[t]*(theta[t+1] - aC1[t]), sH[t])
    }
    return(theta)
}


# Prior Hyperparameters

# theta_01 ~ N(mu_01, sigma2_01)
mu_01 <- 0
sigma2_01 <- 100

# theta_02 ~ N(mu_02, sigma2_02)
mu_02 <- 0
sigma2_02 <- 100

# 1/W1 ~ Gamma(nu_01, eta_01)
nu_01 <- 2
eta_01 <- 0.01

# 1/W2 ~ Gamma(nu_02, eta_02)
nu_02 <- 2
eta_02 <- 0.001


# Initial Values
theta1<- numeric(Tt)
theta2 <- numeric(Tt)
theta_01 <- 0
theta_02 <- 0
W1 <- 0.01
W2 <- 0.01


# Simulation parameters
N <- 10000
burnin <- 1000
M_is <- 3
tol <- 1e-4
M_irls_max <- 20       # maximum iterations for IRLS


# Auxiliary variables
W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <- numeric(N)
theta_02_hist <- numeric(N)
theta1_hist <- matrix(0, N, Tt)
theta2_hist <- matrix(0, N, Tt)
ess_is <- numeric(N)
itr_irls <- numeric(N) # number of iterations of IRLS (for each gibbs step)
theta1_tilde <- numeric(Tt)
Weights <- matrix(0, N, M_is)

#####
start_time <- proc.time()
for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[1]]
        printf("Iteration %d / %d | Elapsed CPU time: %.0f s", n, N, elapsed_time)
    }

    # Sample theta_01 (conjugated normal)
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar * (mu_01 / sigma2_01 +
                                      (theta1[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))

    # Sample theta_02 (conjugated normal)
    sigma2_02_bar <- (1 / sigma2_02 + 1 / W1 + 1 / W2)^(-1)
    mu_02_bar <- sigma2_02_bar * ((theta1[1] - theta_01) / W1 +
                                      theta2[1] / W2 +
                                      mu_02 / sigma2_02)
    theta_02 <- rnorm(1, mean = mu_02_bar, sd = sqrt(sigma2_02_bar))


    # Sample W1 (conjugated invgamma)
    dif1   <- theta1- c(theta_01, theta1[-Tt])
    diffs1 <- dif1 - c(theta_02, theta2[-Tt])
    nu_01_bar <- nu_01 + Tt / 2
    eta_01_bar  <- eta_01 + 0.5 * sum(diffs1^2)
    phi1 <- rgamma(1, shape = nu_01_bar, rate = eta_01_bar)
    W1 <- 1/phi1

    # Sample W2 (conjugated invgamma)
    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    nu_02_bar <- nu_02 + Tt / 2
    eta_02_bar  <- eta_02 + 0.5 * sum(diffs2^2)
    phi2 <- rgamma(1, shape = nu_02_bar, rate = eta_02_bar)
    W2 <- 1/phi2

    sd_W1 <- sqrt(W1)
    sd_W2 <- sqrt(W2)

    #
    # Importance Sampling for  theta_t1
    #

    # Local approximation (IRLS)
    theta1_tilde_old <- theta1_tilde
    for (j in 1:M_irls_max) {

            f_t <- exp(-theta1_tilde)       # observational variance
            z_t <- theta1_tilde + f_t*y - 1 # pseudo-observation
            drift <- c(theta_02, theta2[-Tt])  # pré-computado, elimina o if
            kf <- forward_filter(z_t, theta_01, 0, f_t, W1, drift=drift)
            ks <- kalman_smoother(kf)
            theta1_tilde <- ks$m_s   # update the mode
            if (max(abs(theta1_tilde - theta1_tilde_old)) < tol) break
            theta1_tilde_old <- theta1_tilde
    }
    itr_irls[n] <- j

    # Importance Sampling step
    if (M_is > 1) {
        log_w <- numeric(M_is)
        trajectories <- matrix(0, M_is, Tt)
        for (i in 1:M_is) {
            # sample theta1 proposed
            theta1_prop <- ffbs(kf, ks)
            trajectories[i, ] <- theta1_prop

            # Log-weights: log p(y|theta1) - log g(y|theta1)
            log_p <- sum(y * theta1_prop - exp(theta1_prop))
            log_g <- sum(-0.5 * log(2*pi*f_t) - 0.5 * (theta1_prop - z_t)^2 / f_t)
            log_w[i] <- log_p - log_g
        }

        log_w <- log_w - logsumexp(log_w)
        w <- exp(log_w)
        Weights[n, ] <- w
        ess_is[n] <- 1 / sum(exp(2*log_w))

        idx <- sample(1:M_is, 1, prob=w)   # index for theta1*
        theta1<- trajectories[idx, ] # theta1
    } else {
        theta1<- ffbs(kf, ks)
    }

    #
    # FFBS for  theta_t2 | theta1, W1, W2
    #
    z <- diff(theta1) # pseudo-observation
    kf2 <- forward_filter(z, theta_02, 0, W1, W2)
    ks2 <- kalman_smoother(kf2)
    theta2[1:(Tt-1)] <- ffbs(kf2, ks2)
    theta2[Tt] <- rnorm(1, theta2[Tt-1], sqrt(W2))

    # Store the results
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_hist[n, ] <- theta1
    theta2_hist[n, ] <- theta2

}
# Simulation summary ####
# Execution time
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[1]]
printf("Total elapsed CPU time: %.0f s", elapsed_time)


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
plot(W1_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W1")
plot(W2_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W2")


# Traceplots for theta_t1 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
}


# Traceplots for theta_t2 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta2_hist[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
}


# ACF for theta1 e theta2
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
for (t in t_obs) {
    acf(theta1_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*1]))
    acf(theta2_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*2]))
}

# Prior vs posterior for phi2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
curve(dgamma(x, shape=nu_02, rate=eta_02), from=0, to=max(1/W2_hist[-(1:burnin)]),
      main="phi2 prior vs. posterior", col="red", lwd=2)
lines(density(1/W2_hist[-(1:burnin)]), col="blue", lwd=2)
legend("topright", legend=c("Prior","Posterior"), col=c("red","blue"), lwd=2)


# Weights ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
i_<- min(M_is, 3)
matplot(Weights[,1:i_], type="l")
