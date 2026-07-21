# Poisson - 2nd Order Polynomial Dynamic Model
# Gibbs (PG) Importance Sampling + Component-wise Gibbs for theta_t2
# Strategy: IS for theta_t1 using the Chan Method, theta_t2 sampled via
#           component-wise Gibbs (Normal conjugate full conditionals)
# Author: Cleiton Moya de Almeida

library(invgamma)
library(Matrix)
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
t_observed  <- c(50, 100, 150, 175)
#t_observed  <- c(250, 500, 750, 2000)
#t_observed <- c(250, 500, 750, 1000)

theta1_present <- TRUE
theta2_present <- FALSE

if (theta1_present) {
    theta1_true <- data$theta
    lambda_true <- exp(theta1_true)
}

if (theta2_present)
    theta2_true <- data$theta2

printf <- function(...)
    cat(paste(sprintf(...), "\n"))

logsumexp <- function(x) {
    cc <- max(x)
    return(cc + log(sum(exp(x - cc))))
}

log_p_yt <- function(yt, theta_t1) {
    res <- yt * theta_t1 - exp(theta_t1)
    res[!is.finite(res)] <- -Inf
    return(res)
}

# Chan method
diag_base <- c(rep(2, Tt-1), 1)
off_base  <- rep(-1, Tt-1)
K0 <- bandSparse(Tt, k=c(0,1), diagonals=list(diag_base, off_base), symmetric=TRUE)

chan_build <- function(z, f, theta2, theta_01, theta_02, W1, sample=TRUE) {

    P <- K0/W1 + Diagonal(x=1/f)
    ch <- Cholesky(P, LDL=FALSE, perm=FALSE)
    drift <- (c(theta_02, theta2) - c(theta2, 0)) / W1
    b <- z/f + drift[1:Tt]
    b[1] <- b[1] + (theta_01 + theta_02) / W1
    eta_hat <- as.vector(solve(ch, b))
    return(list(eta_hat=eta_hat, ch=ch, b=b))
}

chan_sample <- function(res) {
    u <- rnorm(Tt)
    x <- as.vector(solve(res$ch, u, system="Lt"))
    return(res$eta_hat + x)
}

# Prior Hyperparameters

# theta_01 ~ N(mu_01, sigma2_01)
mu_01 <- log(y[1] + 0.5)
sigma2_01 <- 10

# theta_02 ~ N(mu_02, sigma2_02)
mu_02 <- 0
sigma2_02 <- 10

# W1 ~ InvGamma(alpha_W1, beta_W1)
alpha_W1 <- 2
beta_W1 <- 0.1

# W2 ~ InvGamma(alpha_W2, beta_W2)
alpha_W2 <- 2
beta_W2 <- 0.1

N <- 10000
burnin <- 1000

W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <- numeric(N)
theta_02_hist <- numeric(N)
theta1_hist <- matrix(0, N, Tt)
theta2_hist <- matrix(0, N, Tt)

# Initial Values
theta1 <- log(y + 0.5)
theta2 <- c(diff(theta1), 0)

theta_01 <- log(y[1] + 0.5)
theta_02 <- y[2]-y[1]
W1 <- 0.01
W2 <- 0.01

ess_smc <- numeric(N)
itr_irls <- numeric(N) # number of iterations of IRLS (for earch gibbs step)
M_irls_max <- 20 # maximum iterations for IRLS

M_is <- 3
tol <- 1e-4
theta1_tilde <- log(y + 0.5)

#####
start_time <- proc.time()
for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[3]]
        printf("Iteration %d / %d | Elapsed time: %.0f s | j mean: %.1f", n, N, elapsed_time, mean(itr_irls[1:(n-1)]))
    }

    # Sample theta_02 (conjugated normal)
    sigma2_02_bar <- (1 / sigma2_02 + 1 / W1 + 1 / W2)^(-1)
    mu_02_bar <- sigma2_02_bar * ((theta1[1] - theta_01) / W1 +
                                      theta2[1] / W2 +
                                      mu_02 / sigma2_02)
    theta_02 <- rnorm(1, mean = mu_02_bar, sd = sqrt(sigma2_02_bar))

    # Sample theta_01 (conjugated normal)
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar * (mu_01 / sigma2_01 +
                                      (theta1[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))

    # Sample W1 (conjugated invgamma)
    dif1   <- theta1 - c(theta_01, theta1[-Tt])
    diffs1 <- dif1 - c(theta_02, theta2[-Tt])
    alpha_W1_bar <- alpha_W1 + Tt / 2
    beta_W1_bar  <- beta_W1 + 0.5 * sum(diffs1^2)
    W1 <- rinvgamma(1, shape = alpha_W1_bar, rate = beta_W1_bar)

    # Sample W2 (conjugated invgamma)
    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    alpha_W2_bar <- alpha_W2 + Tt / 2
    beta_W2_bar  <- beta_W2 + 0.5 * sum(diffs2^2)
    W2 <- rinvgamma(1, shape = alpha_W2_bar, rate = beta_W2_bar)

    sd_W1 <- sqrt(W1)
    sd_W2 <- sqrt(W2)

    # ------------------------------------------------------------------
    # Importance Sampling for  theta_t1
    # ------------------------------------------------------------------

    # Local approximation (IRLS)
    theta1_tilde_old <- theta1_tilde
    for (j in 1:M_irls_max) {

            f_t <- exp(-theta1_tilde)       # observational variance
            z_t <- theta1_tilde + f_t*y - 1 # pseudo-observation

            res <- chan_build(z_t, f_t, theta2, theta_01, theta_02, W1)
            theta1_tilde <- res$eta_hat

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
            theta1_prop  <- chan_sample(res)
            trajectories[i, ] <- theta1_prop

            # Log-weights: log p(y|theta1) - log g(y|theta1)
            log_p <- sum(y * theta1_prop - exp(theta1_prop))
            log_g <- sum(-0.5 * log(2*pi*f_t) - 0.5 * (theta1_prop - z_t)^2 / f_t)
            log_w[i] <- log_p - log_g
        }

        #log_w <- log_w - max(log_w)
        log_w <- log_w - logsumexp(log_w)
        w <- exp(log_w)
        idx <- sample(1:M_is, 1, prob=w)
        theta1 <- trajectories[idx, ]
    } else {
        theta1  <- chan_sample(res)
    }


    # ------------------------------------------------------------------
    # Component-wise Gibbs for  theta_t2 | theta1, W1, W2
    # ------------------------------------------------------------------

    # Conditional variances (constant)
    sigma2_t2_bar_interior <- (1 / W1 + 2 / W2)^(-1)
    sigma2_t2_bar_last <- W2

    # t = 1
    sigma2_bar <- (1 / W1 + 2 / W2)^(-1)
    mu_bar <- sigma2_bar * ((theta1[2] - theta1[1]) / W1 +
                                theta2[2]  / W2 +
                                theta_02 / W2)
    theta2[1] <- rnorm(1, mean = mu_bar, sd = sqrt(sigma2_bar))

    # t = 2, ..., T-1
    # for (t in 2:(Tt - 1)) {
    #     mu_bar <- sigma2_t2_bar_interior * ((theta1[t + 1] - theta1[t]) / W1 +
    #                                             theta2[t + 1] / W2 +
    #                                             theta2[t - 1] / W2)
    #     theta2[t] <- rnorm(1, mean = mu_bar, sd = sqrt(sigma2_t2_bar_interior))
    # }
    sd_t2 <- sqrt(sigma2_t2_bar_interior)
    d1    <- diff(theta1) / W1       # (theta1[t+1]-theta1[t])/W1 para t=1..T-1
    for (t in 2:(Tt-1)) {
        mu_bar  <- sigma2_t2_bar_interior * (d1[t] + (theta2[t+1] + theta2[t-1])/W2)
        theta2[t] <- rnorm(1, mu_bar, sd_t2)
    }

    # t = T
    theta2[Tt] <- rnorm(1, mean = theta2[Tt - 1], sd = sd_W2)

    # Store the results
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_hist[n, ] <- theta1
    theta2_hist[n, ] <- theta2

}

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


# Plots ####
x <- 1:Tt

# lambda ####
lambda_mean <- exp(theta1_mean)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(
    x,
    y,
    type = "l",
    col = "gray",
    xlab = "t",
    ylab = "",
    main = "Poisson 2nd Order Polynomial Model"
)
points(x, y, pch = 20)
lines(x, lambda_mean, col = "red", lwd = 2)
lines(x, lambda_true, col = "blue", lwd = 2)
legend(
    "topright",
    legend = expression(y[t], lambda[t], hat(lambda)[t]),
    col = c("black", "blue", "red"),
    lty = c(NA, 1, 1),
    lwd = c(NA, 2, 2),
    pch = c(20, NA, NA),
    bty = "n"
)

# theta_t1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
ylim_range <- range(theta1_mean, theta1_true)
plot(
    x,
    theta1_mean,
    type = "l",
    col = "red",
    lwd = 2,
    ylim = ylim_range,
    xlab = "t",
    ylab = "",
    main = "theta_t1"
)
lines(x, theta1_true, col = "blue", lwd = 2)
legend(
    "topright",
    legend = expression(hat(theta)[t1], theta[t1]),
    col = c("red", "blue"),
    lwd = 2,
    bty = "n"
)

# theta_t2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(
    x,
    theta2_mean,
    type = "l",
    col = "red",
    lwd = 2,
    xlab = "t",
    ylab = "",
    main = "theta_t2"
)
legend(
    "topright",
    legend = expression(hat(theta)[t2], theta[t2]),
    col = c("red", "blue"),
    lwd = 2,
    bty = "n"
)

# Traceplot W1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(
    W1_hist[-(1:burnin)],
    type = "l",
    xlab = "n",
    ylab = "W1",
    main = "Traceplot of W1"
)

# Traceplot W2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(
    W2_hist[-(1:burnin)],
    type = "l",
    xlab = "n",
    ylab = "W2",
    main = "Traceplot of W2"
)

# Traceplots theta1 ####
par(mfrow = c(2, 2))
for (t in t_observed) {
    plot(
        theta1_hist[, t],
        type = "l",
        main = bquote(theta[list(.(t), 1)]),
        xlab = "n",
        ylab = ""
    )
}

# Traceplots theta1
par(mfrow = c(2, 2))
for (t in t_observed) {
    plot(
        theta2_hist[, t],
        type = "l",
        main = bquote(theta[list(.(t), 2)]),
        xlab = "n",
        ylab = ""
    )
}
