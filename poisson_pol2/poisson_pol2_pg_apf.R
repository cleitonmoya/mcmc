# Poisson - 2nd Order Polynomial Dynamic Model
# Particle Gibbs (PG) with Backward Sampling + Component-wise Gibbs for theta_t2
# Strategy: - APF for theta_t1 (scalar state);
#           - theta_t2 sampled via component-wise Gibbs
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269-342.
# Author: Cleiton Moya de Almeida

library(coda)

#graphics.off()     # close the plots
#cat("\014")        # clear the console
rm(list = ls())     # clear the environment
set.seed(42)
options(error = function() traceback(2))  # more informative traceback

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source_file <- "poisson_pol2_200"
data <- readRDS(paste("../data/", source_file, ".rds", sep=""))
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

# Auxiliary functions
printf <- function(...) cat(paste(sprintf(...), "\n"))

logsumexp <- function(x) {
  cc <- max(x)
  return(cc + log(sum(exp(x - cc))))
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
eta_01  <- 0.01

# phi2 = 1/W2 ~ Gamma(nu_02, eta_02)
nu_02 <- 2
eta_02  <- 0.0001


# Simulation parameters
N <- 10000      # number of Gibbs iterations
K <- 200        # number of particles
burnin <- 1000


# Initial Values
theta1 <- numeric(Tt)
theta2 <- numeric(Tt)
theta_01 <- 0
theta_02 <- 0
W1 <- 0.01
W2 <- 0.01


# Auxiliary variables
W1_hist          <- numeric(N)
W2_hist          <- numeric(N)
theta_01_hist    <- numeric(N)
theta_02_hist    <- numeric(N)
theta1_hist <- matrix(0, N, Tt)
theta2_hist <- matrix(0, N, Tt)

start_time <- proc.time()
ess_smc <- numeric(N)

for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[1]]
        printf("Iteration %d / %d | Elapsed CPU time: %.0f s", n, N, elapsed_time)
    }

    # Sample theta_02
    sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
    mu_02_bar <- sigma2_02_bar * (
        (theta1[1] - theta_01)/W1 + theta2[1]/W2 +
        mu_02/sigma2_02
    )
    theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))

    # Sample theta_01
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar *
        (mu_01 / sigma2_01 + (theta1[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))

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


    #
    # Conditional SMC for theta_t1
    #
    theta_1_k   <- matrix(0, Tt, K)
    log_w_tilde <- matrix(0, Tt, K)

    # t = 1
    theta_1_k[1, ] <- rnorm(K, mean = theta_01 + theta_02, sd = sd_W1)
    theta_1_k[1, K] <- theta1[1]

    log_w_1 <- log_p_yt(y[1], theta_1_k[1, ])
    log_w_tilde[1, ] <- log_w_1 - logsumexp(log_w_1) # normalizing

    # t = 2, ..., T
    for (t in 2:Tt) {
        # Predictor (auxiliary variable)
        theta_hat_t1_k = theta_1_k[t - 1, ] + theta2[t - 1]

        # Auxiliary weights
        log_lambda_k = y[t] * theta_hat_t1_k - exp(theta_hat_t1_k)

        # First stage resampling
        log_aux <- log_w_tilde[t - 1, ] + log_lambda_k
        A <- sample(1:K, K, replace = TRUE, prob = exp(log_aux - max(log_aux)))
        A[K] <- K   # reference path

        # Propagation
        theta_t1_k <- rnorm(K, mean=theta_1_k[t-1, A] + theta2[t-1], sd=sd_W1)
        theta_t1_k[K] <- theta1[t]
        theta_1_k[t, ] <- theta_t1_k

        # Updated weights
        log_w_t <- log_p_yt(y[t], theta_1_k[t, ]) - log_lambda_k[A]
        log_w_tilde[t, ] <- log_w_t - logsumexp(log_w_t)
    }

    ess_smc[n] <- 1 / sum(exp(2 * (log_w_tilde[Tt, ] - logsumexp(log_w_tilde[Tt, ]))))


    # backward sampling for theta1
    # Backward weights: w_t^k * N(theta1[t+1] | theta_t1^k + theta2[t], W1)
    k_final <- sample(1:K, 1, prob = exp(log_w_tilde[Tt, ] - max(log_w_tilde[Tt, ])))
    theta1[Tt] <- theta_1_k[Tt, k_final]

    for (t in (Tt - 1):1) {
        log_bw <- log_w_tilde[t, ] +
            dnorm(theta1[ +1],
                mean = theta_1_k[t, ] + theta2[t],
                sd = sd_W1,
                log = TRUE
            )
        log_bw <- log_bw - max(log_bw)
        bw <- exp(log_bw)
        bw <- bw / sum(bw)
        b  <- sample(1:K, 1, prob = bw)
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
    sigma2_t2_bar_interior <- (1 / W1 + 2 / W2)^(-1)
    sigma2_t2_bar_last <- W2

    # t = 1: theta2[-1] = theta_02 (estado inicial)
    # Tratar separadamente pois theta2[t-1] = theta_02
    sigma2_bar <- (1 / W1 + 2 / W2)^(-1)
    mu_bar <- sigma2_bar * ((theta1[2] - theta1[1]) / W1 +
                                theta2[2]/W2 + theta_02/W2)
    theta2[1] <- rnorm(1, mean = mu_bar, sd = sqrt(sigma2_bar))

    # t = 2, ..., T-1
    for (t in 2:(Tt - 1)) {
        mu_bar <- sigma2_t2_bar_interior * ((theta1[t + 1] - theta1[t]) / W1 +
                                                theta2[t+1]/W2 +
                                                theta2[t-1]/W2)
        theta2[t] <- rnorm(1,
                                mean = mu_bar,
                                sd = sqrt(sigma2_t2_bar_interior))
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

#### Simulation summary

# Execution time
elapsed_time <- (proc.time() - start_time)[[1]]
printf("Execution time: %.0f s", elapsed_time)

# Posterior mean
theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])

if (theta1_present) lambda_true <- exp(theta1_true)
lambda_mean <- exp(theta1_mean)
W1_mean <- mean(W1_hist[-(1:burnin)])
W2_mean <- mean(W2_hist[-(1:burnin)])
printf("W1 posterior mean: %.6f", W1_mean)
printf("W2 posterior mean: %.6f", W2_mean)

# Log-likelihood
loglik <- sum(dpois(y, lambda_mean, log=TRUE))
printf("Log-likelihood: %.2f", loglik)

# Effective sample size
t1  <- theta1_hist[-(1:burnin), ]
t2  <- theta2_hist[-(1:burnin), ]
ess_theta1 <- apply(t1, 2, effectiveSize)
ess_theta2 <- apply(t2, 2, effectiveSize)
ess_W1 <- effectiveSize(W1_hist[-(1:burnin)])
ess_W2 <- effectiveSize(W2_hist[-(1:burnin)])

printf("ESS theta1 (min/median): %.0f / %.0f", min(ess_theta1), median(ess_theta1))
printf("ESS theta2 (min/median): %.0f / %.0f", min(ess_theta2), median(ess_theta2))
printf("ESS W1: %.0f", ess_W1)
printf("ESS W1/sec: %.2f", ess_W1/elapsed_time)
printf("ESS W2: %.0f", ess_W2)
printf("ESS W2/sec: %.2f", ess_W2/elapsed_time)

#### Plots
x <- 1:Tt

# Y and lambda ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(x, y, type="l", col="gray", xlab="t", ylab="",
     main="Poisson 2nd Order Polynomial Model")
points(x, y, pch=20)
lines(x, lambda_mean, col="red",  lwd=2)

if (theta1_present) {
    lines(x, lambda_true, col="blue", lwd=2)
    legend("topright",
           legend=expression(y[t], lambda[t], hat(lambda)[t]),
           col=c("black","blue","red"),
           lty=c(NA,1,1), lwd=c(NA,2,2), pch=c(20,NA,NA), bty="n")
} else {
    legend("topright",
           legend=expression(y[t],  hat(lambda)[t]),
           col=c("black","red"),
           lty=c(NA,1), lwd=c(NA,2), pch=c(20,NA), bty="n")
}

# theta_t1 ####
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

# theta_t2 ####
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

# Traceplot W1 ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(W1_hist, type="l", xlab="n", ylab="W1",
     main="Traceplot of W1")


# Traceplot W2 ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(W2_hist, type="l", xlab="n", ylab="W2",
     main="Traceplot of W2")


# Traceplots theta1 ####
par(mfrow=c(2,2))
for (t in t_obs) {
  plot(theta1_hist[, t], type="l",
       main=bquote(theta[list(.(t),1)]), xlab="n", ylab="")
}

# Traceplots theta2 ####
par(mfrow=c(2,2))
for (t in t_obs) {
  plot(theta2_hist[, t], type="l",
       main=bquote(theta[list(.(t),2)]), xlab="n", ylab="")
}

# Effective sample size ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_smc, type="l")
