# Poisson Local Level Model
#
# Model:
#  y_t ~ Poisson(exp{theta_t})
#  theta_t = theta_{t-1} + omega_t, omega_t ~ N(0, W)
#
# Priors:
#  theta_0 | D_0 ~ N(mu_0, sigma2_0)
#  1/W ~ gamma(shape=nu_0, rate=eta_0)
#
# MCMC:
#  Component-Wise Metropolis within Gibbs for theta_t
#  Metropolis proposal: Random Walking
#
# Reference: Geweke, J., & Tanizaki, H. (2001).
#    Bayesian estimation of state-space models using the
#    Metropolis–Hastings algorithm within Gibbs sampling.
#    Computational statistics & data analysis, 37(2), 151-170
#
# Author: Cleiton Moya de Almeida


#graphics.off()    # close the plots
rm(list = ls())    # clear the environment
#cat("\014")       # clear the console
tp <- Matrix::t    # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback
set.seed(42)
library(coda)      # diagnostics for mcmc

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "doppler" # rds file with data
#t_obs <- c(50, 75, 100, 192)
t_obs <- c(250, 500, 750, 1000)
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
Tt <- length(y) # dimension T

theta_sim_available <- TRUE    # simulated theta is available
if (theta_sim_available) {
    theta_true <- data$theta
    lambda_true <- exp(theta_true)
}

# AUXILIARY FUNCTIONS ####

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))


# Full conditional log-posterior for theta_t, t=1, ..., T-1
# theta_tm1: theta_{t-1}
# theta_tp1: theta_{t+1}
logpost_theta_t <- function(theta_t, theta_tm1, theta_tp1, yt, W) {
    p1 <- yt*theta_t - exp(theta_t) # log-likelihood
    p2 <- -(theta_t - theta_tm1)^2/(2*W)
    p3 <- -(theta_tp1 - theta_t)^2/(2*W)
    logp <- p1+p2+p3
    return(logp)
}


# Full conditional log-posterior for theta_T (t=T)
logpost_theta_T <- function(theta_t, theta_tm1, yt, W) {
    p1 <- yt*theta_t - exp(theta_t) # log-likelihood
    p2 <- -(theta_t - theta_tm1)^2/(2*W)
    logp <- p1+p2
    return(logp)
}


# Sample theta_t ~ logpost_theta_t (Metropolis step)
# final_t: boolean (0: t<T; 1: t=T)
sample_theta_t <- function(theta_t_current, theta_tm1, theta_tp1,
                           yt, W, varsigma2, final_t) {

    # proposed theta
    theta_t_prop <- rnorm(1, mean=theta_t_current, sd=sqrt(varsigma2))

    # acceptance/rejection step
    ac <- 0 # accepted flag
    logu <- log(runif(1))

    if (final_t) {
        logp1 <- logpost_theta_T(theta_t_prop, theta_tm1, yt, W)
        logp2 <- logpost_theta_T(theta_t_current, theta_tm1, yt, W)
    } else {
        logp1 <- logpost_theta_t(theta_t_prop, theta_tm1, theta_tp1, yt, W)
        logp2 <- logpost_theta_t(theta_t_current, theta_tm1, theta_tp1, yt, W)
    }

    logr <- logp1 - logp2

    # acceptance criteria
    if (logu < logr){
        theta_t <- theta_t_prop
        ac <- 1
    } else {
        theta_t <- theta_t_current
    }

    return(list(theta_t=theta_t, ac=ac))
}


# SIMULATION MAIN PARAMETERS ####

# Prior hyperparameters
# theta_0 ~ N(mu_0, sigma2_0)
mu_0 <- 0
sigma2_0 <- 10

# phi = 1/W ~ Gamma(shape=nu_0, rate=eta_0)
nu_0 <- 2
eta_0 <- 0.01

N <- 10000         # number of simulation steps
varsigma2 <- 0.01  # random walkikng variance hyperparameter
burnin <- 1000     # number of burn-in steps

# Auxiliary vectors and matrix to store the results
theta_hist <- matrix(nrow=N, ncol=Tt)
W_hist <- numeric(N)
theta_0_hist <-numeric(N)
ac_theta_hist <- numeric(N)

# Initialization
W <- 0.001
theta <- y
phi <- 1/W


# MAIN LOOP ####

start_time = proc.time() # execution time
for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[3]]
        printf("Iteration %d / %d | Elapsed time: %.0f s", n, N, elapsed_time)
    }

    # 1. Sample theta_0
    sigma2_0_bar <- (1/sigma2_0 +1/W)^(-1)
    mu_0_bar <- sigma2_0_bar*(mu_0/sigma2_0 + theta[1]/W)
    theta_0 <- rnorm(1, mean=mu_0_bar, sd=sqrt(sigma2_0_bar))

    # 2. Sample phi (W^(-1))
    nu_0_bar <- nu_0 + Tt/2
    diffs <- theta - c(theta_0, theta[-Tt])
    eta_0_bar <- eta_0 + 0.5 * sum(diffs^2)
    phi <- rgamma(1, shape=nu_0_bar, rate=eta_0_bar)
    W <- 1/phi

    # 3. Sample \theta (sample \theta_t, t=1,...,T)
    n_ac <- 0 # number of accepçted samples for \theta
    for (t in 1:Tt) {

        # 3.2 sample theta_t (Metropolis)
        if (t < Tt) {
            if (t==1) {
                res <- sample_theta_t(theta[t], theta_0, theta[t+1],
                                      y[t], W, varsigma2, final_t=FALSE)
            } else {
                res <- sample_theta_t(theta[t], theta[t-1], theta[t+1],
                                      y[t], W, varsigma2, final_t=FALSE)
            }

        } else {
            res <- sample_theta_t(theta[t], theta[t-1], NULL,
                                  y[t], W, varsigma2, final_t=TRUE)
        }

        theta[t] <- res$theta_t
        ac <- res$ac # flag: sample accpeted(1) or not (0)
        n_ac <- n_ac + ac
    }

    # Mean acceptance ratio of theta
    ac_theta_hist[n] <- n_ac/Tt

    # Store the sampled values
    theta_0_hist[n] <- theta_0
    W_hist[n] <- W
    theta_hist[n, ] <- theta
}


# RESULTS AND DIAGNOSTICS ####

# Execution time
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)
printf("Mean acception ratio of theta: %.2f", mean(ac_theta_hist))

# Posterior mean
theta_mean <- colMeans(theta_hist[-(1:burnin), ])
lambda_mean <- exp(theta_mean)
W_mean <- mean(W_hist[-(1:burnin)])
W_median <- median(W_hist[-(1:burnin)])
printf("W mean: %.5f", W_mean)
printf("W median: %.5f", W_median)

# Log-likelihood
loglik <- sum(dpois(y, lambda_mean, log=TRUE))
printf("Log-likelihood: %.2f", loglik)

# Effective sample size
printf("Effective Sample Size:")
ess_w <- effectiveSize(mcmc(W_hist[-(1:burnin)]))
printf("\tW: %.0f", ess_w)

ess_theta <- effectiveSize(mcmc(theta_hist[-(1:burnin),]))
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ess_theta, type="l", main=expression("Effective sample size of " * theta[t]), xlab="t")
hist(ess_theta)
printf("\ttheta (mean): %.0f", mean(ess_theta))

# Effective sample size per second
printf("Effective Sample Size / second:")
printf("\tW: %.2f", ess_w/elapsed_time)
ess_s <- ess_theta/elapsed_time
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ess_s, type="l", main=expression("Effective sample size per second of " * theta[t]), xlab="t")

# Geweke diagnostic: Z test for two mean difference
#   H0: segments with different means -> chain has not converged
printf("Geweke convergence diagnostic")
z_w <- unname(geweke.diag(W_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
printf("\tz_w: %.2f", z_w)
z_theta <- unname(geweke.diag(theta_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(z_theta, type="l", main=expression("Geweke diagnostic for " * theta[t]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")


# Plots ####
# y, lambda_true, lambda_estimated ####
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson Local Level Polynomial Model")
points(x, y, pch = 20)
lines(x, lambda_mean, col="red", lwd=2)
if (theta_sim_available) {
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


# theta_true vs theta_estimated ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
if (theta_sim_available) {
    ylim_range <- range(theta_mean, theta_true)
} else {
    ylim_range <- range(theta_mean)
}
plot(x, theta_mean, type="l", col="red", lwd=2, ylim=ylim_range,
     xlab="t", ylab="", main=expression(theta[t]))
if (theta_sim_available) {
    lines(x, theta_true, col="blue", lwd=2)
    legend("topright", legend=expression(hat(theta)[t1], theta[t1]),
           col=c("red","blue"), lwd=2, bty="n")
} else {
    legend("topright", legend=expression(hat(theta)[t1]),
           col="red", lwd=2, bty="n")
}

# Posterior distribution of theta_t #####
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t)]), main = bquote("Posterior of " * theta[.(t)]))
    lines(density(theta_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}

# Posterior distribution of W ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W_hist[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = bquote(theta[.(t)]), main ="Posterior of W")
lines(density(W_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot for W ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W")


# Traceplots for theta_t ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta_hist[, t], type="l", main=bquote(theta[.(t)]), xlab="", ylab="")
}


# Traceplot for the mean acceptance ratio of theta ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_theta_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Mean acceptance ratio of " * theta[t]))
