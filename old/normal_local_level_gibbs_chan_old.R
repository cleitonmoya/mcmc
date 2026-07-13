# Local Level Dynamic Linear Model
#
# Model:
#  y_t ~ theta_t + nu_t,               nu_t ~ N(0, V)
#  theta_t = theta_{t-1} + omega_t, omega_t ~ N(0, W)
#
# Priors:
#  theta_0 | D_0 ~ N(mu_0, sigma2_0)
#  1/V | D_0 ~ gamma(shape=nu_V, rate=eta_V)
#  1/W | D_0 ~ gamma(shape=nu_W, rate=eta_W)
#
# MCMC:
#  Gibbs with Chan Method
#
# Reference: Chan, J. C. C., & Jeliazkov, I. (2009). Efficient simulation
# and integrated likelihood estimation in state space models.
# International Journal of Mathematical Modelling and Numerical Optimisation,
# 1(1/2), 101. https://doi.org/10.1504/IJMMNO.2009.030090
#
# Author: Cleiton Moya de Almeida

library(Matrix)    # sparse matrix manipulation
library(coda)      # diagnostics for mcmc

#graphics.off()    # close the plots
rm(list = ls())    # clear the environment
#cat("\014")       # clear the console
tp <- Matrix::t    # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback
set.seed(42)


# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "normal_local_level_sim_2000" # rds file with data
t_obs <- c(50, 100, 150, 175)
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
y_matr <- Matrix(y, ncol=1)
Tt <- length(y) # dimension T

theta_sim_available <- TRUE    # simulated theta is available
if (theta_sim_available) {
    theta_true <- data$theta
    W_true <- data$W
    V_true <- data$V
}

# AUXILIARY FUNCTIONS ####

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))


# SIMULATION MAIN PARAMETERS ####

# Prior hyperparameters
# theta_0 ~ N(mu_0, sigma2_0)
mu_0 <- 0
sigma2_0 <- 10

# phi_V = 1/V ~ Gamma(shape=nu_V, rate=eta_V)
nu_V <- 2
eta_V <- 0.01

# phi_W = 1/W ~ Gamma(shape=nu_W, rate=eta_W)
nu_W <- 2
eta_W <- 0.01

N <- 10000         # number of simulation steps
burnin <- 1000     # number of burn-in steps

# Initialization
V <- 0.001
W <- 0.001
theta <- y

# Chan prior precision matrix base form (without V and W)
sub_diag_base <- rep(-1, Tt-1)
main_diag_base <- c(rep(2, Tt-1), 1)
K0 <- bandSparse(n=Tt, k=c(0, -1), diagonals=list(main_diag_base, sub_diag_base), symmetric = TRUE)
Ch_factor0 <- Cholesky(K0, perm = FALSE)

# Use the basic form K0
chan_sample <- function(y, V, W, theta_0) {

    # Posterior precision matrix (built uppon K0)
    invVcal <- .sparseDiagonal(n=Tt, x=1/V)

    #sub_diag <- rep(-1/W, Tt-1)
    #main_diag <- c(rep(2/W, Tt-1), 1/W)
    #K <- bandSparse(n=Tt, k=c(0, -1), diagonals=list(main_diag, sub_diag), symmetric = TRUE)
    #P <- forceSymmetric(K + invVcal)
    P <- forceSymmetric(K0/W + invVcal)

    # Cholesky decomposition
    Ch_factor <- update(Ch_factor0, P)

    # Smooothing
    b <- invVcal %*% y_matr + Matrix(data=c(theta_0/W, rep(0, Tt-1)), ncol=1, sparse=TRUE)
    theta_hat <- as.numeric(Matrix::solve(Ch_factor, b, system="A"))

    # Sampling
    d <- Matrix::diag(Ch_factor)
    u <- rnorm(Tt)
    w <- u / sqrt(d)
    x <- as.vector(Matrix::solve(Ch_factor, w, system="Lt"))
    theta <- theta_hat + x
    return(theta)
}


# Auxiliary vectors and matrix to store the results
theta_hist <- matrix(nrow=N, ncol=Tt)
V_hist <- numeric(N)
W_hist <- numeric(N)
theta_0_hist <-numeric(N)


# MAIN LOOP ####

start_time = proc.time() # execution time
for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[3]]
        printf("Iteration %d / %d | Elapsed time: %.0f s", n, N, elapsed_time)
    }

    # Sample theta_0
    sigma2_0_bar <- (1/sigma2_0 +1/W)^(-1)
    mu_0_bar <- sigma2_0_bar*(mu_0/sigma2_0 + theta[1]/W)
    theta_0 <- rnorm(1, mean=mu_0_bar, sd=sqrt(sigma2_0_bar))

    # Sample phi_V
    nu_V_bar <- nu_V + Tt/2
    diffs <- y-theta
    eta_V_bar <- eta_V + 0.5 * sum(diffs^2)
    phi_V <- rgamma(1, shape=nu_V_bar, rate=eta_V_bar)
    V <- 1/phi_V

    # Sample phi_W
    nu_W_bar <- nu_W + Tt/2
    diffs <- theta - c(theta_0, theta[-Tt])
    eta_W_bar <- eta_W + 0.5 * sum(diffs^2)
    phi_W <- rgamma(1, shape=nu_W_bar, rate=eta_W_bar)
    W <- 1/phi_W

    # Sample theta (Chan)
    theta <- chan_sample(y, V, W, theta_0)

    # Store the sampled values
    theta_0_hist[n] <- theta_0
    V_hist[n] <- V
    W_hist[n] <- W
    theta_hist[n, ] <- theta
}


# RESULTS AND DIAGNOSTICS ####

# Execution time
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)


# Posterior mean
theta_mean <- colMeans(theta_hist[-(1:burnin), ])

V_mean <- mean(V_hist[-(1:burnin)])
V_median <- median(V_hist[-(1:burnin)])
if (theta_sim_available) printf("V true: %.5f", V_true)
printf("V mean: %.5f", V_mean)
printf("V median: %.5f", V_median)

W_mean <- mean(W_hist[-(1:burnin)])
W_median <- median(W_hist[-(1:burnin)])
if (theta_sim_available) printf("W true: %.5f", W_true)
printf("W mean: %.5f", W_mean)
printf("W median: %.5f", W_median)


# Log-likelihood
loglik <- sum(dnorm(y, mean=theta_mean, sd = sqrt(V_mean), log=TRUE))
printf("Log-likelihood: %.2f", loglik)


# Effective sample size
printf("Effective Sample Size:")
ess_V <- effectiveSize(mcmc(V_hist[-(1:burnin)]))
ess_W <- effectiveSize(mcmc(W_hist[-(1:burnin)]))
printf("\tV: %.0f", ess_V)
printf("\tW: %.0f", ess_W)

ess_theta <- effectiveSize(mcmc(theta_hist[-(1:burnin),]))
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ess_theta, type="l", main=expression("Effective sample size of " * theta[t]), xlab="t")
hist(ess_theta)
ess_theta_mean <-mean(ess_theta)
printf("\ttheta (mean): %.0f", ess_theta_mean)


# Effective sample size per second
printf("Effective Sample Size / second:")
printf("\tV: %.2f", ess_V/elapsed_time)
printf("\tW: %.2f", ess_W/elapsed_time)
printf("\ttheta (mean): %.2f", ess_theta_mean/elapsed_time)

ess_s_theta <- ess_theta/elapsed_time
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ess_s_theta, type="l", main=expression("Effective sample size per second of " * theta[t]), xlab="t")


# Geweke diagnostic: Z test for two mean difference
#   H0: segments with different means -> chain has not converged
z_V <- unname(geweke.diag(V_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
z_W <- unname(geweke.diag(W_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])

printf("Geweke convergence diagnostic")
printf("\tz_V: %.2f", z_V)
printf("\tz_w: %.2f", z_W)

z_theta <- unname(geweke.diag(theta_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(z_theta, type="l", main=expression("Geweke diagnostic for " * theta[t]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")


# Plots ####
# y, theta_true, theta_mean ####
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Local Level Dynamic Linear Model")
points(x, y, pch = 20)
lines(x, theta_mean, col="red", lwd=2)
if (theta_sim_available) {
    lines(x, theta_true, col="blue", lwd=2)
    legend("topright",
           legend = expression(y[t], theta[t], hat(theta)[t]),
           col = c("black", "blue", "red"),
           lty = c(NA, 1, 1),
           lwd = c(NA, 2, 2),
           pch = c(20, NA, NA),
           bty = "n")
} else {
    legend("topright",
           legend = expression(y[t], hat(theta)[t]),
           col = c("black", "red"),
           lty = c(NA, 1),
           lwd = c(NA, 2),
           pch = c(20, NA),
           bty = "n")
}


# theta_true vs theta_mean ####
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

# Posterior distribution of V ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(V_hist[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = bquote(theta[.(t)]), main ="Posterior of V")
lines(density(V_hist[-(1:burnin)]), col = "blue", lwd = 2)

# Traceplot for V ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(V_hist[-(1:burnin)], type="l", xlab="n", ylab="V", main="Traceplot of V")


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
