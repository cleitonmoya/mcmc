# Poisson - 2nd Order Polynomial Dynamic Model
#
# Model:
#  y_t ~ Poisson(exp{theta_t1})
#  theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1, omega_t1 ~ N(0, W1)
#  theta_t2 =                 theta_{t-1,1} + omega_t2, omega_t2 ~ N(0, W2)
#
# Priors:
#  theta_01 | D_0 ~ N(mu_01, sigma2_01)
#  theta_01 | D_0 ~ N(mu_02, sigma2_02)
#  1/W1 ~ gamma(shape=nu_01, rate=eta_01)
#  1/W2 ~ gamma(shape=nu_02, rate=eta_02)
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


#graphics.off()     # close the plots
rm(list = ls())     # clear the environment
#cat("\014")         # clear the console
tp <- Matrix::t     # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback
set.seed(42)

library(Rfast)      # provide colMedians()
library(coda)

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "poisson_pol2_200"
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
#y <- as.data.frame(Seatbelts)$DriversKilled
Tt <- length(y) # dimension T
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

# Full conditional log-posterior for theta_t1, t=1, ..., T-1
# theta_t1: theta_{t1}
# theta_tm11: theta_{t-1,1}
# theta_tp11: theta_{t+1,1}
# theta_t2: theta_{t2}
# theta_tm12: theta_{t-1,2}
logpost_theta_t1 <- function(theta_t1, theta_tm11, theta_tp11,
                             theta_t2, theta_tm12, yt, W1) {
    sigma2_star <- W1/2
    mu_star <- ((theta_tm11+theta_tm12)+(theta_tp11-theta_t2))/2

    p1 <- yt*theta_t1 - exp(theta_t1) # log-likelihood
    p2 <- -(theta_t1 - mu_star)^2/(2*sigma2_star)
    logp <- p1 + p2
    return(logp)
}

# Full conditional log-posterior for theta_T1 (t=T)
logpost_theta_T1 <- function(theta_t1, theta_tm11, theta_tm12, yt, W1) {
    p1 <- yt*theta_t1 - exp(theta_t1) # log-likelihood
    p2 <- -(theta_t1 - theta_tm11 - theta_tm12)^2/(2*W1)
    logp <- p1+p2
    return(logp)
}

# Sample theta_t1 ~ logpost_theta_t1 (Metropolis step)
# final_t: boolean (0: t<T; 1: t=T)
sample_theta_t1 <- function(theta_t1_current, theta_tm11, theta_tp11,
                            theta_t2, theta_tm12,
                            yt, W1, varsigma2, final_t) {

    # proposed theta
    theta_t1_prop <- rnorm(1, mean=theta_t1_current, sd=sqrt(varsigma2))

    # acceptance/rejection step
    ac <- 0 # accepted flag
    logu <- log(runif(1))

    if (final_t) {
        logp1 <- logpost_theta_T1(theta_t1_prop, theta_tm11, theta_tm12, yt, W1)
        logp2 <- logpost_theta_T1(theta_t1_current, theta_tm11, theta_tm12, yt, W1)

    } else {
        logp1 <- logpost_theta_t1(theta_t1_prop, theta_tm11, theta_tp11,
                                  theta_t2, theta_tm12, yt, W1)
        logp2 <- logpost_theta_t1(theta_t1_current, theta_tm11, theta_tp11,
                                  theta_t2, theta_tm12, yt, W1)
    }
    logr <- logp1 - logp2

    # acceptance criteria
    if (logu < logr){
        theta_t1 <- theta_t1_prop
        ac <- 1
    } else {
        theta_t1 <- theta_t1_current
    }

    return(list(theta_t1=theta_t1, ac=ac))
}

#####
# SIMULATION MAIN PARAMETERS

# Prior hyperparameters
# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- 0
sigma2_01 <- 100

# theta_02 ~ N(mu_2, sigma2_02)
mu_02     <- 0
sigma2_02 <- 100

# phi1 = W1^(-1) ~ Gamma(nu_01, eta_01)
nu_01  <- 2
eta_01 <- 0.01

# phi2 = W2^(-1) ~ Gamma(nu_02, eta_02)
nu_02  <- 2
eta_02 <- 0.0001

# Initialization
W2 <- 0.01
W1 <- 0.01
phi1 <- 1/W1
phi2 <- 1/W2
theta_01 <- 0
theta_02 <- 0
theta1 <- numeric(Tt)
theta2 <- numeric(Tt)

N <- 10000           # Number of steps
burnin <- 1000       # Number of burn-in steps
#varsigma2 <- 0.03   # Random walkikng variance hyperparameter - Doppler
if (Tt == 200) varsigma2 <- 0.05    # Random walkikng variance hyperparameter - Poisson_pol2_200
if (Tt == 2000) varsigma2 <- 0.01    # Random walkikng variance hyperparameter - Poisson_pol2_200

# Auxiliary vectors and matrix to store the results
theta1_hist <- matrix(nrow=N, ncol=Tt)
theta2_hist <- matrix(nrow=N, ncol=Tt)
W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <-numeric(N)
theta_02_hist <-numeric(N)
ac_hist <- numeric(N)

#####
# Gibbs sampling

 # Main loop
start_time = proc.time() # execution time
for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[1]]
        printf("Iteration %d / %d | Elapsed CPU time: %.0f s", n, N, elapsed_time)
    }

    # Sample theta_01
    sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
    mu_01_bar <- sigma2_01_bar*(mu_01/sigma2_01 + (theta1[1]-theta_02)/W1)
    theta_01 <- rnorm(1, mean=mu_01_bar, sd=sqrt(sigma2_01_bar))

    # Sample theta_02
    sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
    mu_02_bar <- sigma2_02_bar*((theta1[1]-theta_01)/W1 +
                                    theta2[1]/W2 + mu_02/sigma2_02)
    theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))

    # Sample phi1
    nu_01_bar <- nu_01 + Tt/2
    dif1 <- theta1 - c(theta_01, theta1[-Tt])
    dif2 <- dif1 - c(theta_02, theta2[-Tt])
    eta_01_bar <- eta_01 + 0.5 * sum(dif2^2)
    phi1 <- rgamma(1, nu_01_bar, eta_01_bar)
    W1 <- 1/phi1

    # Sample phi2
    nu_02_bar <- nu_02 + Tt/2
    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    eta_02_bar <- eta_02 + 0.5 * sum(diffs2^2)
    phi2 <- rgamma(1, nu_02_bar, eta_02_bar)
    W2 <- 1/phi2

    # Sample theta_t1 (random walking Metropolis) and
    #        theta_t2 ~ N()
    n_ac <- 0 # number of accepted samples
    for (t in 1:Tt) {

        if (t < Tt) {

            sigma2_star <- (1/W1 + 2/W2)^(-1) # for theta_t2

            if (t==1) {
                # theta_t11
                res <- sample_theta_t1(theta1[t], theta_01, theta1[t+1],
                                       theta2[t], theta_02,
                                       y[t], W1, varsigma2, final_t=FALSE)
                theta1[t] <- res$theta_t1
                # theta_t12
                mu_star <- sigma2_star*((theta1[t+1] - theta1[t])/W1 +
                                            (theta_02 + theta2[t+1])/W2)
            } else {

                res <- sample_theta_t1(theta1[t], theta1[t-1], theta1[t+1],
                                       theta2[t], theta2[t-1],
                                       y[t], W1, varsigma2, final_t=FALSE)
                theta1[t] <- res$theta_t1
                mu_star <- sigma2_star*((theta1[t+1] - theta1[t])/W1 +
                                            (theta2[t-1] + theta2[t+1])/W2)
            }

        } else {
            res <- sample_theta_t1(theta1[t], theta1[t-1], NULL,
                                   theta2[t], theta2[t-1],
                                   y[t], W1, varsigma2, final_t=TRUE)
            theta1[t] <- res$theta_t1
            mu_star <- theta2[t-1]
            sigma2_star <- W2
        }

        theta2[t] <- rnorm(1, mean=mu_star, sd=sqrt(sigma2_star))
        ac <- res$ac # flag: sample accpeted(1) or not (0)
        n_ac <- n_ac + ac
    }

    # Mean acceptance ratio of theta_t1
    ac_hist[n] <- n_ac/Tt

    # Store the sampled values
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

printf("Mean acception ratio of theta1: %.2f", mean(ac_hist))

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
printf("Effective Sample Size:")
ess_w1 <- effectiveSize(mcmc(W1_hist[-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W2_hist[-(1:burnin)]))
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)

ess_theta1 <- effectiveSize(mcmc(theta1_hist[-(1:burnin),]))
ess_theta2 <- effectiveSize(mcmc(theta2_hist[-(1:burnin),]))
printf("\ttheta1 (mean): %.2f", mean(ess_theta1))
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


# Traceplot of acceptance ratio of theta_t1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Acceptance ratio of " * theta[t*1]))


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
