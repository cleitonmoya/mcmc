# 2nd Order Polynomial Dynamic Linear Model
# MCMC: Gibbs with FFBS
# FFBS: Matrix
# Author: Cleiton Moya de Almeida

library(coda)

#graphics.off()      # close the plots
rm(list = ls())     # clear the environment
#cat("\014")         # clear the console
set.seed(42)
tp <- Matrix::t     # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "normal_pol2_sim_200" # csv file with data
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y

Tt <- length(y) # dimension Tt
if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 2000) t_obs <- c(500, 1000, 1500, 1750)

theta1_true <- data$theta1
theta2_true <- data$theta2
V_true      <- data$V
W1_true     <- data$W1
W2_true     <- data$W2


# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))

inv2x2 <- function(M) {
    det_M <- M[1,1]*M[2,2] - M[1,2]*M[2,1]
    inv <- matrix(c(M[2,2], -M[2,1], -M[1,2], M[1,1]), 2, 2) / det_M
    return(inv)
}


# Sample from a multivariate distribution using the Cholesky decomposition
rmvn_chol <- function(mu, Sigma) {
    L <- chol(Sigma)
    mu + drop(tp(L) %*% rnorm(length(mu)))
}

# Forward Filtering (Kalman Filter)
forward_filter <- function(theta_01, theta_02, Ff, G, W, y, V) {

    m_t <- as.matrix(c(theta_01, theta_02))
    C_t <- matrix(0, nrow=2, ncol=2)

    a <- array(data=NA, dim=c(2,1,Tt))
    R <- array(data=NA, dim=c(2,2,Tt))
    m <- array(data=NA, dim=c(2,1,Tt))
    C <- array(data=NA, dim=c(2,2,Tt))


    for (t in 1:Tt) {

        # Prior: (theta_t | D_{t-1}) ~ N[a_t, R_t]
        a_t <- G %*% m_t
        R_t <- G %*% C_t %*% tp(G) + W

        a[,,t] <- a_t
        R[,,t] <- R_t

        # One-step-ahead forecast
        f_t <- tp(Ff) %*% a_t
        Q_t <- drop(tp(Ff) %*% R_t %*% Ff + V)

        # Posterior: (theta_t | D_t) ~ N[m_t, C_t]
        e_t <- drop(y[t] - f_t)
        m_t <- a_t + (1/Q_t) * R_t %*% Ff * e_t
        C_t <- R_t - (1/Q_t) * R_t %*% Ff %*% tp(Ff) %*% R_t
        m[,,t] <- m_t
        C[,,t] <- C_t
    }

    return(list("a"=a, "R"=R, "m"=m, "C"=C))

}


# FFBS
sample_theta <- function(theta_01, theta_02, Ff, G, W, y, V) {

    # theta = {theta_1, ..., theta_T}
    theta <- array(data=NA, dim=c(2,1,Tt))
    res <- forward_filter(theta_01, theta_02, Ff, G, W, y, V)

    a <- res$a
    R <- res$R
    m <- res$m
    C <- res$C

    # Sample theta_T ~ N[m_T, C_T]
    theta_t <- rmvn_chol(mu=m[,,Tt], Sigma=C[,,Tt])
    theta[,,Tt] <- theta_t

    # Backward sampling: t = Tt-1, ..., 1
    # theta_t | theta_{t+1} ~ N(h_t, H_t) [
    # West and Harrison 1997, eq. 15.8, p. 570]
    B <- array(data=NA, dim=c(2,2,Tt))
    for (t in seq(Tt-1, 1)) {
        B[,,t] <- C[,,t] %*% tp(G) %*% inv2x2(R[,,t+1])
        h_t <- m[,,t] + B[,,t] %*% (theta_t - a[,,t+1])
        H_t <- C[,,t] - B[,,t] %*% R[,,t+1] %*% tp(B[,,t])
        theta_t <- rmvn_chol( mu=h_t, Sigma=H_t)
        theta[,,t] <- theta_t
    }

    return(theta)
}

# SIMULATION MAIN PARAMETERS

# Prior hyperparameters
# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- 0
sigma2_01 <- 100

# theta_02 ~ N(mu_2, sigma2_02)
mu_02     <- 0
sigma2_02 <- 100

# V ~ Gamma(shape=nu_V, rate=eta_V)
nu_V  <- 2
eta_V <- 1000

# phi1 = W1^(-1) ~ Gamma(shape=nu_01, rate=eta_01)
nu_01  <- 2
eta_01 <- 1

# phi2 = W2^(-1) ~ Gamma(shape=nu_02, rate=eta_02)
nu_02  <- 2
eta_02 <- 0.01

# initial values
theta_01 <- 0
theta_02 <- 0
V <- 0.01
W1 <- 0.01
W2 <- 0.01
theta1 <- numeric(Tt)
theta2 <- numeric(Tt)

# DLM main parameters
Ff <- matrix(c(1,0))             # dim = 2x1
G <- rbind(c(1, 1), c(0, 1))    # dim = 2x2
W <- diag(c(W1, W2), nrow=2, ncol=2)

N <- 10000           # Number of steps
burnin <- 1000      # Number of burn-in steps

# Auxiliary vectors and matrix to store the results
theta1_hist <- matrix(nrow=N, ncol=Tt)
theta2_hist <- matrix(nrow=N, ncol=Tt)
V_hist <- numeric(N)
W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <-numeric(N)
theta_02_hist <-numeric(N)

start_time = proc.time() # execution time
for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[3]]
        printf("Iteration %d / %d | Elapsed time: %.0f s", n, N, elapsed_time)
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

    theta <- sample_theta(theta_01, theta_02, Ff, G, W, y, V)
    theta1 <- theta[1,1,]
    theta2 <- theta[2,1,]

    # Sample phi_V
    nu_V_bar <- nu_V + Tt/2
    dif <- y - theta1
    eta_V_bar <- eta_V + 0.5 * sum(dif^2)
    phi_V <- rgamma(1, nu_V_bar, eta_V_bar)
    V <- 1/phi_V

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

    W <- diag(c(W1, W2), nrow=2, ncol=2)


    # Armazenamento
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    V_hist[n] <- V
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_hist[n, ] <- theta1
    theta2_hist[n, ] <- theta2

}

end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]

# Simulation summary ####
# Execution time
printf("Execution time: %.2f s", elapsed_time)


# Posterior mean
theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)

printf("\nV true: %.5f", V_true)
printf("V mean: %.5f", mean(V_hist[-(1:burnin)]))
printf("V median: %.5f", median(V_hist[-(1:burnin)]))

printf("\nW1 true: %.5f", W1_true)
printf("W1 mean: %.5f", mean(W1_hist[-(1:burnin)]))
printf("W1 median: %.5f", median(W1_hist[-(1:burnin)]))

printf("\nW2 true: %.5f", W2_true)
printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))

# Log-likelihood
loglik <- sum(dnorm(y, theta1_mean, log=TRUE))
printf("\nLog-likelihood: %.2f", loglik)

# Effective sample size
printf("\nEffective Sample Size:")
ess_V  <- effectiveSize(mcmc(V_hist[-(1:burnin)]))
ess_w1 <- effectiveSize(mcmc(W1_hist[-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W2_hist[-(1:burnin)]))
ess_theta1 <- effectiveSize(mcmc(theta1_hist[-(1:burnin),]))
ess_theta2 <- effectiveSize(mcmc(theta2_hist[-(1:burnin),]))

printf("\tV: %.0f", ess_V)
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)
printf("\ttheta1 (mean): %.0f", mean(ess_theta1))
printf("\ttheta2 (mean): %.0f", mean(ess_theta2))


# Effective sample size per second
printf("\nEffective Sample Size / second:")
printf("\tV: %.2f", ess_V/elapsed_time)
printf("\tW1: %.2f", ess_w1/elapsed_time)
printf("\tW2: %.2f", ess_w2/elapsed_time)

ess_sec_theta1 <- ess_theta1/elapsed_time
ess_sec_theta2 <- ess_theta2/elapsed_time
printf("\ttheta1 (mean): %.2f", mean(ess_sec_theta1))
printf("\ttheta2 (mean): %.2f", mean(ess_sec_theta2))


# Geweke diagnostic: Z test for two mean difference
#   H0: segments same means -> chain has converged
printf("\nGeweke convergence diagnostic")
z_V <- unname(geweke.diag(V_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
z_W1 <- unname(geweke.diag(W1_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
z_W2 <- unname(geweke.diag(W2_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
printf("\tz_V: %.2f", z_V)
printf("\tz_W1: %.2f", z_W1)
printf("\tz_W2: %.2f", z_W2)

# Percent of instants in the H_0 rejection region:
z_theta1 <- unname(geweke.diag(theta1_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z_theta2 <- unname(geweke.diag(theta2_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z1_out <- sum((z_theta1 < -1.96) | (z_theta1 > 1.96))/Tt
z2_out <- sum((z_theta2 < -1.96) | (z_theta2 > 1.96))/Tt
printf("\tPercent of theta1 out: %.3f", z1_out)
printf("\tPercent of theta2 out: %.3f", z2_out)


# Plots ####
x <- 1:Tt
# theta1_true, theta1_mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_mean, type="l", ylab="", col="blue", lwd=2)
lines(x, theta1_true)
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


# theta2_true, theta2_mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta2_mean, type="l", ylab="", col="blue", lwd=2)
lines(x, theta2_true)
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


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


# Posterior distribution of V ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(V_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of V")
lines(density(V_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W1")
lines(density(W1_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W2")
lines(density(W2_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot of V, W1 and W2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(V_hist, type="l", xlab="n", ylab="", main="Traceplot of V")
abline(v = burnin, col = "red")
plot(V_hist[-(1:burnin)], type="l", xlab="n", ylab="", main="Traceplot of V")

plot(W1_hist, type="l", xlab="n", ylab="", main="Traceplot of W1")
abline(v = burnin, col = "red")
plot(W1_hist[-(1:burnin)], type="l", xlab="n", ylab="", main="Traceplot of W1")

plot(W2_hist, type="l", xlab="n", ylab="", main="Traceplot of W2")
abline(v = burnin, col = "red")
plot(W2_hist[-(1:burnin)], type="l", xlab="n", ylab="", main="Traceplot of W2")


# Traceplots for theta_t1 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
    abline(v = burnin, col = "red")
}
for (t in t_obs) {
    plot(theta1_hist[-(1:burnin), t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
}


# Traceplots for theta_t2 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta2_hist[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
    abline(v = burnin, col = "red")
}
for (t in t_obs) {
    plot(theta2_hist[-(1:burnin), t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
}

# Effective sample size
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ess_theta1, type="l", main=expression("Effective sample of " * theta[t1]), xlab="t")
plot(ess_theta2, type="l", main=expression("Effective sample of " * theta[t2]), xlab="t")

par(mfrow = c(1, 2), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
hist(ess_theta1)
hist(ess_theta2)


# Geweke diagnostic
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(z_theta1, type="l", main=expression("Geweke diagnostic for " * theta[t1]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")

plot(z_theta2, type="l", main=expression("Geweke diagnostic for " * theta[t2]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")
