# 2nd Order Polynomial Dynamic Linear Model
# MCMC: Gibbs with Chan Method
# Author: Cleiton Moya de Almeida

library(Matrix)
library(coda)

graphics.off()      # close the plots
#rm(list = ls())     # clear the environment
#cat("\014")         # clear the console
set.seed(42)
tp <- Matrix::t     # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "normal_pol2_sim1" # csv file with data
df <- read.table(paste("../data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
theta1_true <- df$theta1
theta2_true <- df$theta2

T <- length(y) # dimension T


printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}


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


# Montoril Method
# Matriz H: T x T bidiagonal
H   <- bandSparse(T, T,
                  k         = c(0, -1),
                  diagonals = list(rep(1, T), rep(-1, T - 1)))
B   <- Diagonal(T) - H

HtH <- tp(H) %*% H
BtB <- tp(B) %*% B
BtH <- tp(B) %*% H
HtB <- tp(H) %*% B    # = t(BtH)

e1 <- c(1, rep(0, T - 1))

sample_theta_montoril <- function(y, V, W1, W2, theta_01, theta_02,
                                  vartheta1_prev) {

    # ---- Amostrar vartheta_2 ----
    Phi2_bar <- (1/W1)*BtB + (1/W2)*HtH
    rhs2     <- (1/W1)*(BtH %*% vartheta1_prev) + (theta_02/W2)*e1

    ch2       <- Cholesky(Phi2_bar, LDL = FALSE, perm = FALSE)
    mu2_bar   <- solve(ch2, rhs2)
    u2        <- rnorm(T)
    x2        <- as.vector(solve(ch2, u2, system = "Lt"))
    vartheta2 <- as.vector(mu2_bar) + x2

    # ---- Amostrar vartheta_1 | vartheta_2 ----
    Phi1_bar <- (1/V)*Diagonal(T) + (1/W1)*HtH
    rhs1     <- (1/V)*y +
        ((theta_01 + theta_02)/W1)*e1 +
        (1/W1)*(HtB %*% vartheta2)

    ch1       <- Cholesky(Phi1_bar, LDL = FALSE, perm = FALSE)
    mu1_bar   <- solve(ch1, rhs1)
    u1        <- rnorm(T)
    x1        <- as.vector(solve(ch1, u1, system = "Lt"))
    vartheta1 <- as.vector(mu1_bar) + x1

    return(list(theta1 = vartheta1, theta2 = vartheta2))
}



# SIMULATION MAIN PARAMETERS

# Prior hyperparameters
# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- 0.01
sigma2_01 <- 10
theta_01 <- y[1] # initial value

# theta_02 ~ N(mu_2, sigma2_02)
mu_02     <- 0.01
sigma2_02 <- 10
theta_02 <- y[2]-y[1] # initial value

# V ~ Gamma(nu_V, eta_V)
nu_V  <- 0.01
eta_V <- 0.01
V <- 0.1

# phi1 = W1^(-1) ~ Gamma(nu_01, eta_01)
nu_01  <- 0.01
eta_01 <- 0.01
W1 <- 1

# phi2 = W2^(-1) ~ Gamma(nu_02, eta_02)
nu_02  <- 0.01
eta_02 <- 0.01
W2 <- 1

# initial values for theta_t1 and theta_t2
theta1 <- y
theta2 <- numeric(T)

# DLM main parameters
F <- matrix(c(1,0))             # dim = 2x1
G <- rbind(c(1, 1), c(0, 1))    # dim = 2x2
#W <- diag(c(W1, W2), nrow=2, ncol=2)

N <- 5000           # Number of steps
burnin <- 1000      # Number of burn-in steps

# Auxiliary vectors and matrix to store the results
theta1_samples <- matrix(nrow=N, ncol=T)
theta2_samples <- matrix(nrow=N, ncol=T)
V_samples <- numeric(N)
W1_samples <- numeric(N)
W2_samples <- numeric(N)
theta_01_samples <-numeric(N)
theta_02_samples <-numeric(N)


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


    # Sample theta - Montoril Method
    res    <- sample_theta_montoril(y, V, W1, W2, theta_01, theta_02,
                                    vartheta1_prev = theta1)
    theta1 <- res$theta1
    theta2 <- res$theta2

    # Sample phi_V
    nu_V_bar <- nu_V + T/2
    dif <- y - theta1
    eta_V_bar <- eta_V + 0.5 * sum(dif^2)
    phi_V <- rgamma(1, nu_V_bar, eta_V_bar)
    V <- 1/phi_V

    # Sample phi1
    nu_01_bar <- nu_01 + T/2
    dif1 <- theta1 - c(theta_01, theta1[-T])
    dif2 <- dif1 - c(theta_02, theta2[-T])
    eta_01_bar <- eta_01 + 0.5 * sum(dif2^2)
    phi1 <- rgamma(1, nu_01_bar, eta_01_bar)
    W1 <- 1/phi1

    # Sample phi2
    nu_02_bar <- nu_02 + T/2
    diffs2 <- theta2 - c(theta_02, theta2[-T])
    eta_02_bar <- eta_02 + 0.5 * sum(diffs2^2)
    phi2 <- rgamma(1, nu_02_bar, eta_02_bar)
    W2 <- 1/phi2

    # Store the sampled values
    theta_01_samples[n] <- theta_01
    theta_02_samples[n] <- theta_02
    V_samples[n] <- V
    W1_samples[n] <- W1
    W2_samples[n] <- W2
    theta1_samples[n, ] <- theta1
    theta2_samples[n, ] <- theta2

}
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]

# Simulation summary ####
# Execution time
sink("../summary/pol2_montoril.txt", split = TRUE)
printf("Execution time: %.0f s", elapsed_time)


# Posterior mean
theta1_mean <- colMeans(theta1_samples[-(1:burnin), ])
theta2_mean <- colMeans(theta2_samples[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)
printf("V mean: %.2f", mean(V_samples[-(1:burnin)]))
printf("W1 mean: %.3f", mean(W1_samples[-(1:burnin)]))
printf("W1 median: %.5f", median(W1_samples[-(1:burnin)]))
printf("W2 mean: %.5f", mean(W2_samples[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_samples[-(1:burnin)]))


# Effective sample size
printf("Effective Sample Size:")
ess_V  <- effectiveSize(mcmc(V_samples[-(1:burnin)]))
ess_w1 <- effectiveSize(mcmc(W1_samples[-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W2_samples[-(1:burnin)]))
printf("\tV:  %.0f", ess_V)
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)

observed_times <- c(50, 100, 150, 250)
for (t in observed_times) {
    ess <- effectiveSize(mcmc(theta1_samples[-(1:burnin),t]))
    printf("\ttheta %d,1: %0.f", t, ess)
}

for (t in observed_times) {
    ess <- effectiveSize(mcmc(theta2_samples[-(1:burnin),t]))
    printf("\ttheta %d,2: %0.f", t, ess)
}

# Effective sample size / elapsed time
printf("Effective Sample Size / second:")
printf("\tV1: %.2f", ess_V/elapsed_time)
printf("\tW1: %.2f", ess_w1/elapsed_time)
printf("\tW2: %.2f", ess_w2/elapsed_time)

for (t in observed_times) {
    ess <- effectiveSize(mcmc(theta1_samples[-(1:burnin),t]))
    printf("\ttheta %d,1: %.2f", t, ess/elapsed_time)
}

for (t in observed_times) {
    ess <- effectiveSize(mcmc(theta2_samples[-(1:burnin),t]))
    printf("\ttheta %d,2: %.2f", t, ess/elapsed_time)
}
sink()


# Plots ####

x <- 1:T
# theta1_true, theta1_mean
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_mean, type="l", ylab="", col="blue", lwd=2)
lines(x, theta1_true)
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


# theta2_true, theta2_mean
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta2_mean, type="l", ylab="", col="blue", lwd=2)
lines(x, theta2_true)
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


# Posterior distribution of theta_t1
par(mfrow = c(2, 2))
for (t in observed_times) {
    hist(theta1_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Posterior of " * theta[.(t) * "," * 1]))
    lines(density(theta1_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of theta_t2
par(mfrow = c(2, 2))
for (t in observed_times) {
    hist(theta2_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 2]),
         main = bquote("Posterior of " * theta[.(t) * "," * 2]))
    lines(density(theta2_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}

# Posterior distribution of V
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(V_samples[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of V")
lines(density(V_samples[-(1:burnin)]), col = "blue", lwd = 2)

# Posterior distribution of W1
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_samples[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W1")
lines(density(W1_samples[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_samples[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W2")
lines(density(W2_samples[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot of V, W1 and W2
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(V_samples[-(1:100)], type="l", xlab="n", ylab="", main="Traceplot of V")
plot(W1_samples[-(1:100)], type="l", xlab="n", ylab="", main="Traceplot of W1")
plot(W2_samples[-(1:100)], type="l", xlab="n", ylab="", main="Traceplot of W2")


# Traceplots for theta_t1
par(mfrow = c(2, 2))
for (t in observed_times) {
    plot(theta1_samples[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
}


# Traceplots for theta_t2
par(mfrow = c(2, 2))
for (t in observed_times) {
    plot(theta2_samples[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
}
