# 2nd Order Polynomial Dynamic Linear Model
# MCMC: Gibbs with FFBS
# Author: Cleiton Moya de Almeida

library(coda)

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- Matrix::t     # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "normal_pol2_sim1" # csv file with data
df <- read.table(paste("../data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
T <- length(y) # dimension T
theta1_true <- df$theta1
theta2_true <- df$theta2
V_true      <- df$V[1]
W1_true     <- df$W1[1]
W2_true     <- df$W2[1]


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

# Forward Filtering (Kalman Filter)
forward_filter <- function(mu_01, sigma2_01, mu_02, sigma2_02,
                           F, G, W, y, V) {

    m_t <- as.matrix(c(mu_01, mu_02))
    C_t <- diag(c(sigma2_01, sigma2_02), nrow=2, ncol=2)

    a <- array(data=NA, dim=c(2,1,T))
    R <- array(data=NA, dim=c(2,2,T))
    m <- array(data=NA, dim=c(2,1,T))
    C <- array(data=NA, dim=c(2,2,T))


    for (t in 1:T) {

        # Prior: (theta_t | D_{t-1}) ~ N[a_t, R_t]
        a_t <- G %*% m_t
        R_t <- G %*% C_t %*% tp(G) + W

        a[,,t] <- a_t
        R[,,t] <- R_t

        # One-step-ahead forecast
        f_t <- tp(F) %*% a_t
        Q_t <- drop(tp(F) %*% R_t %*% F + V)

        # Posterior: (theta_t | D_t) ~ N[m_t, C_t]
        e_t <- drop(y[t] - f_t)
        m_t <- a_t + (1/Q_t) * R_t %*% F * e_t
        C_t <- R_t - (1/Q_t) * R_t %*% F %*% tp(F) %*% R_t
        m[,,t] <- m_t
        C[,,t] <- C_t
    }

    return(list("a"=a, "R"=R, "m"=m, "C"=C))

}


# FFBS
sample_theta <- function(mu_01, sigma2_01, mu_02, sigma2_02,
                         F, G, W, y, V) {

    # theta = {theta_1, ..., theta_T}
    theta <- array(data=NA, dim=c(2,1,T))
    res <- forward_filter(mu_01, sigma2_01, mu_02, sigma2_02,
                          F, G, W, y, V)

    a <- res$a
    R <- res$R
    m <- res$m
    C <- res$C

    # Compute B_t
    B <- array(data=NA, dim=c(2,2,T))

    for (t in 1:(T-1)) {
        B[,,t] <- C[,,t] %*% tp(G) %*% inv2x2(R[,,t+1])
    }

    # Sample theta_T ~ N[m_T, C_T]
    theta_t <- rmvn_chol(mu=m[,,T], Sigma=C[,,T])
    theta[,,T] <- theta_t

    # Backward: t = T-1, ..., 1
    # theta_t | theta_{t+1} ~ N(h_t, H_t) [
    # West and Harrison 1997, eq. 15.8, p. 570]
    for (t in seq(T-1, 1)) {
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
mu_01     <- 0.01
sigma2_01 <- 10
theta_01 <- y[1] # initial value

# theta_02 ~ N(mu_2, sigma2_02)
mu_02     <- 0.01
sigma2_02 <- 10
theta_02 <- y[2]-y[1] # initial value

# V ~ Gamma(nu_V, eta_V)
nu_V  <- 1
eta_V <- 1
V <- 0.1 # initial value

# phi1 = W1^(-1) ~ Gamma(nu_01, eta_01)
nu_01  <- 1
eta_01 <- 1
W1 <- 1

# phi2 = W2^(-1) ~ Gamma(nu_02, eta_02)
nu_02  <- 1
eta_02 <- 1
W2 <- 1

# initial values for theta_t1 and theta_t2
theta1 <- y
theta2 <- numeric(T)

# DLM main parameters
F <- matrix(c(1,0))             # dim = 2x1
G <- rbind(c(1, 1), c(0, 1))    # dim = 2x2
W <- diag(c(W1, W2), nrow=2, ncol=2)

N <- 3000           # Number of steps
burnin <- 500       # Number of burn-in steps

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


    theta <- sample_theta(mu_01, sigma2_01, mu_02, sigma2_02, F, G, W, y, V)
    theta1 <- theta[1,1,]
    theta2 <- theta[2,1,]

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

    W <- diag(c(W1, W2))

    # Armazenamento
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
sink("../summary/pol2_ffbs.txt", split = TRUE)
printf("Execution time: %.0f s", elapsed_time)


# Posterior mean
theta1_mean <- colMeans(theta1_samples[-(1:burnin), ])
theta2_mean <- colMeans(theta2_samples[-(1:burnin), ])

printf("V mean: %.2f", mean(V_samples[-(1:burnin)]))
printf("V median: %.2f", median(V_samples[-(1:burnin)]))
printf("V true: %.2f", V_true)
printf("W1 mean: %.3f", mean(W1_samples[-(1:burnin)]))
printf("W1 median: %.5f", median(W1_samples[-(1:burnin)]))
printf("W1 true: %.2f", W1_true)
printf("W2 mean: %.5f", mean(W2_samples[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_samples[-(1:burnin)]))
printf("W2 true: %.5f", W2_true)


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
for (t in observed_times) {
    hist(theta1_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Posterior of " * theta[.(t) * "," * 1]))
    lines(density(theta1_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of theta_t2 ####
par(mfrow = c(2, 2))
for (t in observed_times) {
    hist(theta2_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 2]),
         main = bquote("Posterior of " * theta[.(t) * "," * 2]))
    lines(density(theta2_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of V ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(V_samples[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of V")
lines(density(V_samples[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_samples[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W1")
lines(density(W1_samples[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_samples[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W2")
lines(density(W2_samples[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot of V, W1 and W2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(V_samples[-(1:100)], type="l", xlab="n", ylab="", main="Traceplot of V")
abline(v = burnin, col = "red")
plot(W1_samples[-(1:100)], type="l", xlab="n", ylab="", main="Traceplot of W1")
abline(v = burnin, col = "red")
plot(W2_samples[-(1:100)], type="l", xlab="n", ylab="", main="Traceplot of W2")
abline(v = burnin, col = "red")

# Traceplots for theta_t1 ####
par(mfrow = c(2, 2))
for (t in observed_times) {
    plot(theta1_samples[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
    abline(v = burnin, col = "red")
}


# Traceplots for theta_t2 ####
par(mfrow = c(2, 2))
for (t in observed_times) {
    plot(theta2_samples[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
    abline(v = burnin, col = "red")
}
