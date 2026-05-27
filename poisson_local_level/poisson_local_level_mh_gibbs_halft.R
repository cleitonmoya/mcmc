# Poisson Local Level Model
# MCMC: naive Gibbs Sampling
# Author: Cleiton Moya de Almeida

#graphics.off()      # close the plots
rm(list = ls())     # clear the environment
#cat("\014")         # clear the console
set.seed(42)
tp <- Matrix::t     # matrix transpose alias

library(coda)

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "doppler" # rds file with data
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
theta_true <- data$theta
lambda_true <- exp(theta_true)
t_obs <- c(250, 500, 750, 1000)

T <- length(y) # dimension T

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

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

# SIMULATION MAIN PARAMETERS

# Prior hyperparameters
# theta_0 ~ N(mu_0, sigma2_0)
theta_0 <- log(y + 0.5)  # +0.5 to avoid log(0)
mu_0 <- y[1]
sigma2_0 <- 10

nu  <- 1     # nu=1: half-Cauchy; nu=3: half-t(3)
A_W <- 0.2    # calibrar: limite superior plausível para sigma_W

p_halft <- function(sigma, nu, A) 2*pt(sigma/A, df=nu) - 1
#cat("True sigma_W percentil:",
#    p_halft(sqrt(W_true), nu, A_W), "\n")

N <- 10000           # Number of steps
varsigma2 <- 0.01    # Random walkikng variance hyperparameter
burnin <- 1000       # Number of burn-in steps

# Auxiliary vectors and matrix to store the results
theta_hist <- matrix(nrow=N, ncol=T)
W_hist <- numeric(N)
theta_0_hist <-numeric(N)


# Gibbs sampling ####

# Initialization
W <- 0.001
theta <- as.matrix(theta_0)
phi <- 1/W
a_W <- 1   # variável latente auxiliar da mistura

ac_theta_hist <- numeric(N)

# Main loop
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
    diffs  <- theta - c(theta_0, theta[-T])
    SS_W   <- sum(diffs^2)

    # W | theta, a_W ~ Inv-Gamma((nu+T)/2, nu/a_W + SS/2)
    nu_bar  <- (nu + T) / 2
    eta_bar <- nu/a_W + SS_W/2
    phi <- rgamma(1, shape = nu_bar, rate = eta_bar)
    W   <- 1/phi

    # a_W | W ~ Inv-Gamma(nu, nu*(1/A_W^2 + 1/W))
    a_W <- 1/rgamma(1, shape = nu, rate = nu*(1/A_W^2 + 1/W))

    # 3. Sample \theta (sample \theta_t, t=1,...,T)
    n_ac <- 0 # number of accepçted samples for \theta
    for (t in 1:T) {

        # 3.2 sample theta_t (Metropolis)
        if (t < T) {
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
    ac_theta_hist[n] <- n_ac/T

    # Store the sampled values
    theta_0_hist[n] <- theta_0
    W_hist[n] <- W
    theta_hist[n, ] <- theta
}

#### Simulation summary
# Execution time
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)
printf("Mean acception ratio of theta: %.2f", mean(ac_theta_hist))

# Posterior mean
theta_mean <- colMeans(theta_hist[-(1:burnin), ])
W_mean <- mean(W_hist[-(1:burnin)])
W_median <- median(W_hist[-(1:burnin)])
printf("W mean: %.2f", W_mean)
printf("W median: %.2f", W_median)

# Effective sample size ####
printf("Effective Sample Size:")
ess_w <- effectiveSize(mcmc(W_hist[-(1:burnin)]))
printf("\tW: %.0f", ess_w)
printf("Effective Sample Size / second:")
printf("\tW: %.2f", ess_w/elapsed_time)


# Plots ####
# y, lambda_true, lambda_estimated ####
x <- 1:T
lambda_mean <- exp(theta_mean)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson Local Level Polynomial Model")
points(x, y, pch = 20)
lines(x, lambda_mean, col="red", lwd=2)
lines(x, lambda_true, col="blue", lwd=2)
legend("topright",
       legend = expression(y[t], lambda[t], hat(lambda)[t]),
       col = c("black", "blue", "red"),
       lty = c(NA, 1, 1),
       lwd = c(NA, 2, 2),
       pch = c(20, NA, NA),
       bty = "n")


# theta_true vs theta_estimated ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta_true, type="l", xlab="t", ylab="",
     main="Poisson local level model")
lines(x, theta_mean, col="blue", lwd=2)
legend("topright",
       legend = expression(theta[t], hat(theta)[t]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(2, 2),
       bty = "n")

#####
# Posterior distribution of theta_t
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t)]), main = bquote("Posterior of " * theta[.(t)]))
    lines(density(theta_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}

#####
# Posterior distribution of W
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W_hist[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = bquote(theta[.(t)]), main ="Posterior of W")
lines(density(W_hist[-(1:burnin)]), col = "blue", lwd = 2)

#####
# Traceplot for W
plot(W_hist[-(1:100)], type="l", xlab="n", ylab="W", main="Traceplot of W")

#####
# Traceplots for theta_t
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta_hist[, t], type="l", main=bquote(theta[.(t)]), xlab="", ylab="")
}

#####
# Traceplot of acceptance ratio of theta
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_theta_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Acceptance ratio of " * theta))
