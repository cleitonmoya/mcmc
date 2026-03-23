# Poisson Local Level Model
# MCMC: naive Gibbs Sampling
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- Matrix::t     # matrix transpose alias

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "sin_level" # csv file with data
df <- read.table(paste("data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
theta_true <- df$theta1

Tt <- length(y) # dimension T

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

#####
# SIMULATION MAIN PARAMETERS

# Prior hyperparameters
# theta_0 ~ N(mu_0, sigma2_0)
vartheta_init <- log(y + 0.5)  # +0.5 to avoid log(0)
mu_0 <- y[1] # mean(vartheta_init)     # suggested by Claude
sigma2_0 <- 10

# phi = W^(-1) ~ Gamma(nu_0, eta_0)
nu_0 <- 0.01
eta_0 <- 0.01

N <- 5000           # Number of steps
varsigma2 <- 0.07   # Random walkikng variance hyperparameter
burnin <- 200       # Number of burn-in steps

# Auxiliary vectors and matrix to store the results
vartheta_hist <- matrix(nrow=N, ncol=Tt)
W_hist <- numeric(N)
theta_0_hist <-numeric(N)

#####
# Gibbs sampling

# Initialization
W <- 0.3 #var(diff(vartheta_init)) # suggested by Claude
vartheta <- as.matrix(vartheta_init)
phi <- 1/W


ac_vartheta_hist <- numeric(N)

# Main loop
start_time = proc.time() # execution time
for (n in 1:N) {

    # 1. Sample theta_0
    sigma2_0_bar <- (1/sigma2_0 +1/W)^(-1)
    mu_0_bar <- sigma2_0_bar*(mu_0/sigma2_0 + vartheta[1]/W)
    theta_0 <- rnorm(1, mean=mu_0_bar, sd=sqrt(sigma2_0_bar))

    # 2. Sample phi (W^(-1))
    nu_0_bar <- nu_0 + Tt/2
    diffs <- vartheta - c(theta_0, vartheta[-Tt])
    eta_0_bar <- eta_0 + 0.5 * sum(diffs^2)
    phi <- rgamma(1, nu_0_bar, eta_0_bar)
    W <- 1/phi

    # 3. Sample \vartheta (sample \theta_t, t=1,...,T)
    n_ac <- 0 # number of accepçted samples for \vartheta
    for (t in 1:Tt) {

        # 3.2 sample theta_t (Metropolis)
        if (t < Tt) {
            if (t==1) {
                res <- sample_theta_t(vartheta[t], theta_0, vartheta[t+1],
                                      y[t], W, varsigma2, final_t=FALSE)
            } else {
                res <- sample_theta_t(vartheta[t], vartheta[t-1], vartheta[t+1],
                                      y[t], W, varsigma2, final_t=FALSE)
            }

        } else {
            res <- sample_theta_t(vartheta[t], vartheta[t-1], NULL,
                                  y[t], W, varsigma2, final_t=TRUE)
        }

        vartheta[t] <- res$theta_t
        ac <- res$ac # flag: sample accpeted(1) or not (0)
        n_ac <- n_ac + ac
    }

    # Mean acceptance ratio of vartheta
    ac_vartheta_hist[n] <- n_ac/Tt

    # Store the sampled values
    theta_0_hist[n] <- theta_0
    W_hist[n] <- W
    vartheta_hist[n, ] <- vartheta
}

#### Simulation summary
# Execution time
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)

printf("Mean acception ratio of theta: %.2f", mean(ac_vartheta_hist))

# Posterior mean
theta_mean <- colMeans(vartheta_hist[-(1:burnin), ])
lambda_mean <- exp(theta_mean)
W_mean <- mean(W_hist[-(1:burnin)])
W_median <- median(W_hist[-(1:burnin)])
printf("W mean: %.2f", W_mean)
printf("W median: %.2f", W_median)

#####
# Plots
# y, lambda_true, lambda_estimated
lambda_true = exp(theta_true)
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson local level model")
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

#####
# Posterior distribution of theta_t
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    hist(vartheta_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t)]), main = bquote("Posterior of " * theta[.(t)]))
    lines(density(vartheta_hist[-(1:burnin), t]), col = "blue", lwd = 2)
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
for (t in c(10, 50, 100, 150)) {
    plot(vartheta_hist[, t], type="l", main=bquote(theta[.(t)]), xlab="", ylab="")
}

#####
# Traceplot of acceptance ratio of theta
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_vartheta_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Acceptance ratio of " * theta))
