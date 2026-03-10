# Poisson - 2nd Order Polynomial Dynamic Model
# MCMC: Precision Based & Component Wise Metropolis
# Reference:
# Montoril, Michel H., Leandro T. Correia, e Helio S. Migon.
#     “Bayesian Estimation of Dynamic Weights in Gaussian Mixture Models”.
#      arXiv:2104.03395. Preprint, arXiv, 2022.
#      https://doi.org/10.48550/arXiv.2104.03395.
# Author: Cleiton Moya de Almeida

library(Matrix)     # to deal with band sparse matrix

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(20)
tp <- Matrix::t     # matrix transpose alias

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "pol2_sim1" # csv file with data
df <- read.table(paste("data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
theta1_true <- df$theta1
theta2_true <- df$theta2
Tt <- length(y) # dimension T

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

# H band matrix (TxT)
H <- bandSparse(Tt, Tt, k = c(0, -1),
                diagonals = list(rep(1, Tt), rep(-1, Tt-1)))

# Full conditional posterior for theta_t
logpost_theta_t <- function(theta_t, yt, mu_t_star, tau2_t_star) {
    p1 <- yt*theta_t - exp(theta_t) # log-likelihood
    p2 <- -(theta_t-mu_t_star)^2/(2*tau2_t_star)
    p <- p1 + p2
    return(p)
}


# Sample theta_t ~ logpost_theta_t (Metropolis step)
sample_theta_t <- function(theta_t_current, mu_t_star, tau2_t_star,
                           varsigma2, yt) {

    # proposed theta
    theta_t_prop <- rnorm(1, mean=theta_t_current, sd=sqrt(varsigma2))

    # acceptance/rejection step
    ac <- 0 # accepted flag
    logu <- log(runif(1))
    logp1 <- logpost_theta_t(theta_t_prop, yt, mu_t_star, tau2_t_star)
    logp2 <- logpost_theta_t(theta_t_current, yt, mu_t_star, tau2_t_star)
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
mu_0 <- mean(vartheta_init)    # suggested by Claude
sigma2_0 <- 3

# phi = W^(-1) ~ Gamma(nu_0, eta_0)
nu_0 <- 0.01
eta_0 <- 0.01

N <- 3000           # Number of steps
varsigma2 <- 0.07   # Random walkikng variance hyperparameter
burnin <- 200       # Number of burn-in steps

# Auxiliary vectors and matrix to store the results
vartheta_hist <- matrix(nrow=N, ncol=Tt)
W_hist <- numeric(N)
theta_0_hist <-numeric(N)

#####
# Gibbs sampling

# Initialization
W <- var(diff(vartheta_init)) # suggested by Claude
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
    mu <- as.matrix(theta_0*rep(1, Tt))

    # 2. Sample phi (W)
    nu_0_bar <- nu_0 + Tt/2

    diff <- vartheta - mu
    eta_0_bar <- eta_0 + 0.5 * as.numeric(crossprod(H %*% diff))

    phi <- rgamma(1, nu_0_bar, eta_0_bar)
    W <- 1/phi

    # 3. Sample \vartheta (sample \theta_t, t=1,...,T)
    n_ac <- 0 # number of accepçted samples for \vartheta
    for (t in 1:Tt) {

        # 3.1 update mu_t_star and tau2_t_star
        if (t<Tt) {
            if (t > 1) {
                mu_t_star <- (vartheta[t+1] + vartheta[t-1])/2
            } else {
                mu_t_star <- (vartheta[2] + theta_0)/2
            }
            tau2_t_star <- W/2
        } else {
            mu_t_star <- vartheta[t-1]
            tau2_t_star <- W
        }

        # 3.2 sample theta_t (Metropolis)
        res <- sample_theta_t(vartheta[t], mu_t_star, tau2_t_star, varsigma2, y[t])
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
printf("W mean: %.2f", mean(W_hist))


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
plot(W_hist, type="l", xlab="n", ylab="W", main="Traceplot of W")

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
