# Poisson - 2nd Order Polynomial Dynamic Model
# MCMC: Metropolis within Gibbs
# Reference: Geweke, J., & Tanizaki, H. (2001).
#    Bayesian estimation of state-space models using the
#    Metropolis–Hastings algorithm within Gibbs sampling.
#    Computational statistics & data analysis, 37(2), 151-170
# Author: Cleiton Moya de Almeida

library(Rfast)      # provide colMedians()

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
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
mu_01     <- y[1]
sigma2_01 <- 10

# theta_02 ~ N(mu_2, sigma2_02)
mu_02     <- 0.1
sigma2_02 <- 1

# phi1 = W1^(-1) ~ Gamma(nu_01, eta_01)
nu_01  <- 0.01
eta_01 <- 0.01

# phi2 = W2^(-1) ~ Gamma(nu_02, eta_02)
nu_02  <- 0.01
eta_02 <- 0.01

N <- 10000           # Number of steps
varsigma2 <- 0.07    # Random walkikng variance hyperparameter
burnin <- 1000       # Number of burn-in steps

# Auxiliary vectors and matrix to store the results
vartheta1_hist <- matrix(nrow=N, ncol=Tt)
vartheta2_hist <- matrix(nrow=N, ncol=Tt)
W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <-numeric(N)
theta_02_hist <-numeric(N)
ac_hist <- numeric(N)

#####
# Gibbs sampling

# Initialization
vartheta1_init <- log(y + 0.5)

library(stats)
vartheta1_smooth <- filter(vartheta1_init, rep(1/5, 5), sides=2)
vartheta1_smooth[is.na(vartheta1_smooth)] <- vartheta1_init[is.na(vartheta1_smooth)]
vartheta2_init <- c(diff(vartheta1_smooth), 0)

#vartheta2_init <- c(diff(vartheta1_init), 0)  # diferenças de theta1
W2 <- var(diff(vartheta2_init)) + 1e-6        # escala compatível
W1 <- var(diff(vartheta1_init))
phi1 <- 1/W1
phi2 <- 1/W2
vartheta1 <- as.matrix(vartheta1_init)
vartheta2 <- as.matrix(vartheta2_init)
theta_01 <- 0.01
theta_02 <- 0.01

 # Main loop
start_time = proc.time() # execution time
for (n in 1:N) {

    # Sample theta_01
    sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
    mu_01_bar <- sigma2_01_bar*(mu_01/sigma2_01 + (vartheta1[1]-theta_02)/W1)
    theta_01 <- rnorm(1, mean=mu_01_bar, sd=sqrt(sigma2_01_bar))

    # Sample theta_02
    sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
    mu_02_bar <- sigma2_02_bar*((vartheta1[1]-theta_01)/W1 +
                                    vartheta2[1]/W2 + mu_02/sigma2_02)
    theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))

    # Sample phi1
    nu_01_bar <- nu_01 + Tt/2
    dif1 <- vartheta1 - c(theta_01, vartheta1[-Tt])
    dif2 <- dif1 - c(theta_02, vartheta2[-Tt])
    eta_01_bar <- eta_01 + 0.5 * sum(dif2^2)
    phi1 <- rgamma(1, nu_01_bar, eta_01_bar)
    W1 <- 1/phi1

    # Sample phi2
    nu_02_bar <- nu_02 + Tt/2
    diffs2 <- vartheta2 - c(theta_02, vartheta2[-Tt])
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
                res <- sample_theta_t1(vartheta1[t], theta_01, vartheta1[t+1],
                                       vartheta2[t], theta_02,
                                       y[t], W1, varsigma2, final_t=FALSE)
                vartheta1[t] <- res$theta_t1
                # theta_t12
                mu_star <- sigma2_star*((vartheta1[t+1] - vartheta1[t])/W1 +
                                            (theta_02 + vartheta2[t+1])/W2)
            } else {

                res <- sample_theta_t1(vartheta1[t], vartheta1[t-1], vartheta1[t+1],
                                       vartheta2[t], vartheta2[t-1],
                                       y[t], W1, varsigma2, final_t=FALSE)
                vartheta1[t] <- res$theta_t1
                mu_star <- sigma2_star*((vartheta1[t+1] - vartheta1[t])/W1 +
                                            (vartheta2[t-1] + vartheta2[t+1])/W2)
            }

        } else {
            res <- sample_theta_t1(vartheta1[t], vartheta1[t-1], vartheta1[t+1],
                                   vartheta2[t], vartheta2[t-1],
                                   y[t], W1, varsigma2, final_t=TRUE)
            vartheta1[t] <- res$theta_t1
            mu_star <- vartheta2[t-1]
            sigma2_star <- W2
        }

        vartheta2[t] <- rnorm(1, mean=mu_star, sd=sqrt(sigma2_star))
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
    vartheta1_hist[n, ] <- vartheta1
    vartheta2_hist[n, ] <- vartheta2
}

#### Simulation summary
# Execution time
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)

printf("Mean acception ratio of theta1: %.2f", mean(ac_hist))

# Posterior mean
theta1_mean <- colMeans(vartheta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(vartheta2_hist[-(1:burnin), ])
theta2_median <- colMedians(vartheta2_hist[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)
printf("W1 mean: %.3f", mean(W1_hist))
printf("W2 mean: %.5f", mean(W2_hist))
printf("W2 median: %.5f", median(W2_hist))

#####
# Plots
# y, lambda_true, lambda_estimated
lambda_true <- exp(theta1_true)
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson 2nd Order Polynomial Model")
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


# y, theta2_true, theta2_mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta2_mean, type="l")
lines(x, theta2_true, col="blue", lwd="2")
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


#####
# Posterior distribution of theta_t1
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    hist(vartheta1_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Posterior of " * theta[.(t) * "," * 1]))
    lines(density(vartheta1_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}

# Posterior distribution of theta_t2
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    hist(vartheta2_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 2]),
         main = bquote("Posterior of " * theta[.(t) * "," * 2]))
    lines(density(vartheta2_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


#####
# Posterior distribution of W1
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W1")
lines(density(W1_hist[-(1:burnin)]), col = "blue", lwd = 2)

# Posterior distribution of W2
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W2")
lines(density(W2_hist[-(1:burnin)]), col = "blue", lwd = 2)

#####
# Traceplot for W1 and W2
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W1_hist[-(1:100)], type="l", xlab="n", ylab="W", main="Traceplot of W1")
plot(W2_hist[-(1:100)], type="l", xlab="n", ylab="W", main="Traceplot of W2")

#####
# Traceplots for theta_t1
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    plot(vartheta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
}

# Traceplots for theta_t2
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    plot(vartheta2_hist[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
}


#####
# Traceplot of acceptance ratio of theta_t1
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Acceptance ratio of " * theta[t*1]))
