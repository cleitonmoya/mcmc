# Poisson - 2nd Order Polynomial Dynamic Model
# MCMC: Metropolis within Gibbs
# Reference: Geweke, J., & Tanizaki, H. (2001).
#    Bayesian estimation of state-space models using the
#    Metropolis–Hastings algorithm within Gibbs sampling.
#    Computational statistics & data analysis, 37(2), 151-170
# Author: Cleiton Moya de Almeida

library(Rfast)      # provide colMedians()
library(coda)
library(stats)

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- Matrix::t     # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "poisson_pol2_sim1" # csv file with data
df <- read.table(paste("../data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
theta1_true <- df$theta1
theta2_true <- df$theta2
W1_true <- 0.03
W2_true <- 0.00003

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

p_halft <- function(sigma, nu, A) {
    pt(sigma/A, df = nu) * 2 - 1  # half-t CDF
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

# \sigma_j ~ half-t(nu, A_j)
nu   <- 1      # degree of freedom: nu=1 -> half-Cauchy
A_W1 <- 0.5    # escala para sigma_W1: calibrar via pexp ou análogo
A_W2 <- 0.005  # escala para sigma_W2
cat("Percentil sigma_W1:", p_halft(sqrt(W1_true), nu=nu, A=A_W1), "\n")
cat("Percentil sigma_W2:", p_halft(sqrt(W2_true), nu=nu, A=A_W2), "\n")


N <- 11000           # Number of steps
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


vartheta1_smooth <- filter(vartheta1_init, rep(1/5, 5), sides=2)
vartheta1_smooth[is.na(vartheta1_smooth)] <- vartheta1_init[is.na(vartheta1_smooth)]
vartheta2_init <- c(diff(vartheta1_smooth), 0)

#vartheta2_init <- c(diff(vartheta1_init), 0) # diferenças de theta1
W2 <- var(diff(vartheta2_init)) + 1e-6        # escala compatível
W1 <- var(diff(vartheta1_init))
phi1 <- 1/W1
phi2 <- 1/W2
vartheta1 <- as.matrix(vartheta1_init)
vartheta2 <- as.matrix(vartheta2_init)
theta_01 <- 0.01
theta_02 <- 0.01
a_W1 <- 1
a_W2 <- 1

 # Main loop
start_time = proc.time() # execution time
for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[3]]
        printf("Iteration %d / %d | Elapsed time: %.0f s", n, N, elapsed_time)
    }

    # Sample theta_01
    sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
    mu_01_bar <- sigma2_01_bar*(mu_01/sigma2_01 + (vartheta1[1]-theta_02)/W1)
    theta_01 <- rnorm(1, mean=mu_01_bar, sd=sqrt(sigma2_01_bar))

    # Sample theta_02
    sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
    mu_02_bar <- sigma2_02_bar*((vartheta1[1]-theta_01)/W1 +
                                    vartheta2[1]/W2 + mu_02/sigma2_02)
    theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))

    # Sample W1 | theta, a_W1
    nu_bar_W1   <- (nu + Tt) / 2
    dif1        <- vartheta1 - c(theta_01, vartheta1[-Tt])
    dif2        <- dif1 - c(theta_02, vartheta2[-Tt])
    SS_W1       <- sum(dif2^2)
    eta_bar_W1  <- nu / a_W1 + SS_W1 / 2
    phi1        <- rgamma(1, shape = nu_bar_W1, rate = eta_bar_W1)
    W1          <- 1 / phi1

    # Sample a_W1 | W1
    a_W1 <- 1 / rgamma(1, shape = nu, rate = nu * (1/A_W1^2 + 1/W1))

    # Sample W2 | theta, a_W2
    nu_bar_W2   <- (nu + Tt) / 2
    diffs2      <- vartheta2 - c(theta_02, vartheta2[-Tt])
    SS_W2       <- sum(diffs2^2)
    eta_bar_W2  <- nu / a_W2 + SS_W2 / 2
    phi2        <- rgamma(1, shape = nu_bar_W2, rate = eta_bar_W2)
    W2          <- 1 / phi2

    # Sample a_W2 | W2
    a_W2 <- 1 / rgamma(1, shape = nu, rate = nu * (1/A_W2^2 + 1/W2))

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
            res <- sample_theta_t1(vartheta1[t], vartheta1[t-1], NULL,
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
printf("W1 mean: %.3f", mean(W1_hist[-(1:burnin)]))
printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))


# Effective sample size ####
printf("Effective Sample Size:")
ess_w1 <- effectiveSize(mcmc(W1_hist[-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W2_hist[-(1:burnin)]))
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)

observed_times <- c(50, 100, 200, 250)

for (t in observed_times) {
    ess <- effectiveSize(mcmc(vartheta1_hist[-(1:burnin),t]))
    printf("\ttheta %d,1: %0.f", t, ess)
}

for (t in observed_times) {
    ess <- effectiveSize(mcmc(vartheta2_hist[-(1:burnin),t]))
    printf("\ttheta %d,2: %0.f", t, ess)
}

# Effective sample size / elapsed time ####
printf("Effective Sample Size / second:")
printf("\tW1: %.2f", ess_w1/elapsed_time)
printf("\tW2: %.2f", ess_w2/elapsed_time)


for (t in observed_times) {
    ess <- effectiveSize(mcmc(vartheta1_hist[-(1:burnin),t]))
    printf("\ttheta %d,1: %.2f", t, ess/elapsed_time)
}

for (t in observed_times) {
    ess <- effectiveSize(mcmc(vartheta2_hist[-(1:burnin),t]))
    printf("\ttheta %d,2: %.2f", t, ess/elapsed_time)
}


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


# y, theta1_true, theta1_mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_mean, type="l", ylab="", col="blue", lwd=2)
lines(x, theta1_true)
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


# y, theta2_true, theta2_mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta2_mean, type="l", ylab="", col="blue", lwd=2, ylim=c(-0.05,0.1))
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
    hist(vartheta1_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Posterior of " * theta[.(t) * "," * 1]))
    lines(density(vartheta1_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of theta_t2 ####
par(mfrow = c(2, 2))
for (t in observed_times) {
    hist(vartheta2_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 2]),
         main = bquote("Posterior of " * theta[.(t) * "," * 2]))
    lines(density(vartheta2_hist[-(1:burnin), t]), col = "blue", lwd = 2)
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
plot(W1_hist[-(1:100)], type="l", xlab="n", ylab="W", main="Traceplot of W1")
plot(W2_hist[-(1:100)], type="l", xlab="n", ylab="W", main="Traceplot of W2")


# Traceplots for theta_t1 ####
par(mfrow = c(2, 2))
for (t in observed_times) {
    plot(vartheta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
}


# Traceplots for theta_t2 ####
par(mfrow = c(2, 2))
for (t in observed_times) {
    plot(vartheta2_hist[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
}


# Traceplot of acceptance ratio of theta_t1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Acceptance ratio of " * theta[t*1]))
