# 2nd Order Polynomial Dynamic Linear Model (Local Trend DLM)
#
# Model:
#  y_t = theta_{t1} + nu_t, nu_t ~ N(0,V)
#  theta_{t1} = theta_{t-1,1} + theta_{t-1,2} + omega_{t1}, omega_{t1} ~ N(0, W1)
#  theta_{t2} = theta_{t-1,2} + omega_{t2}, omega_{t2} ~ N(0, W2)
#
# Priors:
#  theta_{01} | D_0 ~ N(mu_{01}, sigma2_{01})
#  theta_{02} | D_0 ~ N(mu_{02}, sigma2_{02})
#  1/V  | D_0 ~ gamma(shape=nu_V, rate=eta_V)
#  1/W1 | D_0 ~ gamma(shape=nu_01, rate=eta_01)
#  1/W2 | D_0 ~ gamma(shape=nu_02, rate=eta_02)
#
# MCMC:
#  FFBS within Gibbs
#
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
V_true <- data$V
W1_true <- data$W1
W2_true <- data$W2


# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))

inv2x2 <- function(M) {
    det_M <- M[1,1]*M[2,2] - M[1,2]*M[2,1]
    inv <- matrix(c(M[2,2], -M[2,1], -M[1,2], M[1,1]), 2, 2) / det_M
    return(inv)
}


# Sample from a bivariate normal distribution
rbvn_cond <- function(mu1, mu2, s11, s12, s22) {
    x1 <- rnorm(1, mu1, sqrt(s11))
    mu2_c <- mu2 + (s12/s11)*(x1 - mu1)
    s22_c <- s22 - s12^2/s11
    x2 <- rnorm(1, mu2_c, sqrt(max(s22_c, 0)))
    c(x1, x2)
}


# Forward Filtering (Kalman Filter)
forward_filter <- function(theta_01, theta_02, W1, W2, y, V) {
    m1 <- theta_01; m2 <- theta_02
    c11 <- 0; c12 <- 0; c22 <- 0

    a1v <- a2v <- R11v <- R12v <- R22v <- m1v <- m2v <- c11v <- c12v <- c22v <- numeric(Tt)

    for (t in 1:Tt) {
        a1 <- m1 + m2
        a2 <- m2
        R11 <- c11 + 2*c12 + c22 + W1
        R12 <- c12 + c22
        R22 <- c22 + W2

        Q <- R11 + V
        e <- y[t] - a1

        m1 <- a1 + (R11/Q)*e
        m2 <- a2 + (R12/Q)*e
        c11 <- R11*V/Q
        c12 <- R12*V/Q
        c22 <- R22 - R12^2/Q

        a1v[t]<-a1; a2v[t]<-a2; R11v[t]<-R11; R12v[t]<-R12; R22v[t]<-R22
        m1v[t]<-m1; m2v[t]<-m2; c11v[t]<-c11; c12v[t]<-c12; c22v[t]<-c22
    }
    list(a1=a1v,a2=a2v,R11=R11v,R12=R12v,R22=R22v,m1=m1v,m2=m2v,c11=c11v,c12=c12v,c22=c22v)
}


# FFBS
sample_theta_ffbs <- function(theta_01, theta_02, W1, W2, y, V) {
    f <- forward_filter(theta_01, theta_02, W1, W2, y, V)
    theta1 <- numeric(Tt); theta2 <- numeric(Tt)

    draw <- rbvn_cond(f$m1[Tt], f$m2[Tt], f$c11[Tt], f$c12[Tt], f$c22[Tt])
    theta1[Tt] <- draw[1]; theta2[Tt] <- draw[2]

    for (t in seq(Tt-1, 1)) {
        c11<-f$c11[t]; c12<-f$c12[t]; c22<-f$c22[t]
        R11p<-f$R11[t+1]; R12p<-f$R12[t+1]; R22p<-f$R22[t+1]
        det <- R11p*R22p - R12p^2

        B11 <- ((c11+c12)*R22p - c12*R12p)/det
        B12 <- (-(c11+c12)*R12p + c12*R11p)/det
        B21 <- ((c12+c22)*R22p - c22*R12p)/det
        B22 <- (-(c12+c22)*R12p + c22*R11p)/det

        d1 <- theta1[t+1] - f$a1[t+1]
        d2 <- theta2[t+1] - f$a2[t+1]

        h1 <- f$m1[t] + B11*d1 + B12*d2
        h2 <- f$m2[t] + B21*d1 + B22*d2

        H11 <- c11 - (B11*R11p+B12*R12p)*B11 - (B11*R12p+B12*R22p)*B12
        H12 <- c12 - (B11*R11p+B12*R12p)*B21 - (B11*R12p+B12*R22p)*B22
        H22 <- c22 - (B21*R11p+B22*R12p)*B21 - (B21*R12p+B22*R22p)*B22

        draw <- rbvn_cond(h1, h2, H11, H12, H22)
        theta1[t] <- draw[1]; theta2[t] <- draw[2]
    }
    list(theta1=theta1, theta2=theta2)
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

N <- 20000           # Number of steps
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

    # Sample theta - FFBS
    theta <- sample_theta_ffbs(theta_01, theta_02, W1, W2, y, V)
    theta1 <- theta$theta1
    theta2 <- theta$theta2

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

    # Store values
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
