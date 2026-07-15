# Poisson - 2nd Order Polynomial Dynamic Model
# Gibbs (PG) Importance Sampling + Component-wise Gibbs for theta_t2
# Strategy: IS for theta_t1, theta_t2 sampled via
#           component-wise Gibbs (Normal conjugate full conditionals)
# Author: Cleiton Moya de Almeida

library(invgamma)
library(coda)

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
#cat("\014")         # clear the console
options(error = function() traceback(2))
tp <- base::t       # alias to transpose function
set.seed(42)

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source_file <- "poisson_pol2_200"
data <- readRDS(paste("../data/", source_file, ".rds", sep = ""))
y <- data$y

Tt <- length(y)
if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 2000) t_obs <- c(500, 1000, 1500, 1750)

theta1_present <- TRUE
theta2_present <- FALSE

if (theta1_present) {
    theta1_true <- data$theta
    lambda_true <- exp(theta1_true)
}

if (theta2_present)
    theta2_true <- data$theta2

printf <- function(...) cat(paste(sprintf(...), "\n"))

logsumexp <- function(x) {
    cc <- max(x)
    return(cc + log(sum(exp(x - cc))))
}

log_p_yt <- function(yt, theta_t1) {
    res <- yt * theta_t1 - exp(theta_t1)
    res[!is.finite(res)] <- -Inf
    return(res)
}


forward_filter_1d <- function(z, f, theta2, theta_01, theta_02, W1) {
    a <- numeric(Tt)
    R <- numeric(Tt)
    m <- numeric(Tt)
    C <- numeric(Tt)

    drift <- c(theta_02, theta2[-Tt])  # pré-computado, elimina o if
    m_t <- theta_01;
    C_t <- 0;

    for (t in 1:Tt) {
        a_t <- m_t + drift[t]
        R_t <- C_t + W1
        Q_t <- R_t + f[t]
        A_t <- R_t / Q_t
        m_t <- a_t + A_t * (z[t] - a_t)
        C_t <- R_t * (1 - A_t)

        a[t] <- a_t
        R[t] <- R_t
        m[t] <- m_t
        C[t] <- C_t
    }
    return(list(a=a, R=R, m=m, C=C))
}


kalman_smoother_1d <- function(kf, W1) {
    m_s <- kf$m
    B   <- kf$C[-Tt] / kf$R[-1]
    for (t in seq(Tt-1, 1)) {
        m_s[t] <- kf$m[t] + B[t] * (m_s[t+1] - kf$a[t+1])
    }
    return(list(m_s=m_s, B=B))
}

ffbs_1d <- function(kf, ks, W1) {
    theta1 <- numeric(Tt)
    theta1[Tt] <- rnorm(1, kf$m[Tt], sqrt(kf$C[Tt]))
    B <- ks$B

    # pré-computar fora do loop
    H   <- kf$C[-Tt] - B^2 * kf$R[-1]
    sH  <- sqrt(H)
    mC  <- kf$m[-Tt]
    aC1 <- kf$a[-1]

    for (t in seq(Tt-1, 1)) {
        theta1[t] <- rnorm(1, mC[t] + B[t]*(theta1[t+1] - aC1[t]), sH[t])
    }
    return(theta1)
}


# Prior Hyperparameters

# theta_01 ~ N(mu_01, sigma2_01)
mu_01 <- log(y[1] + 0.5)
sigma2_01 <- 100

# theta_02 ~ N(mu_02, sigma2_02)
mu_02 <- 0
sigma2_02 <- 100

# W1 ~ InvGamma(alpha_W1, beta_W1)
alpha_W1 <- 2
beta_W1 <- 0.01

# W2 ~ InvGamma(alpha_W2, beta_W2)
alpha_W2 <- 2
beta_W2 <- 0.0001

N <- 10000
burnin <- 1000

W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <- numeric(N)
theta_02_hist <- numeric(N)
theta1_star_hist <- matrix(0, N, Tt)
theta2_hist <- matrix(0, N, Tt)

# Initial Values
theta1_star <- stats::filter(log(y + 0.5), rep(1 / 5, 5), sides = 2)
theta1_star[is.na(theta1_star)] <- log(y[is.na(theta1_star)] + 0.5)
theta1_star <- as.numeric(theta1_star)

theta2 <- c(0, diff(theta1_star))
theta_01 <- log(y[1] + 0.5)

theta_02 <- 0
W1 <- 0.01
W2 <- 0.01

ess_is <- numeric(N)
itr_irls <- numeric(N) # number of iterations of IRLS (for earch gibbs step)
M_irls_max <- 20 # maximum iterations for IRLS

M_is <- 1
tol <- 1e-4
theta1_tilde <- log(y + 0.5)
#theta1_tilde <- numeric(Tt)
Weights <- matrix(0, N, M_is)


#####
start_time <- proc.time()
for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[1]]
        printf("Iteration %d / %d | Elapsed time: %.0f s | j mean: %.1f", n, N, elapsed_time, mean(itr_irls[1:(n-1)]))
    }

    # Sample theta_02 (conjugated normal)
    sigma2_02_bar <- (1 / sigma2_02 + 1 / W1 + 1 / W2)^(-1)
    mu_02_bar <- sigma2_02_bar * ((theta1_star[1] - theta_01) / W1 +
                                      theta2[1] / W2 +
                                      mu_02 / sigma2_02)
    theta_02 <- rnorm(1, mean = mu_02_bar, sd = sqrt(sigma2_02_bar))

    # Sample theta_01 (conjugated normal)
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar * (mu_01 / sigma2_01 +
                                      (theta1_star[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))

    # Sample W1 (conjugated invgamma)
    dif1   <- theta1_star - c(theta_01, theta1_star[-Tt])
    diffs1 <- dif1 - c(theta_02, theta2[-Tt])
    alpha_W1_bar <- alpha_W1 + Tt / 2
    beta_W1_bar  <- beta_W1 + 0.5 * sum(diffs1^2)
    W1 <- rinvgamma(1, shape = alpha_W1_bar, rate = beta_W1_bar)

    # Sample W2 (conjugated invgamma)
    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    alpha_W2_bar <- alpha_W2 + Tt / 2
    beta_W2_bar  <- beta_W2 + 0.5 * sum(diffs2^2)
    W2 <- rinvgamma(1, shape = alpha_W2_bar, rate = beta_W2_bar)

    sd_W1 <- sqrt(W1)
    sd_W2 <- sqrt(W2)

    # ------------------------------------------------------------------
    # Importance Sampling for  theta_t1
    # ------------------------------------------------------------------

    # Local approximation (IRLS)
    theta1_tilde_old <- theta1_tilde
    for (j in 1:M_irls_max) {

            f_t <- exp(-theta1_tilde)       # observational variance
            z_t <- theta1_tilde + f_t*y - 1 # pseudo-observation

            kf  <- forward_filter_1d(z_t, f_t, theta2, theta_01, theta_02, W1)
            ks <- kalman_smoother_1d(kf, W1)
            theta1_tilde <- ks$m_s   # update the mode
            if (max(abs(theta1_tilde - theta1_tilde_old)) < tol) break
            theta1_tilde_old <- theta1_tilde
    }
    itr_irls[n] <- j

    # Importance Sampling step
    if (M_is > 1) {
        log_w <- numeric(M_is)
        trajectories <- matrix(0, M_is, Tt)
        for (i in 1:M_is) {
            # sample theta1 proposed
            theta1_prop <- ffbs_1d(kf, ks, W1=W1)
            trajectories[i, ] <- theta1_prop

            # Log-weights: log p(y|theta1) - log g(y|theta1)
            log_p <- sum(y * theta1_prop - exp(theta1_prop))
            log_g <- sum(-0.5 * log(2*pi*f_t) - 0.5 * (theta1_prop - z_t)^2 / f_t)
            log_w[i] <- log_p - log_g
        }

        log_w <- log_w - logsumexp(log_w)
        w <- exp(log_w)
        Weights[n, ] <- w
        ess_is[n] <- 1 / sum(exp(2*log_w))

        idx <- sample(1:M_is, 1, prob=w)
        theta1_star <- trajectories[idx, ]
    } else {
        theta1_star <- ffbs_1d(kf, ks, W1)
    }


    # ------------------------------------------------------------------
    # Component-wise Gibbs for  theta_t2 | theta1_star, W1, W2
    # ------------------------------------------------------------------

    # Conditional variances (constant)
    sigma2_t2_bar_interior <- (1 / W1 + 2 / W2)^(-1)
    sigma2_t2_bar_last <- W2

    # t = 1
    sigma2_bar <- (1 / W1 + 2 / W2)^(-1)
    mu_bar <- sigma2_bar * ((theta1_star[2] - theta1_star[1]) / W1 +
                                theta2[2]  / W2 +
                                theta_02 / W2)
    theta2[1] <- rnorm(1, mean = mu_bar, sd = sqrt(sigma2_bar))

    # t = 2, ..., T-1
    # for (t in 2:(Tt - 1)) {
    #     mu_bar <- sigma2_t2_bar_interior * ((theta1_star[t + 1] - theta1_star[t]) / W1 +
    #                                             theta2[t + 1] / W2 +
    #                                             theta2[t - 1] / W2)
    #     theta2[t] <- rnorm(1, mean = mu_bar, sd = sqrt(sigma2_t2_bar_interior))
    # }
    sd_t2 <- sqrt(sigma2_t2_bar_interior)
    d1    <- diff(theta1_star) / W1       # (theta1[t+1]-theta1[t])/W1 para t=1..T-1
    for (t in 2:(Tt-1)) {
        mu_bar  <- sigma2_t2_bar_interior * (d1[t] + (theta2[t+1] + theta2[t-1])/W2)
        theta2[t] <- rnorm(1, mu_bar, sd_t2)
    }

    # t = T
    theta2[Tt] <- rnorm(1, mean = theta2[Tt - 1], sd = sd_W2)

    # Store the results
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_star_hist[n, ] <- theta1_star
    theta2_hist[n, ] <- theta2

}

elapsed_time <- (proc.time() - start_time)[[1]]
printf("Execution time: %.0f s", elapsed_time)

# Results
theta1_mean <- colMeans(theta1_star_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
W1_mean <- mean(W1_hist[-(1:burnin)])
W2_mean <- mean(W2_hist[-(1:burnin)])
printf("W1 posterior mean: %.6f", W1_mean)
printf("W2 posterior mean: %.6f", W2_mean)

cat(sprintf("ESS IS (min/mediana/max): %.1f / %.1f / %.1f de M_is=%d\n",
            min(ess_is), median(ess_is), max(ess_is), M_is))
hist(ess_is, breaks=50, main="ESS do IS por iteracao Gibbs")

# Plot
x <- 1:Tt

# lambda ####
lambda_mean <- exp(theta1_mean)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(
    x,
    y,
    type = "l",
    col = "gray",
    xlab = "t",
    ylab = "",
    main = "Poisson 2nd Order Polynomial Model"
)
points(x, y, pch = 20)
lines(x, lambda_mean, col = "red", lwd = 2)
lines(x, lambda_true, col = "blue", lwd = 2)
legend(
    "topright",
    legend = expression(y[t], lambda[t], hat(lambda)[t]),
    col = c("black", "blue", "red"),
    lty = c(NA, 1, 1),
    lwd = c(NA, 2, 2),
    pch = c(20, NA, NA),
    bty = "n"
)

# theta_t1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
ylim_range <- range(theta1_mean, theta1_true)
plot(
    x,
    theta1_mean,
    type = "l",
    col = "red",
    lwd = 2,
    ylim = ylim_range,
    xlab = "t",
    ylab = "",
    main = "theta_t1"
)
lines(x, theta1_true, col = "blue", lwd = 2)
legend(
    "topright",
    legend = expression(hat(theta)[t1], theta[t1]),
    col = c("red", "blue"),
    lwd = 2,
    bty = "n"
)

# theta_t2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(
    x,
    theta2_mean,
    type = "l",
    col = "red",
    lwd = 2,
    xlab = "t",
    ylab = "",
    main = "theta_t2"
)
legend(
    "topright",
    legend = expression(hat(theta)[t2], theta[t2]),
    col = c("red", "blue"),
    lwd = 2,
    bty = "n"
)

# Traceplot W1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(
    W1_hist[-(1:burnin)],
    type = "l",
    xlab = "n",
    ylab = "W1",
    main = "Traceplot of W1"
)

# Traceplot W2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(
    W2_hist[-(1:burnin)],
    type = "l",
    xlab = "n",
    ylab = "W2",
    main = "Traceplot of W2"
)

# Traceplots theta1_star ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(
        theta1_star_hist[, t],
        type = "l",
        main = bquote(theta[list(.(t), 1)]),
        xlab = "n",
        ylab = ""
    )
}

# Traceplots theta2_star
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(
        theta2_hist[, t],
        type = "l",
        main = bquote(theta[list(.(t), 2)]),
        xlab = "n",
        ylab = ""
    )
}

# Weights
i_<- min(M_is,10)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
matplot(Weights[,1:i_], type="l")

# Diagnostics
idx <- (burnin + 1):N
t1  <- theta1_star_hist[idx, ]
t2  <- theta2_hist[idx, ]

rho  <- sapply(1:Tt, function(t)
    cor(t1[, t], t2[, t]))
ess1 <- apply(t1, 2, effectiveSize)
ess2 <- apply(t2, 2, effectiveSize)

cat(sprintf("mediana |rho_t| : %.3f\n", median(abs(rho))))
cat(sprintf(
    "ESS theta1 (min/mediana): %.0f / %.0f\n",
    min(ess1),
    median(ess1)
))
cat(sprintf(
    "ESS theta2 (min/mediana): %.0f / %.0f\n",
    min(ess2),
    median(ess2)
))

# ESS
cat(sprintf(
    "ESS W1: %.0f | ESS W2: %.0f\n",
    +effectiveSize(W1_hist[(burnin + 1):N]),
    +effectiveSize(W2_hist[(burnin + 1):N])
))
