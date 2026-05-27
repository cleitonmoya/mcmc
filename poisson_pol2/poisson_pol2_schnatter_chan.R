# Poisson 2nd Order Polynomial Dynamic Model
# MCMC: Gibbs with data augmentation (Schnatter) +
#       Precision-Based state sampling (Chan)
#
# References:
#    Frühwirth-Schnatter, S., & Wagner, H. (2006).
#    Auxiliary Mixture Sampling for Parameter-Driven Models of
#    Time Series of Counts with Applications to State Space Modelling.
#    Biometrika, 93(4), 827–841.
#    https://doi.org/10.1093/biomet/93.4.827
#
# Author: Cleiton Moya de Almeida

library(Matrix)
library(coda)

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- base::t       # alias to transpose function
options(error = function() traceback(2)) # more informative traceback

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))


# Load the data
source <- "poisson_pol2_sim1" # csv file with data
df <- read.table(paste("../data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
theta1_true <- df$theta1
theta2_true <- df$theta2
T <- length(y)

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

inv2x2 <- function(M) {
    det_M <- M[1,1]*M[2,2] - M[1,2]*M[2,1]
    inv <- matrix(c(M[2,2], -M[2,1], -M[1,2], M[1,1]), 2, 2) / det_M
    return(inv)
}

rmvn_chol <- function(mu, Sigma) {
    L <- chol(Sigma)
    mu + drop(tp(L) %*% rnorm(length(mu)))
}

# Normal mixture for \epsilon_{tj} ~ Gumbel(0,1) approximation (from paper)
mix_params <- matrix(
    c(0.00397, 5.09,   4.50,
      0.00396, 3.29,   2.02,
      0.168,   1.82,   1.10,
      0.147,   1.24,   0.422,
      0.125,   0.764,  0.198,
      0.101,   0.391,  0.107,
      0.104,   0.0431, 0.0778,
      0.116,  -0.306,  0.0766,
      0.107,  -0.673,  0.0947,
      0.088,  -1.06,   0.146),
    nrow = 10, byrow = TRUE,
    dimnames = list(NULL, c("w", "m", "s2"))
)
R_mix <- nrow(mix_params)
w_vec  <- mix_params[, "w"]
w_vec <- w_vec/sum(w_vec)     # normalizing
m_vec  <- mix_params[, "m"]
s2_vec <- mix_params[, "s2"]
inv_s2 <- 1 / s2_vec
sqrt_s2 <- sqrt(s2_vec)


# Sample theta - Chan Method
sample_theta_chan <- function(y_hat, V_hat, G, W1, W2) {

    n    <- 2 * T
    invW <- diag(c(1/W1, 1/W2), nrow=2, ncol=2)

    invD     <- invW
    GtinvWG  <- tp(G) %*% invW %*% G
    blk1     <- GtinvWG + invD
    blk_mid  <- GtinvWG + invW
    blkT     <- invW
    blk_off  <- -invW %*% G

    diag_vals <- c(blk1[1,1], blk1[2,2],
                   rep(c(blk_mid[1,1], blk_mid[2,2]), T-2),
                   blkT[1,1], blkT[2,2])

    sub1_vals <- c(blk1[2,1], rep(blk_mid[2,1], T-2), blkT[2,1])
    i_s1      <- seq(2, n, by=2)

    t_idx <- seq_len(T-1)
    row_b <- c(2*t_idx+1, 2*t_idx+2, 2*t_idx+1, 2*t_idx+2)
    col_b <- c(2*t_idx-1, 2*t_idx-1, 2*t_idx,   2*t_idx  )
    val_b <- c(rep(blk_off[1,1], T-1), rep(blk_off[2,1], T-1),
               rep(blk_off[1,2], T-1), rep(blk_off[2,2], T-1))

    K <- sparseMatrix(
        i         = c(seq_len(n), i_s1,     row_b),
        j         = c(seq_len(n), i_s1-1,   col_b),
        x         = c(diag_vals,  sub1_vals, val_b),
        dims      = c(n, n),
        symmetric = TRUE
    )

    # V_hat é vetor de comprimento T: 1/V_hat[t] na posição (2t-1, 2t-1)
    obs_idx <- seq(1, n, by=2)
    P <- K + sparseMatrix(i=obs_idx, j=obs_idx,
                          x=1/V_hat, dims=c(n,n))   # <-- vetor, não escalar

    ch      <- Cholesky(P, LDL=FALSE, perm=FALSE)

    F_bar_b <- sparseMatrix(i=obs_idx, j=seq_len(T),
                            x=1/V_hat, dims=c(n,T))  # <-- vetor, não escalar
    b       <- as.vector(F_bar_b %*% y_hat)
    eta_hat <- solve(ch, b)
    u       <- rnorm(n)
    x       <- as.vector(solve(ch, u, system="Lt"))

    theta_vec <- as.vector(eta_hat + x)

    # Reformatar: theta[i,t] -> array c(2,1,T)
    theta <- array(NA, dim=c(2,1,T))
    for (t in 1:T) {
        theta[1,,t] <- theta_vec[2*t-1]
        theta[2,,t] <- theta_vec[2*t]
    }
    return(theta)
}


sample_W <- function(theta, nu_01, eta_01, nu_02, eta_02) {

    theta1 <- theta[1,,]
    theta2 <- theta[2,,]
    omega_t1   <- theta1[-1] - theta1[-T] -  theta2[-T]
    omega_t2   <- theta2[-1] - theta2[-T]

    nu_bar_01  <- nu_01 + (T-1)/2    # theta does not include theta_0
    eta_bar_01 <- eta_01 + sum(omega_t1^2)/2
    phi1 <- rgamma(1, shape=nu_bar_01, rate=eta_bar_01)
    W1 <- 1/phi1

    nu_bar_02 <- nu_02 + (T-1)/2
    eta_bar_02 <- eta_02 + sum(omega_t2^2)/2
    phi2 <- rgamma(1, shape=nu_bar_02, rate=eta_bar_02)
    W2 <- 1/phi2
    W <- diag(c(W1, W2), nrow=2, ncol=2)

    return(W)
}


sample_tau <- function(theta, y) {

    # tau = {tau_tj, j = 1, ... ,(y_t + 1), t = 1, ..., T}
    tau <- vector("list", T)

    for (t in 1:T) {
        y_t <- y[t]
        n_t <- y_t + 1

        theta_t1 <- theta[1,,t]
        lambda_t <- exp(theta_t1)
        xi_t <- rexp(1, rate = lambda_t)

        if (y_t > 0) {
            e <- rexp(y_t + 1)
            total_sum <- sum(e)
            tau_t <- e[1:y_t] / total_sum
            tau_t_ytp1 <- 1 - sum(tau_t) + xi_t
            tau_t <- c(tau_t, tau_t_ytp1)
        } else {
            tau_t <- xi_t
        }
        tau[[t]] <- tau_t
    }

    return(tau)
}


sample_S <- function(y, tau, theta) {

    # Pré-computar constantes
    log_w  <- log(mix_params[, "w"])
    log_s  <- 0.5 * log(mix_params[, "s2"])
    s2_inv <- 1 / mix_params[, "s2"]
    m_vec  <- mix_params[, "m"]

    # Empilhar todos os z_{tj} de uma vez: vetor de comprimento sum(y+1)
    n_t_vec <- y + 1L
    z_all <- unlist(lapply(seq_len(T), function(t)
        -log(tau[[t]]) - theta[1,,t]
    ))  # comprimento total = sum(n_t_vec)

    total <- length(z_all)

    # Log-pesos: total x R — uma única operação outer
    dif2  <- outer(z_all, m_vec, function(a, b) (a - b)^2)  # total x R
    log_v <- matrix(log_w - log_s, nrow=total, ncol=R_mix, byrow=TRUE) -
        0.5 * dif2 * matrix(s2_inv, nrow=total, ncol=R_mix, byrow=TRUE)

    # Estabilidade e normalização
    log_v <- log_v - apply(log_v, 1, max)
    v     <- exp(log_v)
    prob  <- v / rowSums(v)                        # total x R

    # Amostragem categórica: total draws, sem nenhum loop em R
    #cump <- matrix(apply(prob, 1, cumsum), nrow=R_mix)  # R x total
    cump <- t(apply(prob, 1, cumsum))
    u    <- runif(total)
    #r    <- colSums(cump < matrix(u, nrow=R_mix, ncol=total, byrow=TRUE)) + 1L
    r    <- rowSums(cump < u) + 1L

    # Reempacotar em lista
    S   <- vector("list", T)
    idx <- 0L
    for (t in seq_len(T)) {
        S[[t]] <- r[idx + seq_len(n_t_vec[t])]
        idx    <- idx + n_t_vec[t]
    }
    return(S)
}

# Simulation parameters

# Model main parameters
F <- matrix(c(1,0))             # dim = 2x1
G <- rbind(c(1, 1), c(0, 1))    # dim = 2x2
W1 <- 1                         # W_01
W2 <- 0.1                       # W_02
W <- diag(c(W1, W2), nrow=2, ncol=2)

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

N <- 1100
burnin <- 100

theta_samples   <- array(data=NA, dim=c(2,1,T,N))
W_samples       <- array(data=NA, dim=c(2,2,N))
tau_samples     <- vector("list", N)
y_tilde_samples <- vector("list", N)
S_samples       <- vector("list", N)

tau     <- vector("list", T)      # tau = {tau_tj, j=1,...,y_t+1, t=1,...,T}
S       <- vector("list", T)      # S =   {r_tj,   j=1,...,y_t+1, t=1,...,T}
y_tilde <- vector("list", T)


update_y_tilde <- function(tau, S) {

    y_tilde <- vector("list", T)

    for (t in 1:T) {
        #m_t <- mix_params[S[[t]], "m"]
        m_t <- m_vec[S[[t]]]
        y_tilde_t <- -(log(tau[[t]]) + m_t)
        y_tilde[[t]] <- y_tilde_t
    }

    return(y_tilde)
}

compute_sufficient <- function(tau, S) {
    y_hat <- numeric(T)
    V_hat <- numeric(T)
    for (t in 1:T) {
        s2_t    <- unname(mix_params[S[[t]], "s2"])
        m_t     <- unname(mix_params[S[[t]], "m"])
        y_tilde_t <- -log(tau[[t]]) - m_t      # tilde_y_tj = -log(tau_tj) - m_{r_tj}
        V_hat[t] <- 1 / sum(1/s2_t)
        y_hat[t] <- V_hat[t] * sum(y_tilde_t / s2_t)
    }
    return(list(y_hat=y_hat, V_hat=V_hat))
}

# Initialization of tau, S and y_tilde
for (t in 1:T) {
    y_t <- y[t]
    n_t <- y_t + 1

    # Initilize tau
    if (y_t > 0) {
        u <- sort.int(runif(y_t), method="shell")
        tau_t <- c(u[1], u[-1]-u[-y_t], 1 - u[y_t] + rexp(1))
    } else {
        tau_t <- 1 + rexp(1)
    }
    tau[[t]] <- tau_t

    # Initialize S
    # sample 1:10 randomly according to component weights probabilit
    S_t <- sample.int(R_mix, size = n_t, replace = TRUE, prob = w_vec)
    S[[t]] <- S_t

    #  Initialize y_tilde
    m_t <- m_vec[S[[t]]]
    y_tilde_t <- -(log(tau[[t]]) + m_t)
    y_tilde[[t]] <- y_tilde_t
}


# Main loop - Gibbs Sampling
start_time = proc.time() # execution time
for (n in 1:N) {

    if (n %% 100 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[3]]
        printf("Iteration %d / %d | Elapsed time: %.0f s", n, N, elapsed_time)
    }

    # theta ~ p(theta | W, tau, S)
    suf   <- compute_sufficient(tau, S)
    theta <- sample_theta_chan(suf$y_hat, suf$V_hat, G, W[1,1], W[2,2])
    theta_samples[,,,n] <- theta

    # W ~ p(W | theta, tau, S)
    W <- sample_W(theta, nu_01, eta_01, nu_02, eta_02)
    W_samples[,,n] <- W

    # tau ~ p(tau | theta, W, S)
    tau <- sample_tau(theta, y)
    tau_samples[[n]] <- tau

    # S ~ p(S | theta, W, tau)
    S <- sample_S(y, tau, theta)
    S_samples[[n]] <- S

    # Update y_tilde
    y_tilde <- update_y_tilde(tau, S)
    y_tilde_samples[[n]] <- y_tilde

}
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]


# Simulation Summary ####
printf("Execution time: %.0f s", elapsed_time)

theta1_hat <- apply(theta_samples[1,1,,-(1:burnin)], 1, mean)
theta2_hat <- apply(theta_samples[2,1,,-(1:burnin)], 1, mean)
W1_hat <- mean(W_samples[1,1,-(1:burnin)])
W2_hat <- mean(W_samples[2,2,-(1:burnin)])
printf("W1 mean: %.3f", W1_hat)
printf("W2 mean: %.5f", W2_hat)

lambda_true <- exp(theta1_true)
lambda_hat <- exp(theta1_hat)


# Effective sample size ####
printf("Effective Sample Size:")
ess_w1 <- effectiveSize(mcmc(W_samples[1,1,-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W_samples[2,2,-(1:burnin)]))
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)

for (t in c(50, 100, 200, 250)) {
    ess <- effectiveSize(mcmc(theta_samples[1,1,t,-(1:burnin)]))
    printf("\ttheta %d,1: %0.f", t, ess)
}

for (t in c(50, 100, 200, 250)) {
    ess <- effectiveSize(mcmc(theta_samples[2,1,t,-(1:burnin)]))
    printf("\ttheta %d,2: %0.f", t, ess)
}


# Effective sample size / elapsed time ####
printf("Effective Sample Size / seconds:")
printf("\tW1: %.2f", ess_w1/elapsed_time)
printf("\tW2: %.2f", ess_w2/elapsed_time)

for (t in c(50, 100, 200, 250)) {
    ess <- effectiveSize(mcmc(theta_samples[1,1,t,-(1:burnin)]))
    printf("\ttheta %d,1: %.2f", t, ess/elapsed_time)
}

for (t in c(50, 100, 200, 250)) {
    ess <- effectiveSize(mcmc(theta_samples[2,1,t,-(1:burnin)]))
    printf("\ttheta %d,2: %.2f", t, ess/elapsed_time)
}


# Plots ####
# y, lambda_true, lambda_estimated #####
x <- 1:T
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson 2nd Order Polynomial Model")
points(x, y, pch = 20)
lines(x, lambda_hat, col="red", lwd=1)
lines(x, lambda_true, col="blue", lwd=1)
legend("topright",
       legend = expression(y[t], lambda[t], hat(lambda)[t]),
       col = c("black", "blue", "red"),
       lty = c(NA, 1, 1),
       lwd = c(NA, 1, 1),
       pch = c(20, NA, NA),
       bty = "n")


# theta1_true, theta2_mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_hat, type="l", ylab="", col="blue", lwd=2)
lines(x, theta1_true)
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


# theta2_true, theta2_mean ####
plot(x, theta2_hat, type="l", ylab="", col="blue", lwd=2)
lines(x, theta2_true)
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


# Posterior distribution of theta_t1 ####
par(mfrow = c(2, 2))
for (t in c(50, 100, 200, 250)) {
    hist(theta_samples[1, 1, t, -(1:burnin)], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Posterior of " * theta[.(t) * "," * 1]))
    lines(density(theta_samples[1, 1, t, -(1:burnin)]), col = "blue", lwd = 2)
}


# Posterior distribution of theta_t2 ####
par(mfrow = c(2, 2))
for (t in c(50, 100, 200, 250)) {
    hist(theta_samples[2, 1, t, -(1:burnin)], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Posterior of " * theta[.(t) * "," * 2]))
    lines(density(theta_samples[2, 1, t, -(1:burnin)]), col = "blue", lwd = 2)
}



# Posterior distribution of W1 ####

par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W_samples[1,1,-(1:burnin)], breaks = 50, freq = FALSE, xlab="",
     main ="Posterior of W1")
lines(density(W_samples[1,1,-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W_samples[2,2,-(1:burnin)], breaks = 50, freq = FALSE, xlab="",
     main ="Posterior of W2")
lines(density(W_samples[2,2,-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot for W1 and W2 ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W_samples[1,1,-(1:100)], type="l", xlab="n", ylab=expression(W[1]),
     main=expression("Traceplot of " * W[1]))
plot(W_samples[2,2,-(1:100)], type="l", xlab="n", ylab=expression(W[2]),
     main=expression("Traceplot of " * W[2]))


# Traceplot theta_t1 ####
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
for (t in c(50,100,200,250)) {
    plot(theta_samples[1,1,50,], type="l", xlab="", ylab="",
         main=bquote(theta[.(t)*","*1]))
}


# Traceplot theta_t2 ####
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
for (t in c(50,100,200,250)) {
    plot(theta_samples[2,1,50,], type="l", xlab="t", ylab="",
         main=bquote(theta[.(t)*","*2]))
}
