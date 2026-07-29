# Poisson 2nd Order Polynomial Dynamic Model
# MCMC: Gibbs with data augmentation
#
# Reference:
#    Frühwirth-Schnatter, S., & Wagner, H. (2006).
#    Auxiliary Mixture Sampling for Parameter-Driven Models of
#    Time Series of Counts with Applications to State Space Modelling.
#    Biometrika, 93(4), 827–841.
#    https://doi.org/10.1093/biomet/93.4.827

# Author: Cleiton Moya de Almeida

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
source <- "poisson_pol2_sim_200" # csv file with data
df <- read.table(paste("../data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
theta1_true <- df$theta1
theta2_true <- df$theta2
Tt <- length(y)

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


# Sample from a bivariate normal distribution
rbvn_cond <- function(mu1, mu2, s11, s12, s22) {
    x1 <- rnorm(1, mu1, sqrt(s11))
    mu2_c <- mu2 + (s12/s11)*(x1 - mu1)
    s22_c <- s22 - s12^2/s11
    x2 <- rnorm(1, mu2_c, sqrt(max(s22_c, 0)))
    c(x1, x2)
}


# Forward Filtering (Kalman Filter)
forward_filter <- function(theta_01, theta_02, W1, W2, y_tilde, S) {
    m1 <- theta_01; m2 <- theta_02
    c11 <- 0; c12 <- 0; c22 <- 0

    a1v <- a2v <- R11v <- R12v <- R22v <- m1v <- m2v <- c11v <- c12v <- c22v <- numeric(Tt)

    for (t in 1:Tt) {

        s2_t    <- s2_vec[S[[t]]]
        V_hat_t <- 1/sum(1/s2_t)
        y_hat_t <- V_hat_t * sum(y_tilde[[t]]/s2_t)

        a1 <- m1 + m2
        a2 <- m2
        R11 <- c11 + 2*c12 + c22 + W1
        R12 <- c12 + c22
        R22 <- c22 + W2

        Q <- R11 + V_hat_t
        e <- y_tilde[t] - a1

        m1 <- a1 + (R11/Q)*e
        m2 <- a2 + (R12/Q)*e
        c11 <- R11*V_hat_t/Q
        c12 <- R12*V_hat_t/Q
        c22 <- R22 - R12^2/Q

        a1v[t]<-a1; a2v[t]<-a2; R11v[t]<-R11; R12v[t]<-R12; R22v[t]<-R22
        m1v[t]<-m1; m2v[t]<-m2; c11v[t]<-c11; c12v[t]<-c12; c22v[t]<-c22
    }
    list(a1=a1v,a2=a2v,R11=R11v,R12=R12v,R22=R22v,m1=m1v,m2=m2v,c11=c11v,c12=c12v,c22=c22v)
}


# FFBS
sample_theta_ffbs <- function(theta_01, theta_02, W1, W2, y_tilde, V, S) {

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

sample_W <- function(theta, nu_01, eta_01, nu_02, eta_02) {

    theta1 <- theta[1,,]
    theta2 <- theta[2,,]
    omega_t1   <- theta1[-1] - theta1[-Tt] -  theta2[-Tt]
    omega_t2   <- theta2[-1] - theta2[-Tt]

    nu_bar_01  <- nu_01 + (Tt-1)/2    # theta does not include theta_0
    eta_bar_01 <- eta_01 + sum(omega_t1^2)/2
    phi1 <- rgamma(1, shape=nu_bar_01, rate=eta_bar_01)
    W1 <- 1/phi1

    nu_bar_02 <- nu_02 + (Tt-1)/2
    eta_bar_02 <- eta_02 + sum(omega_t2^2)/2
    phi2 <- rgamma(1, shape=nu_bar_02, rate=eta_bar_02)
    W2 <- 1/phi2
    W <- diag(c(W1, W2), nrow=2, ncol=2)

    return(W)
}


sample_tau <- function(theta, y) {

    # tau = {tau_tj, j = 1, ... ,(y_t + 1), t = 1, ..., Tt}
    tau <- vector("list", Tt)

    for (t in 1:Tt) {
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

    n_t_vec <- y + 1L
    z_all <- unlist(lapply(seq_len(Tt), function(t)
        -log(tau[[t]]) - theta[1,,t]
    ))  # comprimento total = sum(n_t_vec)

    total <- length(z_all)

    # Log-pesos: total x R
    dif2  <- outer(z_all, m_vec, function(a, b) (a - b)^2)  # total x R
    log_v <- matrix(log_w - log_s, nrow=total, ncol=R_mix, byrow=TRUE) -
        0.5 * dif2 * matrix(s2_inv, nrow=total, ncol=R_mix, byrow=TRUE)

    # Estabilidade e normalização
    log_v <- log_v - apply(log_v, 1, max)
    v     <- exp(log_v)
    prob  <- v / rowSums(v)

    # Amostragem categórica
    cump <- tp(apply(prob, 1, cumsum))
    u    <- runif(total)
    r    <- rowSums(cump < u) + 1L

    # Reempacotar em lista
    S   <- vector("list", Tt)
    idx <- 0L
    for (t in seq_len(Tt)) {
        S[[t]] <- r[idx + seq_len(n_t_vec[t])]
        idx    <- idx + n_t_vec[t]
    }
    return(S)
}


# Simulation parameters

# Model main parameters
F <- matrix(c(1,0))             # dim = 2x1
Ft <- tp(F)
G <- rbind(c(1, 1), c(0, 1))    # dim = 2x2
Gt <- tp(G)
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

N <- 1000
burnin <- 100

theta_samples   <- array(data=NA, dim=c(2,1,Tt,N))
W_samples       <- array(data=NA, dim=c(2,2,N))
tau_samples     <- vector("list", N)
y_tilde_samples <- vector("list", N)
S_samples       <- vector("list", N)

tau     <- vector("list", Tt)      # tau = {tau_tj, j=1,...,y_t+1, t=1,...,Tt}
S       <- vector("list", Tt)      # S =   {r_tj,   j=1,...,y_t+1, t=1,...,Tt}
y_tilde <- vector("list", Tt)

update_y_tilde <- function(tau, S) {

    y_tilde <- vector("list", Tt)

    for (t in 1:Tt) {
        m_t <- m_vec[S[[t]]]
        y_tilde_t <- -(log(tau[[t]]) + m_t)
        y_tilde[[t]] <- y_tilde_t
    }

    return(y_tilde)
}


# Initialization of tau, S and y_tilde
for (t in 1:Tt) {
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
    theta <- sample_theta_ffbs(mu_01, sigma2_01, mu_02, sigma2_02,
                          F, Ft, G, Gt, W, S, y_tilde)
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


# Plots
# y, lambda_true, lambda_estimated #####
x <- 1:Tt
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
