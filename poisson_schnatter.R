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

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- base::t       # alias to transpose function
options(error = function() traceback(2)) # more informative traceback

# Change de directory to the same of the current file
# setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
setwd("C:/Users/cleit/OneDrive/Documentos/Projetos/R/mcmc")

# Load the data
source <- "pol2_sim1" # csv file with data
df <- read.table(paste("data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
theta1_true <- df$theta1
theta2_true <- df$theta2
T <- length(y)

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
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
    dimnames = list(paste0("r=", 1:10), c("w", "m", "s2"))
)
R_mix <- nrow(mix_params)


# Forward Filtering (Kalman Filter)
forward_filter <- function(mu_01, sigma2_01, mu_02, sigma2_02,
                           F, G, W, S, y_tilde) {
    m_t <- as.matrix(c(mu_01, mu_02))
    C_t <- diag(c(sigma2_01, sigma2_02))

    a <- array(data=NA, dim=c(2,1,T))
    R <- array(data=NA, dim=c(2,2,T))
    m <- array(data=NA, dim=c(2,1,T))
    C <- array(data=NA, dim=c(2,2,T))
    B <- array(data=NA, dim=c(2,2,T))

    for (t in 1:T) {

        # Sufficient scalar observations
        s2_t    <- unname(mix_params[S[[t]], "s2"])
        V_hat_t <- 1/sum(1/s2_t)
        y_hat_t <- V_hat_t * sum(y_tilde[[t]]/s2_t)

        # Prior: (theta_t | D_{t-1}) ~ N[a_t, R_t]
        a_t <- G %*% m_t
        R_t <- G %*% C_t %*% tp(G) + W
        a[,,t] <- a_t
        R[,,t] <- R_t

        # One-step-ahead forecast
        f_t <- tp(F) %*% a_t
        Q_t <- drop(tp(F) %*% R_t %*% F + V_hat_t)

        # Posterior: (theta_t | D_t) ~ N[m_t, C_t]
        e_t <- drop(y_hat_t - f_t)
        m_t <- a_t + (1/Q_t) * R_t %*% F * e_t
        C_t <- R_t - (1/Q_t) * R_t %*% F %*% tp(F) %*% R_t
        m[,,t] <- m_t
        C[,,t] <- C_t

    }
    return(list("a"=a, "R"=R, "m"=m, "C"=C))
}

# FFBS
sample_theta <- function(mu_01, sigma2_01, mu_02, sigma2_02,
                         F, G, W, S, y_tilde) {

    # theta = {theta_1, ..., theta_T}
    theta <- array(data=NA, dim=c(2,1,T))

    res <- forward_filter(mu_01, sigma2_01, mu_02, sigma2_02,
                          F, G, W, S, y_tilde)
    a <- res$a
    R <- res$R
    m <- res$m
    C <- res$C

    # Compute B_t
    B <- array(data=NA, dim=c(2,2,T))
    for (t in 1:(T-1)) {
        B[,,t] <- C[,,t] %*% tp(G) %*% solve(R[,,t+1])
    }

    # Sample theta_T ~ N[m_T, C_T]
    theta_t <- as.matrix(MASS::mvrnorm(1, mu=m[,,T], Sigma=C[,,T]))
    theta[,,T] <- theta_t

    # Backward: t = T-1, ..., 1
    # theta_t | theta_{t+1} ~ N(h_t, H_t) [
    # West and Harrison 1997, eq. 15.8, p. 570]
    for (t in seq(T-1, 1)) {
        h_t <- m[,,t] + B[,,t] %*% (theta_t - a[,,t+1])
        H_t <- C[,,t] - B[,,t] %*% R[,,t+1] %*% tp(B[,,t])
        theta_t <- as.matrix(MASS::mvrnorm(1, mu=h_t, Sigma=H_t))
        theta[,,t] <- theta_t
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
    W <- diag(c(W1, W2))

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
            u <- sort(runif(y_t))
            tau_t <- c(u[1], u[-1]-u[-y_t])
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
    S <- vector("list", T)
    for (t in 1:T) {
        n_t  <- y[t] + 1
        log_lambda_t <- theta[1,,t]          # = theta_{t1}
        S_t  <- integer(n_t)
        for (j in 1:n_t) {
            z_tj <- -log(tau[[t]][j]) - log_lambda_t   # = -log(tau_tj) - theta_{t1}
            v <- numeric(R_mix)
            for (k in 1:R_mix) {
                w_k  <- mix_params[k, "w"]
                m_k  <- mix_params[k, "m"]
                s2_k <- mix_params[k, "s2"]
                v[k] <- (w_k/sqrt(s2_k)) * exp(-(z_tj - m_k)^2 / (2*s2_k))
            }
            S_t[j] <- sample(1:R_mix, size=1, prob=v/sum(v))
        }
        S[[t]] <- S_t
    }
    return(S)
}

# Simulation parameters

# Model main parameters
F <- matrix(c(1,0))             # dim = 2x1
G <- rbind(c(1, 1), c(0, 1))    # dim = 2x2
W1 <- 1                         # W_01
W2 <- 0.1                       # W_02
W <- diag(c(W1, W2))            # dim = 2x2

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

theta_samples   <- array(data=NA, dim=c(2,1,T,N))
W_samples       <- array(data=NA, dim=c(2,2,N))
tau_samples     <- vector("list", N)
y_tilde_samples <- vector("list", N)
S_samples       <- vector("list", N)

#theta   <- array(data=NA, dim=c(2,1,T))
tau     <- vector("list", T)      # tau = {tau_tj, j=1,...,y_t+1, t=1,...,T}
S       <- vector("list", T)      # S =   {r_tj,   j=1,...,y_t+1, t=1,...,T}
y_tilde <- vector("list", T)

update_y_tilde <- function(tau, S) {

    y_tilde <- vector("list", T)

    for (t in 1:T) {
        m_t <- unname(mix_params[S[[t]], "m"])
        y_tilde_t <- -(log(tau[[t]]) + m_t)
        y_tilde[[t]] <- y_tilde_t
    }

    return(y_tilde)
}


# Initialization of tau, S and y_tilde
for (t in 1:T) {
    y_t <- y[t]
    n_t <- y_t + 1

    # Initilize tau
    if (y_t > 0) {
        u <- sort(runif(y_t))
        tau_t <- c(u[1], u[-1]-u[-y_t], 1 - u[y_t] + rexp(1))
    } else {
        tau_t <- 1 + rexp(1)
    }
    tau[[t]] <- tau_t

    # Initialize S
    # sample 1:10 randomly according to component weights probabilit
    S_t <- sample(1:R_mix, size = n_t, replace = TRUE, prob = mix_params[, "w"])
    S[[t]] <- S_t

    #  Initialize y_tilde
    m_t <- unname(mix_params[S[[t]], "m"])
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
    theta <- sample_theta(mu_01, sigma2_01, mu_02, sigma2_02,
                          F, G, W, S, y_tilde)
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
printf("Execution time: %.0f s", elapsed_time)

theta1_hat <- apply(theta_samples[1,1,,-(1:burnin)], 1, mean)
theta2_hat <- apply(theta_samples[2,1,,-(1:burnin)], 1, mean)

lambda_true <- exp(theta1_true)
lambda_hat <- exp(theta1_hat)

#####
# Plots
# y, lambda_true, lambda_estimated
x <- 1:T

# theta1_true, theta2_mean
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_hat, type="l", ylab="")
lines(x, theta1_true, col="blue")
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")

# theta2_true, theta2_mean
#par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta2_hat, type="l", ylab="")
lines(x, theta2_true, col="blue")
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")

# Traceplot theta_t2
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
for (t in c(50,100,200,250)) {
    plot(theta_samples[1,1,50,], type="l", xlab="", ylab="")
}

# Traceplot theta_t1
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
for (t in c(50,100,200,250)) {
    plot(theta_samples[2,1,50,], type="l", xlab="t", ylab="")
}
