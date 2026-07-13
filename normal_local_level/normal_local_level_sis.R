# Sequential Importance Sampling

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- base::t       # alias to transpose function
options(error = function() traceback(2)) # more informative traceback


# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))


# Log-sum-exp auxiliary function
logsumexp <- function(x){
    c <- max(x)
    y <- c + log(sum((exp(x-c))))
    return(y)
}

# Load the data
source <- "normal_local_level_sim_200" # csv file with data
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
T <-length(y)
theta_true <- data$theta
W <- data$W
V <- data$V

# initialization
theta_t <- 1

theta <- numeric(T)

# theta_0 ~ N(m_0, C_0)
m0 <- y[1]
C0 <- 10

# Forward/Kalman Filtering (W&H, Seção 4.3, Teorema 4.1)
m_kf <- numeric(T)   # m_t = E[theta_t | D_t]
C_kf <- numeric(T)   # C_t = Var[theta_t | D_t]

# t = 1: prior (theta_1 | D_0) ~ N(a_1, R_1)
a1   <- m0
R1   <- C0 + W
Q1   <- R1 + V
A1   <- R1 / Q1
e1   <- y[1] - a1
m_kf[1] <- a1 + A1 * e1
C_kf[1] <- (1 - A1) * R1

for (t in 2:T) {
    a_t     <- m_kf[t-1]           # prior mean
    R_t     <- C_kf[t-1] + W       # prior variance
    Q_t     <- R_t + V             # forecast variance
    A_t     <- R_t / Q_t           # adaptive coefficient
    e_t     <- y[t] - a_t          # forecast error
    m_kf[t] <- a_t + A_t * e_t     # posterior mean
    C_kf[t] <- (1 - A_t) * R_t     # posterior variance
}

# proposed densitty g(theta_t | theta_{t-1}^(i))
log_p_yt <- function(yt, theta_t, V) {
    res <- -1/(2*V) * (yt - theta_t)^2
    return(res)
}

# Sequential Importance sampling with resampling
N <- 1000
theta_t <- rep(0, N)
log_w_t <- rep(0, N)
N0 <- N/2

Mt <- numeric(T)
Ct <- numeric(T)
Wt <- matrix(0, T, N)
Wt_log <- matrix(0, T, N)
N_eff     <- numeric(T)
idx_w_max <- numeric(T)

for (t in 1:T) {

    # Sample N theta_t
    theta_t <- rnorm(N, mean=theta_t, sd = sqrt(W))

    # Update the weights
    log_w_t <- log_p_yt(y[t], theta_t, V) + log_w_t
    Wt_log[t,] <- log_w_t

    # Normalize the weights
    log_w_tilde_t <- log_w_t - logsumexp(log_w_t)
    w_tilde_t <- exp(log_w_tilde_t)
    Wt[t,] <- w_tilde_t

    # Effective sample size
    N_eff[t]  <- 1 / sum(w_tilde_t^2)
    idx_w_max[t] <- which.max(w_tilde_t)

    # Posterior mean and variance
    mt_hat <- sum(w_tilde_t * theta_t)
    Ct_hat <- sum(w_tilde_t * (theta_t - mt_hat)^2)
    Mt[t] <- mt_hat
    Ct[t] <- Ct_hat
}


# theta filtered ####
x <- 1:T
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta_true, type="l", xlab="t", ylab="")
lines(x, m_kf, col="blue")
lines(x, Mt, col="red")
#polygon(c(x, rev(x)),
#        c(m_kf+ 3*sqrt(C_kf), rev(m_kf - 3*sqrt(C_kf))),
#        col = adjustcolor("steelblue", alpha.f = 0.3),
#        border = NA)


#####
#matplot(Wt_log[, 1:100], type = "l", ylab = "log w_t", xlab = "t",
#        main = "Log-pesos acumulados (sem normalização)")

#####
# matplot(x, Wt[,1:100], type="l", ylim=c(0,1.1))

# Effective sample size ####
plot(x, N_eff, type = "l", ylab = expression(N[eff]), xlab = "t",
     main = "Effective Sample Size")
abline(h = 1, col = "red", lty = 2)

# W index which is max
plot(x, idx_w_max, type="l")
