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
source <- "poisson_sin_200" # csv file with data
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
T <-length(y)
theta_true <- data$theta
W <- data$W

# initialization
theta_t <- 1
theta <- numeric(T)

# likelihood density g(theta_t | theta_{t-1}^(i))
log_p_yt <- function(yt, theta_t) {
    res <- yt*theta_t - exp(theta_t)
    return(res)
}

# Stratified resampling
stratified_resampling <- function(y, p, N) {
    Q <- cumsum(p)
    V <- runif(N)
    U <- (1:N-1+V)/N
    Y <- numeric(N)
    j <- 1
    for (t in 1:N) {
        while(U[t] > Q[j]) {
            j <- j + 1
        }
        Y[t] <- y[j]
    }
    return(Y)
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

    # Sample N theta_t from proposed density
    theta_t <- rnorm(N, mean=theta_t, sd = sqrt(W))

    # Update the weights
    log_w_t <- log_p_yt(y[t], theta_t) + log_w_t
    Wt_log[t,] <- log_w_t

    # Normalize the weights
    log_w_tilde_t <- log_w_t - logsumexp(log_w_t)
    w_tilde_t <- exp(log_w_tilde_t)
    Wt[t,] <- w_tilde_t

    # Effective sample size
    N_eff[t]  <- 1 / sum(w_tilde_t^2)
    idx_w_max[t] <- which.max(w_tilde_t)

    # Resample strategy
    if (N_eff[t] < N0) {
        theta_t <- sample(theta_t, size=N, replace=TRUE, prob=w_tilde_t)
        #theta_t <- stratified_resampling(theta_t, w_tilde_t, N)
        log_w_t  <- rep(0, N)
    }

    # Posterior mean and variance
    mt_hat <- sum(w_tilde_t * theta_t)
    Ct_hat <- sum(w_tilde_t * (theta_t - mt_hat)^2)
    Mt[t] <- mt_hat
    Ct[t] <- Ct_hat
}

N_eff_mean <- mean(N_eff)
print(N_eff_mean)

# theta filtered ####
x <- 1:T
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right

y1 <- theta_true
y2 <- Mt
y3 <- Mt + 3*sqrt(Ct)
y4 <- Mt - 3*sqrt(Ct)
ylim_range <- range(y1, y2, y3, y3, y4)

plot(x, y1, type="l", xlab="t", ylab="", ylim=ylim_range)
lines(x, y2, col="red")
polygon(c(x, rev(x)),
        c(y3, rev(y4)),
        col = adjustcolor("steelblue", alpha.f = 0.3),
        border = NA)

# Log-weigths ####
matplot(Wt_log[, 1:20], type = "l", ylab = "log w_t", xlab = "t",
        main = "Accumulated log-weights")

# Weights ####
matplot(x, Wt[,1:20], type="l")

# Effective sample size ####
plot(x, N_eff, type = "l", ylab = expression(N[eff]), xlab = "t",
     main = "Effective Sample Size")
abline(h = N_eff_mean, col = "red", lty = 2)
