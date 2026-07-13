# Bootstrap Filter
#
# Reference:
#   Gordon, N. J., Salmond, D. J., & Smith, A. F. M. (1993).
#   Novel approach to nonlinear/non-Gaussian Bayesian state estimation.
#   IEE Proceedings F Radar and Signal Processing, 140(2), 107.
#   https://doi.org/10.1049/ip-f-2.1993.0015

# Example 4.2 - Bearings-only tracking

graphics.off()
rm(list = ls())
cat("\014")
options(error = function() traceback(2))
set.seed(42)

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

# Log-sum-exp auxiliary function
logsumexp <- function(x){
    c <- max(x)
    y <- c + log(sum((exp(x-c))))
    return(y)
}

# Load the data
data <- readRDS("angle_normal.rds")
X_true <- data$X
Y_true <- data$Y
Z <- data$Z

# Model
q <- 0.001^2
r <- 0.005^2
T <- length(Z)

# Likelihood function
loglik <- function(z_obs, y, x) {
    zlog <- dnorm(z_obs, mean=atan(y/x), sd=sqrt(r), log=TRUE)
    return(zlog)
}

N <- 4000 # number of particles

# Prior distribution
x_bar <- 0
y_bar <- 0.4
dx_bar <- 0
dy_bar <- -0.05

sig_x <- 0.5
sig_y <- 0.3
sig_dx <- 0.005
sig_dy <- 0.01

# Initial values
x <- rnorm(N, x_bar, sig_x)
y <- rnorm(N, y_bar, sig_y)
dx <- rnorm(N, dx_bar, sig_dx)
dy <- rnorm(N, dy_bar, sig_dy)

W <- matrix(0, T, N)

X <- numeric(T)
Dx <- numeric(T)
Y <- numeric(T)
Dy <- numeric(T)

# t=1
l <- loglik(Z[1], y, x)
w <- exp(l - logsumexp(l))
idx <- sample(1:N, size = N, replace = TRUE, prob = w)

x <- x[idx];
dx <- dx[idx];
y <- y[idx];
dy <- dy[idx]

X[1]  <- mean(x);
Y[1]  <- mean(y)
Dx[1] <- mean(dx);
Dy[1] <- mean(dy)

for (t in 2:T) {

    # Prediction step p(x_{t-1} | D_{t})
    w_x <- rnorm(N, 0, sqrt(q))
    w_y <- rnorm(N, 0, sqrt(q))

    x <- x + dx + 0.5*w_x
    dx <- dx + w_x
    y <- y + dy + 0.5*w_y
    dy <- dy + w_y

    # Update step
    l <- loglik(Z[t], y, x)
    w <- exp(l - logsumexp(l))
    W[t,] <- w

    # resample
    idx <- sample(1:N, size = N, replace = TRUE, prob = w)
    x  <- x[idx]
    y  <- y[idx]
    dx <- dx[idx]
    dy <- dy[idx]

    X[t] <- mean(x)
    Y[t] <- mean(y)
    Dx[t] <- mean(dx)
    Dy[t] <- mean(dy)

}

# Plot ####
t_seq <- 1:T
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(X_true, Y_true, pch=2, cex=0.5, xlim=c(-0.05, 0.01))
lines(X, Y, col="red")
points(0, 0, pch=3)
segments(0, 0, X_true, Y_true, col='gray')
matplot(t_seq, W[,1:10], type="l")
