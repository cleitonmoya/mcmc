# Auxiliary Particle Filter
# Reference:
#    Pitt, M. K., & Shephard, N. (1999).
#    Filtering via Simulation: Auxiliary Particle Filters. Journal of the
#    American Statistical Association, 94(446), 590–599.
#    https://doi.org/10.1080/01621459.1999.10474153
#
# Example 3.4 - Time series of angles (bearings-only tracking)

#####

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

log_dwc <- function(z, mu, rho) {
    # z: escalar (observação)
    # mu: vetor de N ângulos (um por partícula)
    # rho: escalar
    log(1 - rho^2) - log(2 * pi) - log(1 + rho^2 - 2 * rho * cos(z - mu))
}

# Load the data
data <- readRDS("angle_cauchy.rds")
X_true <- data$X
Y_true <- data$Y
Z <- data$Z
T <- length(Z)

# Number of particles
N <- 4000

# Model
q <- 0.001^2
rho <- 1-0.005^2

# Prior distribution
x_bar  <- -0.05;  sig_x  <- sqrt(0.001) * 0.5    # ~ 0.016
y_bar  <-  0.7;   sig_y  <- sqrt(0.001) * 0.3    # ~ 0.0095
dx_bar <-  0.001; sig_dx <- sqrt(0.001) * 0.005  # ~ 0.00016
dy_bar <- -0.055; sig_dy <- sqrt(0.001) * 0.01   # ~ 0.00032

# Initial values
x <- rnorm(N, x_bar, sig_x)
y <- rnorm(N, y_bar, sig_y)
dx <- rnorm(N, dx_bar, sig_dx)
dy <- rnorm(N, dy_bar, sig_dy)

X <- numeric(T)
Dx <- numeric(T)
Y <- numeric(T)
Dy <- numeric(T)
W <- matrix(0, T, N)

#####
for (t in 1:T) {

    # Step 0: particle representative value
    mu_x  <- x
    mu_y  <- y
    mu_bearing <- atan2(mu_y, mu_x)

    # Step 1: First stage weights
    log_l1 <- log_dwc(Z[t], mu = mu_bearing, rho = rho)
    w1 <- exp(log_l1 - logsumexp(log_l1))

    # Step 2: Auxiliary variable smapling (particle choicing)
    k <- sample(1:N, size = N, replace = TRUE, prob = w1)

        # Step 3: Particle propagation
    u1 <- rnorm(N, 0, sqrt(q))
    u2 <- rnorm(N, 0, sqrt(q))
    x  <- x[k] + dx[k] + 0.5*u1
    dx <- dx[k] + u1
    y  <- y[k] + dy[k] + 0.5*u2
    dy <- dy[k] + u2

    # Step 4: Second stage weights
    mu_new <- atan2(y, x)
    mu_rep <- atan2(mu_y[k], mu_x[k])  # representative

    log_l2_num <- log_dwc(Z[t], mu = mu_new,  rho = rho)
    log_l2_den <- log_dwc(Z[t], mu = mu_rep,  rho = rho)
    log_w2 <- log_l2_num - log_l2_den
    w2 <- exp(log_w2 - logsumexp(log_w2))
    W[t,] <- w2

    # Step 5: Second stage resample
    idx <- sample(1:N, size = N, replace = TRUE, prob = w2)
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
plot(X, Y, col="red", type="l")
points(X_true, Y_true, pch=2, cex=0.5)
points(0, 0, pch=3)
segments(0, 0, X_true, Y_true, col='gray')

#####
matplot(t_seq, W[,1:10], type="l")

