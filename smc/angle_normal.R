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

q <- 0.001^2
r <- 0.005^2

# True initial values
x <- -0.05
dx <- 0.001
y <- 0.7
dy <- -0.055

T <- 24

# Simulation ####
X_true <- numeric(T)
Dx_true <- numeric(T)
Y_true <- numeric(T)
Dy_true <- numeric(T)
Z <- numeric(T)

X_true[1] <- x
Dx_true[1] <- dx
Y_true[1] <- y
Dy_true[1] <- dy

Z[1] <- atan(y/x) + rnorm(1, 0, sqrt(r))

for (t in 2:T) {
    w_x <- rnorm(1, 0, sqrt(q))
    w_y <- rnorm(1, 0, sqrt(q))

    x <- x + dx + 0.5*w_x
    dx <- dx + w_x
    y <- y + dy + 0.5*w_y
    dy <- dy + w_y
    z <- atan(y/x) + rnorm(1, 0, sqrt(r))

    X_true[t] <- x
    Dx_true[t] <- dx
    Y_true[t] <- y
    Dy_true[t] <-dy
    Z[t] <- z
}

# Save the data
saveRDS(list(X = X_true, Y=Y_true, Z=Z), "angle_normal.rds")

# Plot ####
t_seq <- 1:T
plot(X_true, Y_true, pch=2, cex=0.5, xlim=c(-0.05, 0.01))
points(0, 0, pch=3)
segments(0, 0, X_true, Y_true, col='gray')
