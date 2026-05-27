# Bootstrap Filter
#
# Reference:
#   Gordon, N. J., Salmond, D. J., & Smith, A. F. M. (1993).
#   Novel approach to nonlinear/non-Gaussian Bayesian state estimation.
#   IEE Proceedings F Radar and Signal Processing, 140(2), 107.
#   https://doi.org/10.1049/ip-f-2.1993.0015

# Example 4.2 - Bearings-only tracking

library(circular)

graphics.off()
rm(list = ls())
cat("\014")
set.seed(42)

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Função auxiliar: mapear qualquer ângulo para (-pi, pi]
wrap <- function(z) ((z + pi) %% (2 * pi)) - pi

# Sample
rwc  <- function(n, mu, rho) {
    u <- runif(n)
    wrap(mu + 2*atan(((1-rho)/(1+rho)) * tan(pi*(u - 0.5))))
}

q <- 0.001^2
r <- 0.005^2
rho <- 1-0.005^2

# True initial values
x <- -0.05
dx <- 0.001
y <- 0.7
dy <- -0.055

T <- 50

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

Z[1] <- rwc(1, atan2(y, x), rho)

for (t in 2:T) {
    w_x <- rnorm(1, 0, sqrt(q))
    w_y <- rnorm(1, 0, sqrt(q))

    x <- x + dx + 0.5*w_x
    dx <- dx + w_x
    y <- y + dy + 0.5*w_y
    dy <- dy + w_y
    mu <- atan2(y, x)
    z <- rwc(1, atan2(y, x), rho)

    X_true[t] <- x
    Dx_true[t] <- dx
    Y_true[t] <- y
    Dy_true[t] <-dy
    Z[t] <- z
}

# Save the data
saveRDS(list(X = X_true, Y=Y_true, Z=Z), "angle_cauchy.rds")

# Plot ####
t_seq <- 1:T
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(X_true, Y_true, pch=2, cex=0.5)
points(0, 0, pch=3)
segments(0, 0, X_true, Y_true, col='gray')
