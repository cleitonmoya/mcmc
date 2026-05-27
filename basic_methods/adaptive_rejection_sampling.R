# Adaptive Rejection Sampling
# (Gilks and Wild, 1992)
# Author: Cleiton Moya de Almeida

library(Deriv)

# Clear the plots, environment and console
graphics.off()
rm(list = ls())
cat("\014")

set.seed(42)

### SIMULATION PARAMETERS
# Target density - Gumbel distribution
mu <- 1
beta <- 2
f <- function(x) {
    z <- (x - mu) / beta
    (1/beta) * exp(-(z + exp(-z)))
}

# Target function (envelope, can be unormalized)
g <- function(x){
    z <- (x - mu)/beta
    exp(-(z + exp(-z)))
}

# Log of envelope function
h <- function(x){
    -((x - mu)/beta) - exp(-(x - mu)/beta)
}

hprime <- Deriv(h)

xmin <- -5       # minimum value to sample
xmax <- 10       # maximum value to sample
X <- c(-1, 0, 1) # (initial) abscissae points (global variable)
k <- length(X)   # (initial) number of points (global variable)
n <- 10000       # number of samples
tol <- 0.01      # tolerance to update the hull

# Upper hull intersection points, j=1..k-1
z_j <- function(j){
    num <- h(X[j+1]) - h(X[j]) - X[j+1]*hprime(X[j+1]) + X[j]*hprime(X[j])
    den <- hprime(X[j]) - hprime(X[j+1])
    z <- num/den
    return(z)
}

# Vector with the interception points (global variable)
Z <- xmin
Z <- c(Z, sapply(seq(1,k-1), z_j))
Z <- c(Z, xmax)

# Segment parameters (global variables)
M <- hprime(X)
B <- h(X) - M * X

# Piecewise line segments (hull for h(x))
u_k <- function(k){
    function(x) B[k]+M[k]*x
}

# Piecewise exp segments (envelope for g(x))
s_k <- function(k){
    function(x) exp(B[k]+M[k]*x)
}

# Computhe the integral exp(b + m*x) from z1 to z2
integral_exp <- function(b, m, z1, z2){
    if (abs(m) < 1e-10) {
        return(exp(b)*(z2 - z1))
    } else {
        return((exp(b)/m)*(exp(m*z2) - exp(m*z1)))
    }
}

# Compute the area and probability of each segment (global variables)
#areas <- numeric(k)
#for (i in 1:k) {
#   areas[i] <- integral_exp(B[k], M[k], Z[i], Z[i+1])
#}
areas <- mapply(integral_exp, B, M, Z[-length(Z)], Z[-1])
probs <- areas / sum(areas)
cumprobs <- cumsum(probs)

# Add a new point to the hull and update global variables
update_hull <- function(x_new) {

    # Add the new point x_new
    X <<- sort(c(X, x_new))
    k <<- length(X)

    # Upate intersections
    Z <<- xmin
    Z <<- c(Z, sapply(seq(1, k-1), z_j))
    Z <<- c(Z, xmax)

    # Update segment parameters
    M <<- hprime(X)
    B <<- h(X) - M*X

    # Update the areas
    areas <<- mapply(integral_exp, B, M, Z[-length(Z)], Z[-1])

    # Update the probability of areas
    probs <<- areas / sum(areas)
    cumprobs <<- cumsum(probs) # cumulative probablities
}


# MAIN LOOP
nit <- 0                # number of iterations
na  <- 0                # number of accepted
samples <- numeric(n)   # accepted samples

while (na < n){

    # STEP 1: Sample x* from envelope s_k(x)

    # Sample a segment
    # j <- sample(1:k, 1, prob = probs)
    j <- findInterval(runif(1), cumprobs) + 1

    # Segment parameter
    mk <- M[k]

    # Sample x*
    u_sample <- runif(1)
    if (abs(mk) < 1e-10) {
        # uniform
        x_star <- Z[j] + u_sample * (Z[j+1] - Z[j])
    } else {
        # exponential
        x_star <- (1/mk)*log(exp(mk*Z[j]) +
                              u_sample*(exp(mk*Z[j+1]) - exp(mk*Z[j])))
    }

    # STEP 2: Compute u_k(x*)
    u_x_star <- u_k(j)(x_star)

    # STEP 3: Sample w ~ Uniform(0,1)
    w <- runif(1)

    # STEP 4: Rejection test
    h_x_star <- h(x_star)
    if (w <= exp(h_x_star - u_x_star)) {
        # accept x*
        na <- na + 1
        samples[na] <- x_star
    }

    # STEP 5: Update the hull
    if (h_x_star < u_x_star - tol) {
        update_hull(x_star)
    }

    nit <- nit + 1
}

ar <- n/nit # acceptance ratio
cat("Total number of evaluations of g(x):", nit, "\n")
cat("Acceptance ratio:", n / nit, "\n")
cat("Final number of abscissae points:", k, "\n")


# %%
# Plots

# hull for h(x) and g(x) (only for plotting purpose)
x_ <- NULL
u_y <- NULL
s_y <- NULL

for (i in 1:k){
    xk <- seq(Z[i], Z[i+1], length.out = 50)
    x_ <- c(x_, xk)
    u_y <- c(u_y, sapply(xk, u_k(i)))
    s_y <- c(s_y, sapply(xk, s_k(i)))
}

#
# target function g(x), targety density f(x), envelope s(x) and abscissae X
hist(samples, breaks=50, freq = FALSE, ylim=c(0, 0.5),
     main="Adaptive Rejection Sampling",
     xlab="x",
     col=rgb(0.7, 0.7, 0.7, 0.5),
     border="white")
lines(x_, g(x_), lwd=2, col="black")
lines(x_, f(x_), lwd=2, col="steelblue")
lines(x_, s_y, col="red3", lty=2)
points(X, g(X), pch=19, col="forestgreen")
grid(col="gray80", lty="dotted")
legend("topright",
       legend=c("g(x) target", "f(x) normalized", "s(x) envelope", "Abscissae"),
       col=c("black", "steelblue", "red3", "forestgreen"),
       lwd=c(2, 2, 2, NA),
       lty=c(1, 1, 2, NA),
       pch=c(NA, NA, NA, 19),
       bty="n",
       cex=0.8)

#####
# Log-envelope h(x) = log(g(x)) and Upper huldl u(x) and abscissae X
plot(x_, h(x_), type="l", lwd=2,
     main="Log-envelope and Upper Hull",
     xlab="x",
     ylab="h(x)")
lines(x_, u_y, col="red3", lty=2)
points(X, h(X), pch=19, col="forestgreen")
grid(col="gray80", lty="dotted")
legend("bottom",
       legend=c("h(x) = log g(x)", "u(x) upper hull", "Abscissae"),
       col=c("black", "red3", "forestgreen"),
       lwd=c(2, 2, NA),
       lty=c(1, 2, NA),
       pch=c(NA, NA, 19),
       bty="n",
       cex=0.8)
