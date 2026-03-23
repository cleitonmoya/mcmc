# Metropolis
# Author: Cleiton Moya de Almeida

# Clear the plots, environment and console
graphics.off()
rm(list = ls())
cat("\014")

set.seed(42)

### SIMULATION PARAMETERS
# Target density - Gumbel distribution

mu_t <- 1
beta_t <- 2
# Target density
f <- function(x) {
    z <- (x - mu_t) / beta_t
    (1/beta_t) * exp(-(z + exp(-z)))
}

# Target function (unormalized)
g <- function(x){
    z <- (x - mu_t)/beta_t
    exp(-(z + exp(-z)))
}

# Sample from proposed distribion
sample_xprop <- function(mu, sigma_p){
    x <- rnorm(1, mean=mu, sd=sigma_p)
    return(x)
}

N <- 10100 # number of steps
burnin <- 100

# Propsosed distribution
sigma_p <- 0.5

x <- 0          # initial value
na <- 0         # number of accepted value
X <- numeric(N) # samples

# Main loop
for (n in 1:N){

    # Draw a proposed value
    x_prop <- sample_xprop(x, sigma_p)

    # Accpeting probability
    p <- g(x_prop) / g(x)
    u <- runif(1)
    if (u < p){
        x <- x_prop
        na <- na + 1
    }
    X[n] <- x
}

ar <- na/N
print(noquote(sprintf("Accepance ratio: %.2f", ar)))

#####
x_ <- seq(-5,10,0.1)
samples <- X[burnin:N]

plot(1:N, X, type="l")

hist(samples, breaks=50, freq=FALSE, xlim=c(-5,15), ylim=c(0,0.5))
lines(x_, g(x_), col="red3")
lines(x_, f(x_), col="steelblue")
