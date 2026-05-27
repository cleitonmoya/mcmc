# Rejection method for sampling from a target disribution
# Envelope function: Uniform
# Author: Cleiton Moya de Almeida

# Clear the plots, environment and console
graphics.off()
rm(list = ls())
cat("\014")

set.seed(42)

# g:    target function  (R function) (does not need to be normalized)
# k:    constant that envelopes the target function
# xmin: lower support of the target function
# xmax: upper support of the target function
# n:    number of samples
g <- function(x) dnorm(x)
k <- 0.45
n <- 10000
xmin <- -5
xmax <- 5

# Auxiliary variables
nit <- 0  # number of iterations
na <- 0   # number accepted observations
X <- NULL # accepted observations vector

while(na < n){

    x <- runif(1, min=xmin, max=xmax)
    u <- runif(1) # u ~ Unif(0,1)

    w <- g(x)/k

    if (u <= w){
        X <- c(X, x)
        na <- na + 1
    }
    nit <- nit + 1
}

ar <- n/nit # acceptance ratio


#####
# Results

print(noquote(sprintf("Accepance ratio: %.2f", ar)))

x_ <- seq(xmin,xmax,0.1)
hist(X, freq=FALSE, ylim=c(0,k))          # sampled observations
lines(x_, rep(k,length(x_)), col="green") # envelope function
lines(x_, sapply(x_, g),col="blue")       # target function
lines(density(X), col="red")              # estimated density
