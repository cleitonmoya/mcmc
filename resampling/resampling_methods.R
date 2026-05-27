# Resampling Methods

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- base::t       # alias to transpose function
options(error = function() traceback(2)) # more informative traceback

set.seed(42)

# Multinomial resampling
y <- c(10, 20, 30, 40, 50, 60)
p <- c(0.2, 0.1, 0.1, 0.2, 0.3, 0.1)
N <- 10 # number of samples


multinomial_resampling <- function(y, p, N) {
    Q <- cumsum(p)
    U <- runif(N)
    Y <- numeric(N)
    for (t in 1:N) {
        j <- 1
        while(U[t] > Q[j]) {
            j <- j + 1
        }
        Y[t] <- y[j]
    }
    return(Y)
}

multinomial_resampling_smart <- function(y, p, N) {
    Q <- cumsum(p)
    U <- sort(runif(N))
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

#Y0 <- sample(y, size=N, replace=TRUE, prob=p)
#Y1 <- multinomial_resampling(y, p, N)
#Y2 <- multinomial_resampling_smart(y, p, N)
Y3 <- stratified_resampling(y, p, N)
print(as.integer(table(Y3))/N)
