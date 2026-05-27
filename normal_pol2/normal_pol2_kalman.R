# Forward Filtering - Backward Sampling

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- base::t       # alias to transpose function
options(error = function() traceback(2)) # more informative traceback


# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "normal_pol2_sim1" # csv file with data
df <- read.table(paste("../data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
theta1_true <- df$theta1
theta2_true <- df$theta2
T <- length(y)


V <- 0.5
W <- matrix(c(1, 0, 0, 1), nrow=2)

theta_t1 <- 1
theta_t2 <- 1
W1 <- 0.02
W2 <- 0.03

theta1 <- numeric(T)
theta2 <- numeric(T)

F <- matrix(c(1, 0), nrow=2)
G <- matrix(c(1, 1, 0, 1), nrow=2, byrow=TRUE)
W <- diag(c(W1, W2), nrow=2, ncol=2)
theta <- cbind(theta1, theta2)


# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- y[1]
sigma2_01 <- 10

# theta_02 ~ N(mu_2, sigma2_02)
mu_02     <- 0.01
sigma2_02 <- 10

# Forward Filtering

a <- array(data=NA, dim=c(2,1,T))
R <- array(data=NA, dim=c(2,2,T))
m <- array(data=NA, dim=c(2,1,T))
C <- array(data=NA, dim=c(2,2,T))
B <- array(data=NA, dim=c(2,2,T))

m_t <- as.matrix(c(mu_01, mu_02))
C_t <- diag(c(sigma2_01, sigma2_02), nrow=2, ncol=2)

for (t in 1:T) {

    y_t <- y[t]

    # Prior: (theta_t | D_{t-1}) ~ N[a_t, R_t]
    a_t <- G %*% m_t
    R_t <- G %*% C_t %*% tp(G) + W
    a[,,t] <- a_t
    R[,,t] <- R_t

    # One-step-ahead forecast
    f_t <- drop(tp(F) %*% a_t)
    Q_t <- drop(tp(F) %*% R_t %*% F + V)

    # Posterior: (theta_t | D_t) ~ N[m_t, C_t]
    e_t <- y_t - f_t
    A_t <-  (R_t %*% F) / Q_t
    m_t <- a_t + A_t * e_t
    C_t <- R_t - (A_t * Q_t) %*% tp(A_t)
    m[,,t] <- m_t
    C[,,t] <- C_t

}

# Compute B_t
for (t in 1:(T-1)) {
    B[,,t] <- C[,,t] %*% tp(G) %*% solve(R[,,t+1])
}

# Smoothing
# (theta_t| D_T) ~ N[a*_t, R*_t], t = 1, ... T
a_star <- array(data=NA, dim=c(2,1,T))
R_star <- array(data=NA, dim=c(2,2,T))
a_star[,,T] <- m[,,T]
R_star[,,T] <- C[,,T]

for (t in seq(T-1,1)) {
    a_star[,,t] <- m[,,t] + B[,,t] %*% (a_star[,,t+1] - a[,,t+1])
    R_star[,,t] <- C[,,t] + B[,,t] %*% (R_star[,,t+1] - R[,,t+1]) %*% tp(B[,,t])
}



x <- 1:T

#####
# theta1 filtered ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_true, type="l", xlab="t", ylab="")
lines(x, m[1,,], col="blue")
polygon(c(x, rev(x)),
        c(m[1,,] + 3*sqrt(C[1,1,]), rev(m[1,,] - 3*sqrt(C[1,1,]))),
        col = adjustcolor("steelblue", alpha.f = 0.3),
        border = NA)


# theta1 smoothed ####
plot(x, theta1_true, type="l", xlab="t", ylab="")
lines(x, a_star[1,,], col="blue")
polygon(c(x, rev(x)),
        c(a_star[1,,] + 3*sqrt(R_star[1,1,]),
          rev(a_star[1,,] - 3*sqrt(R_star[1,1,]))),
        col = adjustcolor("steelblue", alpha.f = 0.3),
        border = NA)


# theta2 filtered #####
plot(x, theta2_true, type="l", xlab="t", ylab="")
lines(x, m[2,,], col="blue")
polygon(c(x, rev(x)),
        c(m[2,,] + 3*sqrt(C[2,2,]), rev(m[2,,] - 3*sqrt(C[2,2,]))),
        col = adjustcolor("steelblue", alpha.f = 0.3),
        border = NA)


# theta2 smoothed ####
plot(x, theta2_true, type="l", xlab="t", ylab="")
lines(x, a_star[2,,], col="blue")
polygon(c(x, rev(x)),
        c(a_star[2,,] + 3*sqrt(R_star[2,2,]),
          rev(a_star[2,,] - 3*sqrt(R_star[2,2,]))),
        col = adjustcolor("steelblue", alpha.f = 0.3),
        border = NA)
