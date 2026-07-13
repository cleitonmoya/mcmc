# Local Level Dynamic Linear Model (DLM)
#
# Model:
#   Y_t = theta_t + \nu_t, nu_t ~ N(O, V)
#   theta_t = theta_{t-1} + \omega_t, \omega_t ~ N(O, W)
#
# Prior:
#   theta_0 | D_0 - N(m_0, C_0)
#
# Algorithms:
#   - Univariate forward (Kalman) filter
#   - Univariate (Kalman) smoothing
#
# Reference:
#   West, M., & Harrison, J. (1997).
#   Bayesian forecasting and dynamic models. 2ed.
#   New York, NY: Springer New York.
#
# Author: Cleiton Moya de Almeida

graphics.off()  # close the plots
rm(list = ls()) # clear the environment
cat("\014")     # clear the console
tp <- base::t   # alias to transpose function
options(error = function() traceback(2)) # more informative traceback
set.seed(42)

# Change the directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))

# Load the data
source <- "normal_local_level_sim"
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
Tt <-length(y)
theta_true <- data$theta

# Fixed parameters
W <- data$W
V <- data$V

# Prior: theta_0 ~ N(m_0, C_0)
m0 <- 0
C0 <- 10

#####
# Forward Filtering (W&H Session 4.3, 4.2)
#
# Posterior at t-1:
#   theta_{t-1} | D_{t-1} ~ N(m_{t-1}, C_{t-1})
#
# Prior at t:
#   theta_t | D_{t-1})    ~ N(a_t, R_t)
#     a_t = m_{t-1}
#     R_t = C_{t-1} + W_t
#
# One-step forecast:
#   Y_t | D_{t-1}         ~ N(f_t, Q_t)
#     f_t = m_{t-1}
#     Q_t = R_t + V_t

# Posterior at t:
#   theta_t | D_t         ~ N(m_t, C_t)
#     m_t = a_t + A_t*e_t
#     C_t = A_t*V_t = (1-A_t)*R_t
#       e_t = Y_t - a_t
#       A_t = R_t / Q_t

# Auxiliary auxiiary vectors
theta <- numeric(Tt)
a <- numeric(Tt)
m <- numeric(Tt)   # m_t = E[theta_t | D_t]
C <- numeric(Tt)   # C_t = Var[theta_t | D_t]
R <- numeric(Tt)
B <- numeric(Tt)

start_time = proc.time() # execution time
# For t=1
a1 <- m0
R1 <- C0 + W
Q1 <- R1 + V
A1 <- R1 / Q1
e1 <- y[1] - a1
m[1] <- a1 + A1 * e1
C[1] <- (1 - A1) * R1
R[1] <- R1
a[1] <- a1

for (t in 2:Tt) {
    at <- m[t-1]           # prior mean
    Rt <- C[t-1] + W       # prior variance

    Qt <- Rt + V           # forecast variance
    At <- Rt / Qt          # adaptive coefficient
    et <- y[t] - at        # forecast error

    m[t] <- at + At * et   # posterior mean
    C[t] <- (1 - At) * Rt  # posterior variance
    R[t] <- Rt
    a[t] <- at

    # Backward matrix (for smoothing alghorithm)
    B[t-1] <- C[t-1] / R[t]
}


#####
# Smoothing (W&H Session 4.7 Theorem 4.4)
#

# Compute the smoothed marginal distributions
#   (theta_t| D_T) ~ N[s_t, S_t], t = 1,...T
s <- numeric(Tt)
S <- numeric(Tt)

# for t=T
s[Tt] <- m[Tt]
S[Tt] <- C[Tt]

for (t in seq(Tt-1,1)) {
    s[t] <- m[t] + B[t] * (s[t+1] - a[t+1])
    S[t] <- C[t] + B[t]**2 * (S[t+1] - R[t+1])
}

end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.2f ms", elapsed_time*1000)


#####
# Plots
x <- 1:Tt

# Filtering
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="n", xlab="t", ylab="", col="gray",
     main="Local Level DLM - Forward (Kalman) Filtering")
lines(x, y)
polygon(c(x, rev(x)),
        c(m+ 3*sqrt(C), rev(m - 3*sqrt(C))),
        col = adjustcolor("gray", alpha.f = 0.5),
        border = NA)
points(x, y, pch = 20)
lines(x, theta_true, type="l", xlab="t", ylab="", col="blue")
lines(x, m, col="red")
legend("topright",
       legend = expression(y[t], theta[t], hat(theta)[t], C[t]),
       col = c("black", "blue", "red", adjustcolor("gray", alpha.f = 0.5)),
       lty = c(NA, 1, 1, NA),
       lwd = c(NA, 1, 1, NA),
       pch = c(20, NA, NA, 15),
       pt.cex = c(1, NA, NA, 3),
       bty = "n")


# Smoothing
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="n", xlab="t", ylab="", col="gray",
     main="Local Level DLM - Kalman Smoothing")
lines(x, y)
polygon(c(x, rev(x)),
        c(s+ 3*sqrt(S), rev(s - 3*sqrt(S))),
        col = adjustcolor("gray", alpha.f = 0.5),
        border = NA)
points(x, y, pch = 20)
lines(x, theta_true, type="l", xlab="t", ylab="", col="blue")
lines(x, s, col="red")
legend("topright",
       legend = expression(y[t], theta[t], hat(theta)[t], S[t]^"*"),
       col = c("black", "blue", "red", adjustcolor("gray", alpha.f = 0.5)),
       lty = c(NA, 1, 1, NA),
       lwd = c(NA, 1, 1, NA),
       pch = c(20, NA, NA, 15),
       pt.cex = c(1, NA, NA, 3),
       bty = "n")
