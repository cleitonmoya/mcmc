# Local Level Dynamic Linear Model - Chan Smoother
#
# Model:
#  y_t ~ theta_t + nu_t,               nu_t ~ N(0, V)
#  theta_t = theta_{t-1} + omega_t, omega_t ~ N(0, W)
#
# Prior:
#  theta_0 | D_0 ~ N(mu_0, sigma2_0)
#
# Posterior:
#  theta_t | D_T ~ N(theta_hat, P^(-1))
#    P: Precision matrix
#
# Reference: Chan, J. C. C., & Jeliazkov, I. (2009). Efficient simulation
# and integrated likelihood estimation in state space models.
# International Journal of Mathematical Modelling and Numerical Optimisation,
# 1(1/2), 101. https://doi.org/10.1504/IJMMNO.2009.030090
#
# Author: Cleiton Moya de Almeida

library(Matrix)    # sparse matrix manipulation
#library(sparseinv) # efficient inverse of sparce matrix

graphics.off()    # close the plots
rm(list = ls())    # clear the environment
cat("\014")       # clear the console
tp <- Matrix::t    # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback
set.seed(42)


# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "normal_local_level_sim" # rds file with data
data <- readRDS(paste("../data/", source, ".rds", sep=""))

y <- data$y
Tt <- length(y) # dimension T
theta_true <- data$theta
W <- data$W
V <- data$V

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))

# Prior hyperparameters
# theta_0 ~ N(m0, C0)
m0 <- 0
C0 <- 10

start_time = proc.time() # execution time

# CHAN METHOD ####
# 1. Build the Prior Precision Matrix K
# H \vect{theta} = \vect{gamma} + \vect{\omega}
y_matr <- Matrix(y, ncol=1)

# Prior precision matrix K - symmetric sparse banded (dsCMatrix)
D <- C0 + W
sub_diag <- rep(-1/W, Tt-1)
main_diag <- c(1/D+1/W, rep(2/W, Tt-2), 1/W)
K <- bandSparse(n=Tt, k=c(0, -1), diagonals=list(main_diag, sub_diag), symmetric = TRUE)

# 2. Compute the Posterior Precision Matrix P
# Observational variance matrix - Vcal: \mathcal{V} (TxT)
invVcal <- .sparseDiagonal(n=Tt, x=1/V)

# Posterior precision matrix
P <- forceSymmetric(K + invVcal)


# 3. Smoothing through Cholesky decomposition
# P * hat_theta = invVcal*y + gamma
gamma_over_D_matr <- Matrix(data=c(m0/D, rep(0, Tt-1)), ncol=1, sparse=TRUE)
b <- invVcal %*% y_matr + gamma_over_D_matr

# Cholesky factorization considering sparse symmetric matrix (dsCMatrix)
# By default, perm = TRUE, LDL = !super, super = FALSE, Imult = 0.
# R consider the simplicial factorization, returning object of class dCHMsimpl.
# Permutation is not necessary for tridiagonal blocked matrix
Ch_factor <- Cholesky(P, perm=FALSE)

# 3.1 forward-substitution
z <- Matrix::solve(Ch_factor, b, system="LD")

# 3.2 back-substitution
theta_hat <- as.numeric(Matrix::solve(Ch_factor, z, system="Lt"))

# Equivalent to:
# theta_hat <- as.numeric(Matrix::solve(Ch_factor, b, system="A"))


# Variance estimation - O(T^2)
invP <- Matrix::solve(Ch_factor, Diagonal(Tt), system = "A")
S <- as.numeric(diag(invP))

# Alterantive: Use Takahashi inverse, O(Tp^2) instead of O(T^2)
# But would need to use Cholesky with perm=TRUE, LDL=FALSE
# P_inv_partial <- sparseinv::Takahashi_Davis(P, cholQp = Ch_factor)
# S <- (diag(P_inv_partial)


# RESULTS ####

# Execution time
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.2f ms", elapsed_time*1000)


# Plot - Smoothing ####
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="n", xlab="t", ylab="", col="gray",
     main="Local Level DLM - Chan Smoothing")
lines(x, y)
polygon(c(x, rev(x)),
       c(theta_hat+ 3*sqrt(S), rev(theta_hat - 3*sqrt(S))),
       col = adjustcolor("gray", alpha.f = 0.5),
       border = NA)
points(x, y, pch = 20)
lines(x, theta_true, type="l", xlab="t", ylab="", col="blue")
lines(x, theta_hat, col="red")
legend("topright",
       legend = expression(y[t], theta[t], hat(theta)[t], S[t]^"*"),
       col = c("black", "blue", "red", adjustcolor("gray", alpha.f = 0.5)),
       lty = c(NA, 1, 1, NA),
       lwd = c(NA, 1, 1, NA),
       pch = c(20, NA, NA, 15),
       pt.cex = c(1, NA, NA, 3),
       bty = "n")
