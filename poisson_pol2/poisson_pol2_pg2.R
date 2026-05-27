# Poisson - 2nd Order Polynomial Dynamic Model
# Particle Gibbs (PG)
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269–342.
#    https://doi.org/10.1111/j.1467-9868.2009.00736.x
# Author: Cleiton Moya de Almeida


#graphics.off()      # close the plots
rm(list = ls())     # clear the environment
#cat("\014")         # clear the console
set.seed(42)
tp <- Matrix::t     # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "poisson_model2" # rds file with data
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
theta1_true <- data$theta
Tt <- length(y) # dimension T
t_observed <- c(50, 75, 100, 150)

# Print auxiliary function
printf <- function(...) {
  x = paste(sprintf(...),"\n")
  return(cat(x))
}


logsumexp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(-Inf)
  c <- max(x)
  return(c + log(sum(exp(x - c))))
}

log_p_yt <- function(yt, theta_t1) {
  res <- yt * theta_t1 - exp(theta_t1)
  res[!is.finite(res)] <- -Inf
  return(res)
}


# Prior hyperparameters
# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- log(y[1] + 0.5)
sigma2_01 <- 10

# theta_02 ~ N(mu_2, sigma2_02)
mu_02     <- 1
sigma2_02 <- 10

# phi1 = W1^(-1) ~ Gamma(nu_01, eta_01)
nu_01  <- 0.1
eta_01 <- 0.1

# phi2 = W2^(-1) ~ Gamma(nu_02, eta_02)
nu_02  <- 0.1
eta_02 <- 0.1

N <- 100     # Number of steps
K <- 100        # number of particles
burnin <- 1000  # Number of burn-in steps


# Auxiliary vectors and matrix to store the results
W1_hist <- numeric(N)
W2_hist <- numeric(N)

theta_01_hist <- numeric(N)
theta_02_hist <- numeric(N)
theta1_star_hist <-matrix(0, N, Tt)
theta2_star_hist <-matrix(0, N, Tt)


# 1) Initial values (t=0)
theta1_star <- stats::filter(log(y + 0.5), rep(1/5, 5), sides=2)
theta1_star[is.na(theta1_star)] <- log(y[is.na(theta1_star)] + 0.5)  # bordas
theta1_star <- as.numeric(theta1_star)
theta2_star <- c(0, diff(theta1_star))  # diferença primeira como tendência inicial

theta_01 <- 0.1
theta_02 <- 0.1

W1 <- 0.001
W2 <- 0.001

k_star <- sample(1:K, 1)


start_time = proc.time() # execution time
for (n in 1:N) {

  # Sample theta_02
  sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
  mu_02_bar <- sigma2_02_bar*((theta1_star[1]-theta_01)/W1 +
                                theta2_star[1]/W2 + mu_02/sigma2_02)
  theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))

  # Sample theta_01
  sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
  mu_01_bar <- sigma2_01_bar*(mu_01/sigma2_01 + (theta1_star[1]-theta_02)/W1)
  theta_01 <- rnorm(1, mean=mu_01_bar, sd=sqrt(sigma2_01_bar))


  # 2) Sample W | theta_1, theta_2
  # Sample phi1
  nu_01_bar <- nu_01 + Tt/2
  dif1 <- theta1_star - c(theta_01, theta1_star[-Tt])
  dif2 <- dif1 - c(theta_02, theta2_star[-Tt])
  eta_01_bar <- eta_01 + 0.5 * sum(dif2^2)
  phi1 <- rgamma(1, nu_01_bar, eta_01_bar)
  W1 <- 1/phi1


  # Sample phi2
  nu_02_bar <- nu_02 + Tt/2
  diffs2 <- theta2_star - c(theta_02, theta2_star[-Tt])
  eta_02_bar <- eta_02 + 0.5 * sum(diffs2^2)
  phi2 <- rgamma(1, nu_02_bar, eta_02_bar)
  W2 <- 1/phi2

  # 3) SMC step

  # t = 1
  theta_1_k <- matrix(0, Tt, K)
  theta_2_k <- matrix(0, Tt, K)
  log_w_tilde <- matrix(0, Tt, K)


  theta_1_k[1, ] <- rnorm(K, mean=theta_01 + theta_02, sd=sqrt(W1))
  theta_1_k[1, k_star] <- theta1_star[1]
  theta_2_k[1, ] <- rnorm(K, mean=theta_02, sd=sqrt(W2))
  theta_2_k[1, k_star] <- theta2_star[1]

  log_w_1 <- log_p_yt(y[1], theta_1_k[1, ])
  log_w_tilde[1, ] <- log_w_1 - logsumexp(log_w_1)


  for (t in 2:Tt) {

    if (n %% 1000 == 0) {
      time <- proc.time()
      elapsed_time <- (time - start_time)[[3]]
      printf("Iteration %d / %d | Elapsed time: %.0f s", n, N, elapsed_time)
    }

    # particles resampling
    A <- sample(1:K, K, replace=TRUE, prob=exp(log_w_tilde[t-1,] - max(log_w_tilde[t-1,])))
    A[k_star] <- k_star

    theta_t1_k <- theta_1_k[t-1, ]
    theta_t2_k <- theta_2_k[t-1, ]

    theta_t1_k <- rnorm(K, mean=theta_t1_k[A] + theta_t2_k[A], sd=sqrt(W1))
    theta_t1_k[k_star] <- theta1_star[t]

    theta_t2_k <- rnorm(K, theta_t2_k[A], sd=sqrt(W2))
    theta_t2_k[k_star] <- theta2_star[t]

    theta_1_k[t, ] <- theta_t1_k
    theta_2_k[t, ] <- theta_t2_k

    # weights updating
    log_w_t <- log_p_yt(y[t], theta_1_k[t, ])
    log_w_tilde[t, ] <- log_w_t - logsumexp(log_w_t)

  }

  # Choose the star particle
  k_star <- sample(1:K, 1, prob=exp(log_w_tilde[Tt,] - max(log_w_tilde[Tt,])))
  theta1_star <- theta_1_k[, k_star]
  theta2_star <- theta_2_k[, k_star]

  theta_01_hist[n] <- theta_01
  theta_02_hist[n] <- theta_02
  W1_hist[n] <- W1
  W2_hist[n] <- W2
  theta1_star_hist[n, ] <- theta1_star
  theta2_star_hist[n, ] <- theta2_star

}
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]

# Results ####
printf("Execution time: %.0f s", elapsed_time)

theta1_mean <- colMeans(theta1_star_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_star_hist[-(1:burnin), ])


# Plots ####
# theta_t1 ####
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_true, type="l", xlab="t", ylab="",
     main="Poisson local trend  model", ylim=c(-1.5, 1.5))
lines(x, theta1_mean, col="blue", lwd=2)
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(2, 2),
       bty = "n")


# theta_t2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta2_mean, type="l", xlab="t", ylab="",
     main="Poisson local trend  model", ylim=c(-1.5, 1.5))
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(2, 2),
       bty = "n")


# Traceplots for W1, W2 #####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W1_hist, type="l", xlab="n", ylab="W", main="Traceplot of W1")
plot(W2_hist, type="l", xlab="n", ylab="W", main="Traceplot of W2")

# Traceplot for theta1_star_hist #####
par(mfrow = c(2, 2))
for (t in t_observed) {
  plot(theta1_star_hist[, t], type="l",
       main=bquote(theta[list(.(t), 1)]), xlab="n", ylab="")
}

# Traceplot for theta2_star_hist #####
par(mfrow = c(2, 2))
for (t in t_observed) {
  plot(theta2_star_hist[, t], type="l",
       main=bquote(theta[list(.(t), 2)]), xlab="n", ylab="")
}

