# Poisson - Local Level Dynamic Model
# Particle Gibbs (PG)
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269–342.
#    https://doi.org/10.1111/j.1467-9868.2009.00736.x
# Author: Cleiton Moya de Almeida

library(invgamma)

#graphics.off()     # close the plots
rm(list = ls())     # clear the environment
#cat("\014")        # clear the console
set.seed(42)
tp <- Matrix::t     # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "poisson_pol2_200" # rds file with data
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
theta_true <- data$theta
Tt <- length(y) # dimension Tt

# Print auxiliary function
printf <- function(...) {
    args <- list(...)
    args[[1]] <- gsub("\n\\s+", " ", args[[1]])
    x <- paste(do.call(sprintf, args), "\n")
    return(cat(x))
}
# Log-sum-exp auxiliary function
logsumexp <- function(x){
  c <- max(x)
  y <- c + log(sum((exp(x-c))))
  return(y)
}

# likelihood density g(theta_t | theta_{t-1})
log_p_yt <- function(yt, theta_t) {
  res <- yt * theta_t - exp(theta_t)
  return(res)
}


#####
# SIMULATION MAIN PARAMETERS

# Prior hyperparameters
# theta_0 ~ N(mu_0, sigma2_0)
mu_0     <- log(y[1] + 0.5)
sigma2_0 <- 10

# W ~ InvGamma(alpha_W, beta_W)
alpha_W  <- 2
beta_W <- 0.1
x_ <- seq(1, 10, 0.1)
y_ <- dinvgamma(x_, alpha_W, beta_W)
plot(x_, y_, type="l")

N <- 10000      # Number of steps
K <- 200        # number of particles
burnin <- 1000  # Number of burn-in steps

# Auxiliary vectors and matrix to store the results
W_hist <- numeric(N)

theta_0_hist <- numeric(N)
theta_star_hist <-matrix(0, N, Tt)

# 1) Initial values (t=0)
theta_star <- numeric(Tt)
theta_star[] <- log(y + 0.5)
theta_0 <- 0
W <- 0.0001

#####
start_time = proc.time() # execution time
for (n in 1:N) {

  # 1. Sample theta_0
  sigma2_0_bar <- (1/sigma2_0 +1/W)^(-1)
  mu_0_bar <- sigma2_0_bar*(mu_0/sigma2_0 + theta_star[1]/W)
  theta_0 <- rnorm(1, mean=mu_0_bar, sd=sqrt(sigma2_0_bar))

  # 2. Sample W
  alpha_W_bar <- alpha_W + Tt/2
  diffs <- theta_star - c(theta_0, theta_star[-Tt])
  beta_W_bar <- beta_W + 0.5 * sum(diffs^2)
  W <- rinvgamma(1, shape=alpha_W_bar, rate=beta_W_bar)
  sd_W <- sqrt(W)

  # 3) SMC step
  # t = 1
  theta_k <- matrix(0, Tt, K)
  log_w_tilde <- matrix(0, Tt, K)

  theta_k[1, ] <- rnorm(K, mean=theta_0, sd=sd_W)
  theta_k[1, K] <- theta_star[1]

  log_w <- log_p_yt(y[1], theta_k[1, ])
  log_w_tilde[1, ] <- log_w - logsumexp(log_w)

  for (t in 2:Tt) {

    # particles resampling
    A <- sample(1:K, K, replace=TRUE, prob=exp(log_w_tilde[t-1,] - max(log_w_tilde[t-1,])))
    #A[K] <- K

    theta_t_k <- theta_k[t-1, ]
    theta_t_k <- rnorm(K, mean=theta_t_k[A], sd=sd_W)
    theta_t_k[K] <- theta_star[t]
    theta_k[t, ] <- theta_t_k

    # weights updating
    log_w_t <- log_p_yt(y[t], theta_k[t, ])
    log_w_tilde[t, ] <- log_w_t - logsumexp(log_w_t)

  }

  # Choose the star particle
  k_final <- sample(1:K, 1, prob=exp(log_w_tilde[Tt,] - max(log_w_tilde[Tt,])))
  theta_star[Tt] <- theta_k[Tt, k_final]

  # Backward pass
  for (t in (Tt-1):1) {
      log_bw <- log_w_tilde[t, ] +
          dnorm(theta_star[t+1], mean=theta_k[t, ], sd=sd_W, log=TRUE)
      log_bw <- log_bw - max(log_bw)
      bw <- exp(log_bw)
      bw <- bw / sum(bw)
      b <- sample(1:K, 1, prob=bw)
      theta_star[t] <- theta_k[t, b]
  }

  theta_0_hist[n] <- theta_0
  W_hist[n] <- W
  theta_star_hist[n, ] <- theta_star

  if (n %% 100 == 0) {
      time <- proc.time()
      elapsed_time <- (time - start_time)[[3]]
      printf("Iteration %d / %d |
           Elapsed time: %.0f s |
           W=%.6f |
           sum_dif2=%.4f |
           range_theta_star=[%.3f, %.3f]",
        n, N, elapsed_time, W, sum(diffs^2), min(theta_star), max(theta_star))
  }
}

end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]

# Results ####
printf("Execution time: %.0f s", elapsed_time)
theta_mean <- colMeans(theta_star_hist[-(1:burnin), ])

#####
printf("Effective Sample Size:")
ess_w <- effectiveSize(mcmc(W_hist[-(1:burnin)]))
printf("\tW: %.0f", ess_w)
printf("Effective Sample Size / second:")
printf("\tW: %.2f", ess_w/elapsed_time)


# Plots ####
# theta_t1 ####
t_observed <- c(50, 100, 150, 200)
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta_mean, type="l", xlab="t", ylab="",
     main="Poisson local level  model")
lines(x, theta_true, col="blue", lwd=2)
legend("topright",
       legend = expression(hat(theta)[t], theta[t]),
       col = c( "black", "blue"),
       lwd = c(1, 2),
       bty = "n")


# Traceplots for W #####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W")

# Traceplot for theta_star_hist #####
par(mfrow = c(2, 2))
for (t in t_observed) {
  plot(theta_star_hist[, t], type="l",
       main=bquote(theta[.(t)]), xlab="n", ylab="")
}
