# Poisson - 2nd Order Polynomial Dynamic Model
# Particle Gibbs (PG) with Backward Sampling
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269-342.
#    https://doi.org/10.1111/j.1467-9868.2009.00736.x
# Author: Cleiton Moya de Almeida

library(invgamma)

#graphics.off()      # close the plots
rm(list = ls())     # clear the environment
#cat("\014")         # clear the console
set.seed(42)
options(error = function() traceback(2))

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "poisson_model2"
data   <- readRDS(paste("../data/", source, ".rds", sep=""))
y           <- data$y
theta1_true <- data$theta
Tt          <- length(y)
t_observed  <- c(50, 75, 100, 150)

# Print auxiliary function
printf <- function(...) {
  x <- paste(sprintf(...), "\n")
  return(cat(x))
}

# Log-sum-exp
logsumexp <- function(x) {
  cc <- max(x)
  return(cc + log(sum(exp(x - cc))))
}

# Log-likelihood: log g(y_t | theta_t1)
log_p_yt <- function(yt, theta_t1) {
  res <- yt * theta_t1 - exp(theta_t1)
  res[!is.finite(res)] <- -Inf
  return(res)
}

#####
# PRIOR HYPERPARAMETERS

# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- log(y[1] + 0.5)
sigma2_01 <- 10

# theta_02 ~ N(mu_02, sigma2_02)
mu_02     <- 0        # BUG 1 CORRIGIDO: era mu_02=1, tendência inicial deve ser 0
sigma2_02 <- 1

# W1 ~ InvGamma(alpha_W1, beta_W1)
alpha_W1 <- 2
beta_W1  <- 0.01

# W2 ~ InvGamma(alpha_W2, beta_W2)
alpha_W2 <- 2
beta_W2  <- 0.001

N      <- 1000   # BUG 2 CORRIGIDO: era N=100, insuficiente
K      <- 200    # número de partículas
burnin <- 200    # BUG 3 CORRIGIDO: era burnin=1000 > N=100 -> colMeans retornava NaN

# Armazenamento
W1_hist          <- numeric(N)
W2_hist          <- numeric(N)
theta_01_hist    <- numeric(N)
theta_02_hist    <- numeric(N)
theta1_star_hist <- matrix(0, N, Tt)
theta2_star_hist <- matrix(0, N, Tt)

#####
# INITIAL VALUES

theta1_star <- stats::filter(log(y + 0.5), rep(1/5, 5), sides=2)
theta1_star[is.na(theta1_star)] <- log(y[is.na(theta1_star)] + 0.5)
theta1_star <- as.numeric(theta1_star)
theta2_star <- c(0, diff(theta1_star))

theta_01 <- log(y[1] + 0.5)
theta_02 <- 0
W1 <- 0.01
W2 <- 0.001

#####
start_time <- proc.time()

for (n in 1:N) {

  # ------------------------------------------------------------------
  # STEP 1: Sample theta_02 | theta_01, theta1_star, theta2_star, W1, W2
  # ------------------------------------------------------------------
  sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
  mu_02_bar     <- sigma2_02_bar * (
    (theta1_star[1] - theta_01) / W1 +
    theta2_star[1]              / W2 +
    mu_02 / sigma2_02
  )
  theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))

  # ------------------------------------------------------------------
  # STEP 2: Sample theta_01 | theta_02, theta1_star, W1
  # ------------------------------------------------------------------
  sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
  mu_01_bar     <- sigma2_01_bar * (
    mu_01 / sigma2_01 +
    (theta1_star[1] - theta_02) / W1
  )
  theta_01 <- rnorm(1, mean=mu_01_bar, sd=sqrt(sigma2_01_bar))

  # ------------------------------------------------------------------
  # STEP 3: Sample W1 | theta_01, theta_02, theta1_star, theta2_star
  # omega_t1 = theta_t1* - theta_{t-1,1}* - theta_{t-1,2}*
  # BUG 4 CORRIGIDO: priori Gamma sobre phi substituída por InvGamma sobre W
  # ------------------------------------------------------------------
  dif1   <- theta1_star - c(theta_01, theta1_star[-Tt])
  diffs1 <- dif1 - c(theta_02, theta2_star[-Tt])
  alpha_W1_bar <- alpha_W1 + Tt/2
  beta_W1_bar  <- beta_W1  + 0.5 * sum(diffs1^2)
  W1 <- rinvgamma(1, shape=alpha_W1_bar, rate=beta_W1_bar)

  # ------------------------------------------------------------------
  # STEP 4: Sample W2 | theta_02, theta2_star
  # omega_t2 = theta_t2* - theta_{t-1,2}*
  # BUG 4 CORRIGIDO: priori Gamma sobre phi substituída por InvGamma sobre W
  # ------------------------------------------------------------------
  diffs2 <- theta2_star - c(theta_02, theta2_star[-Tt])
  alpha_W2_bar <- alpha_W2 + Tt/2
  beta_W2_bar  <- beta_W2  + 0.5 * sum(diffs2^2)
  W2 <- rinvgamma(1, shape=alpha_W2_bar, rate=beta_W2_bar)

  sd_W1 <- sqrt(W1)
  sd_W2 <- sqrt(W2)

  # ------------------------------------------------------------------
  # STEP 5: Conditional SMC (forward pass)
  # Estado: (theta_t1, theta_t2) — trajetória condicionada na coluna K
  # BUG 5 CORRIGIDO: k_star removido — com backward sampling usa-se
  #   coluna fixa K; A[k_star]<-k_star também removido
  # ------------------------------------------------------------------
  theta_1_k   <- matrix(0, Tt, K)
  theta_2_k   <- matrix(0, Tt, K)
  log_w_tilde <- matrix(0, Tt, K)

  # t = 1
  theta_1_k[1, ] <- rnorm(K, mean=theta_01 + theta_02, sd=sd_W1)
  theta_2_k[1, ] <- rnorm(K, mean=theta_02,            sd=sd_W2)

  # Fixar trajetória condicionada na coluna K
  theta_1_k[1, K] <- theta1_star[1]
  theta_2_k[1, K] <- theta2_star[1]

  log_w_1 <- log_p_yt(y[1], theta_1_k[1, ])
  log_w_tilde[1, ] <- log_w_1 - logsumexp(log_w_1)

  # t = 2, ..., T
  for (t in 2:Tt) {

    # Reamostragem — sem forçar A[K]=K (backward sampling dispensa)
    A <- sample(1:K, K, replace=TRUE,
                prob=exp(log_w_tilde[t-1, ] - max(log_w_tilde[t-1, ])))

    # Proposta bootstrap
    theta_t1_k <- rnorm(K,
                        mean=theta_1_k[t-1, A] + theta_2_k[t-1, A],
                        sd=sd_W1)
    theta_t2_k <- rnorm(K,
                        mean=theta_2_k[t-1, A],
                        sd=sd_W2)

    # Fixar trajetória condicionada na coluna K
    theta_t1_k[K] <- theta1_star[t]
    theta_t2_k[K] <- theta2_star[t]

    theta_1_k[t, ] <- theta_t1_k
    theta_2_k[t, ] <- theta_t2_k

    # Atualizar pesos
    log_w_t <- log_p_yt(y[t], theta_1_k[t, ])
    log_w_tilde[t, ] <- log_w_t - logsumexp(log_w_t)
  }

  # ------------------------------------------------------------------
  # STEP 6: Backward sampling
  # Pesos backward: w_t^k * p(theta_{t+1}* | theta_t^k)
  #   = w_t^k * N(theta_{t+1,1}* | theta_t1^k + theta_t2^k, W1)
  #           * N(theta_{t+1,2}* | theta_t2^k, W2)
  # ------------------------------------------------------------------

  # Selecionar theta_T
  k_final <- sample(1:K, 1,
                    prob=exp(log_w_tilde[Tt, ] - max(log_w_tilde[Tt, ])))
  theta1_star[Tt] <- theta_1_k[Tt, k_final]
  theta2_star[Tt] <- theta_2_k[Tt, k_final]

  # Backward pass: t = T-1, ..., 1
  for (t in (Tt-1):1) {
    log_bw <- log_w_tilde[t, ] +
              dnorm(theta1_star[t+1],
                    mean=theta_1_k[t, ] + theta_2_k[t, ],
                    sd=sd_W1, log=TRUE) +
              dnorm(theta2_star[t+1],
                    mean=theta_2_k[t, ],
                    sd=sd_W2, log=TRUE)
    log_bw <- log_bw - max(log_bw)
    bw     <- exp(log_bw)
    bw     <- bw / sum(bw)
    b      <- sample(1:K, 1, prob=bw)
    theta1_star[t] <- theta_1_k[t, b]
    theta2_star[t] <- theta_2_k[t, b]
  }

  # ------------------------------------------------------------------
  # Armazenar resultados
  # ------------------------------------------------------------------
  theta_01_hist[n]      <- theta_01
  theta_02_hist[n]      <- theta_02
  W1_hist[n]            <- W1
  W2_hist[n]            <- W2
  theta1_star_hist[n, ] <- theta1_star
  theta2_star_hist[n, ] <- theta2_star

  if (n %% 100 == 0) {
    elapsed <- (proc.time() - start_time)[[3]]
    printf("n=%d/%d | t=%.0fs | W1=%.6f | W2=%.6f | sum_diffs1^2=%.4f | sum_diffs2^2=%.6f",
           n, N, elapsed, W1, W2, sum(diffs1^2), sum(diffs2^2))
  }
}

elapsed_time <- (proc.time() - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)

# ------------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------------
theta1_mean <- colMeans(theta1_star_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_star_hist[-(1:burnin), ])
W1_mean <- mean(W1_hist[-(1:burnin)])
W2_mean <- mean(W2_hist[-(1:burnin)])
printf("W1 posterior mean: %.6f", W1_mean)
printf("W2 posterior mean: %.6f", W2_mean)

# ------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------
x <- 1:Tt

# theta_t1
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(x, theta1_true, type="l", xlab="t", ylab="",
     main="Poisson local trend model - theta_t1", col="blue", lwd=2)
lines(x, theta1_mean, col="black")
legend("topright",
       legend=expression(theta[t1], hat(theta)[t1]),
       col=c("blue","black"), lwd=c(2,1), bty="n")

# theta_t2
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(x, theta2_mean, type="l", xlab="t", ylab="",
     main="Poisson local trend model - theta_t2")
legend("topright",
       legend=expression(hat(theta)[t2]),
       col="black", lwd=1, bty="n")

# Traceplot W1
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(W1_hist[-(1:burnin)], type="l", xlab="n", ylab="W1",
     main="Traceplot of W1")

# Traceplot W2
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(W2_hist[-(1:burnin)], type="l", xlab="n", ylab="W2",
     main="Traceplot of W2")

# Traceplots theta1_star
par(mfrow=c(2,2))
for (t in t_observed) {
  plot(theta1_star_hist[, t], type="l",
       main=bquote(theta[list(.(t),1)]), xlab="n", ylab="")
}

# Traceplots theta2_star
par(mfrow=c(2,2))
for (t in t_observed) {
  plot(theta2_star_hist[, t], type="l",
       main=bquote(theta[list(.(t),2)]), xlab="n", ylab="")
}
