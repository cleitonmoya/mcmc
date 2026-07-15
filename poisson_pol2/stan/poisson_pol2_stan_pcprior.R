# Poisson - 2nd Order Polynomial Dynamic Model
# MCMC: NUTS (Stan)
# Author: Cleiton Moya de Almeida

graphics.off()    # close the plots
rm(list = ls())   # clear the environment
cat("\014")       # clear the console
library(rstan)

# compile flag
compile <- FALSE

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

# RStan configuration
options(mc.cores = parallel::detectCores())  # paralelizar chains
rstan_options(auto_write = TRUE)             # recompila só se o .stan mudar

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "poisson_pol2_sim1" # csv file with data
df <- read.table(paste("../data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
Tt <- length(y)
theta1_true <- df$theta1
theta2_true <- df$theta2
lambda_true <- exp(theta1_true)
W1_true <- 0.03     # state variance W1
W2_true <- 0.00003  # state variance W2

# Compile or load the model
if (compile) {
    model <- stan_model(file = "poisson_pol2_pcprior.stan")
    saveRDS(model, file = "poisson_pol2_pcprior.rds")
} else {
    model <- readRDS("poisson_pol2_pcprior.rds")
}

# sigma_W1 ####
alpha <- 0.05
U_W1 <- 0.5
U_W2 <- 0.02
lambda_W1 <- -log(alpha) / U_W1
lambda_W2 <- -log(alpha) / U_W2

sigmas_W1 <- seq(0, 3*U_W1, length.out = 500)
sigmas_W2 <- seq(0, 3*U_W2, length.out = 500)


printf("True percentil: %.2f", pexp(sqrt(W1_true), rate = lambda_W1))
plot(sigmas_W1, dexp(sigmas_W1, rate = lambda_W1), type = "l", lwd = 2,
     main = expression(sigma[W1]), xlab = expression(sigma), ylab = "densidade")
abline(v = sqrt(W1_true), col = "red",  lty = 2, lwd = 2)
abline(v = U_W1,          col = "blue", lty = 3, lwd = 2)
legend("topright", legend = c("verdadeiro", "U"),
       col = c("red", "blue"), lty = c(2, 3), bty = "n")


printf("True percentil: %.2f", pexp(sqrt(W2_true), rate = lambda_W2))
plot(sigmas_W2, dexp(sigmas_W2, rate = lambda_W2), type = "l", lwd = 2,
     main = expression(sigma[W2]), xlab = expression(sigma), ylab = "densidade")
abline(v = sqrt(W2_true), col = "red",  lty = 2, lwd = 2)
abline(v = U_W2,          col = "blue", lty = 3, lwd = 2)
legend("topright", legend = c("verdadeiro", "U"),
       col = c("red", "blue"), lty = c(2, 3), bty = "n")

#####
stan_data <- list(
    T         = Tt,
    y         = y,

    mu_01     = log(y[1]),  # theta_01 ~ N(mu_01, sigma2_01)
    sigma2_01 = 10,

    mu_02     = 0.1,        # theta_02 ~ N(mu_02, sigma2_02)
    sigma2_02 = 1,

    # PC priors
    alpha     = alpha,
    U_W1      = U_W1,
    U_W2      = U_W2
)

burnin <- 2000

fit <- sampling(
    object = model,
    data   = stan_data,
    chains = 1,
    iter   = 10000,
    warmup = burnin,
    thin   = 1,
    seed   = 42,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
)
warmup_time <- get_elapsed_time(fit)[1, 1]
sample_time <- get_elapsed_time(fit)[1, 2]
elapsed_time <- warmup_time + sample_time


# Results ####
samples <- extract(fit, inc_warmup = TRUE)
sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
ac_hist <- sampler_params[[1]][, "accept_stat__"]
ac_mean <- mean(ac_hist)

W1_samples   <- samples$W1
W2_samples   <- samples$W2

theta1_samples <- samples$theta1   # matrix: n_iter x T
theta1_mean  <- apply(theta1_samples[-(1:burnin), ], 2, mean)

theta2_samples <- samples$theta2   # matrix: n_iter x T
theta2_mean  <- apply(theta2_samples[-(1:burnin), ], 2, mean)

lambda_samples <- samples$lambda_hat
lambda_mean  <- apply(lambda_samples[-(1:burnin), ], 2, mean)

W1_mean  <- mean(W1_samples[-(1:burnin)])
W2_mean  <- mean(W2_samples[-(1:burnin)])

printf("Execution time: %.0f s", elapsed_time)
printf("Mean acception ratio: %.2f", ac_mean)
print(fit, pars = c("W1", "W2",
                    "theta1[10]", "theta1[50]", "theta1[100]", "theta1[150]",
                    "theta2[10]", "theta2[50]", "theta2[100]", "theta2[150]"),
      probs = c(0.25, 0.5, 0.75))
printf("W1_mean: %.5f", W1_mean)
printf("W2 mean: %.5f", W2_mean)

# Effective sample size / elapsed time ####
s <- summary(fit)$summary
ess_bulk <- s[, "n_eff"]
ess_per_sec <- ess_bulk / elapsed_time

ess_per_sec_w1 <- ess_per_sec[["W1"]]

printf("Effective Sample Size / second:")
printf("\tW1: %.2f", ess_per_sec[["W1"]])
printf("\tW2: %.2f", ess_per_sec[["W2"]])

for (t in c(50, 100, 200, 250)) {
    printf("\ttheta1[%d]: %.2f", t, ess_per_sec[[sprintf("theta1[%d]", t)]])
}

for (t in c(50, 100, 200, 250)) {
    printf("\ttheta2[%d]: %.2f", t, ess_per_sec[[sprintf("theta2[%d]", t)]])
}


# Plots #####
# y, lambda_true, lambda_mean
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson local level model")
points(x, y, pch = 20)
lines(x, lambda_mean, col="red", lwd=2)
lines(x, lambda_true, col="blue", lwd=2)
legend("topright",
       legend = expression(y[t], lambda[t], hat(lambda)[t]),
       col = c("black", "blue", "red"),
       lty = c(NA, 1, 1),
       lwd = c(NA, 2, 2),
       pch = c(20, NA, NA),
       bty = "n")


# theta1_true, theta1_mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_mean, type="l", col="blue", lwd="2")
lines(x, theta1_true)
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")

# theta2_true, theta2_mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta2_mean, type="l", col="blue", lwd="2")
lines(x, theta2_true)
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


# Posterior distribution of theta_t1 #####
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    hist(theta1_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[1][.(t)]), main = bquote("Posterior of " * theta[1][.(t)]))
    lines(density(theta1_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}

# Posterior distribution of theta_t2 #####
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    hist(theta2_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[2][.(t)]), main = bquote("Posterior of " * theta[2][.(t)]))
    lines(density(theta2_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of W1 #####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_samples[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = bquote(theta[1][.(t)]), main ="Posterior of W1")
lines(density(W1_samples[-(1:burnin)]), col = "blue", lwd = 2)

# Posterior distribution of W2 #####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_samples[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = bquote(theta[2][.(t)]), main ="Posterior of W2")
lines(density(W2_samples[-(1:burnin)]), col = "blue", lwd = 2)

# Traceplot for W1 and W2 #####
plot(W1_samples, type="l", xlab="n", ylab="W1", main="Traceplot of W1")
plot(W2_samples, type="l", xlab="n", ylab="W2", main="Traceplot of W2")

##### # Traceplots for theta_t1
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    plot(theta1_samples[, t], type="l", main=bquote(theta[1][.(t)]), xlab="", ylab="")
}

##### # Traceplots for theta_t2
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    plot(theta2_samples[, t], type="l", main=bquote(theta[2][.(t)]), xlab="", ylab="")
}


##### # Traceplot of acceptance ratio
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Acceptance ratio of " * theta))
