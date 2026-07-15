# Poisson - 2nd Order Polynomial Dynamic Model
# MCMC: NUTS (Stan)
# Author: Cleiton Moya de Almeida

graphics.off()    # close the plots
rm(list = ls())   # clear the environment
cat("\014")       # clear the console
library(rstan)

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
source <- "normal_pol2_sim1" # csv file with data
df <- read.table(paste("../data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
Tt <- length(y)
theta1_true <- df$theta1
theta2_true <- df$theta2

#####
# Complile once and store the object
# model <- stan_model(file = "normal_pol2.stan")
# saveRDS(model, file = "normal_pol2.rds")

#####
model <- readRDS("normal_pol2.rds")

stan_data <- list(
    y         = y,
    T         = Tt,

    mu_01     = 0.01,       # theta_01 ~ N(mu_01, sigma2_01)
    sigma2_01 = 10,

    mu_02     = 0.01,       # theta_02 ~ N(mu_02, sigma2_02)
    sigma2_02 = 10,

    nu_V      = 0.01,       # phi_V    ~ Gamma(nu_V, eta_V)
    eta_V     = 0.01,

    nu_01     = 0.01,       # phi_01   ~ Gamma(nu_01, eta_01)
    eta_01    = 0.01,

    nu_02     = 0.01,       # phi_02   ~ Gamma(nu_02, eta_02)
    eta_02    = 0.01
)

burnin <- 1000

fit <- sampling(
    object = model,
    data   = stan_data,
    chains = 1,
    iter   = 5000,
    warmup = burnin,
    thin   = 1,
    seed   = 42
)
warmup_time <- get_elapsed_time(fit)[1, 1]
sample_time <- get_elapsed_time(fit)[1, 2]
elapsed_time <- warmup_time + sample_time


# Results ####
samples <- extract(fit, inc_warmup = TRUE)
sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
ac_hist <- sampler_params[[1]][, "accept_stat__"]
ac_mean <- mean(ac_hist)

V_samples  <- samples$V
W1_samples <- samples$W1
W2_samples <- samples$W2
theta1_samples <- samples$theta1   # matrix: n_iter x T
theta2_samples <- samples$theta2

theta1_mean  <- apply(theta1_samples[-(1:burnin), ], 2, mean)
theta2_mean  <- apply(theta2_samples[-(1:burnin), ], 2, mean)
V_mean   <- mean(V_samples[-(1:burnin)])
W1_mean  <- mean(W1_samples[-(1:burnin)])
W2_mean  <- mean(W2_samples[-(1:burnin)])

printf("Execution time: %.0f s", elapsed_time)
printf("Mean acception ratio: %.2f", ac_mean)
print(fit, pars = c("V", "W1", "W2",
                    "theta1[50]", "theta1[100]", "theta1[150]", "theta1[250]",
                    "theta2[50]", "theta2[100]", "theta2[150]", "theta2[250]"),
      probs = c(0.25, 0.5, 0.75))
printf("V mean: %.5f", V_mean)
printf("W1 mean: %.5f", W1_mean)
printf("W2 mean: %.5f", W2_mean)


# Effective sample size / elapsed time ####
s <- summary(fit)$summary
ess_bulk <- s[, "n_eff"]
ess_per_sec <- ess_bulk / elapsed_time

printf("Effective Sample Size / second:")
printf("\tV: %.2f", ess_per_sec[["V"]])
printf("\tW1: %.2f", ess_per_sec[["W1"]])
printf("\tW2: %.2f", ess_per_sec[["W2"]])

observed_times <- c(50, 100, 150, 250)
for (t in observed_times) {
    printf("\ttheta1[%d]: %.2f", t, ess_per_sec[[sprintf("theta1[%d]", t)]])
}

for (t in observed_times) {
    printf("\ttheta2[%d]: %.2f", t, ess_per_sec[[sprintf("theta2[%d]", t)]])
}

#####

# theta1_true, theta1_mean
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_mean, type="l", col="blue", lwd="2")
lines(x, theta1_true)
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")

# theta2_true, theta2_mean
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta2_mean, type="l", col="blue", lwd="2")
lines(x, theta2_true)
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")


# Posterior distribution of theta_t1
par(mfrow = c(2, 2))
for (t in observed_times) {
    hist(theta1_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[1][.(t)]),
         main = bquote("Posterior of " * theta[.(t) * "," * 1]))
    lines(density(theta1_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}

# Posterior distribution of theta_t2
par(mfrow = c(2, 2))
for (t in observed_times) {
    hist(theta2_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[2][.(t)]),
         main = bquote("Posterior of " * theta[.(t) * "," * 2]))
    lines(density(theta2_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}

# Posterior distribution of V
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(V_samples[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = "V", main ="Posterior of V")
lines(density(V_samples[-(1:burnin)]), col = "blue", lwd = 2)

# Posterior distribution of W1
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_samples[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = "W1", main ="Posterior of W1")
lines(density(W1_samples[-(1:burnin)]), col = "blue", lwd = 2)

# Posterior distribution of W2
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_samples[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = "W2", main ="Posterior of W2")
lines(density(W2_samples[-(1:burnin)]), col = "blue", lwd = 2)

# Traceplot for V, W1 and W2
plot(V_samples, type="l", xlab="n", ylab="V", main="Traceplot of V")
plot(W1_samples, type="l", xlab="n", ylab="W1", main="Traceplot of W1")
plot(W2_samples, type="l", xlab="n", ylab="W2", main="Traceplot of W2")

# Traceplots for theta_t1
par(mfrow = c(2, 2))
for (t in observed_times) {
    plot(theta1_samples[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
}

# Traceplots for theta_t2
par(mfrow = c(2, 2))
for (t in observed_times) {
    plot(theta2_samples[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
}

# Traceplot of acceptance ratio
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_hist, type="l", xlab="n", ylab="ratio", main="Acceptance ratio")
