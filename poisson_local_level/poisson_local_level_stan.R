# ============================================================
# Poisson Local Level Model via RStan
# ============================================================

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

source <- "sin_level" # csv file with data
df <- read.table(paste("data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
Tt <- length(y)
theta_true <- df$theta1
lambda_true <- exp(theta_true)

# Complile once and store the object
# model <- stan_model(file = "poisson_local_level.stan")
# saveRDS(model, file = "poisson_local_level_compilado.rds")

#####
model <- readRDS("poisson_local_level_compilado.rds")
stan_data <- list(
    T        = Tt,
    y        = y,
    mu_0     = y[1],       # theta_0 ~ N(mu_0, sigma2_0)
    sigma2_0 = 10,
    nu_0     = 0.01,       # phi_0 ~ Gamma(nu_0, eta_0)
    eta_0    = 0.01
)

fit <- sampling(
    object = model,
    data   = stan_data,
    chains = 1,
    iter   = 5000,
    warmup = 200,
    thin   = 1,
    seed   = 42
)
elapsed_time <- get_elapsed_time(fit)

# Results ####
burnin <- 200
samples <- extract(fit, inc_warmup = TRUE)
sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
ac_hist <- sampler_params[[1]][, "accept_stat__"]
ac_mean <- mean(ac_hist)

W_samples   <- samples$W
theta_samples <- samples$theta   # matrix: n_iter x T
theta_mean  <- apply(theta_samples[-(1:burnin)], 2, mean)
theta_lower <- apply(theta_samples[-(1:burnin)], 2, quantile, 0.025)
theta_upper <- apply(theta_samples[-(1:burnin)], 2, quantile, 0.975)
lambda_samples <- samples$lambda_hat
lambda_mean  <- apply(lambda_samples[-(1:burnin)], 2, mean)

printf("Execution time: %.0f s", elapsed_time)
printf("Mean acception ratio of theta: %.2f", ac_mean)
print(fit, pars = c("W", "theta[10]", "theta[50]", "theta[100]", "theta[150]"),
      probs = c(0.25, 0.5, 0.75))

# Plots #####
# y, lambda_true, lambda_estimated
lambda_true = exp(theta_true)
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


# Posterior distribution of theta_t #####
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    hist(theta_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t)]), main = bquote("Posterior of " * theta[.(t)]))
    lines(density(theta_samples[, t]), col = "blue", lwd = 2)
}


# Posterior distribution of W #####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W_samples[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = bquote(theta[.(t)]), main ="Posterior of W")
lines(density(W_samples[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot for W #####
plot(W_samples, type="l", xlab="n", ylab="W", main="Traceplot of W")


##### # Traceplots for theta_t
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    plot(theta_samples[-(1:burnin), t], type="l", main=bquote(theta[.(t)]), xlab="", ylab="")
}


##### # Traceplot of acceptance ratio of theta
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Acceptance ratio of " * theta))
