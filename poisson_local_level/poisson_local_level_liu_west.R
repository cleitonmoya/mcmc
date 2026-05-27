# Poisson Local Trend Model
# Liu and West Filter
# Model:
#   y_t ~ Poisson(exp(theta_t))
#   theta_t = theta_{t-1} + omega_t,  omega_t ~ N(0, W)

graphics.off()
rm(list = ls())
cat("\014")
set.seed(42)
tp <- base::t
options(error = function() traceback(2))

library(pomp)

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

source <- "poisson_sin_200"
data   <- readRDS(paste("../data/", source, ".rds", sep = ""))
theta_true <- data$theta
y           <- data$y
T           <- length(y)

df <- data.frame(y = y, time = seq_len(T))

# Model
modelo <- pomp(
    data  = df,
    times = "time",
    t0    = 0,

    rprocess = discrete_time(
        step.fun = Csnippet("
            theta += rnorm(0, exp(0.5 * log_W));
        "), delta.t = 1
    ),

    rmeasure = Csnippet("
        y = rpois(exp(theta));
    "),

    rprior = Csnippet("
        log_W = log(rgamma(nu_0, 1.0/eta_0));
    "),

    dmeasure = Csnippet("
        lik = dpois(y, exp(theta), give_log);
    "),

    statenames = c("theta"),

    paramnames = c("log_W", "m0", "C0", "nu_0", "eta_0"),

    rinit = Csnippet("
        theta = rnorm(m0, sqrt(C0));
    ")
)

# Parameters
params_base <- c(
    log_W = 0,
    m0 = 0,
    C0 = 1,
    nu_0 = 0.1,
    eta_0 = 0.1
)

# Filter
start_time = proc.time() # execution time
lw_fit <- bsmc2(modelo, params = params_base, Np = 5000, smooth = 0.1)
end_time <- proc.time()

# Execution time
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)


log_W <- mean(lw_fit@post["log_W", ])
printf("W estimado: %.5f", exp(log_W))

params_lw <- params_base
params_lw["log_W"] <- log_W
pf_lw  <- pfilter(modelo, params = params_lw, Np = 5000, filter.mean = TRUE)
theta_mean <- as.numeric(pf_lw@filter.mean["theta", ])

# Effective sample size
ess_lw <- lw_fit@eff.sample.size
printf("ESS/s: %.2f", ess_lw / elapsed_time)

# Plots ####
x <- seq_len(T)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(x, theta_mean, type="l", ylab="", col="blue", lwd=2)
lines(x, theta_true)
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")

#####
plot(lw_fit@eff.sample.size, type = "l", ylab = "ESS", xlab = "t")
