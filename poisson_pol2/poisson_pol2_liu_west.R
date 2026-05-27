# Liu and West Filter — Poisson Polinomial de Segunda Ordem
# Modelo:
#   y_t ~ Poisson(exp(theta1_t))
#   theta1_t = theta1_{t-1} + theta2_{t-1} + omega1_t,  omega1_t ~ N(0, W1)
#   theta2_t = theta2_{t-1}                + omega2_t,  omega2_t ~ N(0, W2)

#graphics.off()
rm(list = ls())
#cat("\014")
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

source <- "poisson_sincos_2000"
data   <- readRDS(paste("../data/", source, ".rds", sep = ""))
theta1_true <- data$theta   # apenas theta1 é observável como referência
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
            double w1 = rnorm(0, exp(0.5 * log_W1));
            double w2 = rnorm(0, exp(0.5 * log_W2));
            theta1 += theta2 + w1;
            theta2 += w2;
        "), delta.t = 1
    ),

    rmeasure = Csnippet("
        y = rpois(exp(theta1));
    "),

    rprior = Csnippet("
        log_W1 = log(rgamma(nu_01, eta_01));
        log_W2 = log(rgamma(nu_02, eta_02));
    "),

    dmeasure = Csnippet("
        lik = dpois(y, exp(theta1), give_log);
    "),


    statenames = c("theta1", "theta2"),

    paramnames = c("log_W1", "log_W2",
                   "m01", "m02",
                   "C01", "C02",
                   "nu_01", "nu_02",
                   "eta_01", "eta_02"
                   ),

    rinit = Csnippet("
        theta1 = rnorm(m01, C01);
        theta2 = rnorm(m02, C02);
    ")
)

# Parameters
params_base <- c(
    log_W1 = 0,
    log_W2 = 0,
    m01 = 0,
    C01 = 1,
    m02 = 0,
    C02 = 1,
    nu_01 = 0.1,
    eta_01 = 0.1,
    nu_02 = 0.1,
    eta_02 = 0.1
)

# Filter
start_time = proc.time() # execution time
lw_fit <- bsmc2(modelo, params = params_base, Np = 5000, smooth = 0.1)
end_time <- proc.time()

# Execution time
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.0f s", elapsed_time)


log_W1 <- mean(lw_fit@post["log_W1", ])
log_W2 <- mean(lw_fit@post["log_W2", ])
cat("log_W1 estimado:", log_W1, "\n")
cat("log_W2 estimado:", log_W2, "\n")

params_lw <- params_base
params_lw["log_W1"] <- log_W1
params_lw["log_W2"] <- log_W2
pf_lw  <- pfilter(modelo, params = params_lw, Np = 5000, filter.mean = TRUE)
theta1_mean <- as.numeric(pf_lw@filter.mean["theta1", ])
theta2_mean <- as.numeric(pf_lw@filter.mean["theta2", ])

# Effective sample size
ess_lw <- tail(lw_fit@eff.sample.size, 1)
printf("ESS/s: %.2f", ess_lw / elapsed_time)


# Plots ####
x <- seq_len(T)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(x, theta1_mean, type="l", ylab="", col="blue", lwd=1)
lines(x, theta1_true)
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 1),
       bty = "n")

#####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(theta2_mean[-(1:100)], type="l", ylab="", col="blue", lwd=2)
legend("topright",
       legend = expression(theta[t2], hat(theta)[t2]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")

#####
plot(lw_fit@eff.sample.size, type = "l", ylab = "ESS", xlab = "t")
