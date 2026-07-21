# Liu and West Filter — Poisson Polinomial de Segunda Ordem
# Modelo:
#   y_t ~ Poisson(exp(theta1_t))
#   theta1_t = theta1_{t-1} + theta2_{t-1} + omega1_t,  omega1_t ~ N(0, W1)
#   theta2_t = theta2_{t-1}               + omega2_t,  omega2_t ~ N(0, W2)

graphics.off()
rm(list = ls())
cat("\014")
set.seed(42)
tp <- base::t
options(error = function() traceback(2))
library(pomp)

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

source <- "poisson_pol2_200"
data   <- readRDS(paste("../data/", source, ".rds", sep = ""))

theta1_true <- data$theta   # apenas theta1 é observável como referência
y           <- data$y
T           <- length(y)

df <- data.frame(y = y, time = seq_len(T))

# -----------------------------------------------------------------------
# Definição do modelo
# -----------------------------------------------------------------------
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
    dmeasure = Csnippet("
        lik = dpois(y, exp(theta1), give_log);
    "),

    rprior = Csnippet("
        log_W1 = rnorm(log_W1_mean, log_W1_sd);
        log_W2 = rnorm(log_W2_mean, log_W2_sd);
    "),

    statenames = c("theta1", "theta2"),
    paramnames = c("log_W1", "log_W2",
                   "m01", "m02", "C01", "C02",
                   "log_W1_mean", "log_W1_sd",
                   "log_W2_mean", "log_W2_sd"),

    rinit = Csnippet("
        theta1 = rnorm(m01, C01);
        theta2 = rnorm(m02, C02);
    ")
)

# -----------------------------------------------------------------------
# Parâmetros
# -----------------------------------------------------------------------
params_base <- c(
    log_W1      = 0,
    log_W2      = 0,
    m01         = 0,   m02 = 0,
    C01         = 1,   C02 = 1,
    log_W1_mean = 0,  log_W1_sd = 1,
    log_W2_mean = 0,  log_W2_sd = 1
)

# -----------------------------------------------------------------------
# 1. Liu & West (estima log_W1 e log_W2)
# -----------------------------------------------------------------------
lw_fit    <- bsmc2(modelo, params = params_base, Np = 2000, smooth = 0.1)
log_W1_lw <- mean(lw_fit@post["log_W1", ])
log_W2_lw <- mean(lw_fit@post["log_W2", ])
cat("log_W1 estimado (L&W) :", log_W1_lw, "\n")
cat("log_W2 estimado (L&W) :", log_W2_lw, "\n")

params_lw           <- params_base
params_lw["log_W1"] <- log_W1_lw
params_lw["log_W2"] <- log_W2_lw
pf_lw  <- pfilter(modelo, params = params_lw, Np = 2000, filter.mean = TRUE)
theta1_lw <- as.numeric(pf_lw@filter.mean["theta1", ])

# -----------------------------------------------------------------------
# 2. Particle MCMC (PMMH)
# -----------------------------------------------------------------------
pmcmc_fit <- pmcmc(
    modelo,
    params   = params_base,
    Np       = 500,
    Nmcmc    = 2000,
    proposal = mvn_diag_rw(c(log_W1 = 0.3, log_W2 = 0.3))
)

chain        <- traces(pmcmc_fit)
burnin       <- 500
log_W1_pmcmc <- mean(chain[-(1:burnin), "log_W1"])
log_W2_pmcmc <- mean(chain[-(1:burnin), "log_W2"])
cat("log_W1 estimado (PMCMC):", log_W1_pmcmc, "\n")
cat("log_W2 estimado (PMCMC):", log_W2_pmcmc, "\n")

params_pmcmc           <- params_base
params_pmcmc["log_W1"] <- log_W1_pmcmc
params_pmcmc["log_W2"] <- log_W2_pmcmc
pf_pmcmc <- pfilter(modelo, params = params_pmcmc, Np = 2000, filter.mean = TRUE)
theta1_pmcmc <- as.numeric(pf_pmcmc@filter.mean["theta1", ])

# -----------------------------------------------------------------------
# Plot: apenas theta1
# -----------------------------------------------------------------------
t_seq <- seq_len(T)

par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(t_seq, theta1_true, type = "l", col = "black", lwd = 1.5,
     xlab = "t", ylab = expression(theta[t1]),
     main = expression(paste(theta[t1], ": verdadeiro vs. filtrado")))
lines(t_seq, theta1_lw,    col = "firebrick", lwd = 1.5, lty = 2)
lines(t_seq, theta1_pmcmc, col = "darkgreen", lwd = 1.5, lty = 3)
legend("topright",
       legend = c("Verdadeiro",
                  paste0("L&W (log_W1=",   round(log_W1_lw,    2), ")"),
                  paste0("PMCMC (log_W1=", round(log_W1_pmcmc, 2), ")")),
       col = c("black", "firebrick", "darkgreen"),
       lty = c(1
