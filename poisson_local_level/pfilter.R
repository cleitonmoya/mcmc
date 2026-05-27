# Liu and West Filter

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- base::t       # alias to transpose function
options(error = function() traceback(2)) # more informative traceback
library(pomp)

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "poisson_local_level_sim1" # csv file with data
data <- readRDS(paste("../data/", source, ".rds", sep=""))

theta_true <- data$theta
W_true     <- data$W
y          <- data$y
T          <- length(y)

df <- data.frame(y = y, time = seq_len(T))


# Defina o modelo uma vez
modelo <- pomp(
    data        = df,
    times       = "time",
    t0          = 0,
    rprocess    = discrete_time(
        step.fun = Csnippet("
            theta = theta + rnorm(0, exp(log_W));
        "), delta.t = 1
    ),
    rmeasure = Csnippet("
        y = rpois(exp(theta));
    "),
    dmeasure = Csnippet("
        lik = dpois(y, exp(theta), give_log);
    "),
    rprior = Csnippet("
        log_W = rnorm(log_W_mean, log_W_sd);
    "),
    statenames = c("theta"),
    paramnames = c("log_W", "m0", "C0", "log_W_mean", "log_W_sd"),
    rinit      = Csnippet("
        theta = rnorm(m0, C0);
    ")
)

# Parâmetros: log_W fixado no valor verdadeiro, m0 e C0 do prior inicial
log_W_true <- log(W_true)
params_pf  <- c(log_W = log_W_true,
                m0 = 0,
                C0 = 1,
                log_W_mean = log_W_true,  # centro do prior em log_W
                log_W_sd   = 1)          # incerteza do prior

# Particle Filter with adaptive resampling
pf_fit <- pfilter(
    modelo,
    params      = params_pf,
    Np          = 2000,
    filter.mean = TRUE
)

lw_fit <- bsmc2(
    modelo,
    params = params_pf,
    Np     = 2000,
    smooth = 0.1
)

# Extrair média posterior de log_W estimada pelo L&W
log_W_post <- lw_fit@post["log_W", ]          # partículas da posteriori
log_W_lw   <- mean(log_W_post)                # média posterior
cat("log_W verdadeiro:", log_W_true, "\n")
cat("Log-verossimilhança (wpfilter):", logLik(pf_fit), "\n")
cat("log_W estimado (L&W):", log_W_lw, "\n")

# Rodar pfilter com log_W estimado pelo L&W para obter estados filtrados
params_lw <- params_pf
params_lw["log_W"] <- log_W_lw

pf_lw_est <- pfilter(
    modelo,
    params      = params_lw,
    Np          = 2000,
    filter.mean = TRUE
)

# Extrair estados filtrados
theta_filtrado <- as.numeric(filter_mean(pf_fit)["theta", ])
theta_lw       <- as.numeric(pf_lw_est@filter.mean["theta", ])

# Como o modelo está em escala log (exp(theta) = lambda),
# a média filtrada de theta é E[log(lambda_t) | y_{1:t}]
# Para obter E[lambda_t | y_{1:t}] use exp(theta_filtrado) como aproximação
lambda_filtrado <- exp(theta_filtrado)
lambda_true     <- exp(theta_true)   # se theta_true está em log; ajuste se não

# Plot
t_seq <- seq_len(T)

#####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(t_seq, theta_true, type = "l", col = "black", lwd = 1.5,
     xlab = "t", ylab = expression(theta[t]),
     main = "Estado latente: verdadeiro vs. filtrado")
lines(t_seq, theta_filtrado, col = "steelblue", lwd = 1.5, lty = 2)
lines(t_seq, theta_lw,       col = "firebrick",  lwd = 1.5, lty = 3)
legend("topright",
       legend = c("Verdadeiro",
                  paste0("PF (log_W=", round(log_W_true, 2), " verdadeiro)"),
                  paste0("L&W (log_W=", round(log_W_lw, 2), " estimado)")),
       col    = c("black", "steelblue", "firebrick"),
       lty    = c(1, 2, 3), lwd = 1.5, bty = "n")
