# ==============================================================================
# PARTICLE GIBBS WITH ANCESTRAL SAMPLING (PGAS)
# Local Level Model
#
# ==============================================================================

library(coda)

#graphics.off()      # close the plots
rm(list = ls())     # clear the environment
#cat("\014")        # clear the console
set.seed(42)
options(error = function() traceback(2)) # more informative traceback

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))

# Load the data
# source <- "poisson_pol2_2000" # rds file with data
# data <- readRDS(paste("../data/", source, ".rds", sep=""))
# y <- data$y
# theta_true <- data$theta
# Tt <- length(y)


# Load the data
functions_grid <- c("constant","linear", "quadratic", "sinusoidal")

f <- 3
Tt <- 1600
replica <- 1

source <- sprintf("%s_%s_%s", functions_grid[f], Tt, replica)
data <- readRDS(paste("../../cobalebeb2027/data/simulated/", source, ".rds", sep=""))
y <- data$y
theta_true <- data$theta
Tt <- length(y)

#
# Set the seed according to the task_grid
#
Tt_grid <- c(200, 400, 800, 1600)
method <- 2    # grid: (montoril, pg_apf, sir_laplace, sir_collapsed, stan)

tau <- match(Tt, Tt_grid)
seed = method*1e5 + f*1e4 + tau*1e3 + replica

set.seed(seed)

printf("Executing %s, seed=%d" , source, seed)


# Log-verossimilhança da observação Poisson: p(y_t | x_t)
log_obs_density <- function(y_t, x_t) {
	lambda <- exp(x_t)
	dpois(y_t, lambda = lambda, log = TRUE)
}

# Log-densidade de transição de estado: p(x_t | x_{t-1}, sigma_v)
log_trans_density <- function(x_t, x_prev, sigma_v) {
	dnorm(x_t, mean = x_prev, sd = sigma_v, log = TRUE)
}

# 4. AMOSTRAGEM DE SIGMA_V VIA GIBBS (INVERSA-GAMA) ---------------------------
sample_W <- function(theta, nu_0, eta_0) {
	Tt <- length(theta)

	# Diferenças do passeio aleatório: dx_t = x_t - x_{t-1}
	diffs <- diff(theta)

	# Atualização dos hiperparâmetros conjugados para Inversa-Gama
	a_post <- nu_0 + (Tt - 1) / 2
	b_post <- eta_0 + sum(diffs^2) / 2

	# Amostragem da variância sigma_v^2 ~ Inv-Gamma(a_post, b_post)
	W <- 1 / rgamma(1, shape = a_post, rate = b_post)
	return(W)
}



# 6. EXECUÇÃO -----------------------------------------------------------------
N <- 1000
K <- 30  # Com Ancestral Sampling, 30 partículas costumam ser suficientes!
burnin <- 200

# Hyperparameters
nu_0 <- 2
eta_0 <- 0.1

mu_0 <- 0
sigma2_0 <- 100

Tt <- length(y)
if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 400) t_obs <- c(75, 100, 200, 300)
if (Tt == 800) t_obs <- c(200, 300, 500, 700)
if (Tt == 1600) t_obs <- c(400, 800, 1200, 1600)
if (Tt == 2000) t_obs <- c(500, 1000, 1500, 1750)

# Initial values
W <- 0.01
theta <- numeric(Tt)
theta_0 <- 0



# Matrizes/Vetores para guardar amostras da a posteriori
theta_hist <- matrix(0, nrow = N, ncol = Tt)
W_hist <- numeric(N)
theta_0_hist <- numeric(N)

for (n in 1:N) {

    # STEP 1: Update theta_{1:T} |  W
    particles <- matrix(0, nrow = K, ncol = Tt)
    ancestors <- matrix(0, nrow = K, ncol = Tt)
    weights   <- matrix(0, nrow = K, ncol = Tt)

	# --- Tempo t = 1 ---
	# Propõe K-1 partículas a partir da distribuição a priori do estado inicial
	particles[1:(K-1), 1] <- theta[1] + rnorm(K - 1, 0, sqrt(W))
	particles[K, 1]       <- theta[1]  # Partícula de referência

	# Ponderação no tempo t = 1
	log_w <- log_obs_density(y[1], particles[, 1])
	max_w <- max(log_w)
	weights[, 1] <- exp(log_w - max_w) / sum(exp(log_w - max_w))


	for (t in 2:Tt) {

		# a) Reamostragem para as K-1 partículas
		a_indices <- sample(1:K, size = K - 1, replace = TRUE, prob = weights[, t-1])
		ancestors[1:(K-1), t] <- a_indices

		# Propagação via Passeio Aleatório: x_t = x_{t-1} + v_t
		particles[1:(K-1), t] <- particles[a_indices, t-1] + rnorm(K - 1, 0, sqrt(W))

		# b) ANCESTRAL SAMPLING para a partícula K (Referência)
		# Avalia probabilidade de x_ref[t] ter vindo de cada partícula em t-1
		log_p_trans <- log_trans_density(theta[t], particles[, t-1], sqrt(W))
		log_w_ancestral <- log(weights[, t-1] + 1e-300) + log_p_trans

		max_wa <- max(log_w_ancestral)
		w_ancestral <- exp(log_w_ancestral - max_wa) / sum(exp(log_w_ancestral - max_wa))

		ancestors[K, t] <- sample(1:K, size = 1, prob = w_ancestral)
		particles[K, t] <- theta[t]

		# c) Ponderação pela verossimilhança Poisson
		log_w <- log_obs_density(y[t], particles[, t])
		max_w <- max(log_w)
		weights[, t] <- exp(log_w - max_w) / sum(exp(log_w - max_w))
	}

	# --- RECONSTRUÇÃO DA TRAJETÓRIA (BACKWARD TRACKING) ---
	k <- sample(1:K, size = 1, prob = weights[, Tt])

	theta_new <- numeric(Tt)
	theta_new[Tt] <- particles[k, Tt]

	for (t in Tt:2) {
		k <- ancestors[k, t]
		theta_new[t-1] <- particles[k, t-1]
	}

	theta <- theta_new

	# Sample W
	W <- sample_W(theta, nu_0, eta_0)

	# Sample theta_0
	# 1. Sample theta_0
	sigma2_0_bar <- (1/sigma2_0 +1/W)^(-1)
	mu_0_bar <- sigma2_0_bar*(mu_0/sigma2_0 + theta[1]/W)
	theta_0 <- rnorm(1, mean=mu_0_bar, sd=sqrt(sigma2_0_bar))

	# Store samples
	theta_hist[n, ] <- theta
	W_hist[n] <- W
	theta_0_hist[n] <- theta_0

	if (n %% 100 == 0) {
		printf("Iteration %d/%d | W_hat = %.4f", n, N, mean(W_hist[(burnin+1):n]))
	}
}


# Results
theta_mean <- colMeans(theta_hist[(burnin + 1):N, ])
lambda <- exp(theta_mean)


# Effective sample size ####
ess_theta_0 <- effectiveSize(mcmc(theta_0_hist[-(1:burnin)]))
ess_W <- effectiveSize(mcmc(W_hist[-(1:burnin)]))
ess_theta <- effectiveSize(mcmc(theta_hist[-(1:burnin),]))

printf("Effective Sample Size:")
printf("\ttheta_0: %.0f", ess_theta_0)
printf("\tW: %.0f", ess_W)
printf("\ttheta (mean): %.2f", mean(ess_theta))


#####
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", col="gray", xlab="t", ylab=expression(y[t]))
points(x, y, pch = 20)
lines(x, exp(theta_true), col="blue", lwd=1)
lines(x, lambda, col="red", lwd=1)

# Traceplot for theta0 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(theta_0_hist, type="l", xlab="n", ylab="W", main="Traceplot of theta_0")

# Traceplot for W ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W_hist[-(1:100)], type="l", xlab="n", ylab="W", main="Traceplot of W")

# Traceplot for theta #####
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), cex = 0.8)
for (t in t_obs) {
	plot(theta_hist[, t], type="l", main=bquote(theta[.(t)]), xlab="n", ylab="")
}



