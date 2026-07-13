# Approximating Posterior Distributions by Mixtures
# Reference: West, Mike. “Approximating Posterior Distributions by Mixtures”.
#              Journal of the Royal Statistical Society Series B: Statistical
#              Methodology 55, n. 2 (1993): 409–22.
#              https://doi.org/10.1111/j.2517-6161.1993.tb01911.x.
#
# Example 2.4 - Bimodal posterior
#
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
tp <- Matrix::t     # matrix transpose alias
options(error = function() traceback(2)) # more informative traceback

library(mvtnorm)    # multivariate normal and t distributions

set.seed(42)

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))

# Log-sum-exp auxiliary function
logsumexp <- function(x){
    c <- max(x)
    y <- c + log(sum((exp(x-c))))
    return(y)
}

# Posterior density (exact)
dpost <- function(pi1, pi2, n, r, s) {
	y <- (pi1*(1-pi2))^r *
		 (pi2*(1-pi1))^s *
		 (1 - pi1*(1-pi2) - pi2*(1-pi1))^(n-r-s)
	return(y)
}

log_dpost <- function(pi1, pi2, n, r, s) {
	y <- r*log(pi1*(1-pi2)) +
		 s*log(pi2*(1-pi1)) +
		 (n-r-s)*log(1 - pi1*(1-pi2) - pi2*(1-pi1))
	return(y)
}

# logit function for theta
theta_i <- function(pi_i) {
	log(pi_i/(1-pi_i))
}

# logistic function for pi (inverse-logit)
pi_i <- function(theta_i) {
	exp(theta_i) / (1 + exp(theta_i))
}

# Jacobian determinant of the transformation theta_i = logit(pi_i)
jacobian <- function(pi1, pi2) {
	1 / (pi1*(1-pi1) * pi2*(1-pi2))
}

# Initial ISF g0 (Multivariate T)
g0 <- Vectorize(
	function(pi1, pi2, delta, sigma, df) {
		theta1 <- theta_i(pi1)
		theta2 <- theta_i(pi2)
		f_theta <- dmvt(c(theta1, theta2), delta=delta, sigma=sigma, df=df, log = FALSE)
		f_pi <- f_theta*jacobian(pi1, pi2)
		return(f_pi)
	},
	vectorize.args = c("pi1", "pi2")
)


# Figure 2.4a - True posterior
nn <- 45
r <- 5
s <- 3
p <- 2

pi1 <- seq(0.001, 0.999, length.out=200)
pi2 <- seq(0.001, 0.999, length.out=200)

Z1 <- outer(pi1, pi2, dpost, n = nn, r = r, s = s)
levels1 <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9) * max(Z1)
contour(pi1, pi2, Z1, levels = levels1, drawlabels = FALSE, main="Figure 2.4a")


# Figure 2.4b - Initial ISF g0
nu <- 9
mode <- c(0, 0)
Sigma <- diag(c(1.8, 1.8))

Z2 <- outer(pi1, pi2, g0, delta = mode, sigma = Sigma, df = nu)
levels2 <- c(0.01, 0.1, 0.25, 0.5, 0.6, 0.75, 0.9) * max(Z2)
contour(pi1, pi2, Z2, levels = levels2, drawlabels = FALSE, main = "Figure 2.4b")



# Initial mixture
# 1) Draw from g0
N <- 1000
Theta0 <- rmvt(N, sigma=Sigma, df=nu, delta=mode) # dim: N x p (p=2)

# Replicate the points due to the simetry
Theta0_ref <- cbind(-Theta0[,2], -Theta0[,1])
Theta0 <- rbind(Theta0, Theta0_ref)
N <- 2*N

# 2) Weights
log_g0 <- dmvt(Theta0, delta = mode, sigma = Sigma, df = nu, log = TRUE)
pii0 <- pi_i(Theta0)
log_f0 <- log_dpost(pii0[,1], pii0[, 2], nn, r ,s)
log_w0 <- log_f0 - log_g0
w0 <- matrix(exp(log_w0 - logsumexp(log_w0)), ncol=1) # dim: N x 1


# 3) Variance matrix (MC estimation)
n <- 1
c <- (4/(1+2*p))^(1/(1+4*p))
h <- c/(n^(1/(1+4*p)))
x <- sqrt(1-h^2)

theta_bar0 <- t(Theta0) %*% w0 # dim:
M0 <- sweep(x*Theta0, 2, (1-x)*theta_bar0, "+") # Shrinkage matrix

D0 <- sweep(Theta0, 2, theta_bar0, "-") # D[j,] = Theta[j,] - theta_bar
D_w0 <- sqrt(drop(w0)) * D0
V0 <- crossprod(D_w0)

# First refinement

# a) Sample from current mixture
j <- sample(1:N, N, replace=TRUE, prob=w0)
M_j0 <- M0[j,]

Sigma1 <- V0 * h^2
Theta1 <- t(sapply(1:N, function(i) {
    rmvt(1, sigma =  Sigma1, df = nu, delta = M_j0[i, ])
}))

# b) Update the weights
log_g1_theta_i  <- function(theta_i, w, M_j, Sigma1) {
    den_j <- function(j) {
        log(w[j]) + dmvt(theta_i, delta = M_j[j, ], sigma = Sigma1, df = nu, log = TRUE)
    }
    vec <- sapply(1:N, den_j)
    y <- logsumexp(vec)
    return(y)
}

pii1 <- pi_i(Theta1)
log_f1 <- log_dpost(pii1[,1], pii1[, 2], nn, r ,s)
log_g1 <- apply(Theta1, 1, log_g1_theta_i, w = w0, M_j = M_j0, Sigma1 = Sigma1)
log_w1 <- log_f1 - log_g1
w1 <- matrix(exp(log_w1 - logsumexp(log_w1)), ncol=1) # dim: N x 1


# c) Collapsed mixture
k <- 500
r <- nrow(M_j0)

while (r > k) {
    # (b) ordenar por peso crescente
    ord <- order(w1)
    M <- M[ord, ]
    w <- w[ord]
    # (c) vizinho mais próximo do componente 1
    dists <- rowSums(sweep(M[-1, ], 2, M[1, ])^2)
    i <- which.min(dists) + 1  # +1 pois removemos componente 1
    # (d) fundir
    w_new <- w[1] + w[i]
    m_new <- (w[1]*M[1,] + w[i]*M[i,]) / w_new
    # remover 1 e i, inserir fusionado
    M <- rbind(M[-c(1,i), ], m_new)
    w <- c(w[-c(1,i)], w_new)
    r <- r - 1
}

# Figure 2.4c
cex_scale <- w1 / max(w1) * 3
plot(pii1, cex = cex_scale, pch = 1, xlim = c(0,1), ylim = c(0,1),
     main = "Figure 2.4c")
contour(pi1, pi2, Z1, levels = 0.01 * max(Z1),
        drawlabels = FALSE, add = TRUE)

# Figure 2.4d

