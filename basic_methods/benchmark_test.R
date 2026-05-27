library(Matrix)     # to deal with band sparse matrix
library(microbenchmark)

# Parâmetros
Tt <- 200

# Matriz H
H <- bandSparse(Tt, Tt, k = c(0, -1),
                diagonals = list(rep(1, Tt), rep(-1, Tt-1)))
B <- Diagonal(Tt) - H

# Valores numéricos
set.seed(42)
vartheta_1 <- rnorm(Tt, mean = log(3), sd = 0.1)
vartheta_2 <- rnorm(Tt, mean = 0,      sd = 0.05)
theta_01   <- log(3)
eta_0      <- 0.01

# mu_1 (necessário para forma matricial)
mu_1 <- theta_01 * rep(1, Tt) + B %*% vartheta_2

# Benchmark
resultado <- microbenchmark(
    matricial = {
        diff_1 <- vartheta_1 - mu_1
        eta_0 + 0.5 * as.numeric(crossprod(H %*% diff_1))
    },
    vetorial = {
        diffs <- vartheta_1 - c(theta_01, vartheta_1[-Tt] + vartheta_2[-Tt])
        eta_0 + 0.5 * sum(diffs^2)
    },
    times = 1000
)
print(resultado)
