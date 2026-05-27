set.seed(42)

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_teste_200" # csv file to save the data

Tt <- 200       # série mais curta — menos tempo para deriva acumular
W1 <- 0.003     # ruído de nível
W2 <- 5e-5      # slope quase constante — muda muitíssimo devagar

theta1 <- numeric(Tt)
theta2 <- numeric(Tt)
theta1[1] <- 2.0
theta2[1] <- 0.015   # slope inicial pequeno mas positivo

for (t in 2:Tt) {
	theta2[t] <- theta2[t-1] + rnorm(1, 0, sqrt(W2))
	theta1[t] <- theta1[t-1] + theta2[t-1] + rnorm(1, 0, sqrt(W1))
}

# Verifica deriva de theta1
cat("range theta1:", round(range(theta1), 2), "\n")

# Se deriva for grande, re-centra apenas para visualização — não para ajuste
lambda <- exp(theta1)
y      <- rpois(Tt, lambda)
cat("range lambda:", round(range(lambda), 1), "\n")

par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(y, pch=20)
lines(lambda, type="l", col="red")

saveRDS(list(y = y, theta=theta1), paste(filename, ".rds", sep=""))