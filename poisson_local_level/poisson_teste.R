set.seed(42)

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_teste_300" # csv file to save the data

Tt <- 300

# Slope determinístico suave: AR(1) com delta forte, simulado sem acumulação
delta <- 0.85
W2    <- 0.004
W1    <- 0.001

# Gera slope estacionário (media zero)
theta2 <- numeric(Tt)
theta2[1] <- 0
for (t in 2:Tt)
	theta2[t] <- delta * theta2[t-1] + rnorm(1, 0, sqrt(W2))

# theta1 = integral do slope + ruído pequeno
# Para evitar drift: usa soma cumulativa re-centrada a cada bloco
theta1 <- 2.0 + cumsum(theta2) / sqrt(Tt) * 1.2 + cumsum(rnorm(Tt, 0, sqrt(W1)))

# Garante range controlado: normaliza amplitude
theta1 <- (theta1 - mean(theta1)) / sd(theta1) * 0.7 + 2.0

lambda <- exp(theta1)
y      <- rpois(Tt, lambda)

cat("range lambda:", round(range(lambda), 2), "\n")
cat("range y:", range(y), "\n")

par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(y, pch=20)
lines(lambda, type="l", col="red")

saveRDS(list(y = y, theta=theta1), paste(filename, ".rds", sep=""))