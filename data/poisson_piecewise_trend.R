set.seed(7)

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_piecewisetrend_300" # csv file to save the data

Tt <- 300

# Slope muda de sinal em breakpoints — tendência local por partes
breakpoints <- c(0, 75, 150, 225, 300)
slopes      <- c(0.025, -0.02, 0.018, -0.015)

theta_true <- numeric(Tt)
theta_true[1] <- 2.0
for (t in 2:Tt) {
	seg <- findInterval(t-1, breakpoints, rightmost.closed = TRUE)
	theta_true[t] <- theta_true[t-1] + slopes[seg]
}

# Adiciona ruído de nível pequeno
theta1 <- theta_true + cumsum(rnorm(Tt, 0, 0.04))

# Centraliza
theta1 <- theta1 - mean(theta1) + 2.5

lambda <- exp(theta1)
y <- rpois(Tt, lambda)

plot(y, pch=20)
lines(exp(theta1), type="l", col="red")

# Save the data
saveRDS(list(y = y, theta=theta1), paste(filename, ".rds", sep=""))