graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)

filename1 <- "poisson_pol2_200_teste"
filename2 <- "poisson_pol2_2000_teste"

gera_theta <- function(T) {

	# Breakpoints e valores definidos na escala [0, 1] (proporções)
	props  <- c(0, 0.10, 0.23, 0.25, 0.42, 0.45, 0.65, 0.70, 0.82, 0.86, 1.00)
	values <- c(2.2, 3.2, 2.7, 0.5, 3.0, 1.5, 2.5, 3.3, 3.2, 2.2, 2.1)
	bulges <- c(0.4, -0.3, -0.6, 1.0, -0.6, 0.5, 0.5, -0.2, -0.4, 0.1)

	breakpoints <- round(props * T)

	theta_true <- numeric(T)
	theta_true[1] <- values[1]

	for (i in seq_along(bulges)) {
		t0 <- breakpoints[i]
		t1 <- breakpoints[i + 1]
		v0 <- values[i]
		v1 <- values[i + 1]
		b  <- bulges[i]

		tt <- (t0 + 1):t1
		s  <- (tt - t0) / (t1 - t0)

		theta_true[tt] <- v0 + (v1 - v0) * s + b * s * (1 - s)
	}

	return(theta_true)
}


theta200 <- gera_theta(200)
theta2000 <- gera_theta(2000)


y200 <- rpois(200, exp(theta200))
y2000 <- rpois(2000, exp(theta2000))
x200 <- 1:200
x2000 <- 1:2000

# Plot ####
plot(x200, y200, type="l", xlab="t", ylab="", col="gray",
	 main="Poisson Local Level Polynomial Model")
points(x200, y200, pch = 20)
lines(x200, exp(theta200), col="red", lwd=2)

plot(x2000, y2000, type="l", xlab="t", ylab="", col="gray",
	 main="Poisson Local Level Polynomial Model")
points(x2000, y2000, pch = 20)
lines(x2000, exp(theta2000), col="red", lwd=2)

# Save the data ####
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
saveRDS(list(y = y200, theta=theta200),
		paste("../data/", filename1, ".rds", sep=""))

saveRDS(list(y = y2000, theta=theta2000),
		paste("../data/", filename2, ".rds", sep=""))
