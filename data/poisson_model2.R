rm(list = ls())     # clear the environment
set.seed(42)

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_model2" # csv file to save the data

Tt <- 200
theta <- numeric(Tt)
y <- numeric(Tt)

gamma <- 0.85
alpha <- 10
beta <- 5
theta_t <- 10

for (t in 1:Tt) {
    e <- rbeta(1, gamma*alpha, (1-gamma)*alpha)
    theta_t <- theta_t*e/gamma
    theta[t] <- theta_t
    y_t <- rpois(1, theta_t)
    y[t] <- y_t
    alpha <- gamma*alpha + y[t]
    beta <- gamma*beta +1

}

# Plot
plot(y, pch=20)
lines(theta, type="l", col="red")

# Save the data
saveRDS(list(y = y, theta=theta), paste(filename, ".rds", sep=""))
