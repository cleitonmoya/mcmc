# Poisson 2nd Order Polynomial Model simulation
# Author: Cleiton Moya de Almeida

#####
graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(6)

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_pol2_100_teste"

theta1 <- 0.01  # theta_01
theta2 <- 0.05 # theta_02
W1 <- 0.03     # state variance W1
W2 <- 0.03  # state variance W2
Tt <- 100      # number of observations

y_sim <- numeric(Tt)
theta1_sim <- numeric(Tt)
theta2_sim <- numeric(Tt)

for (t in 1:Tt) {
    theta1 <- theta1 + theta2 + rnorm(1, mean=0, sd=sqrt(W1))
    theta2 <- 0.75*theta2 + rnorm(1, mean=0, sd=sqrt(W2))

    theta1_sim[t] <- theta1
    theta2_sim[t] <- theta2
    y_sim[t] <- rpois(1, exp(theta1))
}


# Save the data
saveRDS(list(y = y_sim, theta1=theta1_sim, theta2=theta2_sim),
        paste(filename, ".rds", sep=""))

x <- 1:Tt
par(mar = c(4, 4, 2, 4)) # bottom left, top, right
plot(x, y_sim, type="l", col="gray", xlab="t", ylab=expression(y[t]))
points(x, y_sim, pch = 20)
lines(x, exp(theta1_sim), col="red")

plot(x, theta1_sim, type="l")
plot(x, theta2_sim, type="l")
