# Poisson Local Level Model simulation
# Author: Cleiton Moya de Almeida

#####
graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_local_level_sim1" # rds file to save the data

set.seed(42)

theta <- log(3) # theta_0
W <- 0.03      # state variance W
Tt <- 200       # number of observations

y_sim <- numeric(Tt)
theta_sim <- numeric(Tt)

for (t in 1:Tt) {
    theta <- theta + rnorm(1, mean=0, sd=sqrt(W))
    y <- rpois(1, exp(theta))

    theta_sim[t] <- theta
    y_sim[t] <- y
}

# Save the data
saveRDS(list(y = y_sim,
             theta = theta_sim,
             W = W),
        paste("../data/", filename, ".rds", sep=""))

# Plot ####
x <- 1:Tt
par(mar = c(4, 4, 2, 4)) # bottom left, top, right
plot(x, y_sim, type="l", col="gray", xlab="t", ylab=expression(y[t]))
points(x, y_sim, pch = 20)
lines(x, exp(theta_sim), col="red", lwd=2)

plot.ts(theta_sim)

