# 2nd Order Polynomial DLM
# Data simulation
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(41)

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "normal_local_level_sim3" # csv file to save the data

T <- 10000
V <- 0.5       # tau_v = 2
W <- 1     # tau_w1 = 50
theta_t <- 1

theta <- numeric(T)
y <- numeric(T)

# Simulation
for (t in 1:T) {

    theta_t <- theta_t + rnorm(1, mean=0, sd=sqrt(W))
    y_t <- theta_t + rnorm(1, mean=0, sd=sqrt(V))

    y[t] <- y_t
    theta[t] <- theta_t
}

# Save the data ####
saveRDS(list(y = y,
             theta = theta,
             V = V,
             W = W),
        paste("../data/", filename, ".rds", sep=""))


# Plot ####
x <- 1:T
par(mar = c(4, 4, 2, 4)) # bottom left, top, right
plot(x, y, type="l", col="gray", xlab="t", ylab=expression(y[t]))
points(x, y, pch = 20, cex=0.5)
lines(x, theta, col="red")

#####
plot(x, theta, type="l")
