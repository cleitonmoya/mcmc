# 2nd Order Polynomial DLM
# Data simulation
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "normal_local_level_sim_2000" # csv file to save the data

Tt <- 2000
V <- 10    # 0.1
W <- 0.1   # 0.01

theta_t <- 1
theta <- numeric(Tt)
y <- numeric(Tt)

# Simulation
for (t in 1:Tt) {
    theta_t <- theta_t + rnorm(1, mean=0, sd=sqrt(W))
    y[t]  <- theta_t + rnorm(1, mean=0, sd=sqrt(V))

    theta[t] <- theta_t
}

# Save the data ####
saveRDS(list(y = y,
             theta = theta,
             V = V,
             W = W),
        paste("../data/", filename, ".rds", sep=""))

# Plot ####
x <- 1:Tt
par(mar = c(4, 4, 2, 4)) # bottom left, top, right
plot(x, y, type="l", col="gray", xlab="t", ylab=expression(y[t]))
points(x, y, pch = 20, cex=0.5)
lines(x, theta, col="blue")
