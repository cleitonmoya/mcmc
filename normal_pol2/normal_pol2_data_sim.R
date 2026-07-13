# 2nd Order Polynomial DLM
# Data simulation
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "normal_pol2_sim_2000" # csv file to save the data

Tt <- 2000
V <- 1000
theta_t1 <- 1
theta_t2 <- 1
W1 <- 1
W2 <- 0.05

theta1 <- numeric(Tt)
theta2 <- numeric(Tt)
y <- numeric(Tt)

# Simulation
for (t in 1:Tt) {

    theta_t1 <- theta_t1 + theta_t2 + rnorm(1, mean=0, sd=sqrt(W1))
    theta_t2 <-            theta_t2 + rnorm(1, mean=0, sd=sqrt(W2))
    y_t <- theta_t1 + rnorm(1, mean=0, sd=sqrt(V))

    y[t] <- y_t
    theta1[t] <- theta_t1
    theta2[t] <- theta_t2
}

# Save the data ####
saveRDS(list(y = y,
             theta1 = theta1,
             theta2 = theta2,
             V = V,
             W1 = W1,
             W2 = W2),
        paste("../data/", filename, ".rds", sep=""))

# Plot ####
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 4), cex=0.8) # bottom left, top, right
plot(x, y, type="l", col="gray", xlab="t", ylab=expression(y[t]))
points(x, y, pch = 20, cex=0.5)
lines(x, theta1, col="red")

plot(x, theta1, type="l")
plot(x, theta2, type="l")
