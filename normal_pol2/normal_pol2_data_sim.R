# 2nd Order Polynomial DLM
# Data simulation
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "normal_pol2_sim1" # csv file to save the data

T <- 300
V <- 0.5       # tau_v = 2
theta_t1 <- 1
theta_t2 <- 1
W1 <- 0.02     # tau_w1 = 50
W2 <- 0.005    # tau1 = 200

theta1 <- numeric(T)
theta2 <- numeric(T)
y <- numeric(T)

# Simulation
for (t in 1:T) {

    theta_t1 <- theta_t1 + theta_t2 + rnorm(1, mean=0, sd=sqrt(W1))
    theta_t2 <-            theta_t2 + rnorm(1, mean=0, sd=sqrt(W2))
    y_t <- theta_t1 + rnorm(1, mean=0, sd=sqrt(V))

    y[t] <- y_t
    theta1[t] <- theta_t1
    theta2[t] <- theta_t2
}

# Save the data ####
df <- data.frame(y, theta1, theta2, V, W1, W2)
write.table(df, file = paste("../data/", filename, ".csv", sep=""),
            row.names=FALSE, col.names=c("y", "theta1", "theta2",
                                         "V", "W1", "W2"))

# Plot ####
x <- 1:T
par(mar = c(4, 4, 2, 4)) # bottom left, top, right
plot(x, y, type="l", col="gray", xlab="t", ylab=expression(y[t]))
points(x, y, pch = 20, cex=0.5)
lines(x, theta1, col="red")

#####
plot(x, theta1, type="l")
plot(x, theta2, type="l")
