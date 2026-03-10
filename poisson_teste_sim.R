# Poisson 2nd Order Polynomial Model simulation
# Author: Cleiton Moya de Almeida

#####
graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "pol2_sim1" # csv file to save the data

set.seed(42)

theta1 <- 2.5  # theta_01
theta2 <- 0.05 # theta_02
W1 <- 0.3     # state variance W1
W2 <- 0.01   # state variance W2
Tt <- 1000      # number of observations

y_sim <- numeric(Tt)
theta1_sim <- numeric(Tt)
theta2_sim <- numeric(Tt)

for (t in 1:Tt) {
    theta1 <- theta1 + theta2 + rnorm(1, mean=0, sd=sqrt(W1))
    theta2 <- theta2 + rnorm(1, mean=0, sd=sqrt(W2))
    y <- thetat1 
    theta1_sim[t] <- theta1
    theta2_sim[t] <- theta2


}

#####
# Save the data
# df <- data.frame(y_sim, theta1_sim, theta2_sim)
# write.table(df, file = paste("data/", filename, ".csv", sep=""),
#            row.names=FALSE, col.names=c("y", "theta1", "theta2"))

x <- 1:Tt
par(mar = c(4, 4, 2, 4)) # bottom left, top, right
plot(x, theta1_sim, type="l")
lines(x, theta2_sim, lty=2)
