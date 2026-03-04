# Poisson Local Level Model simulation
# Author: Cleiton Moya de Almeida

#####
graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "sim1" # csv file to save the data

set.seed(42)

theta <- log(3) # theta_0
W2 <- 0.03      # state variance W^2
Tt <- 200       # number of observations

y_sim <- numeric(Tt)
theta_sim <- numeric(Tt)

for (t in 1:Tt) {
    theta <- theta + rnorm(1, mean=0, sd=sqrt(W2))
    y <- rpois(1, exp(theta))

    theta_sim[t] <- theta
    y_sim[t] <- y
}

#####
# Save the data
df <- data.frame(y_sim, theta_sim)
write.table(df, file = paste("data/", filename, ".csv", sep=""),
            row.names=FALSE, col.names=c("y", "mu"))

x <- 1:Tt
par(mar = c(4, 4, 2, 4)) # bottom left, top, right
plot(x, y_sim, type="l", xlab="t", ylab="y_t")
points(x, y_sim, pch = 20)
lines(x, exp(theta_sim), col="red")
axis(side = 4)
mtext("mu_t", side = 4, line = 2.2)
