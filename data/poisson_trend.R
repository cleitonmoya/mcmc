# Poisson with Sin Level Model simulation
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)


# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_trend_1000" # csv file to save the data
Tt <- 1000

W1_true <- 0.005
W2_true <- 0.0005

theta2 <- numeric(Tt); theta <- numeric(Tt)
theta[1] <- 3; theta2[1] <- 0.003

for(t in 2:Tt){
    theta2[t] <- theta2[t-1] + rnorm(1, 0, sqrt(W2_true))
    theta[t]  <- theta[t-1] + theta2[t-1] + rnorm(1, 0, sqrt(W1_true))
    # reflective boundary em [1.5, 5.5]
    if(theta[t] > 5.5){ theta[t] <- 2*5.5 - theta[t]; theta2[t] <- -theta2[t] }
    if(theta[t] < 1.5){ theta[t] <- 2*1.5 - theta[t]; theta2[t] <- -theta2[t] }
}

y <- rpois(Tt, exp(theta))

plot(y, pch=20)
lines(exp(theta), type="l", col="red")

# Save the data
saveRDS(list(y = y, theta=theta), paste("../data/", filename, ".rds", sep=""))
