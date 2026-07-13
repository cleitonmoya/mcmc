# Poisson with Sin Level Model simulation
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)


# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_sin_noise_200" # csv file to save the data
W1 <- 0.1

Tt <- 200       # number of observations
theta <- sin(2*pi*seq_len(Tt)/Tt) + rnorm(Tt, 0, sqrt(W1))

y <- rpois(Tt, exp(theta))

plot(theta, type="l")
plot(y)

# Save the data
saveRDS(list(y = y, theta=theta), paste("../data/", filename, ".rds", sep=""))
