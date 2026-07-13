# Poisson with Sin Level Model simulation
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)


# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_triang_1000" # csv file to save the data
W1 <- 0.1

Tt <- 1000       # number of observations
theta <- c(seq(2, 5, length.out=Tt/2), seq(5, 2, length.out=Tt/2))

y <- rpois(Tt, exp(theta))

plot(y, pch=20)
lines(exp(theta), type="l", col="red")

# Save the data
saveRDS(list(y = y, theta=theta), paste("../data/", filename, ".rds", sep=""))
