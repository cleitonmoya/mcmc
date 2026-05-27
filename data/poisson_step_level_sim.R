# Poisson with Sin Level Model simulation
# Author: Cleiton Moya de Almeida

#####
graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename1 <- "poisson_step_200" # csv file to save the data
filename2 <- "poisson_step_2000" # csv file to save the data

set.seed(42)

alpha1 <- function(t) {
    0.2 * (t >= 0 & t < 0.3) +
        0.8 * (t >= 0.3 & t < 0.7) +
        0.3 * (t >= 0.7 & t <= 1.0)
}


T1 <- 200       # number of observations
T2 <- 2000
theta1 <- sapply(seq_len(T1)/T1, alpha1)
theta2 <- sapply(seq_len(T2)/T2, alpha1)

y1 <- rpois(T1, exp(theta1))
y2 <- rpois(T2, exp(theta2))

plot(theta1, type="l")
plot(y1, type="l")

plot(theta2, type="l")
plot(y2, type="l")

# Save the data
saveRDS(list(y = y1, theta = theta1),
        paste("../data/", filename1, ".rds", sep=""))

saveRDS(list(y = y2, theta = theta2),
        paste("../data/", filename2, ".rds", sep=""))
