# Poisson with Sin Level Model simulation
# Author: Cleiton Moya de Almeida

#####
graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)


# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename1 <- "poisson_sin_200" # csv file to save the data
filename2 <- "poisson_sin_2000" # csv file to save the data

T1 <- 200       # number of observations
T2 <- 2000
theta1 <- sin(2*pi*seq_len(T1)/T1)
theta2 <- sin(2*pi*seq_len(T2)/T2)

y1 <- rpois(T1, exp(theta1))
y2 <- rpois(T2, exp(theta2))

plot(theta1, type="l")
plot(theta2, type="l")
plot(y1)
plot(y2)

# Save the data
saveRDS(list(y = y1, theta=theta1),
        paste("../data/", filename1, ".rds", sep=""))

saveRDS(list(y = y2, theta=theta2),
        paste("../data/", filename2, ".rds", sep=""))
