# Poisson with Sin Level Model simulation
# Author: Cleiton Moya de Almeida

#####
graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console

# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "sin_level" # csv file to save the data

set.seed(42)

Tt <- 200       # number of observations
theta_sim <- sin(2*pi*seq_len(Tt)/(Tt))
y_sim <- rpois(Tt, exp(theta_sim))

#####
# Save the data
df <- data.frame(y_sim, theta_sim)
write.table(df, file = paste("data/", filename, ".csv", sep=""),
            row.names=FALSE, col.names=c("y", "theta1"))

#####
x <- 1:Tt
par(mar = c(4, 4, 2, 4)) # bottom left, top, right
plot(x, y_sim, type="l", col="gray", xlab="t", ylab=expression(y[t]))
points(x, y_sim, pch = 20)
lines(x, exp(theta_sim), col="red")
