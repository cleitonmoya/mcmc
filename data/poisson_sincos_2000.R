#####
graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)


# change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
filename <- "poisson_sincos_2000" # csv file to save the data

Tt <- 2000
t <- seq(1, Tt, 1)
theta <- sin(4*pi*t/Tt) + cos(3*pi*t/(Tt/5))+2
y <- rpois(Tt, exp(theta))

plot(t, theta, type="l")
plot(1:Tt, y, pch=20, ylim=c(0,60))
lines(t, exp(theta), type="l", col="red")


# Save the data
saveRDS(list(y = y, theta=theta), paste(filename, ".rds", sep=""))
