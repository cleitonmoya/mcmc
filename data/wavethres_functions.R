library(wavethresh)
set.seed(42)

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

data <- DJ.EX(plotfn = FALSE)

theta1 <- log(data$blocks + 9)
theta2 <- log(data$heavi + 15)
theta3 <- log(data$doppler + 13)

par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)

Tt <- 1024
y1 <- rpois(Tt, exp(theta1))
y2 <- rpois(Tt, exp(theta2))
y3 <- rpois(Tt, exp(theta3))

plot(y1, pch=20)
lines(exp(theta1), col="red", lw=2)

plot(y2, pch=20)
lines(exp(theta2), col="red", lw=2)

plot(y3, pch=20)
lines(exp(theta3), col="red", lw=2)


# Save the data
saveRDS(list(y = y1, theta=theta1), "blocks.rds")
saveRDS(list(y = y2, theta=theta2), "heavi.rds")
saveRDS(list(y = y3, theta=theta3), "doppler.rds")
