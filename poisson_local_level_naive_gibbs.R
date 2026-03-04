# Poisson Local Level Model
# MCMC: naive Gibbs Sampling
# Author: Cleiton Moya de Almeida

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- Matrix::t     # matrix transpose alias

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "sim1" # csv file with data
df <- read.table(paste("data/", source, ".csv", sep=""), header = TRUE)
y <- df$y
theta_true <- df$mu

Tt <- length(y) # dimension T

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

# fot t=1<Tt <- length(y) # dimension T

dTt <- length(y) # dime


Tt <- length(y) # dimension T
Tt <- length(y) # dimension T


dsfsdfdTt <- length(y) # dimension T





