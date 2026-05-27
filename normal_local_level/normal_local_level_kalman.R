# Forward Filtering (Kalman Filter)

graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)
tp <- base::t       # alias to transpose function
options(error = function() traceback(2)) # more informative traceback


# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source <- "normal_local_level_sim" # csv file with data
data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
T <-length(y)
theta_true <- data$theta
W_true <- data$W
V_true <- data$V

# initialization
V <- 0.1
W <- 0.01
theta_t <- 1

theta <- numeric(T)

# theta_0 ~ N(m_0, C_0)
m0     <- y[1]
C0 <- 10

# Forward Filtering (W&H, Seção 4.3, Teorema 4.1)
m_kf <- numeric(T)   # m_t = E[theta_t | D_t]
C_kf <- numeric(T)   # C_t = Var[theta_t | D_t]

# t = 1: prior (theta_1 | D_0) ~ N(a_1, R_1)
a1   <- m0
R1   <- C0 + W
Q1   <- R1 + V
A1   <- R1 / Q1
e1   <- y[1] - a1
m_kf[1] <- a1 + A1 * e1
C_kf[1] <- (1 - A1) * R1

for (t in 2:T) {
    a_t     <- m_kf[t-1]           # prior mean
    R_t     <- C_kf[t-1] + W       # prior variance
    Q_t     <- R_t + V             # forecast variance
    A_t     <- R_t / Q_t           # adaptive coefficient
    e_t     <- y[t] - a_t          # forecast error
    m_kf[t] <- a_t + A_t * e_t     # posterior mean
    C_kf[t] <- (1 - A_t) * R_t     # posterior variance
}


x <- 1:T

# theta filtered ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta_true, type="l", xlab="t", ylab="")
lines(x, m_kf, col="blue")
polygon(c(x, rev(x)),
        c(m_kf+ 3*sqrt(C_kf), rev(m_kf - 3*sqrt(C_kf))),
        col = adjustcolor("steelblue", alpha.f = 0.3),
        border = NA)


