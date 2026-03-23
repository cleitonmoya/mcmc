# Calculation of Posterior Distributions by Data Augmentation
# Tanner and Wong, 1984
# Linkage Model example
# Author: Cleiton Moya de Almeida

#####
graphics.off()      # close the plots
rm(list = ls())     # clear the environment
cat("\014")         # clear the console
set.seed(42)

# Observed likelihood:
# y|theta ~ Multinomial(1/2+theta/4, (1-theta)/4, (1-theta)/4, theta/4)

# Linkage model - section 2,  p. 531
# y1 <- 125
# y2 <- 18
# y3 <- 20
# y4 <- 34

# section 7, p. 540 (Practical implementation example)
y1 <- 12
y2 <- 2
y3 <- 2
y4 <- 3

# Augmented likelihood:
# x|theta ~ Beta(x2+x1+1, x3+x4+1)
# x1 + x2 = y1
x3 <- y2
x4 <- y3
x5 <- y4

# Prior:
#    theta ~ Uniform(0, 1) \equiv Beta(alpha_star=1, beta_star=1)
# Approximated density function to p(theta | y)
#    p(theta | y) \approx p(theta | x) ~ 1/m \sum_{i=1}^m beta(alpha_star_i, beta_star_i)
g_params <- list(m=1, alpha_star_vect = 1, beta_star = 1)


# Sample g(theta | x)
sample_theta <- function(g) {
    m <- g$m
    alpha_star_vect <- g$alpha_star_vect
    beta_star <- g$beta_star
    j <- sample(1:m, 1)
    theta <- rbeta(1, alpha_star_vect[j], beta_star)
    return(theta)
}


# "Adaptive" m (section 7)
m_vect <- c(rep(1600,40), rep(400,20), rep(1600, 10))
nit <- length(m_vect)
theta25_hist <- numeric(nit)
theta50_hist <- numeric(nit)
theta75_hist <- numeric(nit)
g_list <- vector("list", nit)

#####
for (t in 1:nit) {

    m <- m_vect[t]
    g_list[[t]] <- g_params
    alpha_star_vect <- numeric(m)
    theta_vect <- numeric(m)
    x2_vect <- numeric(m)

    # Imputation Step
    for (i in 1:m) {

        # sampling theta
        theta <- sample_theta(g_params)
        theta_vect[i] <- theta

        # sampling x2
        x2 <- rbinom(1, y1, theta/(theta+2))
        x2_vect[i] <- x2
    }

    # Posterior step
    # update the g function
    alpha_star_vect <- x2_vect + x5 + 1
    beta_star <- x3 + x4 + 1
    g_params <- list(m=m, alpha_star_vect = alpha_star_vect, beta_star = beta_star)

    # Percentiles monitoring
    theta25_hist[t] <- quantile(theta_vect, 0.25)[[1]]
    theta50_hist[t] <- median(theta_vect)
    theta75_hist[t] <- quantile(theta_vect, 0.75)[[1]]
}


#####
# Percentiles monitoring
plot(theta50_hist, type="l", ylim=c(0, 1),
     xlab="Iteration", ylab="Percentile")
lines(theta25_hist, lty=3)
lines(theta75_hist, lty=2)
legend(x="top", legend=c("25%","50%","75%"),
       lty=c(3,1,2), bty="n", horiz=TRUE, cex=0.8)

# Estimated posterior density
