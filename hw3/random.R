#
library(MASS)

rm(list = ls())

x <- rnorm(1000000)

# read in the code
data <- read.csv("/Users/nitish/MSMF/CompFin/DesignCPP/hw3/parkmiller.dat", header = TRUE)

par(mfrow=c(3, 1))

truehist(data$Inverse, h = 0.05, prob = TRUE, xlim = c(-4, 4), col="white", xlab="Inverse Transform Method")
lines(density(x, kernel="gaussian"), col = "red")

truehist(data$BoxMuller, h = 0.05, prob = TRUE, xlim = c(-4, 4), col="white", xlab="Box-Muller Method")
lines(density(x, kernel="gaussian"), col = "red")

truehist(data$Fishman, h = 0.05, prob = TRUE, xlim = c(-4, 4), col="white", xlab="Fishman Method")
lines(density(x, kernel="gaussian"), col = "red", xlim=c(-4, 4))

