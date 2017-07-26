## ------------
## example code
## ------------
source("simulation.R")

sigList <- seq(0.01, 1, length = 50)
loadList <- genLoadings(50, 5)
dataList <- matrix(0.0, nrow = 500, ncol = 50)

for (ii in 1:500) {
    dataList[ii, ] <- loadList %*% rnorm(ncol(loadList))
    for (jj in 1:50) {
        dataList[ii, jj] <- dataList[ii, jj] + rnorm(1, sd = sqrt(sigList[jj]))
    }
}

res <- xfa (dataList, maxFactor = 6, nalpha = 3, neta = 2)
