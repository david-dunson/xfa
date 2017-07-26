## Run this file to store all the data sets. The data sets will be used in the simulation.

rm(list=ls())
## Simulation set up of Kneip and Sarda (2011). This not a general code 
## and limited to generating the loadings described in the paper. 
## Kneip & Sarda (2011). Factor models and variable selection in
## high-dimensional regression analysis. The Annals of Statistics 39, 2410–2447.
genLoadings <- function(ndim, nfactor = 5) {
    lambdaMat <- matrix(0.0, nrow = ndim, ncol = nfactor)
    nzeros <- ndim / nfactor  ## we assume the nfactor = 5 currently
    loads <- rev(seq(2, 10, length = 5))

    for (ff in 1:nfactor) {
        idx <- ((ff - 1) * nzeros + 1):(ff * nzeros)
        lambdaMat[idx, ff] <- loads[ff]
    }
    for (ff in 2:nfactor) {
        idx <- ((ff - 1) * nzeros - (0.02 * ndim)):((ff - 1) * nzeros)
        lambdaMat[idx, ff] <- -loads[ff]
    }
    for (ff in 1:(nfactor - 1)) {
        idx <- (ff * nzeros + 1):(ff * nzeros + 1 + (0.02 * ndim))
        lambdaMat[idx, ff] <- loads[ff]
    }
    lambdaMat
}

ndim <- c(50, 100, 250, 500, 2000)
loadList <- vector("list", length(ndim))
sigList <- vector("list", length(ndim))
names(loadList) <- paste0("p", ndim)
names(sigList) <- paste0("p", ndim)
for (dd in seq_along(ndim)) {
    loadList[[dd]] <- genLoadings(ndim[dd])
    sigList[[dd]] <- seq(0.01, 1, length = ndim[dd])
}

nobs <- c(50, 100, 250, 500, 5000)
dataList <- vector("list", length(nobs))
names(dataList) <- paste0("n", nobs)
for (nn in seq_along(nobs)) {
    dataList[[nn]] <- vector("list", length(ndim))
    names(dataList[[nn]]) <- paste0("p", ndim)
}

set.seed(12345)

## simulate data for all combinations of n and p
cvdata <- vector("list", 10)
names(cvdata) <- paste0("cv", 1:length(cvdata))
for (cc in 1:10) {
    for (nn in seq_along(nobs)) {
        for (dd in seq_along(ndim)) {
            dataList[[nn]][[dd]] <- matrix(0.0, nrow = nobs[nn], ncol = ndim[dd])
            for (ii in 1:nobs[nn]) {
               dataList[[nn]][[dd]][ii, ] <- loadList[[dd]] %*% rnorm(ncol(loadList[[dd]]))
                for (jj in 1:ndim[dd]) {
                    dataList[[nn]][[dd]][ii, jj] <- dataList[[nn]][[dd]][ii, jj] + rnorm(1, sd = sqrt(sigList[[dd]][jj]))
                }
            }
        }
    }
    cvdata[[cc]] <- dataList
}
saveRDS(cvdata, "../data/simulated_data.rds")

## When n > p, store separate cases
## n = 5000
dat5000 <- list("vector", 10)
for (ii in 1:10) {
    dat5000[[ii]] <- cvdata[[ii]][[5]]
}
saveRDS(dat5000, "../data/n5000.rds")

## n = 500
dat500 <- list("vector", 10)
for (ii in 1:10) {
    dat500[[ii]] <- cvdata[[ii]][[4]][1:3]
}
saveRDS(dat500, "../data/n500.rds")

## n = 250
dat250 <- list("vector", 10)
for (ii in 1:10) {
    dat250[[ii]] <- cvdata[[ii]][[3]][1:2]
}
saveRDS(dat250, "../data/n250.rds")

## n = 100
dat100 <- list("vector", 10)
for (ii in 1:10) {
    dat100[[ii]] <- cvdata[[ii]][[2]][[1]]
}
saveRDS(dat100, "../data/n100.rds")


## When n <= p, store separate cases
## n = 50
dat50 <- list("vector", 10)
for (ii in 1:10) {
    dat50[[ii]] <- cvdata[[ii]][[1]]
}
saveRDS(dat50, "../data/p50.rds")

## n = 100
dat100 <- list("vector", 10)
for (ii in 1:10) {
    dat100[[ii]] <- cvdata[[ii]][[2]][2:5]
}
saveRDS(dat100, "../data/p100.rds")

## n = 250
dat250 <- list("vector", 10)
for (ii in 1:10) {
    dat250[[ii]] <- cvdata[[ii]][[3]][3:5]
}
saveRDS(dat250, "../data/p250.rds")

## n = 500
dat500 <- list("vector", 10)
for (ii in 1:10) {
    dat500[[ii]] <- cvdata[[ii]][[4]][4:5]
}
saveRDS(dat500, "../data/p500.rds")
