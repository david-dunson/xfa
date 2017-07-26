## decide combination of (n, p) when n <= p using 'mtd' and 'id'
## mtd = 5 => n = 50   & p = (50, 100, 250, 500, 2000) depends on id = 1 ... 50
##     = 4 => n = 100  & p = (100, 250, 500, 2000) depends on id = 1 ... 40
##     = 3 => n = 250  & p = (250, 500, 2000) depends on id = 1 ... 30
##     = 2 => n = 500  & p = (500, 2000) depends on id = 1 ... 20
## id  determines which replication and combination of (n, p)

## provide the (mtd, id) combination via command line
## for example, replace mtd and id below by 5 and 1 to run the analysis for n = 50, p = 50, and rep = 1
## R CMD BATCH --no-save --no-restore "--args mtd id" competitors_nlp.R comp_nlp_mtd_id.rout
cmdArgs <- commandArgs(trailingOnly = TRUE)
mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

## 1. run the simulation.R file to store the *.rds data set in data directory
## 2. change the 'mtd' and 'id' appropriately to load the  correct data set
## 3. change the fname* to store the result at the correct location
if (mtd == 5) {
    ndim <- c(50, 100, 250, 500, 2000)
    cvs <- rep(1:10, each = 5)
    dims <- rep(1:5, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p50.rds")
    train <- cvtrain[[cid]][[did]]

    fname1 <- paste0("../result/nlp/fanc_n_50_cv_", cid, "_p_", ndim[did], ".rds")
    fname2 <- paste0("../result/nlp/spca_n_50_cv_", cid, "_p_", ndim[did], ".rds")
    fname3 <- paste0("../result/nlp/rg_n_50_cv_", cid, "_p_", ndim[did], ".rds")
    fname4 <- paste0("../result/nlp/rgp_n_50_cv_", cid, "_p_", ndim[did], ".rds")
} else if (mtd == 4) {
    ndim <- c(100, 250, 500, 2000)
    cvs <- rep(1:10, each = 4)
    dims <- rep(1:4, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p100.rds")
    train <- cvtrain[[cid]][[did]]

    fname1 <- paste0("../result/nlp/fanc_n_100_cv_", cid, "_p_", ndim[did], ".rds")
    fname2 <- paste0("../result/nlp/spca_n_100_cv_", cid, "_p_", ndim[did], ".rds")
    fname3 <- paste0("../result/nlp/rg_n_100_cv_", cid, "_p_", ndim[did], ".rds")
    fname4 <- paste0("../result/nlp/rgp_n_100_cv_", cid, "_p_", ndim[did], ".rds")
} else if (mtd == 3) {
    ndim <- c(250, 500, 2000)
    cvs <- rep(1:10, each = 3)
    dims <- rep(1:3, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p250.rds")
    train <- cvtrain[[cid]][[did]]

    fname1 <- paste0("../result/nlp/fanc_n_250_cv_", cid, "_p_", ndim[did], ".rds")
    fname2 <- paste0("../result/nlp/spca_n_250_cv_", cid, "_p_", ndim[did], ".rds")
    fname3 <- paste0("../result/nlp/rg_n_250_cv_", cid, "_p_", ndim[did], ".rds")
    fname4 <- paste0("../result/nlp/rgp_n_250_cv_", cid, "_p_", ndim[did], ".rds")
} else {
    ndim <- c(500, 2000)
    cvs <- rep(1:10, each = 2)
    dims <- rep(1:2, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p500.rds")
    train <- cvtrain[[cid]][[did]]

    fname1 <- paste0("../result/nlp/fanc_n_500_cv_", cid, "_p_", ndim[did], ".rds")
    fname2 <- paste0("../result/nlp/spca_n_500_cv_", cid, "_p_", ndim[did], ".rds")
    fname3 <- paste0("../result/nlp/rg_n_500_cv_", cid, "_p_", ndim[did], ".rds")
    fname4 <- paste0("../result/nlp/rgp_n_500_cv_", cid, "_p_", ndim[did], ".rds")
}

nobs <- nrow(train)
ndim <- ncol(train)

## data matrix; overfitted number of factors; what is small
dataMat <- train; nfactor <- 20; tol = 1e-5

## FANC
library(fanc)
ctrl <- list(length.rho = 20, length.gamma = 20, maxit.em = 1000, maxit.cd = 1000, tol.cd = tol, tol.em = tol)
startTime <- proc.time()
fancfit <- fanc(train, factors = nfactor, control = ctrl, normalize = FALSE)
## rows represent rhos and columns represent gammas; notice that BIC
## does not change along the columns; only rhos lead to difference in
## the BIC due to desriable property of MC+ penalty.
idx <- which(fancfit$BIC == min(fancfit$BIC), arr.ind = TRUE)[1, ]
fancout <- out(fancfit, rho = fancfit$rho[idx[1]], gamma = fancfit$gamma[idx[2]])
fancLoad <- as.matrix(fancout$loadings)
fancIdx <- which(colSums(abs(fancLoad)) > 0)
fancRank <- length(fancIdx)
fancLoad <- as.matrix(fancLoad[ , fancIdx])
fancUniq <- as.numeric(fancout$uniquenesses)
endTime <- proc.time()
saveRDS(list(fanc = fancout, time = endTime - startTime, rank = fancRank, load = fancLoad, uniq = fancUniq), fname1)

## SPCA
library(PMA)
startTime <- proc.time()
spccv <- SPC.cv(train)
spcfit <- SPC(train, K = nfactor, sumabsv = spccv$bestsumabs)
spcRank <- which(diff(spcfit$prop.var.explained) < 0.05)[1] ## Algorithm 3 of Shen and Huang (2008)
if (is.na(spcRank)) {
  spcRank <- length(spcfit$prop.var.explained)
}
spcLoad <- spcfit$v * matrix(sqrt(spcfit$d), nrow = ndim, ncol = length(spcfit$d), byrow = TRUE)
endTime <- proc.time()
saveRDS(list(spca = spcfit, time = endTime - startTime, rank = spcRank, load = spcLoad), fname2)

## Rockova & George (no varimax)
source("rockova_george.R")

library(mvtnorm)
library(partitions)
library(nloptr)

## see the rockova_george.R file for appropriate definition of Y, n, G, p, K
Y <- train
n <- nrow(train)
G <- ncol(train)
p <- 10
K <- 20

startB <-matrix(rnorm(G*K),G,K)
alpha <- 1/G

lambda1<-0.001
epsilon<-0.05

start<-list(B=startB,sigma=rep(1,p),theta=rep(0.5,K))
startTime <- proc.time()
lambda0<-5
result_5<-FACTOR_ROTATE(Y,lambda0,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100, 0)

lambda0<-10
result_10<-FACTOR_ROTATE(Y,lambda0,lambda1,result_5,K,epsilon,alpha,TRUE,TRUE,100, 0)

lambda0<-20
result_20<-FACTOR_ROTATE(Y,lambda0,lambda1,result_10,K,epsilon,alpha,TRUE,TRUE,100, 0)

lambda0<-30
result_30<-FACTOR_ROTATE(Y,lambda0,lambda1,result_20,K,epsilon,alpha,TRUE,TRUE,100, 0)
endTime <- proc.time()

saveRDS(list("lam5" = result_5, "lam10" = result_10, "lam20" = result_20, "lam30" = result_30, time = endTime - startTime), fname3)

## Rockova & George (with varimax)
source("rockova_george.R")

library(mvtnorm)
library(partitions)
library(nloptr)

## see the rockova_george.R file for appropriate definition of Y, n, G, p, K
Y <- train
n <- nrow(train)
G <- ncol(train)
p <- 10
K <- 20

startB <-matrix(rnorm(G*K),G,K)
alpha <- 1/G

lambda1<-0.001
epsilon<-0.05

start<-list(B=startB,sigma=rep(1,p),theta=rep(0.5,K))
startTime <- proc.time()
lambda0<-5
result_5<-FACTOR_ROTATE(Y,lambda0,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100, 1)

lambda0<-10
result_10<-FACTOR_ROTATE(Y,lambda0,lambda1,result_5,K,epsilon,alpha,TRUE,TRUE,100, 1)

lambda0<-20
result_20<-FACTOR_ROTATE(Y,lambda0,lambda1,result_10,K,epsilon,alpha,TRUE,TRUE,100, 1)

lambda0<-30
result_30<-FACTOR_ROTATE(Y,lambda0,lambda1,result_20,K,epsilon,alpha,TRUE,TRUE,100, 1)
endTime <- proc.time()

saveRDS(list("lam5" = result_5, "lam10" = result_10, "lam20" = result_20, "lam30" = result_30, time = endTime - startTime), fname4)
