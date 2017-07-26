## decide combination of (n, p) when n > p using 'mtd' and 'id'
## mtd = 5 => n = 50   & p = (50, 100, 250, 500, 2000) depends on id = 1 ... 50
##     = 4 => n = 100  & p = (100, 250, 500, 2000) depends on id = 1 ... 40
##     = 3 => n = 250  & p = (250, 500, 2000) depends on id = 1 ... 30
##     = 2 => n = 500  & p = (500, 2000) depends on id = 1 ... 20
## id  determines which replication and combination of (n, p)

## provide the (mtd, id) combination via command line
## for example, replace mtd and id below by 5 and 1 to run the analysis for n = 50, p = 50, and rep = 1
## R CMD BATCH --no-save --no-restore "--args mtd id" mdl_wts_nlp.R mdl_wts_nlp_mtd_id.rout
cmdArgs <- commandArgs(trailingOnly = TRUE)
mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

## 1. change the 'mtd' and 'id' appropriately to load the  correct data set
## 2. change the fname to store the result at the correct location
if (mtd == 5) {
    ## n = 50
    ndim <- c(50, 100, 250, 500, 2000)
    cvs <- rep(1:10, each = 5)
    dims <- rep(1:5, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p50.rds")
    train <- cvtrain[[cid]][[did]]

    obj <- paste0("../result/nlp/xfa_n_50_cv_", cid, "_p_", ndim[did], ".rds")
    fname <- paste0("../result/nlp/mdl_wts_n_50_cv_", cid, "_p_", ndim[did], ".rds")
} else if (mtd == 4) {
    ## n = 100
    ndim <- c(100, 250, 500, 2000)
    cvs <- rep(1:10, each = 4)
    dims <- rep(1:4, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p100.rds")
    train <- cvtrain[[cid]][[did]]

    obj <- paste0("../result/nlp/xfa_n_100_cv_", cid, "_p_", ndim[did], ".rds")
    fname <- paste0("../result/nlp/mdl_wts_n_100_cv_", cid, "_p_", ndim[did], ".rds")
} else if (mtd == 3) {
    ## n = 250
    ndim <- c(250, 500, 2000)
    cvs <- rep(1:10, each = 3)
    dims <- rep(1:3, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p250.rds")
    train <- cvtrain[[cid]][[did]]

    obj <- paste0("../result/nlp/xfa_n_250_cv_", cid, "_p_", ndim[did], ".rds")
    fname <- paste0("../result/nlp/mdl_wts_n_250_cv_", cid, "_p_", ndim[did], ".rds")
} else {
    ## n = 500
    ndim <- c(500, 2000)
    cvs <- rep(1:10, each = 2)
    dims <- rep(1:2, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p500.rds")
    train <- cvtrain[[cid]][[did]]

    obj <- paste0("../result/nlp/xfa_n_500_cv_", cid, "_p_", ndim[did], ".rds")
    fname <- paste0("../result/nlp/mdl_wts_n_500_cv_", cid, "_p_", ndim[did], ".rds")
}

woodMatInv <- function (Ainv, U, Cinv, V) {
    tmp <- chol(Cinv + V %*% Ainv %*% U)
    Ainv - Ainv %*% U %*% chol2inv(tmp) %*% V %*% Ainv
}

calcEBIC <- function(dataMat, lambdaMat, sigma2) {
    library(Matrix)
    xfactor <- ncol(lambdaMat)
    dat <- dataMat - matrix(colMeans(dataMat), nrow(dataMat), ncol(dataMat), byrow = TRUE)
    sYY <- crossprod(dat) / nrow(dat)
    omegaMatInv  <- woodMatInv(diag(1 / sigma2), lambdaMat, diag(1, xfactor), t(lambdaMat))
    gammaMat     <- omegaMatInv %*% lambdaMat
    psiMat       <- diag(1.0, xfactor) - crossprod(lambdaMat, gammaMat) + crossprod(gammaMat, sYY %*% gammaMat)
    lambdaHatMat <- sYY %*% gammaMat
    xmat         <- chol(psiMat)
    xmatTransInv <- forwardsolve(t(xmat), diag(1, xfactor))
    ## set up for ebic
    xx <- bdiag(rep(list(xmat), nrow(lambdaMat)))
    yy <- as.numeric(xmatTransInv %*% t(lambdaHatMat))
    wts <- 1 / rep(2 * sigma2, each = xfactor)

    nzeros <- as.numeric(abs(t(lambdaMat))) > 0
    wxx <- xx[ , nzeros, drop = FALSE] / sqrt(wts)
    wyy <- yy / sqrt(wts)
    chList <- list()
    for (pp in 1:nrow(lambdaMat)) {
        if (sum(abs(lambdaMat[pp, ]) > 0)) {
            tmpwxx <- xmat[ , abs(lambdaMat[pp, ]) > 0, drop = FALSE] * sqrt(2 * sigma2[pp])
            chList[[pp]] <- chol(crossprod(tmpwxx))
        } else {
            chList[[pp]] <- NULL
        }
    }
    ch0 <- bdiag(chList[!sapply(chList, is.null)])
    lscoef0 <- backsolve(ch0, forwardsolve(ch0, crossprod(wxx, wyy), upper = TRUE, trans = TRUE))
    resids0 <- yy - xx[ , nzeros, drop = FALSE] %*% lscoef0
    dev0 <- sum(resids0^2 * wts)
    ebic0 <- length(yy) * log(dev0 / length(yy)) + sum(nzeros) * (log(nrow(dataMat)) + 2 * log(ncol(dataMat) * log(ncol(dataMat))))
    ebic0
}

mdlMat <- matrix(0.0, 50, 20) ## nalpha <- 20; neta = 50 in xfa_nlp.R
for (aa in 1:20) {
    for (ee in 1:50) {
        warmLambda <- obj$xfa[[aa]][[ee]]$lambdaMat
        if (all(warmLambda == 0) || is.null(warmLambda)) {
            for (eee in ee:50) {
                mdlMat[eee, aa] <- Inf
            }
            break
        } else {
            warmSigma <- obj$xfa[[aa]][[ee]]$sigma2
            mdlMat[ee, aa] <- calcEBIC(train, warmLambda, warmSigma)
        }
    }
}

idx <- which(mdlMat == min(mdlMat[is.finite(mdlMat)]), arr.ind = TRUE)
colnames(idx) <- c("ee", "aa")
saveRDS(list(wts = mdlMat, idx = idx), fname)
