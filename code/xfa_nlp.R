## decide combination of (n, p) when n > p using 'mtd' and 'id'
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
    ## n = 50
    ndim <- c(50, 100, 250, 500, 2000)
    cvs <- rep(1:10, each = 5)
    dims <- rep(1:5, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p50.rds")
    train <- cvtrain[[cid]][[did]]

    fname <- paste0("../result/nlp/xfa_n_50_cv_", cid, "_p_", ndim[did], ".rds")
} else if (mtd == 4) {
    ## n = 100
    ndim <- c(100, 250, 500, 2000)
    cvs <- rep(1:10, each = 4)
    dims <- rep(1:4, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p100.rds")
    train <- cvtrain[[cid]][[did]]

    fname <- paste0("../result/nlp/xfa_n_100_cv_", cid, "_p_", ndim[did], ".rds")
} else if (mtd == 3) {
    ## n = 250
    ndim <- c(250, 500, 2000)
    cvs <- rep(1:10, each = 3)
    dims <- rep(1:3, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p250.rds")
    train <- cvtrain[[cid]][[did]]

    fname <- paste0("../result/nlp/xfa_n_250_cv_", cid, "_p_", ndim[did], ".rds")
} else {
    ## n = 500
    ndim <- c(500, 2000)
    cvs <- rep(1:10, each = 2)
    dims <- rep(1:2, times = 10)
    cid <- cvs[id]
    did <- dims[id]

    cvtrain <- readRDS("../data/p500.rds")
    train <- cvtrain[[cid]][[did]]

    fname <- paste0("../result/nlp/xfa_n_500_cv_", cid, "_p_", ndim[did], ".rds")
}

nobs <- nrow(train)
ndim <- ncol(train)

## center data
dataMat <- train - matrix(colMeans(train), nobs, ndim, byrow = TRUE)

woodMatInv <- function (Ainv, U, Cinv, V) {
    tmp <- chol(Cinv + V %*% Ainv %*% U)
    Ainv - Ainv %*% U %*% chol2inv(tmp) %*% V %*% Ainv
}

objVal <- function (lambdaMatOld, lambdaMatNew, sigma2, psiMat, lambdaHatMat, alphas, etas, nobs) {
    nfactor <- ncol(lambdaMatOld)
    ndim <- nrow(lambdaMatOld)

    psiMatSqrt <- chol(psiMat)
    alphasMat <- matrix(alphas + 1, nrow = ndim, ncol = nfactor, byrow = TRUE)
    etasMat   <- matrix(etas, nrow = ndim, ncol = nfactor, byrow = TRUE)

    loss <- 0.5 * nobs * (rowSums((lambdaHatMat %*% chol2inv(psiMatSqrt)) * lambdaHatMat) - 2 * rowSums(lambdaMatNew * lambdaHatMat) + rowSums((lambdaMatNew %*% psiMat) * lambdaMatNew)) / sigma2

    reg  <- rowSums(alphasMat / (etasMat + abs(lambdaMatOld)) * abs(lambdaMatNew))

    sum(reg) + sum(loss)
}

library(glmnet)

nobs <- nrow(dataMat)
ndim <- ncol(dataMat)

sYY <- crossprod(dataMat) / nobs
eig0 <- eigen(sYY)

## overfitted number of factors; no. of one-step polishes; grid size of alpha and eta; what is small
nfactor <- min(sum(eig0$value > 0), 20)
npolish <- 50; nalpha <- 20; neta <- 50; tol <- 1e-6
## root-n consistent estimates of lambda and sigma
lambdaMat <- eig0$vectors[ , 1:nfactor] * matrix(sqrt(eig0$values[1:nfactor]), nrow = ndim, ncol = nfactor, byrow = TRUE)
sigma2 <- diag(sYY - tcrossprod(lambdaMat, lambdaMat))

## define the alpha and eta grid
alphaMax <- 10
alphaMin <- 2
alphaSeq <- seq(log10(alphaMin), log10(alphaMax), length = nalpha)
etaSeq <- rev(seq(log10(0.01), log10(10000), length = neta))

## save result for every grid point
xfaList <- vector("list", length(alphaSeq))
names(xfaList) <- paste0("alpha.", round(alphaSeq, 4))
for (aa in 1:length(alphaSeq)) {
    xfaList[[aa]] <- vector("list", length(etaSeq))
    names(xfaList[[aa]]) <- paste0("eta.", etaSeq)
}

## ------------------------------------------------------------
## being xfa iterations ...
## ------------------------------------------------------------
startTime <- proc.time()
for (aa in 1:length(alphaSeq)) {
    for (ee in 1:length(etaSeq)) {
        cat("alpha: ", aa, "eta: ", ee, "\n")
        if (aa == 1 && ee == 1) {
            warmLambda <- lambdaMat
            warmSigma <- sigma2
        } else if (aa > 1 && ee == 1) {
            warmLambda <- xfaList[[aa - 1]][[1]]$lambdaMat
            warmSigma <- xfaList[[aa - 1]][[1]]$sigma2
            if (all(warmLambda == 0)) {
                break
            }
        } else {
            warmLambda <- xfaList[[aa]][[ee - 1]]$lambdaMat
            warmSigma <- xfaList[[aa]][[ee - 1]]$sigma2
            if (all(warmLambda == 0)) {
                break
            }
        }
        xfactor <- ncol(warmLambda)  # xfa factor <= overfitted factor
        sig1 <- rep(0, ndim)
        lam1 <- matrix(0.0, nrow = ndim, ncol = xfactor)
        obj <- rep(1e10, npolish + 1)
        for (pp in 1:npolish) {
            ## set up for calculations involving EM
            omegaMatInv  <- woodMatInv(diag(1 / warmSigma), warmLambda, diag(1, xfactor), t(warmLambda))
            gammaMat     <- omegaMatInv %*% warmLambda
            psiMat       <- diag(1.0, xfactor) - crossprod(warmLambda, gammaMat) + crossprod(gammaMat, sYY %*% gammaMat)
            lambdaHatMat <- sYY %*% gammaMat
            xmat         <- chol(psiMat)
            xmatTransInv <- forwardsolve(t(xmat), diag(1, xfactor))
            ## set up for glmnet
            xx <- bdiag(rep(list(xmat), ndim))
            yy <- as.numeric(xmatTransInv %*% t(lambdaHatMat))
            aseq <- rep((10^(alphaSeq[aa]))^(1:xfactor), ndim)
            penFac <- ((aseq + 1) / as.numeric((10^(etaSeq[ee]) * sqrt(ndim)) + t(abs(warmLambda)))) * 1 / (nobs * ndim * xfactor)
            wts <- 1 / rep(warmSigma, each = xfactor)
            tryfit <- tryCatch(
                cvres <- glmnet(xx, yy, weights = wts, intercept = FALSE, standardize = FALSE,
                                penalty.factor = penFac / sum(penFac), family = "gaussian",
                                thres = 1e-20)
               ,
                error = function(e) e
            )
            if (class(tryfit)[2] == "glmnet") {
                ## obtain glmnet-based estimates
                lam1 <- matrix(coef(cvres, s = sum(penFac), exact = TRUE)[-1, ], nrow = ndim, ncol = xfactor, byrow = TRUE)
                for (ss in 1:ndim) {
                    sig1[ss] <- (sYY[ss, ss] + sum(lam1[ss, ] * (psiMat %*% lam1[ss, ])) - 2 * sum(lam1[ss, ] * lambdaHatMat[ss, ])) * (nobs / (nobs + 2))
                }
                sig1 <- pmax(sig1, tol^2) ## avoids numerical errors when k is large
                if (pp %% 25 == 0) cat("*** done with polish no", pp, "***\n")
                ## have we converged?
                obj[pp + 1] <- objVal(warmLambda, lam1, sig1, psiMat, lambdaHatMat, (10^(alphaSeq[aa]))^(1:xfactor), rep(10^(etaSeq[ee]), xfactor), nobs)
                if (abs(obj[pp + 1] - obj[pp]) < tol) {
                    cat("done with grid: (", aa, ", ", ee, ")", " on polish ", pp, "\n")
                    break
                } else {
                    warmLambda <- lam1
                    warmSigma <- sig1
                }
            } else {
                break
            }
        }
        if (class(tryfit)[2] == "glmnet" && any(lam1 != 0)) {
            xfaList[[aa]][[ee]] <- list(lambdaMat = lam1[ , which(colSums(abs(lam1)) > 0), drop = FALSE],
                                        sigma2 = sig1
                                        )
        } else {
            for (eee in ee:length(etaSeq)) {
                xfaList[[aa]][[eee]] <- list(lambdaMat = matrix(0.0, ncol = 1, nrow = ndim),
                                             sigma2 = sig1
                                             )
            }
            break
        }
    }
}
endTime <- proc.time()

saveRDS(list(xfa = xfaList, time = endTime - startTime, agrid = alphaSeq, egrid = etaSeq), fname)
