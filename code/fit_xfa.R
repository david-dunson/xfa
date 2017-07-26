## Implementation of matrix inversion lemma
## (A + UCV)^{-1} = A^{-1} - A^{-1} U (C^{-1} + V A^{-1} U) V A^{-1}, where
## Ainv = A^{-1} and Cinv = C^{-1}.
woodMatInv <- function (Ainv, U, Cinv, V) {
    tmp <- chol(Cinv + V %*% Ainv %*% U)
    Ainv - Ainv %*% U %*% chol2inv(tmp) %*% V %*% Ainv
}

## value of the objective function in the expandable factor analysis
## model; see equation (6) in Srivastava et al. (2017).
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

##' Fit expandable factor analysis model across a specified alpha-eta grid.
##'
##' This function implements Algorithm in Srivastava et
##'     al. (2017) for fitting an expandable factor analysis model.
##'     The rank of the loadings matrix is assumed to be unknown but
##'     the rank is O(log p), where p is the dimension of the observations.
##' @title Fit expandable factor analysis model
##' @param y data matrix. Rows represent samples and columns index samples.
##' @param maxFactor upper bound on the true number of factors.
##' @param nalpha size of the $\alpha$ grid.
##' @param neta size of the $\eta$ grid.
##' @param alphaRange largest and smallest values of alpha for which
##'     the size 'nalpha' grid is constructed.
##' @param etaRange largest and smallest values of eta for which the
##'     size 'neta' grid is constructed.
##' @param npolish number of one-step iterations to be performed.
##' @param tol our definition of small.
##' @return a list containing the estimates of loadings matrices and
##'     residual variances across the alpha-eta grid, values of the
##'     alpha and eta across the grid, and total time required in
##'     estimating all the parameters across the grid.
##' @author Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
fitXfa <- function (y, maxFactor = 20, nalpha = 20, neta = 50, alphaRange = c(2, 10), etaRange = c(0.01, 1e5), npolish = 50, tol = 1e-6) {
    library(matrixStats)
    library(glmnet)

    ## dimensions and no. of obs.
    nobs <- nrow(y)
    ndim <- ncol(y)

    ## center data
    dataMat <- y - matrix(colMeans(y), nobs, ndim, byrow = TRUE)
    sYY <- crossprod(dataMat) / nobs
    eig0 <- eigen(sYY)

    ## overfitted number of factors
    if (nobs > ndim) {
        nfactor <- maxFactor
    } else {
        nfactor <- min(sum(eig0$value > 0), maxFactor)
    }

    ## root-n consistent estimates of Lambda and Sigma
    lambdaMat <- eig0$vectors[ , 1:nfactor] * matrix(sqrt(eig0$values[1:nfactor]), nrow = ndim, ncol = nfactor, byrow = TRUE)
    sigma2 <- diag(sYY - tcrossprod(lambdaMat, lambdaMat))

    ## define the alpha and eta grid
    alphaMax <- alphaRange[2]; alphaMin <- alphaRange[1]
    alphaSeq <- seq(log10(alphaMin), log10(alphaMax), length = nalpha)
    etaMax <- etaRange[2]; etaMin <- etaRange[1]
    etaSeq <- rev(seq(log10(etaMin), log10(etaMax), length = neta))

    ## save result for every grid point
    xfaList <- vector("list", length(alphaSeq))
    names(xfaList) <- paste0("alpha.", round(alphaSeq, 4))
    for (aa in 1:length(alphaSeq)) {
        xfaList[[aa]] <- vector("list", length(etaSeq))
        names(xfaList[[aa]]) <- paste0("eta.", etaSeq)
    }

    ## being xfa iterations ...
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
            sig1 <- rep(0.0, ndim)
            lam1 <- matrix(0.0, nrow = ndim, ncol = xfactor)
            obj <- rep(1e10, npolish + 1)
            for (pp in 1:npolish) {
                ## set up statistics for calculations involving EM
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
                    cvres <- glmnet(x = xx, y = yy, weights = wts, intercept = FALSE, standardize = FALSE,
                                    penalty.factor = penFac / sum(penFac), family = "gaussian",
                                    thres = 1e-20)
                   ,
                    error = function(e) e
                )
                if (class(tryfit)[2] == "glmnet") {
                    ## obtain glmnet-based estimates
                    gres <- predict(object = cvres, newx = xx, s = sum(penFac), type = "coefficients", exact = TRUE, offset = NULL,
                                    x = xx, y = yy, weights = wts, intercept = FALSE, standardize = FALSE, penalty.factor = penFac / sum(penFac))
                    lam1 <- matrix(gres[-1, ], nrow = ndim, ncol = xfactor, byrow = TRUE)
                    for (ss in 1:ndim) {
                        sig1[ss] <- (sYY[ss, ss] + sum(lam1[ss, ] * as.numeric(psiMat %*% lam1[ss, ])) - 2 * sum(lam1[ss, ] * lambdaHatMat[ss, ])) * (nobs / (nobs + 2))
                    }
                    sig1 <- pmax(sig1, tol^2) ## avoids numerical errors when k is large
                    if (pp %% 50 == 0) cat("*** done with polish no", pp, "***\n")
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
                                                 sigma2 = rep(0.0, ndim)
                                                 )
                }
                break
            }
        }
    }
    endTime <- proc.time()

    res <- list(xfa = xfaList, time = endTime - startTime, agrid = alphaSeq, egrid = etaSeq)

    res
}

##' Caluclate the value of EBIC criterion in a factor model
##'
##' This function calculates the EBIC criterion of Chen and Chen
##'     (2008) in factor models for a given data, loadings matrix, and
##'     residual variances.
##' @title EBIC criterion in a factor model
##' @param dataMat data matrix. Rows represent samples and columns index samples.
##' @param lambdaMat loadings matrix. This is a p x k matrix, where 'p'
##'     is the column dimension of 'dataMat' and 'k' is the number of
##'     factors.
##' @param sigma2 residual variances. This is p x 1 vector.
##' @return value of the EBIC criterion.
##' @author  Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
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

##' Fit exandable factor analysis model
##'
##' This function estimates the loadings matrix in a factor model. It
##'     first uses the 'fitXfa' to estimate loadings and residual
##'     variance across the specified alpha-eta grid. The estimates of
##'     loadings and residual variances are used to calculate the EBIC
##'     criterion for factor analysis using the 'calcEBIC'. The 
##'     selected loadings matrix corresponds to the alpha-eta grid
##'     index with mininum value of EBIC.
##' @title Expandable factor analysis
##' @param y data matrix. Rows represent samples and columns index samples.
##' @param maxFactor upper bound number of factors.
##' @param nalpha size of the $\alpha$ grid.
##' @param neta size of the $\eta$ grid.
##' @param alphaRange largest and smallest values of alpha for which
##'     the size 'nalpha' grid is constructed.
##' @param etaRange largest and smallest values of alpha for which the
##'     size 'neta' grid is constructed.
##' @param npolish number of one-step iterations to be performed.
##' @param tol our definition of small.
##' @return return an object of class 'xfa' consisting of a list of
##'     selected loadings and residual variances corresponding to this
##'     loadings matrix. The estimates of loadings and residual
##'     variances across the alpha-eta are also returned.
##' @author Sanvesh Srivastava (sanvesh-srivastava@uiowa.edu)
xfa <- function (y, maxFactor = 20, nalpha = 20, neta = 50, alphaRange = c(2, 10), etaRange = c(0.01, 1e5), npolish = 50, tol = 1e-6) {

    xfaFit <- fitXfa(y, maxFactor = maxFactor, nalpha = nalpha, neta = neta, alphaRange = alphaRange, etaRange = etaRange, npolish = npolish, tol = tol)

    mdlMat <- matrix(0.0, neta, nalpha)
    startTime <- proc.time()
    for (aa in 1:nalpha) {
        for (ee in 1:neta) {
            warmLambda <- xfaFit$xfa[[aa]][[ee]]$lambdaMat
            if (all(warmLambda == 0)) {
                for (eee in ee:neta) {
                    mdlMat[eee, aa] <- Inf
                }
                break
            } else {
                warmSigma <- xfaFit$xfa[[aa]][[ee]]$sigma2
                mdlMat[ee, aa] <- calcEBIC(y, warmLambda, warmSigma)
            }
        }
    }
    endTime <- proc.time()

    idx <- which(mdlMat == min(mdlMat[is.finite(mdlMat)]), arr.ind = TRUE)[1, ]
    names(idx) <- c("ee", "aa")

    res <- list(loadings = xfaFit$xfa[[idx["aa"]]][[idx["ee"]]]$lambdaMat,
                resid.var = xfaFit$xfa[[idx["aa"]]][[idx["ee"]]]$sigma2,
                mdl = idx,
                grid = list(fit = xfaFit, wts = mdlMat)
                )

    attr(res, "max.factor") <- maxFactor
    attr(res, "log10.alpha.grid") <- xfaFit$agrid
    attr(res, "log10.eta.grid") <- xfaFit$egrid
    attr(res, "time") <- xfaFit$time[3] + endTime[3] - startTime[3]

    class(res) <- "xfa"

    res
}
