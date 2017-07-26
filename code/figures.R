rm(list = ls())

## load the saved files from 'ngp' (n > p) and 'nlp' (n <= p) analyses
rnks.ngp <- readRDS( "../result/ngp_rnks.rds")
ms.ngp <- readRDS( "../result/ngp_mse.rds")
ts.ngp <- readRDS( "../result/ngp_tpr.rds")
fs.ngp <- readRDS( "../result/ngp_fdr.rds")

rnks.nlp <- readRDS( "../result/nlp_rnks.rds")
ms.nlp <- readRDS( "../result/nlp_mse.rds")
ts.nlp <- readRDS( "../result/nlp_tpr.rds")
fs.nlp <- readRDS( "../result/nlp_fdr.rds")

res.rnks <- rbind.data.frame(rnks.nlp, rnks.ngp)
res.ms <- rbind.data.frame(ms.nlp, ms.ngp)
res.ts <- rbind.data.frame(ts.nlp, ts.ngp)
res.fs <- rbind.data.frame(fs.nlp, fs.ngp)

nsamp <- c(50, 100, 250, 500, 5000)
ndims <- c(50, 100, 250, 500, 2000)

res.ts$p[res.ts$p == 50] <- 1
res.ts$p[res.ts$p == 100] <- 2
res.ts$p[res.ts$p == 250] <- 3
res.ts$p[res.ts$p == 500] <- 4
res.ts$p[res.ts$p == 2000] <- 5

colors <- rep("black", 6)

pdf("../result/img/tpr.pdf", width = 45, height = 20)
par(mfrow = c(2, 3))
par(cex = 1)
par(mar = c(7, 0, 3, 0), oma = c(1, 12, 0, 0.4))
par(tcl = -0.1)
par(mgp = c(2, 0.6, 0))
for (nn in 1:5) {
    dat <- res.ts[res.ts$n == nsamp[nn], ]
    xx <- dat$p
    yy <- dat$mn
    yy1 <- pmax(dat$mn - dat$sd, 0)
    yy2 <- pmin(dat$mn + dat$sd, 1)
    plot(xx[dat$mtd == "fanc"], yy[dat$mtd == "fanc"], pch = 15, cex = 6, ylim = c(0, 1.1), xlim = c(0.5, 5.5), col = colors[1],
         ylab = NA, xlab = NA, axes = FALSE)
    yylab <- axTicks(2)
    points(xx[dat$mtd == "spc"] - 0.2, yy[dat$mtd == "spc"], pch = 18, cex = 7, col = colors[2])
    points(xx[dat$mtd == "xfa"] - 0.1, yy[dat$mtd == "xfa"], pch = 17, cex = 6, col = colors[3])
    points(xx[dat$mtd == "rock"] + 0.1, yy[dat$mtd == "rock"], pch = 21, cex = 6, col = colors[4], lwd = 10)
    points(xx[dat$mtd == "vari"] + 0.2, yy[dat$mtd == "vari"], pch = 19, cex = 6, col = colors[5])
    arrows(xx[dat$mtd == "fanc"], yy1[dat$mtd == "fanc"], xx[dat$mtd == "fanc"], yy2[dat$mtd == "fanc"], col = colors[1],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    arrows(xx[dat$mtd == "spc"] - 0.2, yy1[dat$mtd == "spc"], xx[dat$mtd == "spc"] - 0.2, yy2[dat$mtd == "spc"], col = colors[2],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    arrows(xx[dat$mtd == "xfa"] - 0.1, yy1[dat$mtd == "xfa"], xx[dat$mtd == "xfa"] - 0.1, yy2[dat$mtd == "xfa"], col = colors[3],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    arrows(xx[dat$mtd == "rock"] + 0.1, yy1[dat$mtd == "rock"], xx[dat$mtd == "rock"] + 0.1, yy2[dat$mtd == "rock"], col = colors[4],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    arrows(xx[dat$mtd == "vari"] + 0.2, yy1[dat$mtd == "vari"], xx[dat$mtd == "vari"] + 0.2, yy2[dat$mtd == "vari"], col = colors[5],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    mtext("p", at = 3, side = 1, line = 6, cex = 5)
    if (nn == 1 | nn == 4) {
        axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 5, las = 2)
        mtext("True positive rate", side = 2, line = 8, cex = 5)
    }
    axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
    mtext(c(50, 100, 250, 500, 2000), at = 1:5, side = 1, line = 2.5, cex = 5)
    box(col = "grey40", lwd = 5)
    mtext(bquote(n == .(nsamp[nn])), side = 3, line = -3.5, adj = 0.1, cex = 5)
}
dev.off()

res.fs$p[res.fs$p == 50] <- 1
res.fs$p[res.fs$p == 100] <- 2
res.fs$p[res.fs$p == 250] <- 3
res.fs$p[res.fs$p == 500] <- 4
res.fs$p[res.fs$p == 2000] <- 5
pdf("../result/img/fdr.pdf", width = 45, height = 20)
par(mfrow = c(2, 3))
par(cex = 1)
par(mar = c(7, 0, 3, 0), oma = c(1, 12, 0, 0.4))
par(tcl = -0.1)
par(mgp = c(2, 0.6, 0))
for (nn in 1:5) {
    dat <- res.fs[res.fs$n == nsamp[nn], ]
    xx <- dat$p
    yy <- dat$mn
    yy1 <- pmax(dat$mn - dat$sd, 0)
    yy2 <- pmin(dat$mn + dat$sd, 1)
    plot(xx[dat$mtd == "fanc"], yy[dat$mtd == "fanc"], pch = 15, cex = 6, ylim = c(-0.1, 1.1), xlim = c(0.5, 5.5), col = colors[1],
         ylab = NA, xlab = NA, axes = FALSE)
    yylab <- axTicks(2)
    points(xx[dat$mtd == "spc"] - 0.2, yy[dat$mtd == "spc"], pch = 18, cex = 7, col = colors[2])
    points(xx[dat$mtd == "xfa"] - 0.1, yy[dat$mtd == "xfa"], pch = 17, cex = 6, col = colors[3])
    points(xx[dat$mtd == "rock"] + 0.1, yy[dat$mtd == "rock"], pch = 21, cex = 6, col = colors[4], lwd = 10)
    points(xx[dat$mtd == "vari"] + 0.2, yy[dat$mtd == "vari"], pch = 19, cex = 6, col = colors[5])
    arrows(xx[dat$mtd == "fanc"], yy1[dat$mtd == "fanc"], xx[dat$mtd == "fanc"], yy2[dat$mtd == "fanc"], col = colors[1],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    arrows(xx[dat$mtd == "spc"] - 0.2, yy1[dat$mtd == "spc"], xx[dat$mtd == "spc"] - 0.2, yy2[dat$mtd == "spc"], col = colors[2],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    arrows(xx[dat$mtd == "xfa"] - 0.1, yy1[dat$mtd == "xfa"], xx[dat$mtd == "xfa"] - 0.1, yy2[dat$mtd == "xfa"], col = colors[3],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    arrows(xx[dat$mtd == "rock"] + 0.1, yy1[dat$mtd == "rock"], xx[dat$mtd == "rock"] + 0.1, yy2[dat$mtd == "rock"], col = colors[4],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    arrows(xx[dat$mtd == "vari"] + 0.2, yy1[dat$mtd == "vari"], xx[dat$mtd == "vari"] + 0.2, yy2[dat$mtd == "vari"], col = colors[5],
           angle = 90, length = 0.1, code = 3, lwd = 6)
    mtext("p", at = 3, side = 1, line = 6, cex = 5)
    if (nn == 1 | nn == 4) {
        axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 5, las = 2)
        mtext("False discovery rate", side = 2, line = 8, cex = 5)
    }
    axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
    mtext(c(50, 100, 250, 500, 2000), at = 1:5, side = 1, line = 2.5, cex = 5)
    box(col = "grey40", lwd = 5)
    mtext(bquote(n == .(nsamp[nn])), side = 3, line = -3.5, adj = 0.1, cex = 5)
}
dev.off()

res.ms$p[res.ms$p == 50] <- 1
res.ms$p[res.ms$p == 100] <- 2
res.ms$p[res.ms$p == 250] <- 3
res.ms$p[res.ms$p == 500] <- 4
res.ms$p[res.ms$p == 2000] <- 5
pdf("../result/img/mse.pdf", width = 45, height = 20)
par(mfrow = c(2, 3))
par(cex = 1)
par(mar = c(7, 0, 3, 0), oma = c(1, 10, 0, 0.4))
par(tcl = -0.1)
par(mgp = c(2, 0.6, 0))
for (nn in 1:5) {
    dat <- res.ms[res.ms$n == nsamp[nn], ]
    xx <- dat$p
    yy <- dat$mn
    yy1 <- pmax(dat$mn - dat$sd, 0)
    yy2 <- pmin(dat$mn + dat$sd)
    plot(xx[dat$mtd == "fanc"], yy[dat$mtd == "fanc"], pch = 15, cex = 6, ylim = c(0, 5), xlim = c(0.5, 5.5),
         col = colors[1], ylab = NA, axes = FALSE, xlab = NA)
    yylab <- axTicks(2)
    points(xx[dat$mtd == "spc"] - 0.2, yy[dat$mtd == "spc"], pch = 18, cex = 7, col = colors[2])
    points(xx[dat$mtd == "xfa"] - 0.1, yy[dat$mtd == "xfa"], pch = 17, cex = 6, col = colors[3])
    points(xx[dat$mtd == "rock"] + 0.1, yy[dat$mtd == "rock"], pch = 21, cex = 6, col = colors[4], lwd = 10)
    points(xx[dat$mtd == "vari"] + 0.2, yy[dat$mtd == "vari"], pch = 19, cex = 6, col = colors[5])
    arrows(xx[dat$mtd == "fanc"], yy1[dat$mtd == "fanc"], xx[dat$mtd == "fanc"], yy2[dat$mtd == "fanc"], col = colors[1],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    arrows(xx[dat$mtd == "spc"] - 0.2, yy1[dat$mtd == "spc"], xx[dat$mtd == "spc"] - 0.2, yy2[dat$mtd == "spc"], col = colors[2],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    arrows(xx[dat$mtd == "xfa"] - 0.1, yy1[dat$mtd == "xfa"], xx[dat$mtd == "xfa"] - 0.1, yy2[dat$mtd == "xfa"], col = colors[3],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    arrows(xx[dat$mtd == "rock"] + 0.1, yy1[dat$mtd == "rock"], xx[dat$mtd == "rock"] + 0.1, yy2[dat$mtd == "rock"], col = colors[4],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    arrows(xx[dat$mtd == "vari"] + 0.2, yy1[dat$mtd == "vari"], xx[dat$mtd == "vari"] + 0.2, yy2[dat$mtd == "vari"], col = colors[5],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    mtext("p", at = 3, side = 1, line = 6, cex = 5)
    if (nn == 1 | nn == 4) {
        axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 5, las = 2)
        mtext(expression("Root mean square error"), side = 2, line = 4.5, cex = 5)
    }
    axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
    mtext(c(50, 100, 250, 500, 2000), at = 1:5, side = 1, line = 2.5, cex = 5)
    box(col = "grey40", lwd = 5)
    mtext(bquote(n == .(nsamp[nn])), side = 3, line = -3.5, adj = 0.1, cex = 5)
}
dev.off()

res.rnks$p[res.rnks$p == 50] <- 1
res.rnks$p[res.rnks$p == 100] <- 2
res.rnks$p[res.rnks$p == 250] <- 3
res.rnks$p[res.rnks$p == 500] <- 4
res.rnks$p[res.rnks$p == 2000] <- 5
pdf("../result/img/rnks.pdf", width = 45, height = 20)
par(mfrow = c(2, 3))
par(cex = 1)
par(mar = c(7, 0, 3, 0), oma = c(1, 10, 0, 0.4))
par(tcl = -0.1)
par(mgp = c(2, 0.6, 0))
for (nn in 1:5) {
    dat <- res.rnks[res.rnks$n == nsamp[nn], ]
    xx <- dat$p
    yy <- dat$mn
    yy1 <- pmax((dat$mn - dat$sd), 0)
    yy2 <- pmin((dat$mn + dat$sd), 20)
    plot(xx[dat$mtd == "fanc"], yy[dat$mtd == "fanc"], pch = 15, cex = 6, ylim = c(-0.1, 21.5), xlim = c(0.5, 5.5),
         col = colors[1], ylab = NA, axes = FALSE, xlab = NA)
    abline(h = 5, col = "black", lwd = 10, lty = "dashed")
    yylab <- axTicks(2)
    points(xx[dat$mtd == "spc"] - 0.2, yy[dat$mtd == "spc"], pch = 18, cex = 7, col = colors[2])
    points(xx[dat$mtd == "xfa"] - 0.1, yy[dat$mtd == "xfa"], pch = 17, cex = 6, col = colors[3])
    points(xx[dat$mtd == "rock"] + 0.1, yy[dat$mtd == "rock"], pch = 21, cex = 6, col = colors[4], lwd = 10)
    points(xx[dat$mtd == "vari"] + 0.2, yy[dat$mtd == "vari"], pch = 19, cex = 6, col = colors[5])
    points(xx[dat$mtd == "bridge"], yy[dat$mtd == "bridge"], pch = 4, cex = 6, col = "black", lwd = 10)
    arrows(xx[dat$mtd == "fanc"], yy1[dat$mtd == "fanc"], xx[dat$mtd == "fanc"], yy2[dat$mtd == "fanc"], col = colors[1],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    arrows(xx[dat$mtd == "spc"] - 0.2, yy1[dat$mtd == "spc"], xx[dat$mtd == "spc"] - 0.2, yy2[dat$mtd == "spc"], col = colors[2],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    arrows(xx[dat$mtd == "xfa"] - 0.1, yy1[dat$mtd == "xfa"], xx[dat$mtd == "xfa"] - 0.1, yy2[dat$mtd == "xfa"], col = colors[3],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    arrows(xx[dat$mtd == "rock"] + 0.1, yy1[dat$mtd == "rock"], xx[dat$mtd == "rock"] + 0.1, yy2[dat$mtd == "rock"], col = colors[4],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    arrows(xx[dat$mtd == "vari"] + 0.2, yy1[dat$mtd == "vari"], xx[dat$mtd == "vari"] + 0.2, yy2[dat$mtd == "vari"], col = colors[5],
           angle = 90, length = 0.1, code = 3, lwd = 10)
    arrows(xx[dat$mtd == "bridge"], yy1[dat$mtd == "bridge"], xx[dat$mtd == "bridge"], yy2[dat$mtd == "bridge"], col = "black",
           angle = 90, length = 0.1, code = 3, lwd = 10)
    mtext("p", at = 3, side = 1, line = 6, cex = 5)
    if (nn == 1 | nn == 4) {
        axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 5, las = 2)
        mtext("Rank", side = 2, line = 6, cex = 5)
    }
    axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
    mtext(c(50, 100, 250, 500, 2000), at = 1:5, side = 1, line = 2.5, cex = 5)
    box(col = "grey40", lwd = 5)
    mtext(bquote(n == .(nsamp[nn])), side = 3, line = -3.5, adj = 0.1, cex = 5)
}
dev.off()
