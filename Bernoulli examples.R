require(SemiPar)
require(mgcv)
require(ggplot2)

data(bpd)
bpd <- bpd[order(bpd$birthweight),]
x <- bpd$birthweight
y <- bpd$BPD
plot(x,y)
fit <- dpd(x, y, family = "b")
fit.gam <- mgcv::gam(y~s(x, k = 40), family = binomial())

lines(x, fitbpd$fitted, lwd = 3, col = "blue", type = "l")
lines(x, fit.gambpd$fitted.values, lwd = 3, col = "red")
fitbpd$alpha


