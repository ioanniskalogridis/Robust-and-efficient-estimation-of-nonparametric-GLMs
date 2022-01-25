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

outliers <- ifelse(abs(fit$a.resids)>2.6, 1, 0) # Anscombe residuals
sum(outliers)

df <- data.frame(xb = as.numeric(x), yb = y, fit.dpd =  c(fit$fitted), fit.gam = fit.gam$fitted.values, out = outliers )
p <- ggplot(data = df, aes(x = xb, y = yb, colour = as.factor(out), shape = as.factor(out) )) + labs(x = "Birthweight", y = "Probability of dysplasia")
p <- p + geom_jitter(data = df, aes(y = yb), height = 0.04, size = 4) + theme_bw(base_size = 35) + guides(color = "none") + guides(shape = "none")
p <- p + geom_line(data = df, aes(y = fit.dpd, x = xb, group = 1), colour = "blue", linetype = 1, size = 1.4) # + geom_line(aes(x = x, y = fit.gam), colour = "red", linetype = 2, size = 1.4)
p <- p + scale_color_manual(values=c("gray", "purple"))
p

