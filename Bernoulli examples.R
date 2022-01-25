require(SemiPar)
require(mgcv)
require(ggplot2)
require(mlbench)

############################################# Birthweight and dysplasia #############################################################

data(bpd)
bpd <- bpd[order(bpd$birthweight),]
x <- bpd$birthweight
y <- bpd$BPD
plot(x,y)
fit <- dpd(x, y, family = "b")
fit.gam <- mgcv::gam(y~s(x, k = 40), family = binomial)

lines(x, fit$fitted, lwd = 3, col = "blue", type = "l")
lines(x, fit.gam$fitted.values, lwd = 3, col = "red")
fit$alpha

outliers <- ifelse(abs(fit$a.resids)>2.6, 1, 0) # Anscombe residuals
sum(outliers)

df <- data.frame(x = x, y = y, fit.dpd =  fit$fitted, fit.gam = fit.gam$fitted.values, out = outliers )
p <- ggplot(data = df, aes(x = x, y = y, colour = as.factor(out), shape = as.factor(out) )) + labs(x = "Birthweight", y = "Probability of dysplasia")
p <- p + geom_jitter(data = df, aes(y = y), height = 0.04, size = 4) + theme_bw(base_size = 35) + guides(color = "none") + guides(shape = "none")
p <- p + geom_line(data = df, aes(y = fit.dpd, x = x, group = 1), colour = "blue", linetype = 1, size = 1.4) + geom_line(aes(x = x, y = fit.gam, group = 1), colour = "red", linetype = 2, size = 1.4)
p <- p + scale_color_manual(values=c("gray", "purple"))
p

############################################ Pima Indian Diabetes ###################################################################

data(PimaIndiansDiabetes)
pima <- PimaIndiansDiabetes
pima$diabetes <- ifelse(pima$diabetes == "pos", 1, 0)
pima <- pima[order(pima$glucose),]
x <- pima$glucose
y <- pima$diabetes
plot(x, y)
fit <- dpd(x, y, family = "b")
fit.gam <- mgcv::gam(y~s(x, k = 40), family = binomial)

lines(x, fit$fitted, lwd = 3, col = "blue", type = "l")
lines(x, fit.gam$fitted.values, lwd = 3, col = "red")
fit$alpha

hist(fit$a.resids, xlab = "Anscombe residuals") # Anscombe residuals 

outliers <- ifelse(abs(fit$a.resids)>2.6, 1, 0) 
sum(outliers)
df <- data.frame(x = x, y = y, fit.dpd =  fit$fitted, fit.gam = fit.gam$fitted.values, out = outliers )
p <- ggplot(df, aes(x = x, y = y, colour = as.factor(out), shape = as.factor(out) ))  + labs(x = "Plasma glucose", y = "Probability of Diabetes")
p <- p + geom_jitter(data = df, aes(y = y), height = 0.04, size = 4) + theme_bw(base_size = 35)+ guides(color = "none") + guides(shape = "none")
p <- p + geom_line(aes(x = x, y = fit.dpd, group = 1), colour = "blue", linetype = 1, size = 1.4) + geom_line(aes(x = x, y = fit.gam, group = 1), colour = "red", linetype = 2, size = 1.4)
p <- p + scale_color_manual(values=c("gray", "purple"))
p
