require(ggplot2)
require(mgcv)

# Import data
zf <- read.table(file = file.choose(), header = TRUE, sep = ";")
# select the file "LOSdata.csv"
Pop  = data.frame(LOS=zf$LOS,Age=zf$Age,Sex=zf$Sexe,DRG=zf$DRG,MDC=zf$MDC,Nbdg=zf$NbDiag,Nbtt=zf$NbAct,CW=zf$CWeff)
table(Pop$MDC) # Number of observations per DRG

########################################### DRG: Diseases and disorders of the eye #################################################
pop = Pop[Pop$MDC=="02",]
pop <- pop[order(pop$Age),]
x <- pop$Age
y <- pop$LOS-1
plot(x, y, cex = 1.2, pch = 19)
pr <- dpd(x, y, family = "p")
pr.gam <- gam(y~s(x, bs = "cr", k = 40), family = poisson)
lines(x, pr$fitted, lwd = 3, col = "blue")
lines(x, pr.gam$fitted.values, lwd = 3, col = "red", lty = 2)
pr$alpha

## Plot in the log scake
outliers <- ifelse(abs(pr$a.resids) > 2.6, 1, 0)
sum(outliers)
library(ggplot2)
df <- data.frame(x = x, y = log(y+1), pr =  log(pr$fitted+1), fit.gam = log(pr.gam$fitted.values+1), out = as.vector(outliers ))
p1 <- ggplot(df, aes(x = x, y = y, colour = as.factor(out), size = 4)) + geom_point( aes(shape = as.factor(out)))  + labs(x = "Age", y = "log(LOS+1)")
p1 <- p1 + theme_bw(base_size = 35) + guides(color= "none", shape = "none", size = "none")
p1 <- p1 + geom_line(aes(x = x, y = pr), data = df, colour = "blue", linetype = 1, size = 1.4) + geom_line(aes(x = x,y = fit.gam), colour = "red", linetype = 5, size = 1.4)
p1 <- p1 + scale_color_manual(values=c("gray", "purple"))
p1

########################################### DRG: Diseases and disorders of the male reproductive system ############################
pop = Pop[Pop$MDC=="12",]
pop <- pop[order(pop$Age),]
x <- pop$Age
y <- pop$LOS-1
plot(x, y, cex = 1.2, pch = 19)
pr <- dpd(x, y, family = "p")
pr.gam <- gam(y~s(x, bs = "cr", k = 40), family = poisson)
lines(x, pr$fitted, lwd = 3, col = "blue")
lines(x, pr.gam$fitted.values, lwd = 3, col = "red", lty = 2)
pr$alpha

## Plot in the log scake
outliers <- ifelse(abs(pr$a.resids) > 2.6, 1, 0)
sum(outliers)
library(ggplot2)
df <- data.frame(x = x, y = log(y+1), pr =  log(pr$fitted+1), fit.gam = log(pr.gam$fitted.values+1), out = as.vector(outliers ))
p2 <- ggplot(df, aes(x = x, y = y, colour = as.factor(out), size = 4)) + geom_point( aes(shape = as.factor(out)))  + labs(x = "Age", y = "log(LOS+1)")
p2 <- p2 + theme_bw(base_size = 35) + guides(color= "none", shape = "none", size = "none")
p2 <- p2 + geom_line(aes(x = x, y = pr), data = df, colour = "blue", linetype = 1, size = 1.4) + geom_line(aes(x = x,y = fit.gam), colour = "red", linetype = 5, size = 1.4)
p2 <- p2 + scale_color_manual(values=c("gray", "purple"))
p2

