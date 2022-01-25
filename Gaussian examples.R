require(SemiPar)
require(quantreg)
require(MASS)


########################################################## LIDAR ###############################################################
data(lidar)
y <- lidar$logratio
x <- lidar$range
plot(x, y, cex = 1.2, pch = 19) ; grid()
fit.AIC <- dpd(x, y, family = "g", sel = "AIC")
fit.GIC <- dpd(x, y, family = "g", sel = "GIC")
lines(x, fit.AIC$fitted, lwd = 4, col = "blue")
lines(x, fit.GIC$fitted, lwd = 4, col = "darkgray",)
legend("bottomleft", legend = c("DPD with AIC", "DPD with GIC"), lwd = c(4,4), col = c("blue", "darkgray"), cex = 1.4)
fit.AIC$alpha

######################################################### Motorcycle impact ####################################################

y <- mcycle$accel
x <- mcycle$times
plot(x, y, cex = 1.2, pch = 19) ; grid()
fit.mcycle <- dpd(x, y, family = "g")
lines(x, fit.mcycle$fitted, lwd = 4, col = "blue")
fit.mcycle$alpha

######################################################## Age and income ########################################################

data("age.income")
age.income <- age.income[order(age.income$age),]
x <- age.income$age
y <- age.income$log.income
plot(x, y, cex = 1.2, pch = 19, xlab = "Income", ylab = "Age") ; grid()
fit.ai <- dpd(x, y, family= "g")
fit.ai$alpha
lines(x, fit.ai$fitted, lwd = 4, col = "blue")

######################################################## Mammals weight and speed ###############################################

data("Mammals")
Mammals <- Mammals[order(Mammals$weight), ]
Mammals$weight <- log(Mammals$weight)
Mammals$speed <- log(Mammals$speed)
x <- Mammals$weight
y <- Mammals$speed
plot(x, y, cex = 1.2, pch = 19, xlab = "Weight (log)", ylab = "Speed (log)")
fit.mamAIC <- dpd(x, y, sel = "AIC")
fit.mamGIC <- dpd(x, y, sel = "GIC")
fit.gam <- gam(y~s(x, k = 40, bs = "cr"))
lines(x, fit.mamAIC$fitted, lwd = 4, col = "blue")
lines(x, fit.mamGIC$fitted, lwd = 4, col = "darkgray")
lines(x, fit.gam$fitted.values, lwd = 4, col = "red", lty = 2)
legend("bottomright", legend = c("DPD with AIC", "DPD with GIC", "GAM"), lwd = c(4,4, 4), lty = c(1, 1, 2),
       col = c("blue", "darkgray", "red"), cex = 1.4)
fit.mamGIC$alpha


#hist(fit.mamGIC$p.resids, col = "gray") 
# Numerous outliers

