require(EnvStats)
require(mgcv)
remotes::install_github("ilapros/DoubleRobGam")
require(DoubleRobGam)

#################################################### Gaussian data ################################################################


nrep <- 20
n <- 200

mse.dpd <- rep(NA, nrep)
alpha.dpd <- rep(NA, nrep)
mse.dpd1 <- rep(NA, nrep)
mse.gam <- rep(NA, nrep)
mse.il <- rep(NA, nrep)
# mse.lo <- rep(NA, nrep)

for(i in 1:nrep){
  print(i)
  x <- seq(1/(n+1), n/(n+1), len = n)
  
  # f1 <- -sin(5*x/1.2)/0.8-1
  f1 <- 1.8*sin(3.4*x^2)
  
  #y <- f1 + rnormMix(n, 0, 1, 0, 9, p.mix = 0)
  y <- f1 + rnormMix(n, 0, 1, 0, 9, p.mix = 0.05)
  #y <- f1 + rnormMix(n, 0, 1, 0, 9, p.mix = 0.1)
  
  fit <- dpd(x, y, family = "g")
  alpha.dpd[i] <- fit$alpha
  fit1 <- dpd(x, y, family = "g", alpha.cand = c(1))
  fit.gam <- mgcv::gam(y~s(x, k = 60, bs = "cr"), method = "GCV.Cp")
  fit.il <- DoubleRobGam(y ~ bsp(x, nknots = 36, order = 2, p = 3), family="gaussian", selection = "RAIC",
                         scale = scaleM(diff(y), delta = 3/4, tuning.chi = sqrt(2)*0.70417 ))
  # fit.lo <- loess_wrapper_extrapolate(x, y)

  mse.dpd[i] <- mean( (f1-fit$fitted )^2)
  mse.gam[i] <- mean((f1-predict.gam(fit.gam))^2)
  mse.dpd1[i] <- mean((f1 - fit1$fitted)^2)
  mse.il[i] <- mean( (f1-fit.il$fitted.values )^2)
  # mse.lo[i] <- mean((f1 - predict(fit.lo))^2)
}

mean(mse.dpd)*100 ; median(mse.dpd)*100
mean(mse.dpd1)*100 ; median(mse.dpd1)*100
mean(mse.gam)*100 ; median(mse.dpd)*100
mean(mse.il)*100 ; median(mse.il)*100
hist(alpha.dpd)


nrep <- 1000
n <- 200
mse.dpe <- rep(NA, nrep)
mse.dpe1 <- rep(NA, nrep)
mse.gam <- rep(NA, nrep)
alpha.dpe <- rep(NA, nrep)
mse.rgam <- rep(NA, nrep)
mse.il <- rep(NA, nrep)

for(j in 616:nrep){
  print(j)
  x <- seq(1/(n+1), n/(n+1), len = n)
  # x <- runif(n)
  # x <- sort(x)
  y <- rep(NA, n)
  # f1 <- -sin(5*x/6)/0.8-1
  # f1 <- cos(2*pi*x)
  # f1 <- 3*atan(10*(x-0.5))
  f1 <- -sin(5*x/1.2)/0.8-1
  # f3 <- 4*cos(2*pi*(1-x)^2)
  # y <- rep(NA, n)
  # for(i in 1:n){
  #   # y[i] <- rbinom(n = 1, size = 1, prob = inv.logit( cos(2*pi*x[i])))
  #   # y[i] <- rbinom(n = 1, size = 1, prob = inv.logit( 1.3*x[i]^2))
  #   y[i] <- rbinom(n = 1, size = 1, prob = inv.logit(f1[i]) )
  #   # y[i] <- rbinom(n = 1, size = 1, prob = inv.logit( 4*cos(2*pi*(1-x[i])^2) ))
  # }
  y <- rbinom(n, size = 1, prob = inv.logit(f1))
  # set.seed(1)
  # smpl <- sample(1:n, 0.05*n)
  # y[smpl] <- 1-y[smpl]
  fit <- dpe(x, y, family = "b")
  # fit1 <- dpe(x, y, family = "b", alpha.cand = c(1))
  # fit.gam <- gam(y~s(x, k = 60), family = binomial)
  # fit.rgam <- rgam(x,y, family = "binomial", ni = rep(1, n))
  # fit.il <- DoubleRobGam(y ~ bsp(x, nknots = 36, order = 2, p = 3), family="binomial", selection = "RAIC")
  # curve(3*atan(10*(x-0.5)), lwd = 3, col = "black")
  # curve(4*cos(2*pi*x), lwd = 3, col = "black")
  # lines(x, fit$est, lwd = 3, col = "blue")
  # lines(x, predict.gam(fit.gam), lwd = 3, col = "red")
  
  mean.t <- inv.logit(f1)
  # mean.rgam <- fit.rgam$fitted.values
  mean.dpe <- inv.logit(fit$est)
  # mean.gam <- inv.logit(predict.gam(fit.gam))
  # mean.t <- inv.logit( 3*atan(10*(x-0.5)))
  # mean.dpe1 <- inv.logit(fit1$est)
  # mean.gam <- inv.logit(predict.gam(fit.gam))
  # mean.il <- fit.il$fitted.values
  
  # mse.rgam[j] <- mean(( mean.t - mean.rgam)^2 )
  mse.dpe[j] <- mean((mean.dpe-mean.t)^2)
  # mse.dpe1[j] <- mean((mean.dpe1-mean.t)^2)
  alpha.dpe[j] <- fit$alpha
  # mse.gam[j] <- mean((mean.gam-mean.t)^2)
  # mse.il[j] <- mean((mean.il-mean.t)^2)
  
  # plot(x, y)
  # lines(x, inv.logit(f2), type = "l", lwd = 3)
  # lines(x, inv.logit(fit$est), lwd = 3, col = "blue")
  # lines(x, inv.logit(predict.gam(fit.gam)), lwd = 3, col = "red")
}
# mean(mse.rgam, na.rm = TRUE)*100 ; median(mse.rgam, na.rm = TRUE)*100
mean(mse.dpe, na.rm = TRUE)*100 ; median(mse.dpe, na.rm = TRUE)*100
hist(alpha.dpe)
mean(mse.dpe1, na.rm = TRUE)*100 ; median(mse.dpe1, na.rm = TRUE)*100
mean(mse.gam, na.rm = TRUE)*100 ; median(mse.gam, na.rm = TRUE)*100
mean(mse.il, na.rm = TRUE)*100 ; median(mse.il, na.rm = TRUE)*100
mean(mse.rgam, na.rm = TRUE)*100 ; median(mse.rgam, na.rm = TRUE)*100
