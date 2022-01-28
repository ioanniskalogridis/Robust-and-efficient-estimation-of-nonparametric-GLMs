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
# mean(mse.lo)*100 ; median(mse.lo)*100
hist(alpha.dpd)

###################################################### Bernoulli data ############################################################

nrep <- 10
n <- 200

mse.dpd <- rep(NA, nrep)
alpha.dpd <- rep(NA, nrep)
mse.dpd1 <- rep(NA, nrep)
mse.gam <- rep(NA, nrep)
mse.il <- rep(NA, nrep)
# mse.lo <- rep(NA, nrep)

for(j in 1:nrep){
  print(j)
  x <- seq(1/(n+1), n/(n+1), len = n)
  
  f1 <- -sin(5*x/1.2)/0.8-1
  # f1 <- 1.8*sin(3.4*x^2)
  
  inv.logit <- function(x) 1/(1+exp(-x))
  y <- rbinom(n, size = 1, prob = inv.logit(f1))
  # smpl <- sample(1:n, 0.05*n)
  # y[smpl] <- 1-y[smpl]
  
  fit <- dpd(x, y, family = "b")
  fit1 <- dpd(x, y, family = "b", alpha.cand = c(1))
  fit.gam <- gam(y~s(x, k = 60, bs = "cr"), family = binomial)
  # fit.rgam <- rgam(x,y, family = "binomial", ni = rep(1, n))
  fit.il <- DoubleRobGam(y ~ bsp(x, nknots = 36, order = 2, p = 3), family="binomial", selection = "RAIC")
 
  mean.t <- inv.logit(f1)
  # mean.rgam <- fit.rgam$fitted.values
  mean.dpd <- fit$fitted
  mean.dpd1 <- fit1$fitted
  mean.gam <- fit.gam$fitted.values
  mean.il <- fit.il$fitted.values
  
  # mse.rgam[j] <- mean(( mean.t - mean.rgam)^2 )
  mse.dpd[j] <- mean((mean.dpd-mean.t)^2)
  mse.dpd1[j] <- mean((mean.dpd1-mean.t)^2)
  alpha.dpd[j] <- fit$alpha
  mse.gam[j] <- mean((mean.gam-mean.t)^2)
  mse.il[j] <- mean((mean.il-mean.t)^2)
}

mean(mse.dpd)*100 ; median(mse.dpd)*100
mean(mse.dpd1)*100 ; median(mse.dpd1)*100
mean(mse.gam)*100 ; median(mse.dpd)*100
mean(mse.il)*100 ; median(mse.il)*100

# mean(mse.rgam, na.rm = TRUE)*100 ; median(mse.rgam, na.rm = TRUE)*100
#hist(alpha.dpd)

###################################################### Poisson data ############################################################

nrep <- 10
n <- 200

mse.dpd <- rep(NA, nrep)
alpha.dpd <- rep(NA, nrep)
mse.dpd1 <- rep(NA, nrep)
mse.gam <- rep(NA, nrep)
mse.il <- rep(NA, nrep)
# mse.lo <- rep(NA, nrep)

for(j in 1:nrep){
  print(j)
  x <- seq(1/(n+1), n/(n+1), len = n)
  
  f1 <- -sin(5*x/1.2)/0.8-1
  # f1 <- 1.8*sin(3.4*x^2)
  
  y <- rpois(n, lambda = exp(f1))
  # smpl <- sample(1:n, 0.1*n)
  # y[smpl] <- rpois(length(smpl), lambda = 3*exp(f1[smpl]))
  
  fit <- dpd(x, y, family = "p")
  fit1 <- dpd(x, y, family = "p", alpha.cand = c(1))
  fit.gam <- gam(y~s(x, k = 60, bs = "cr"), family = poisson)
  # fit.rgam <- rgam(x,y, family = "binomial", ni = rep(1, n))
  fit.il <- DoubleRobGam(y ~ bsp(x, nknots = 36, order = 2, p = 3), family="poisson", selection = "RAIC")
  
  mean.t <- exp(f1)
  # mean.rgam <- fit.rgam$fitted.values
  mean.dpd <- fit$fitted
  mean.dpd1 <- fit1$fitted
  mean.gam <- fit.gam$fitted.values
  mean.il <- fit.il$fitted.values
  
  # mse.rgam[j] <- mean(( mean.t - mean.rgam)^2 )
  mse.dpd[j] <- mean((mean.dpd-mean.t)^2)
  mse.dpd1[j] <- mean((mean.dpd1-mean.t)^2)
  alpha.dpd[j] <- fit$alpha
  mse.gam[j] <- mean((mean.gam-mean.t)^2)
  mse.il[j] <- mean((mean.il-mean.t)^2)
}

mean(mse.dpd)*100 ; median(mse.dpd)*100
mean(mse.dpd1)*100 ; median(mse.dpd1)*100
mean(mse.gam)*100 ; median(mse.dpd)*100
mean(mse.il)*100 ; median(mse.il)*100

# mean(mse.rgam, na.rm = TRUE)*100 ; median(mse.rgam, na.rm = TRUE)*100
hist(alpha.dpd)