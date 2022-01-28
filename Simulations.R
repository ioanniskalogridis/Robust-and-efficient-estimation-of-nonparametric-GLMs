require(EnvStats)
require(mgcv)
remotes::install_github("ilapros/DoubleRobGam")
require(DoubleRobGam)

nrep <- 10
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
  y <- f1 + rnormMix(n, 0, 1, 0, 9, p.mix = 0.05)
  fit <- dpd(x, y, family = "g")
  alpha.dpe[i] <- fit$alpha
  fit1 <- dpd(x, y, fam = "g", alpha.cand = c(1))
  fit.gam <- mgcv::gam(y~s(x, k = 60, bs = "cr"), method = "GCV.Cp")
  fit.il <- DoubleRobGam(y ~ bsp(x, nknots = 36, order = 2, p = 3), family="gaussian", selection = "RAIC",
                         scale = scaleM(diff(y), delta = 3/4, tuning.chi = sqrt(2)*0.70417 ))
  fit.lo <- loess_wrapper_extrapolate(x, y)

  mse.dpd[i] <- mean( (f1-fit$est )^2)
  mse.gam[i] <- mean((f1-predict.gam(fit.gam))^2)
  # mse.dpe1[i] <- mean((f1 - fit1$est)^2)
  # mse.gam[i] <- mean((f1 - fit.gam$est)^2)
  # mse.il[i] <- mean( (f1-fit.il$fitted.values )^2)
  # mse.lo[i] <- mean((f1 - predict(fit.lo))^2)
  # mse.aeb[i] <- mean( (f1-fit.aeb$fit$eta1 )^2)
}
