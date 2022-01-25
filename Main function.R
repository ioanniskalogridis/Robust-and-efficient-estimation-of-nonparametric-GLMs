###################################################################################################################################
############## Robust density power divergence estimators for GLMs (Claeskens, Kalogridis and Van Aelst (2022)) ###################
require(fda)
require(MASS)
require(RobStatTM)
require(robustbase)
require(GJRM)
require(remotes)
remotes::install_github("ilapros/DoubleRobGam")
require(DoubleRobGam)
require(maxLik)

dpd <- function(x, y, family = 'g', m = 2, toler = 1e-10, maxiter = 500, alpha.cand = seq(1e-02, 1, len = 10), nsteps = 25,
                sel = "AIC"){
  # Main function
  # x is the one-dimensional predictor
  # y is the response variable
  # family is the probability distribution of the response
  # valid values currently: "g" for Gaussian, "b" for Bernoulli, "p" for Poisson
  # m is the order of the penalty, by default m=2 corresponding to cubic splines
  # toler is the tolerance level of the updating algorithm
  # maxiter is the naximum number of allowed iterations
  # alpha.cand is a vector of candidate values for the tuning parameter (alpha in the notation of the paper)
  # nsteps denotes the number of allowed iterations for the selection of the tuning parameter
  # sele refers to the information criterion used in the determination of the penalty parameter
  # acceptable values are: AIC for the AIC of Hastie and Tibshirani (1990) and GIC for the criterion of Konichi and Kitagawa (2008)
  
  # The estimator is a smoothing-spline for n<50 and a penalized spline for n>=50
  
  
  x.t <- (x-min(x))/(max(x)-min(x))
  n.t <- length(unique(x.t))
  x = x.t
  
  n.knot.f <- function (n)
  {
    if (n < 50L)
      n
    else trunc({
      a1 <- log2(40)
      a2 <- log2(60)
      a3 <- log2(70)
      a4 <- log2(80)
      if (n < 200L) 2^(a1 + (a2 - a1) * (n - 50)/150) else if (n <
                                                               800L) 2^(a2 + (a3 - a2) * (n - 200)/600) else if (n <
                                                                                                                 3200L) 2^(a3 + (a4 - a3) * (n - 800)/2400) else 79 +
        (n - 3200)^(1/(2*m+1))
    })
  }
  n <- length(y)
  if(n.t <= 50){
    bsb <- create.bspline.basis(rangeval = c(min(x),max(x)), breaks = unique(x),
                                norder = 2*m )
    bsbe <-  eval.basis(bsb, x)
  } else if(n.t > 50){
    nknots <- n.knot.f(length(unique(x)))
    bsb <- create.bspline.basis(rangeval = c(min(x), max(x)), nbasis = (nknots + 2*m), norder = 2*m )
    bsbe <-  eval.basis(bsb, x)
  }
  
  if(m==1){
    P <- diff(diag(dim(bsbe)[2]), differences = 1)
    Pen.matrix <- t(P)%*%diag(1/(diff( seq(0, 1, len = (nbasis + 2*m)), differences = 1)))%*%P
  } else {
    Pen.matrix <- bsplinepen(bsb, Lfdobj = m) }
  
  if(family == "g"){
    rho.f <- function(x, tun ) -exp(-tun*x^2/2)/tun
    psi.f <- function(x, tun ) x*exp(-tun*x^2/2)
    psi.p.f <- function(x, tun ) exp(-tun*x^2/2) - x^2*tun*exp(-tun*x^2/2)
    weight.f <-function(x, tun ) exp(-tun*x^2/2)
    
    fit.in <-  DoubleRobGam(y ~ bsp(x, nknots = 35, order = m, p = 2*m-1), family="gaussian",
                            control = DoubleRobGamControl(tccM = 1.345), selection = "RAIC")
    resids.in <- as.vector(fit.in$residuals)
    
    scale = scaleM(diff(y), delta = 3/4, tuning.chi = sqrt(2)*0.70417 )
    
    irls <- function(lambda, alpha, maxit = maxiter, tol = toler){
      ic = 0
      istop = 0
      while(istop == 0 & ic <= maxit){
        ic = ic + 1
        weights.prelim <- weight.f(resids.in/scale, tun = alpha)
        M1 <- scale(t(bsbe), center = FALSE, scale = c(1/weights.prelim))%*%bsbe + 2*scale^2*lambda*Pen.matrix
        M2 <- t(bsbe)%*%(y*c(weights.prelim))
        v1 = SparseM::solve(M1, M2)
        resids1 <- as.vector(y-bsbe%*%v1)
        check = max( abs(resids1-resids.in)/scale) 
        if(check < tol){istop =1}
        resids.in <- resids1
      }
      weights1 = as.vector(weight.f(resids1/scale, tun = alpha) )
      resids <- c(y - bsbe%*%v1)
      return(list(beta = v1, weights = weights1, resids = resids, ic = ic, check = check ) )
    }
    if(sel == "AIC"){
      Pen.cr <- function(lambda, alpha){
        solv.irls <- irls(lambda, alpha)
        beta <- solv.irls$beta
        weights <- psi.p.f(solv.irls$resids/scale, tun = alpha)/scale^2
        NP <- scale(t(bsbe), center = FALSE, scale = c(1/weights))%*%bsbe
        hessian.m <-  NP + 2*lambda*Pen.matrix
        hat.tr <- sum(diag(solve(hessian.m)%*%NP))
        AIC <- sum(rho.f(solv.irls$resids/scale, tun = alpha)) + 2*hat.tr
        return(AIC)
      }
    } else if(sel == "GIC"){
      Pen.cr <- function(lambda, alpha){
        solv.irls <- irls(lambda, alpha)
        beta <- solv.irls$beta
        weights <- psi.p.f(solv.irls$resids/scale, tun = alpha)/scale^2
        NP <- scale(t(bsbe), center = FALSE, scale = c(1/weights))%*%bsbe
        hessian.m <-  NP + 2*lambda*Pen.matrix
        gr <- -psi.f(solv.irls$resids/scale, tun = alpha)/scale
        Q.m <-  scale(t(bsbe), center = FALSE, scale = 1/(gr*c(solv.irls$resids/scale^2)) )%*%bsbe
        hat.tr <- sum(diag(  -solve(hessian.m)%*%Q.m   ))
        GIC <- sum(rho.f(solv.irls$resids/scale, tun = alpha)) + 2*hat.tr
        return(GIC)
      }
    }
    lambda.cand <- c(1e-09, 3e-08, 8e-08, 3e-07, 8e-07, 3e-06, 8e-06,
                     3e-05, 8e-05, 3e-04, 8e-04, 3e-03, 8e-03, 3e-02, 8e-02, 3)
    B.m <-  bsplinepen(bsb, Lfdobj = 0)
    
    comp.alpha <-  function(alpha){
      lambda.e <- vapply(lambda.cand, FUN = Pen.cr, alpha = alpha, FUN.VALUE = numeric(1))
      wm <- which.min(lambda.e)
      if(wm== 1){wm <- 2}
      if(wm==length(lambda.cand)){ wm<- (length(lambda.cand)-1)}
      lambda.opt  <- optimize(f = Pen.cr, interval = c(lambda.cand[wm-1], lambda.cand[wm+1]), alpha)$minimum
      fit.opt <- irls(lambda.opt, maxit = maxiter, tol = toler, alpha)
      gr <- -psi.f(fit.opt$resids/scale, tun = alpha)/scale
      weights <- psi.p.f(fit.opt$resids/scale, tun = alpha)/scale^2
      hessian.m <- scale(t(bsbe), center = FALSE, scale = c(1/weights))%*%bsbe + 2*lambda.opt*Pen.matrix
      K.m <- scale(t(bsbe), center = FALSE, scale = c(1/gr^2))%*%bsbe
      msq1 <- sum( diag( solve(hessian.m, K.m%*%solve(hessian.m, B.m))   ) )
      return(list(msq1 = msq1,  beta.opt = fit.opt$beta, lambda.opt = lambda.opt))
    }
    compts <-  vapply(alpha.cand, FUN = comp.alpha, FUN.VALUE = rep(list(1), 3) )
    
    ic <- 0
    istop <- 0
    alpha.in <- 1
    while(ic <= nsteps & istop ==0){
      ic = ic + 1 
      msqs <- rep(NA, length(alpha.cand))
      for(j in 1:length(alpha.cand)){
        msqs[j] <- compts[[1, j]] + t(compts[[2, j]]-compts[[2, which(alpha.cand == alpha.in)]])%*%B.m%*%(compts[[2, j]]-compts[[2, which(alpha.cand == alpha.in)]])
      }
      alpha.up <- alpha.cand[which.min(msqs)]
      check <- abs(alpha.up-alpha.in)
      if(check< 1e-04){istop = 1}
      alpha.in <- alpha.up
    }
    
    alpha.opt <- alpha.up
    beta.opt.f <- compts[[2, which(alpha.cand == alpha.up)]]
    lambda.opt.f <- compts[[3, which(alpha.cand == alpha.up)]]
    
    est <- bsbe%*%beta.opt.f
    fitted <- est
    resids <- as.vector(y - est)
    p.resids <- as.vector(resids/mad(resids))
    a.resids <- p.resids
    
    ######################################################################################################################
    ######################################################################################################################
    
  } else if(family == "b") {
    inv.logit <- function(x) 1/(1+exp(-x))
    beta.il <- rep(0, dim(bsbe)[2])
    
    obj.f <- function(beta, pen, tun){
      lambda = pen
      alpha = tun
      t1 <- (1-inv.logit(bsbe%*%beta))^{1+alpha} + inv.logit(bsbe%*%beta)^{1+alpha}
      t2 <- (1+1/alpha)*dbinom(y, size = 1, prob = inv.logit(bsbe%*%beta))^{alpha}
      contr <- t1 - t2
      obj <- as.numeric(sum(contr)+ lambda*t(beta)%*%(Pen.matrix%*%beta))
      return(-obj)
    }
    
    gradient <- function(beta, pen, tun){
      lambda = pen
      alpha = tun
      gr1 <- (0-inv.logit(bsbe%*%beta))*(1-inv.logit(bsbe%*%beta))^{1+alpha}+ (1-inv.logit(bsbe%*%beta))*(inv.logit(bsbe%*%beta))^{1+alpha}
      gr2 <-  (y-inv.logit(bsbe%*%beta))*dbinom(y, size = 1 , prob = inv.logit( bsbe%*%beta ))^{alpha}
      gr.t <- gr1 - gr2
      gr <- (1+alpha)*diag(c(gr.t))%*%bsbe
      gradient <- as.vector(colSums(gr) + 2*lambda*(Pen.matrix%*%beta))
      return(-gradient)
    }
    
    hessian <-  function(beta, pen, tun){
      lambda = pen
      alpha = tun
      wt <- inv.logit(bsbe%*%beta)^2*(1-inv.logit(bsbe%*%beta))^{1+alpha} + (1-inv.logit(bsbe%*%beta ) )^2*(inv.logit(bsbe%*%beta))^{1+alpha}
      wt <- (1+alpha)*wt
      hs <- scale(t(bsbe), center = FALSE, scale = c(1/wt))%*%bsbe + 2*lambda*Pen.matrix
      return(-hs)
    }
    
    if(sel == "AIC"){
      newraph <- function(lambda, beta.in, maxit = maxiter, tol = toler, alpha){
        mnr <- maxNR(fn = obj.f, grad = gradient, hess = hessian, start = beta.in, pen =  lambda, tun =  alpha, 
                     control = list(tol = tol, iterlim = maxit))
        beta = mnr$estimate
        weights <-  inv.logit(bsbe%*%beta)^2*(1-inv.logit(bsbe%*%beta))^{1+alpha} + (1-inv.logit(bsbe%*%beta ) )^2*(inv.logit(bsbe%*%beta))^{1+alpha}
        weights <- c((1+alpha)*weights)
        gr1 <- (0-inv.logit(bsbe%*%beta))*(1-inv.logit(bsbe%*%beta))^{1+alpha}+ (1-inv.logit(bsbe%*%beta))*(inv.logit(bsbe%*%beta))^{1+alpha}
        gr2 <-  (y-inv.logit(bsbe%*%beta))*dbinom(y, size = 1 , prob = inv.logit( bsbe%*%beta ))^{alpha}
        gr.t <- gr1 - gr2
        gr <- c((1+alpha)*gr.t)
        hess.f <- mnr$hessian
        hat.tr <-  sum(diag(  solve(-hess.f, scale(t(bsbe), center = FALSE, scale = c(1/weights))%*%bsbe    )))
        return(list(beta = beta, hat.tr = hat.tr, gr = gr, hess.f = hess.f))
      }
    } else if(sel == "GIC"){
      newraph <- function(lambda, beta.in, maxit = maxiter, tol = toler, alpha){
        mnr <- maxNR(fn = obj.f, grad = gradient, hess = hessian, start = beta.in, pen =  lambda, tun =  alpha, 
                     control = list(tol = tol, iterlim = maxit))
        beta = mnr$estimate
        gr1 <- (0-inv.logit(bsbe%*%beta))*(1-inv.logit(bsbe%*%beta))^{1+alpha}+ (1-inv.logit(bsbe%*%beta))*(inv.logit(bsbe%*%beta))^{1+alpha}
        gr2 <-  (y-inv.logit(bsbe%*%beta))*dbinom(y, size = 1 , prob = inv.logit( bsbe%*%beta ))^{alpha}
        gr.t <- gr1 - gr2
        gr <- c((1+alpha)*gr.t)
        hess.f <- mnr$hessian
        Q.m <-  scale(t(bsbe), center = FALSE, scale = 1/(gr*c(y-inv.logit(bsbe%*%beta))) )%*%bsbe ##+ 2*lambda*Pen.matrix%*%(beta%*%t(colSums(c(y-inv.logit(bsbe%*%beta))*bsbe)))/n      #GIC
        hat.tr <- sum(diag(  solve(hess.f)%*%Q.m   ))
        return(list(beta = beta, hat.tr = hat.tr, gr = gr, hess.f = hess.f))
      }
    }
  
    obj.f.np <- function(beta, alpha){
      t1 <- (1-inv.logit(bsbe%*%beta))^{1+alpha} + inv.logit(bsbe%*%beta)^{1+alpha}
      t2 <- (1+1/alpha)*dbinom(y, size = 1, prob = inv.logit(bsbe%*%beta))^{alpha}
      contr <- t1 - t2
      obj <- as.numeric(sum(contr))
      return(obj)
    }
    
    Pen.cr <- function(lambda, beta.in, alpha){
      nr <- newraph(lambda, beta.in = beta.in, maxit = maxiter, alpha = alpha)
      beta <- nr$beta
      hat.tr <- nr$hat.tr
      IC <- 2*(obj.f.np(beta, alpha)) + 2*hat.tr 
      return(IC)
    }
    lambda.cand <- c(2e-06, 7e-06, 2e-05, 7e-05, 2e-04, 7e-04, 2e-03, 7e-03, 2e-02, 7e-02, 2e-01, 7e-01, 3)
    lambda.e.in <-  vapply(lambda.cand, FUN = Pen.cr, beta.in = beta.il, alpha = 1, FUN.VALUE = numeric(1))
    wm <- which.min(lambda.e.in)
    if(wm== 1){wm <- 2}
    if(wm==length(lambda.cand)){ wm<- (length(lambda.cand)-1)}
    lambda.opt.in  <- optimize(f = Pen.cr, interval = c(lambda.cand[wm-1], lambda.cand[wm+1]), beta.in = beta.il, alpha = 1)$minimum
    beta.opt.in <- newraph(lambda.opt.in, beta.in = beta.il, maxit = maxiter, alpha = 1)$beta
    B.m <-  bsplinepen(bsb, Lfdobj = 0)
    
    comp.alpha <-  function(alpha){
      lambda.e <- vapply(lambda.cand, FUN = Pen.cr, beta.in = beta.opt.in, alpha = alpha, FUN.VALUE = numeric(1))
      wm <- which.min(lambda.e)
      if(wm== 1){wm <- 2}
      if(wm==length(lambda.cand)){ wm<- (length(lambda.cand)-1)}
      lambda.opt  <- optimize(f = Pen.cr, interval = c(lambda.cand[wm-1], lambda.cand[wm+1]), alpha = alpha, beta.in = beta.opt.in)$minimum
      opt. <- newraph(lambda.opt, beta.in = beta.opt.in, maxit = maxiter, alpha = alpha)
      beta.opt  <- opt.$beta
      gr <- opt.$gr
      hessian.m <- -opt.$hess.f
      K.m <- scale(t(bsbe), center = FALSE, scale = c(1/gr^2))%*%bsbe
      msq1 <- sum( diag( solve(hessian.m, K.m%*%solve(hessian.m, B.m))   ) ) 
      return(list(msq1 = msq1,  beta.opt = beta.opt, lambda.opt = lambda.opt))
    }
    
    compts <-  vapply(alpha.cand, FUN = comp.alpha, FUN.VALUE = rep(list(1), 3) )
    
    ic <- 0
    istop <- 0
    alpha.in <- 1
    while(ic <= nsteps & istop ==0){
      ic = ic + 1 
      msqs <- rep(NA, length(alpha.cand))
      for(j in 1:length(alpha.cand)){
        msqs[j] <- compts[[1, j]] + t(compts[[2, j]]-compts[[2, which(alpha.cand == alpha.in)]])%*%B.m%*%(compts[[2, j]]-compts[[2, which(alpha.cand == alpha.in)]])
      }
      alpha.up <- alpha.cand[which.min(msqs)]
      check <- abs(alpha.up-alpha.in)
      if(check< 1e-04){istop = 1}
      alpha.in <- alpha.up
    }
    alpha.opt <- alpha.up
    beta.opt.f <- compts[[2, which(alpha.cand == alpha.up)]]
    lambda.opt.f <- compts[[3, which(alpha.cand == alpha.up)]]
    
    est <- bsbe%*%beta.opt.f
    fitted <- inv.logit(est)
    resids <- as.vector(y - inv.logit(est))
    p.resids <- as.vector(resids/sqrt(inv.logit(est)*(1-inv.logit(est))))
    
    ibeta <- function(x,a,b){ pbeta(x,a,b)*beta(a,b)}
    ansch.r <- function(y, mu){
      ansch.r <- (ibeta(y, 2/3, 2/3)-ibeta(mu, 2/3, 2/3))/(mu*(1-mu))^{1/6}
      return(ansch.r)
    }
    a.resids <- ansch.r(y, fitted)
    
    
    #####################################################################################################################
    #####################################################################################################################
    
  } else if(family =="p"){
    
    beta.il <- rep(0, dim(bsbe)[2])
    
    obj.f <- function(beta, pen, tun){
      lambda = pen
      alpha = tun
      means <- c(exp(bsbe%*%beta))
      s <- seq(0, min(3*max(y),round(max(means) + 4*(max(means))^{1/3}) ) )
      t1.f <- function(s) dpois(s, lambda = as.vector(exp(bsbe%*%beta)))^{1+alpha}
      t1 <- sapply(s, FUN = t1.f)
      t1 <- rowSums(t1)
      t2 <- (1+1/alpha)*dpois(y, lambda = exp(bsbe%*%beta))^{alpha}
      contr <- t1 - t2
      obj <- as.numeric(sum(contr)+ lambda*t(beta)%*%Pen.matrix%*%beta)
      return(-obj)
    }
    
    gradient <- function(beta, pen, tun){
      lambda = pen
      alpha = tun
      means <- c(exp(bsbe%*%beta))
      s <- seq(0, min(3*max(y),round(max(means) + 4*(max(means))^{1/3}) ) )
      gr1.f <-  function(s) (s-exp(bsbe%*%beta))*dpois(s, lambda = as.vector(exp(bsbe%*%beta)))^{1+alpha}
      gr1 <- sapply(s, FUN = gr1.f)
      gr1 <- rowSums(gr1)
      gr2 <- (y-exp(bsbe%*%beta))*dpois(y, lambda = exp( bsbe%*%beta ))^{alpha}
      gr.t <- gr1 - gr2
      gr <- (1+alpha)*t(scale(t(bsbe), center = FALSE, scale = c(1/gr.t) ))# *diag(c(gr.t))%*%bsbe
      gradient <- as.vector(colSums(gr) + 2*lambda*(Pen.matrix%*%beta))
      return(-gradient)
    }
    
    hessian <-  function(beta, pen, tun){
      lambda = pen
      alpha = tun
      means <- c(exp(bsbe%*%beta))
      s <- seq(0, min(3*max(y),round(max(means) + 4*(max(means))^{1/3}) ) )
      w.f <- function(s) (s-exp(bsbe%*%beta))^2*dpois(s, lambda = exp(bsbe%*%beta))^{1+alpha}
      wt <- sapply(s, FUN = w.f) 
      wt <- (1+alpha)*rowSums(wt)
      hs <- scale(t(bsbe), center = FALSE, scale = c(1/wt))%*%bsbe + 2*lambda*Pen.matrix
      return(-hs)
    }
    
    if(sel == "AIC"){
      newraph <- function(lambda, beta.in, maxit = maxiter, tol = toler, alpha){
        mnr <- maxNR(fn = obj.f, grad = gradient, hess = hessian, start = beta.in, pen =  lambda, tun =  alpha, 
                     control = list(tol = tol, iterlim = maxit))      
        beta = mnr$estimate
        means <- c(exp(bsbe%*%beta))
        s <- seq(0, min(3*max(y),round(max(means) + 4*(max(means))^{1/3}) ) )
        w.f <- function(s) (s-exp(bsbe%*%beta))^2*dpois(s, lambda = exp(bsbe%*%beta))^{1+alpha}
        wt <- sapply(s, FUN = w.f)
        weights <- (1+alpha)*rowSums(wt)
        gr1.f <-  function(s) (s-exp(bsbe%*%beta))*dpois(s, lambda = as.vector(exp(bsbe%*%beta)))^{1+alpha}
        gr1 <- sapply(s, FUN = gr1.f)
        gr1 <- rowSums(gr1)
        gr2 <- (y-exp(bsbe%*%beta))*dpois(y, lambda = exp( bsbe%*%beta ))^{alpha}
        gr.t <- gr1 - gr2
        gr <- as.vector((1+alpha)*gr.t)
        hess.f <- mnr$hessian
        hat.tr <- sum(diag(  solve(-hess.f)%*%scale(t(bsbe), center = FALSE, scale = c(1/weights))%*%bsbe    ))
        return(list(beta = beta, hat.tr = hat.tr, weights = weights, gr = gr, hess.f = hess.f))
      }
    } else if(sel == "GIC"){
      newraph <- function(lambda, beta.in, maxit = maxiter, tol = toler, alpha){
        mnr <- maxNR(fn = obj.f, grad = gradient, hess = hessian, start = beta.in, pen =  lambda, tun =  alpha, 
                     control = list(tol = tol, iterlim = maxit))      
        beta = mnr$estimate
        means <- c(exp(bsbe%*%beta))
        s <- seq(0, min(3*max(y),round(max(means) + 4*(max(means))^{1/3}) ) )
        gr1.f <-  function(s) (s-exp(bsbe%*%beta))*dpois(s, lambda = as.vector(exp(bsbe%*%beta)))^{1+alpha}
        gr1 <- sapply(s, FUN = gr1.f)
        gr1 <- rowSums(gr1)
        gr2 <- (y-exp(bsbe%*%beta))*dpois(y, lambda = exp( bsbe%*%beta ))^{alpha}
        gr.t <- gr1 - gr2
        gr <- as.vector((1+alpha)*gr.t)
        hess.f <- mnr$hessian
        Q.m <-  scale(t(bsbe), center = FALSE, scale = 1/(gr*c(y-exp(bsbe%*%beta))) )%*%bsbe  # + 2*lambda*Pen.matrix%*%(beta%*%t(colSums(c(y-exp(bsbe%*%beta))*bsbe)))/n      #GIC
        hat.tr <- sum(diag(  solve(hess.f, Q.m   )))  #GIC
        return(list(beta = beta, hat.tr = hat.tr, gr = gr, hess.f = hess.f))
      }
    }
    
    obj.f.np <- function(beta, alpha){
      means <- c(exp(bsbe%*%beta))
      s <- seq(0, min(3*max(y),round(max(means) + 4*(max(means))^{1/3}) ) )
      t1.f <- function(s) dpois(s, lambda = as.vector(exp(bsbe%*%beta)))^{1+alpha}
      t1 <- sapply(s, FUN = t1.f)
      t1 <- rowSums(t1)
      t2 <- (1+1/alpha)*dpois(y, lambda = exp(bsbe%*%beta))^{alpha}
      contr <- t1 - t2
      obj <- as.numeric(sum(contr))
      return(obj)
    }
    
    Pen.cr <- function(lambda, beta.in, alpha){
      nr <- newraph(lambda, beta.in = beta.in, maxit = maxiter, alpha = alpha)
      beta <- nr$beta
      hat.tr <- nr$hat.tr
      IC <- 2*(obj.f.np(beta, alpha)) + 2*hat.tr
      return(IC)
    }
    lambda.cand <- c(2e-06, 7e-06, 2e-05, 7e-05, 2e-04, 7e-04, 2e-03, 7e-03, 2e-02, 7e-02, 2e-01, 7e-01, 3)
    lambda.e.in <-  vapply(lambda.cand, FUN = AIC.f, beta.in = beta.il, alpha = 1, FUN.VALUE = numeric(1))
    
    wm <- which.min(lambda.e.in)
    if(wm== 1){wm <- 2}
    if(wm==length(lambda.cand)){ wm<- (length(lambda.cand)-1)}
    lambda.opt.in  <- optimize(f = Pen.cr, interval = c(lambda.cand[wm-1], lambda.cand[wm+1]), beta.in = beta.il, alpha = 1)$minimum
    beta.opt.in <- newraph(lambda.opt.in, beta.in = beta.il, maxit = maxiter, alpha = 1)$beta
    
    B.m <-  bsplinepen(bsb, Lfdobj = 0)
    
    comp.alpha <-  function(alpha){
      lambda.e <- vapply(lambda.cand, FUN = AIC.f, beta.in = beta.opt.in, alpha = alpha, FUN.VALUE = numeric(1))
      wm <- which.min(lambda.e)
      if(wm== 1){wm <- 2}
      if(wm==length(lambda.cand)){ wm<- (length(lambda.cand)-1)}
      lambda.opt  <- optimize(f = AIC.f, interval = c(lambda.cand[wm-1], lambda.cand[wm+1]), alpha = alpha, beta.in = beta.opt.in)$minimum
      opt. <- newraph(lambda.opt, beta.in = beta.opt.in, maxit = maxiter, alpha = alpha)
      beta.opt  <- opt.$beta
      gr <- opt.$gr
      hessian.m <- -opt.$hess.f
      K.m <- scale(t(bsbe), center = FALSE, scale = c(1/gr^2))%*%bsbe
      msq1 <- sum( diag( solve(hessian.m, K.m%*%solve(hessian.m, B.m))   ) ) 
      return(list(msq1 = msq1,  beta.opt = beta.opt, lambda.opt = lambda.opt))
    }
    
    compts <-  vapply(alpha.cand, FUN = comp.alpha, FUN.VALUE = rep(list(1), 3) )
    
    ic <- 0
    istop <- 0
    alpha.in <- 1
    while(ic <= nsteps & istop ==0){
      ic = ic + 1 
      msqs <- rep(NA, length(alpha.cand))
      for(j in 1:length(alpha.cand)){
        msqs[j] <- compts[[1, j]] + t(compts[[2, j]]-compts[[2, which(alpha.cand == alpha.in)]])%*%B.m%*%(compts[[2, j]]-compts[[2, which(alpha.cand == alpha.in)]])
      }
      alpha.up <- alpha.cand[which.min(msqs)]
      check <- abs(alpha.up-alpha.in)
      if(check< 1e-04){istop = 1}
      alpha.in <- alpha.up
    }
    alpha.opt <- alpha.up
    beta.opt.f <- compts[[2, which(alpha.cand == alpha.up)]]
    lambda.opt.f <- compts[[3, which(alpha.cand == alpha.up)]]
    
    est <- bsbe%*%beta.opt.f
    fitted <- exp(est)
    resids <- as.vector(y - exp(est))
    p.resids <- as.vector(resids/sqrt(exp(est)))
    ansch.r <- function(y, mu){
      ansch.r <- 3*(y^{2/3}-mu^{2/3})/(2*mu^{1/6})
      return(ansch.r)
    }
    a.resids <- ansch.r(y, fitted)
    
  }
  return(list(est = est, resids = resids, fitted = fitted, a.resids = a.resids,
              nsteps = ic, p.resids = p.resids, lambda = lambda.opt.f, alpha = alpha.opt, family = family))
}
