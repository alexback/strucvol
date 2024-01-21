#'@title Fit a standard stochastic volatility model.
#'@param y a numeric vector or time series containing log returns.
#'@param start starting parameters for the optimization. 
#'@param N number of importance samples to draw for the Monte Carlo ll evaluation.
#'@import Rsolnp 
#'@return A list containing the output from the solver (model) and the outputs from 
#'the Kalman filter and monte carlo evaluation routine (fit).
#'@export
fitsv <- function(y, N = 5, start = c(0.95, 0.3)){
  mlmodel <- gosolnp(pars = start, fun = svloglik, LB = c(0,-100), UB = c(1,1),
                   yt = y, Ht = rep(pi^2 / 2, length(y)), fit = T, n.sim = 100)
  smlmodel <- solnp(pars = mlmodel$pars, fun = skloglik, LB = c(0,-100), UB = c(1,1),
                    yt = y, N = N, fit = T, control = list(tol = 1e-12))
  errors <- sqrt(diag(solve(smlmodel$hessian)))
  
  fit <- skloglik(param = smlmodel$pars, yt = y, N = N, fit = F)
  
  return(list(model = smlmodel, fit = fit, errors = errors))
  
}

#'@title Fit a structural stochastic volatility model.
#'@param y a numeric vector or time series containing log returns.
#'@param x an explanatory variable, presumably the log of a leverage multiplier.
#'@param N number of importance samples to draw for the Monte Carlo ll evaluation.
#'@param start starting parameters for the optimization.
#'@return A list containing the output from the solver ("model") and the outputs from 
#'the Kalman filter and monte carlo evaluation routine ("fit").
#'@export
fitssv <- function(y, x, N = 5, start = c(0.95, 0.3)){
  mlmodel <- gosolnp(pars = start, fun = ssvloglik, LB = c(0,0), UB = c(1,1),
                   yt = y, xt = x, Ht = rep(pi^2 / 2, length(y)), fit = T, n.sim = 100)
  
  smlmodel <- solnp(pars = mlmodel$pars, fun = sskloglik, LB = c(0,0), UB = c(1,1),
                    yt = y, xt = x, N = N, fit = T, control = list(tol = 1e-12))
  
  errors <- sqrt(diag(solve(smlmodel$hessian)))
  
  fit <- sskloglik(param = smlmodel$pars, yt = y, xt = x, N = N, fit = F)
  
  return(list(model = smlmodel, fit = fit, errors = errors))
}

#'@title Fit a bivariate structural stochastic volatility model.
#'@param y a bivariate numeric or time series containing log returns. The first column
#'should contain the "structural" series, while the second corresponds to the market.
#'@param x an explanatory variable, presumably the log of a leverage multiplier.
#'@param start starting parameters for the optimization.
#'@import FKF
#'@return A list containing the output from the solver ("model") and the outputs from 
#'the Kalman filter ("fit").
#'@export
fitmssv <- function(y, x, start = c(0.95, 0.95, 0.3, 0.3, 0.02)){

 model <- solnp(pars = start, fun = mssvloglik, LB = c(0,0,0,0,0), UB = c(1,1,1,1,1),
        yt = y, xt = x, fit = T)  
  
 errors <- mssverrors(model = model, y = y, x = x)
 fit <- mssvloglik(param = model$pars, yt = y, x = x, fit = F)
 
 return(list(model = model, fit = fit, errors = errors))
 
}
#'@title Fit an asymmetric stochastic volatility model (QML).
#'@param y a numeric or time series containing log returns.
#'@return Mainly the fitted parameters ("pars") and the standard errors (toterrors) 
#'@export
fitasysv <- function(y){

  
  sfit <- gosolnp(pars = c(0.9, 0.3, -0.3), fun = asysvloglik, LB = c(0,0,-1), UB = c(1,1,1),
                  yt = y, n.sim = 100, control = list(trace = F))
  
  finalsp <- asyspmaker(sfit$pars, yt = y)
  finalkf <- asysvkf(finalsp$Tt, finalsp$Qt, finalsp$Zt,
                     finalsp$Ht, finalsp$Gt, finalsp$ct,
                     finalsp$a0, finalsp$P0, finalsp$yt)
  
  finalpars <- c( finalkf$at[2,length(y)], sfit$pars)
  omegaerror <- sqrt(finalkf$Pt[2,2,length(y)])
  
  l <- length(y) - 1
  pars <- sfit$pars
  p <- length(pars)
  
  step = 1e-7* pars
  gr <- matrix(data = rep(0, l * p), ncol = p, nrow = l, byrow = T)
  
  for (i in 1:p){
    
    h = step[i]
    delta = rep(0,p)
    delta[i] = h
    
    logliksplus <- asysvloglikt(param = pars + delta, yt = y)
    
    logliksminus <- asysvloglikt(param = pars - delta, yt = y)
    
    
    
    gr[,i] <- (logliksplus - logliksminus) / (2 * h) 
    
    
  }
  
  
  
  A <- t(gr)%*%gr
  
  
  B <- solve(sfit$hessian)
  
  errors <- sqrt(diag(B))
  
  roberrors <- sqrt(diag((B%*%A%*%B)))
  
  
  
  return(list(pars = finalpars, toterror = c(omegaerror, roberrors), 
              sfit = sfit, A = A, B = B,
              errors = errors, roberrors = roberrors, kf = finalkf))
  
}


#'@title Fit an asymmetric stochastic volatility model with a covariate (QML).
#'@param y a numeric or time series containing log returns.
#'@param x an explanatory variable, presumably the log of a leverage multiplier.
#'@return Mainly the fitted parameters ("pars") and the standard errors (toterrors) 
#'@export
fitasysvx <- function(y, x){
  
  
  loglm2 <- log(x^2)
  
  
  sfit <- gosolnp(pars = c(0.9, 0.3, -0.3), fun = asysvloglikx, LB = c(0,0,-1), UB = c(1,1,1),
                  yt = y, x = loglm2, n.sim = 100, control = list(trace = F))
  
  finalsp <- asyspmakerx(sfit$pars, yt = y, x = loglm2)
  finalkf <- asysvkfx(finalsp$Tt, finalsp$Qt, finalsp$Zt,
                      finalsp$Ht, finalsp$Gt, finalsp$ct,
                      finalsp$a0, finalsp$P0, finalsp$yt)
  
  finalpars <- c( finalkf$at[2,length(y)], sfit$pars, finalkf$at[3,length(y)])
  omegaerror <- sqrt(finalkf$Pt[2,2,length(y)])
  xerror <- sqrt(finalkf$Pt[3,3,length(y)])
  
  l <- length(y) - 2
  pars <- sfit$pars
  p <- length(pars)
  
  step = 1e-7 * pars
  gr <- matrix(data = rep(0, l * p), ncol = p, nrow = l, byrow = T)
  
  for (i in 1:p){
    
    h = step[i]
    delta = rep(0,p)
    delta[i] = h
    
    logliksplus <- asysvlogliktx(param = pars + delta, yt = y,
                                 x = loglm2)
    
    logliksminus <- asysvlogliktx(param = pars - delta, yt = y,
                                  x = loglm2)
    
    
    
    gr[,i] <- (logliksplus - logliksminus) / (2 * h) 
    
    
  }
  
  
  
  A <- t(gr)%*%gr
  
  
  B <- solve(sfit$hessian)
  
  errors <- sqrt(diag(B))
  
  roberrors <- sqrt(diag((B%*%A%*%B)))
  
  
  
  return(list(pars = finalpars, toterror = c(omegaerror, roberrors, xerror), 
              sfit = sfit, loglm2 = loglm2, A = A, B = B,
              errors = errors, roberrors = roberrors, kf = finalkf))
  
}



