#'@title Fit a standard stochastic volatility model.
#'@param y a numeric vector or time series containing log returns.
#'@param start starting parameters for the optimization. 
#'@param N number of importance samples to draw for the Monte Carlo ll evaluation.
#'@import Rsolnp 
#'@return A list containing the output from the solver (model) and the outputs from 
#'the Kalman filter and monte carlo evaluation routine (fit).
#'@export
fitsv <- function(y, N = 5, start = c(0.95, 0.3)){
  mlmodel <- solnp(pars = start, fun = svloglik, LB = c(0,0), UB = c(1,1),
                   yt = y, Ht = rep(pi^2 / 2, length(y)), fit = T)
  smlmodel <- solnp(pars = mlmodel$pars, fun = skloglik, LB = c(0,0), UB = c(1,1),
                    yt = y, fit = T, control = list(tol = 1e-8))
  
  fit <- skloglik(param = smlmodel$pars, yt = y, N = N, fit = F)
  
  return(list(model = mlmodel, fit = fit))
  
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
  mlmodel <- solnp(pars = start, fun = ssvloglik, LB = c(0,0), UB = c(1,1),
                   yt = y, xt = x, Ht = rep(pi^2 / 2, length(y)), fit = T)
  
  smlmodel <- solnp(pars = mlmodel$pars, fun = sskloglik, LB = c(0.001,0.001), UB = c(0.99,0.99),
                    yt = y, xt = x, N = N, fit = T, control = list(tol = 1e-8))
  
  fit <- sskloglik(param = smlmodel$pars, yt = y, xt = x, N = N, fit = F)
  
  return(list(model = smlmodel, fit = fit))
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
  
 fit <- mssvloglik(param <- model$pars, yt = y, x = x, fit = F)
 
 return(list(model = model, fit = fit))
 
}