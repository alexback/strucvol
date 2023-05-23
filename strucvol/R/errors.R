#'@title MSSV Errors
#'@description Calculate numerical "sandwich" variance-covariance errors 
#'of a multivariate stochastic volatility model.
#'@param model Output "model" from function "fitmssv".
#'@param y The bivariate time series of log returns that the model has been fitted to.
#'@param x The explanatory variable for the model of the first column in y, presumably 
#'a log leverage multiplier.

mssverrors <- function(model, y, x){
  
  
  lg <- length(y[,1])
  pars <- model$pars
  p <- length(pars)
  step = 1e-7 * pars
  gr <- matrix(data = rep(0, lg * p), ncol = p, nrow = lg, byrow = T)
  
  for (i in 1:p){
    
    h = step[i]
    delta = rep(0,p)
    delta[i] = h
    
    spplus <- msspmaker(param = pars + delta, yt = y, xt = x)
    spminus <- msspmaker(param = pars - delta, yt = y, xt = x)
    # kf <- svkf(Tt = sp$Tt, Qt = sp$Qt, Zt = sp$Zt, Ht = sp$Ht, a0 = sp$a0,
    #  P0 = sp$P0, yt = sp$yt)
    
    kfplus <- fkf(a0 = spplus$a0, P0 = spplus$P0, dt = c(0,0,0,0,0), ct = c(0,0),
              Tt = spplus$Tt, Zt = spplus$Zt, HHt = spplus$Qt, GGt = spplus$Ht, 
              yt = t(spplus$yt))
    kfminus <- fkf(a0 = spminus$a0, P0 = spminus$P0, dt = c(0,0,0,0,0), ct = c(0,0),
                  Tt = spminus$Tt, Zt = spminus$Zt, HHt = spminus$Qt, GGt = spminus$Ht, 
                  yt = t(spminus$yt))
    
    
    
    logliksplus <- rep(0, lg)
                            
    logliksminus <- rep(0, lg)
    
    for (j in 1:lg){
      
    logliksplus[j] <- -0.5 * log(det(kfplus$Ft[,,j])) - 0.5 * t(kfplus$vt[,j]) %*% kfplus$Ftinv[,,j] %*% kfplus$vt[,j]  
    logliksminus[j] <- -0.5 * log(det(kfminus$Ft[,,j])) - 0.5 * t(kfminus$vt[,j]) %*% kfminus$Ftinv[,,j] %*% kfminus$vt[,j]
      
    }
    
    
    
    gr[,i] <- (logliksplus - logliksminus) / (2 * h) 
    
    
  }
  
  
  
  A <- t(gr) %*% gr
  
  
  B <- solve(model$hessian)
  
  seerrors <- sqrt(diag((B%*%A%*%B)))
  

  
  return(list(A = A, B = B, errors = seerrors))
  
}
