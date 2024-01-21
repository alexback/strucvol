########################## Test functions ######################################

lgarch_loglik <- function(param, y, init, fit = T){
  t <- length(y)
  ## Initialize parameters from vector
  
  alpha0 <- param[1]
  beta1 <- param[2]
  alpha1 <- param[3]
  
  logsigma2 <- rep(0, t)
  
  
  logsigma2[1] <- init
  
  logsigma2 <- lcrec(y = y, logsigma2 = logsigma2,
                     alpha0 = alpha0, alpha1 = alpha1, beta1 = beta1)
  
  sigma2 <- exp(logsigma2)
  
  
  if(fit == T){
    
    return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
  }
  
  else{
    return(list(lls = dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE), sigma2 = sigma2, logsigma2 = logsigma2))
  
  }
}


asylgarch_loglik <- function(param, y, init, fit = T){
  t <- length(y)
  ## Initialize parameters from vector
  
  alpha0 <- param[1]
  beta1 <- param[2]
  alpha1 <- param[3]
  alpha2 <- param[4]
  
  logsigma2 <- rep(0, t)
  
  
  logsigma2[1] <- init
  
  ind <- as.numeric(y<0) * log(y^2)
  
  logsigma2 <- asylcrec(y = y, ind = ind, logsigma2 = logsigma2,
                     alpha0 = alpha0, alpha1 = alpha1,
                     alpha2 = alpha2, beta1 = beta1)
  
  sigma2 <- exp(logsigma2)
  
  
  if(fit == T){
    
    return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
  }
  
  else{
    return(list(lls = dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE), sigma2 = sigma2, logsigma2 = logsigma2))
    
  }
}


################ Convert between stochastic volatility and lgarch ##############

svol2lgarch <- function(sigma0, beta1, sigma1){
  
  ceta <- digamma(1/2) - log(1/2)
  theta <- (1 / (2 * beta1)) * (1 + beta1^2 + 2 * sigma1^2 / pi^2 - sqrt((1 - beta1^2 + (2 * sigma1^2) / pi^2)^2 + (8 * beta1^2 * sigma1^2) / pi^2))
  
  sum <- sum(log(gamma((theta - beta1) * theta^(1:2000) + 0.5)) - log(gamma(1/2)))
  
  czeta <- exp(sigma1^2 / (4 * (1 - theta^2)) - ceta/2 - digamma(1/2) * (theta - beta1) * theta / (2 * (1 - theta)) +
                 0.5 * sum)
  
  alpha0 <- sigma0 + (1 - beta1) * ceta + (1 - theta) * log(czeta^2)
  alpha1 <- beta1 - theta
  lbeta1 <- theta
  
  return(c(alpha0, lbeta1, alpha1))
  
  
}

lgarchderiv <- function(pars, y){
  
  
  
  l <- length(y)
  p <- length(pars)
  step = 1e-12 * pars
  
  pars2 <- svol2lgarch(pars[1], pars[2], pars[3])
  gr <- matrix(data = rep(0, l * p), ncol = p, nrow = l, byrow = T)
  
  for (i in 1:p){
    
    h = step[i]
    delta = rep(0,p)
    delta[i] = h
    
    plus <- lgarch_loglik(param = pars2 + delta, y = y, init = mean(log(y^2)), fit = F)$logsigma2
    
    minus <- lgarch_loglik(param = pars2 - delta, y = y, init = mean(log(y^2)), fit = F)$logsigma2
    
    
    
    gr[,i] <- (plus - minus) / (2 * h) 
    
    
  }
  
  return(gr)
  
}

asylgarchderiv <- function(pars, y){
  
  
  
  l <- length(y)

  
  pars2 <- c(svol2lgarch(pars[1], pars[2], pars[3]),0)
  
  p <- length(pars2)
  step = 1e-12 * pars2
  step[4] <- step[3]
  gr <- matrix(data = rep(0, l * p), ncol = p, nrow = l, byrow = T)
  
  for (i in 1:p){
    
    h = step[i]
    delta = rep(0,p)
    delta[i] = h
    
    plus <- asylgarch_loglik(param = pars2 + delta, y = y, init = mean(log(y^2)), fit = F)$logsigma2
    
    minus <- asylgarch_loglik(param = pars2 - delta, y = y, init = mean(log(y^2)), fit = F)$logsigma2
    
    
    
    gr[,i] <- (plus - minus) / (2 * h) 
    
    
  }
  
  return(gr)
  
}


#'@title Test for a leverage effect in the data.
#'@param data The data to be tested.
#'@param model The null model.
#'@import Rcpp
#'@return A test statistic with an asymptotic chi^2(1) distribution under the null model.
#'@export
llevtest <- function(data, model)
{
  
  
  svfit <- svloglik2(param = model$pars, yt = data, Ht = rep(pi^2 / 2, length(data)), fit = F)
  
  
  svpars <- c((1 - model$pars[1]) * svfit$kf$at[2,length(data)], model$pars)
  
  
  
  
  lgarchpars <- svol2lgarch(sigma0 = svpars[1], beta1 = svpars[2], sigma1 = svpars[3])
  
  
  fit <- lgarch_loglik(param = lgarchpars, y = data, init = mean(log(data^2)), fit = F)
  
  sigma2 <- fit$sigma2
  lsigma2 <- fit$logsigma2
  
  leta2 <- log(data^2)
  
  l <- length(data)
  
  signind <- as.numeric(data < 0)

  
  chi <- data / sqrt(sigma2)
  
  chisquaredx <- chi^2 - 1
  
  # Find the derivatives recursively
  
  betahat <- lgarchpars[2]
  
  
  # Starts at time t = 2 gradients.
  derivdf <- derivreclev(leta2 = leta2, logsigma2 = lsigma2,
                        beta1 = rep(betahat, length(leta2)))
  #derivdf <- asylgarchderiv(svpars, data)
  
  gr1 <- derivdf[,1] 
  gr2 <- derivdf[,2]
  gr3 <- derivdf[,3]

  reg1 <- lm(chisquaredx[-1] ~ gr1[-length(gr1)] + gr2[-length(gr2)]  + gr3[-length(gr3)] + 0)
   #reg1 <- lm(chisquaredx[-1] ~ gr1[-1] + gr2[-1]  + gr3[-1] + 0)
  chisquaredx <- residuals(reg1) 
  
  ## Linear model
  
  # We get rid of the first observation in chisquaredx to be able to lag the gradient.
  #  This is in effect the time t gradient, which is a function of time t - 1 data.
  # To make the vectors have equal lengths, we remove the last observation of the gradient. 
  wmodel1 <- lm(signind[-length(signind)] ~ gr1[-length(gr1)] + gr2[-length(gr2)] + gr3[-length(gr3)] + 0)
  ### W vector components
  
  
  w1 <- as.vector(residuals(wmodel1))
  
  ###################################
  
  # Here, chisquaredx starts at t = 2 because the time t = 2 gradient is a function of time t = 1 data. 
  
  w1x <- chisquaredx * w1


  
  ###################################
  
  ones <- rep(1,length(w1x))
  
  robustauxmodel <- lm(ones ~ w1x + 0)
  
  # Obtain robust stat
  
  LMR <- (length(w1x)) * summary(robustauxmodel)$r.squared
  #t <- summary(robustauxmodel)$coefficients[3]
  return(LMR)
  
}

#'@title Test for a leverage effect in the data based on the log garch representation.
#'@param data The data to be tested.
#'@param model The null model.
#'@import Rcpp
#'@return A test statistic with an asymptotic chi^2(1) distribution under the null model.
#'@export
llevtest2 <- function(data)
{
  
  
  model <- solnp(pars = c(0, 0.8, 0.1), fun = lgarch_loglik, y = data,
                 init = mean(log(data^2)),
                 LB = c(-100, 0,0), UB = c(100,1,1))

  init <- mean(log(data^2))
  print(model$pars)
  
  fit <- lgarch_loglik(param = model$pars, y = data, init = init,
                       fit = F)

  sigma2 <- fit$sigma2
  lsigma2 <- fit$logsigma2
  
  leta2 <- log(data^2)
  
  l <- length(data)
  
  signind <- as.numeric(data < 0)
  sizeind <- signind * log(data^2)
  
  ones <- rep(1,l - 1)
  
  chi <- data / sqrt(sigma2)
  
  chisquaredx <- chi^2 - 1
  
  # Find the derivatives recursively
  
  betahat <- model$pars[2]
  
  
  # Starts at time t = 2 gradients.
  derivdf <- derivreclev(leta2 = leta2, logsigma2 = lsigma2,
                         beta1 = rep(betahat, length(leta2)))
  
  gr1 <- derivdf[,1] 
  gr2 <- derivdf[,2]
  gr3 <- derivdf[,3] 
  
  reg1 <- lm(chisquaredx[-1] ~ gr1[-length(gr1)] + gr2[-length(gr2)]  + gr3[-length(gr3)] + 0)
  
  chisquaredx <- residuals(reg1)
  
  ## Linear model
  
  # We get rid of the first observation in chisquaredx to be able to lag the gradient.
  #  This is in effect the time t gradient, which is a function of time t - 1 data.
  # To make the vectors have equal lengths, we remove the last observation of the gradient. 
  wmodel1 <- lm(signind[-length(signind)] ~ gr1[-length(gr1)] + gr2[-length(gr2)] + gr3[-length(gr3)] + 0)
  
  wmodel2 <- lm(sizeind[-length(sizeind)] ~ gr1[-length(gr1)] + gr2[-length(gr2)] + gr3[-length(gr3)] + 0)
  ### W vector components
  
  
  w1 <- as.vector(residuals(wmodel1))
  w2 <- as.vector(residuals(wmodel2))
  
  ###################################
  
  # Here, chisquaredx starts at t = 2 because the time t = 2 gradient is a function of time t = 1 data. 
  
  w1x <- chisquaredx * w1
  w2x <- chisquaredx * w2
  
  ###################################
  
  robustauxmodel <- lm(ones ~ w1x + w2x + 0)
  
  # Obtain robust stat
  
  LMR <- (length(w1x)) * summary(robustauxmodel)$r.squared
  
  
  return(LMR)
  
}




################################################################################
#'@title Test for misspecification in the form of an excluded leverage multiplier.
#'@param data The data to be tested.
#'@param model The null model.
#'@return A test statistic with an asymptotic chi^2(1) distribution under the null model.
#'@export 
levmulttest <- function(data, lmt, model)
{
  
  
  svfit <- svloglik(param = model$pars, yt = data, Ht = rep(pi^2 / 2, length(data)), fit = F)
  
  svpars <- c((1 - model$pars[1]) * svfit$kf$at[2,length(data)], model$pars)
  
  # Convert the fitted stochastic volatility pars to their log-garch counterparts
  
  lgarchpars <- svol2lgarch(sigma0 = svpars[1], beta1 = svpars[2], sigma1 = svpars[3])
  
  
  fit <- lgarch_loglik(param = lgarchpars, y = data, init = mean(log(data^2)), fit = F)
  
  sigma2 <- fit$sigma2
  lsigma2 <- fit$logsigma2
  
  leta2 <- log(data^2)
  
  l <- length(data)
  
  signind <- as.numeric(data < 0)
  
  ones <- rep(1,l - 1)
  
  chi <- data / sqrt(sigma2)
  
  chisquaredx <- chi^2 - 1
  
  # Find the derivatives recursively
  
  betahat <- lgarchpars[2]
  
  derivdf <- derivreclm(leta2 = leta2, logsigma2 = lsigma2,
                        beta1 = rep(betahat, length(leta2)), lm2 = c(lmt[1], lmt[-length(lmt)]))
  
  gr1 <- derivdf[,1] 
  gr2 <- derivdf[,2]
  gr3 <- derivdf[,3] 
  gr4 <- lmt
  
  reg1 <- lm(chisquaredx[-1] ~ gr1[-length(gr1)] + gr2[-length(gr2)]  + gr3[-length(gr3)] + 0)
  
  chisquaredx <- residuals(reg1)
  
  ## Linear model
  
  # We get rid of the first observation in chisquaredx to be able to lag the gradient.
  #  This is in effect the time t gradient, which is a function of time t - 1 data.
  # To make the vectors have equal lengths, we remove the last observation of the gradient. 
  
  wmodel1 <- lm(gr4[-length(gr4)] ~ gr1[-length(gr1)] + gr2[-length(gr2)] + gr3[-length(gr3)] + 0)
  
  ### W vector components
  
  
  w1 <- as.vector(residuals(wmodel1))
  
  ###################################
  
  # Here, chisquaredx starts at t = 2 because the time t = 2 gradient is a function of time t = 1 data. 
  
  w1x <- chisquaredx * w1
  
  
  ###################################
  
  robustauxmodel <- lm(ones ~ w1x + 0)
  
  
  # Obtain robust stat
  
  LMR <- (length(w1x)) * summary(robustauxmodel)$r.squared
  
  
  return(LMR)
  
}
################################################################################
#'@title Likelihood ratio test for two competing stochastic volatility models.
#'@param model0 The null model
#'@param model1 The alternative model
#'@return The likelihood ratio test statistic, here presumably asymptotically
#'distributed as chi^2(1) under the null.
#'@export
llratiotest <- function(model0 , model1){
  
  ll0 <- model0$values[length(model0$values)]
  ll1 <- model1$values[length(model1$values)]
  
  return(2 * (ll0 - ll1))
}
