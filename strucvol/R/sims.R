#'@title simstrucsystem
#'@description Simulate a bivariate structural volatility model.
#'@importFrom MASS mvrnorm 
#'@param len Integer, the length of the simulated time period.
#'@param pars A vector with the parameter values of the structural multivariate stochastic volatility model.
#'@param Ain Numeric, initial asset value.
#'@param Ein Numeric, initial equity value.
#'@param K Numeric, debt value over the time period.
#'@param r Numeric, the risk-free rate of interest
#'@param uv Numeric, the daily unconditional variance of the assets over the period.
#'@param ttv Numeric, time to maturity of the debt over the period (measured in years).
#'@param crisis_sim Boolean, should the function return the leverage ratio conditional on a crisis?
#'@param crisisret Numeric, return that defines the upper limit of a crisis  during the period.
#'@return A list containing relevant simulated quantities.
#'@export

simstrucsystem <- function(len = 30, pars = c(-10, 0.95, 0.3, -12, 0.9, 0.2, 0.6, 1),
                           Ain = 100, Ein = 80, K = 20, r = 0.001, uv = 0.0005, ttv = 5,
                           crisis_sim = F, crisisret = -0.1)
  {
   
    mu1 <- pars[1]
    beta1 <- pars[2]
    sigma1 <- pars[3]
    
    mu2 <- pars[4]
    beta2 <- pars[5]
    sigma2 <- pars[6]
  
    rho <- pars[7]
    phi <- pars[8]

    len <- as.integer(len)
    
    h_1 <- rep_len(0, length.out = len)
    h_2 <- rep_len(0, length.out = len)
    
    h0_1 <- rnorm(1, mean = mu1, sd = sigma1 / sqrt(1 - beta1^2))
    h0_2 <- rnorm(1, mean = mu2, sd = sigma2 / sqrt(1 - beta2^2))
    
    
    eta1 <- rnorm(len)
    eta2 <- rnorm(len)
    
    
    sigmaeps <- matrix(c(1, rho, rho, 1), byrow = T, ncol = 2, nrow = 2) 
    eps <- mvrnorm(n = len, mu = c(0,0), Sigma = sigmaeps)
    
    eps1 <- eps[,1]
    eps2 <- eps[,2]
    
    
    h_1[1] <- mu1 + beta1 * (h0_1 - mu1) + sigma1 * eta1[1]  # Initialize
    h_2[1] <- mu2 + beta2 * (h0_2 - mu2) + sigma2 * eta2[1]  # Initialize
    
    for (i in 2:length(h_1)){
    h_1[i] <- mu1 + beta1 * (h_1[i-1] - mu1) + sigma1 * eta1[i]
    h_2[i] <- mu2 + beta2 * (h_2[i-1] - mu2) + sigma2 * eta2[i]
    
    }
    
    ra <- exp(h_1 / 2) * eps1
    rm <- exp(h_2 / 2) * eps2
    
    
    A <- cumprod(c(Ain, exp(ra)))

    AK <- A / K

    d1 <- (log(AK) + (r + 0.5 * sqrt(uv * 252)^2 ) * ttv) / (sqrt(uv * 252 * ttv))
    
    
    E <- rep(0, len)
    E[1] <- Ein
    lm <- rep(0, len-1)
    
    re <- rep(0, len)
    re[1] <- ra[1] * pnorm(d1[1]) * A[1] / E[1] 
    
    for (j in 2:len){
    lm[j -1] <- pnorm(d1[j - 1]) * A[j - 1] / E[j - 1]
    re[j] <- (lm[j - 1]^phi) * ra[j]
    
    E[j] <- E[j - 1] * exp(re[j])
    }
    
    EK = E / K
    
    levrat <- 1 - (E / (E + K))
   
      
    cind <- ifelse(sum(rm) < crisisret, T, F)
        
    if (crisis_sim == T){
      
      return(ifelse(cind, levrat[len], 999))
    }
    
    else{
    return(list(retseries = cbind(ra, re, rm), h = cbind(h_1, h_2),
                lm = lm, E = E, A = A, EK = EK, levrat = levrat, AK = AK, crisis = cind))
    }
}

#'@title simdevent
#'@description Simulate the minimum leverage ratio during some period given inputs.
#'@importFrom MASS mvrnorm
#'@param N Integer, the number of time series to use for the simulation. 
#'@param len Integer, the length of the simulated time period.
#'@param pars A vector with the parameter values of the structural multivariate stochastic volatility model.
#' Should be in the order: c(mu_1, beta_1, sigma_1, mu_2, beta_2, sigma_2, rho, phi),
#' where subscript 1 corresponds to the "structural" series that has an explanatory variable in the state equation.
#'@param Ain Numeric, initial asset value.
#'@param Ein Numeric, initial equity value.
#'@param K Numeric, debt value over the time period.
#'@param r Numeric, the risk-free rate of interest
#'@param uv Numeric, the daily unconditional variance of the assets over the period.
#'@param ttv Numeric, time to maturity of the debt over the period (measured in years).
#'@param thd Numeric, specifies a threshold leverage ratio that the user wants to monitor. Used to draw a vertical abline in the histogram.
#'@param plot Boolean, should a plot be created (T) or not (F)?
#'@param trim Boolean, should the simulated terminal leverage ratios be trimmed? 
#'@param trimquants Vector indicating the quantiles to trim at if argument "trim" is set to TRUE. 
#'@return Plots the histogram with a vertical abline for a user-specified threshold and returns the simulated minimum leverage ratios.  
#'@export

simdevent <- function(N = 1000, len = 500, pars = c(-10, 0.95, 0.3, -12, 0.9, 0.2, 0.9, 1),
                      Ain = 80, Ein = 60, K = 60, r = 0.001, uv = 0.0005, ttv = 5, thd = 0.8, plot = T, trim = T, 
                      trimquants = c(0.01, 0.99), crisisret = -0.1, seed = NULL){
  
  
  set.seed(seed)
  
  sims <- replicate(n = N, simstrucsystem(len = len, pars = pars,
                                          Ain = Ain, Ein = Ein, K = K, r = r, uv = uv, ttv = ttv,
                                          crisis_sim = T, crisisret = crisisret))
  set.seed(NULL)
  
    sims <- sims[sims!=999]
    
    if (trim == T) {

      sims <-  sims[(sims > quantile(sims, trimquants[1])) & (sims < quantile(sims, trimquants[2]))]
    }
    
    
  if (plot == T){
  hist(sims, main = "Distribution of the terminal crisis leverage ratio.", xlab = "Leverage ratio", col = 5,
       density = 20,
       angle = 20)
    
  abline(v = thd, col ="red")
  legend("topleft", legend = "Threshold", pch = "|", col = "red")
  }
  
  
 return(sims)
  
} 


#'@title simintprob
#'@description Simulate the probability that leverage ratio lies in some range at the end of the period.
#'@importFrom MASS mvrnorm 
#'@param len Integer, the length of the simulated time period.
#'@param pars A vector with the parameter values of the structural multivariate stochastic volatility model.
#'@param Ain Numeric, initial asset value.
#'@param Ein Numeric, initial equity value.
#'@param K Numeric, debt value over the time period.
#'@param r Numeric, the risk-free rate of interest
#'@param uv Numeric, the daily unconditional variance of the assets over the period.
#'@param ttv Numeric, time to maturity of the debt over the period (measured in years).
#'@param crisis Boolean, only keeps paths where the market return is below some threshold during the period.
#'@param crisisret The return that defines a crisis and sets Boolean crisis to TRUE. 
#'@return The probability that the leverage ratio lies in the specified range at the end of the period.
#'@export

simintprob <- function(N = 1000, len = 500, pars = c(-10, 0.95, 0.3, -12, 0.9, 0.2, 0.9, 1),
                       Ain = 80, Ein = 60, K = 60, r = 0.001, uv = 0.0005, ttv = 5, lower = 0.4, upper = 0.6, seed = NULL,
                       crisis = F, crisisret = -0.1){
  
  set.seed(seed) 
  if(crisis == F){
  
  termlrs <- replicate(N, simstrucsystem(len = len, pars = pars,
                        Ain = Ain, Ein = Ein, K = K, r = r, uv = uv, ttv = ttv,
                        crisis_sim = F)$levrat[len])
  
  print(head(termlrs))
  }
  
  else{
    
    termlrs <- replicate(N, simstrucsystem(len = len, pars = pars,
                                           Ain = Ain, Ein = Ein, K = K, r = r, uv = uv, ttv = ttv,
                                           crisis_sim = T, crisisret = crisisret))
    
    print(head(termlrs))
    termlrs <- termlrs[termlrs!=999]
    
    N <- length(termlrs)
  }
  
  set.seed(NULL)
  
  prob <- sum(termlrs >= lower & termlrs <= upper) / N
  
  return(list(prob = prob, termlrs = termlrs))
}


